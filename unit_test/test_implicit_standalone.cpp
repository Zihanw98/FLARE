#include "test_support.h"
#include "time/implicit.h"
#include "hydro/riemann/exact.h"
#include "hydro/wave_timestep.h"
#include "user/user_api.h"
#include <fstream>
#include <limits>

namespace
{
class ZeroFlux : public hydro::RiemannSolver
{
    hd_core::Conserved flux(const hd_core::Primitives&, const hd_core::Primitives&,
                           int, double, double) const override
    {
        hydro::wave_timestep.face_speed = 0;
        return {};
    }
};
hd_core::MeshBlock mesh(double scale = 1.0)
{
    hd_core::MeshBlock result;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        result.nactive[d] = 3;
        result.nghost[d] = 2;
        result.bc_lower[d] = result.bc_upper[d] = boundary::BoundaryCondition::Periodic;
        result.set_linear_uniform_faces(d, 1.0, 2.0);
    }
    result.allocate(); result.build_geometry();
    for (auto& c : result.cells)
    {
        hd_core::Primitives w{};
        w[hd_core::RHO] = scale;
        w[hd_core::PRE] = 2*scale;
        c = hd_core::Cell(w, 1.4);
    }
    boundary::FillBoundary(result, 0, {});
    return result;
}
void source_tests()
{
    for (double scale : {1.0, 1e-18})
    {
        auto state = mesh(scale), old = state;
        std::size_t count = 0;
        double initial_mean = 0;
        for (std::size_t i = 0; i < state.cells.size(); ++i) if (state.is_active(i))
        {
            state.cells[i].U[hd_core::ENERGY] *= 1+0.05*count++;
            initial_mean += state.cells[i].U[hd_core::ENERGY];
        }
        initial_mean /= count;
        old = state;
        hd_source::SourceList sources{{"nonlocal cooling", hd_source::SourceKind::Thermal,
            [](const hd_core::MeshBlock& s, double time, const config::SimulationConfig&,
               std::vector<hd_core::Conserved>& rhs, double&) {
                CHECK(near(time, 0.1));
                double mean = 0; std::size_t n = 0;
                for (std::size_t i = 0; i < s.cells.size(); ++i) if (s.is_active(i))
                { mean += s.cells[i].U[hd_core::ENERGY]; ++n; }
                mean /= n;
                for (std::size_t i = 0; i < s.cells.size(); ++i) if (s.is_active(i))
                    rhs[i][hd_core::ENERGY] -= 30*s.cells[i].U[hd_core::ENERGY]+mean;
            }}};
        config::SimulationConfig cfg; cfg.time_integrator = "implicit";
        cfg.dt_implicit = 0.1; cfg.implicit_tolerance = 1e-9; cfg.implicit_adaptive = false;
        hd_implicit::Integrator integrator(state, cfg);
        double next = std::numeric_limits<double>::infinity();
        const auto advance = [&] {
            return hd_time::advance_one_step(state, 0, 0.1, hd_time::parse_time_integrator(cfg.time_integrator),
                hydro::SpatialScheme::Godunov, ZeroFlux{}, sources, cfg, {}, cfg.cfl, next,
                nullptr, nullptr, &integrator);
        };
        CHECK(advance() == 0.1);
        CHECK(next == 0.1);
        const double mean = initial_mean/4.1;
        for (std::size_t i = 0; i < state.cells.size(); ++i) if (state.is_active(i))
        {
            const double expected = (old.cells[i].U[hd_core::ENERGY]-0.1*mean)/4.0;
            CHECK(std::abs(state.cells[i].U[hd_core::ENERGY]/expected-1) < 1e-8);
        }
        // No positive solution exists; even modified trial ghosts must roll back.
        old = state;
        sources[0].func = [scale](const hd_core::MeshBlock& s, double,
            const config::SimulationConfig&, std::vector<hd_core::Conserved>& rhs, double&) {
            for (std::size_t i = 0; i < s.cells.size(); ++i)
                rhs[i][hd_core::RHO] -= 100*scale;
        };
        expect_throw<std::runtime_error>([&] {
            advance();
        });
        for (std::size_t i = 0; i < state.cells.size(); ++i)
            CHECK(state.cells[i].U.data == old.cells[i].U.data);
    }
}

void hydro_test(bool spherical)
{
    auto state = mesh();
    config::SimulationConfig cfg; cfg.reconstructor = "PLM_corrected"; cfg.limiter = "MC";
    if (spherical)
    {
        state.geometry = hd_core::Geometry::Spherical;
        state.build_geometry();
        cfg.spherical_gravity_enabled = true;
        cfg.spherical_self_gravity_enabled = true;
        cfg.spherical_central_mass = 0.1;
    }
    for (std::size_t i = 0; i < state.cells.size(); ++i) if (state.is_active(i))
        state.cells[i].U[hd_core::ENERGY] *= 1+0.03*state.coordinates_active(i)[0];
    // An intentionally nonlocal user boundary: every ghost reads the far end.
    state.bc_lower[0] = state.bc_upper[0] = boundary::BoundaryCondition::User;
    hd_user::UserHooks hooks;
    hooks.userbc = [](hd_core::MeshBlock& s, double, int d, boundary::BoundarySide side) {
        auto source = hd_core::Index{};
        source[d] = side == boundary::BoundarySide::Lower ? s.nactive[d]-1 : 0;
        const auto value = s.at_active(source);
        for (std::size_t i = 0; i < s.cells.size(); ++i)
        {
            const int x = s.coordinates_active(i)[d];
            if ((side == boundary::BoundarySide::Lower && x < 0)
                || (side == boundary::BoundarySide::Upper && x >= s.nactive[d])) s.cells[i] = value;
        }
    };
    boundary::FillBoundary(state, 0, hooks);
    const auto old = state;
    const auto sources = hd_source::make_sources(cfg, state, hooks);
    hydro::ExactSolver solver;
    cfg.time_integrator = "implicit"; cfg.dt_implicit = 0.001; cfg.implicit_tolerance = 1e-8;
    hd_implicit::Integrator integrator(state, cfg);
    double next = std::numeric_limits<double>::infinity();
    const auto dt = hd_time::advance_one_step(state, 0, cfg.dt_implicit,
        hd_time::parse_time_integrator(cfg.time_integrator), hydro::SpatialScheme::Godunov, solver, sources,
        cfg, hooks, cfg.cfl, next, nullptr, nullptr, &integrator);
    double rate, limit;
    const auto evaluated = hd_time::evaluate_RHS(state, dt, hydro::SpatialScheme::Godunov,
                                                sources, solver, cfg, cfg.cfl, rate, limit);
    for (std::size_t i = 0; i < state.cells.size(); ++i) if (state.is_active(i))
        for (int k = 0; k < hd_core::NVAR; ++k)
            CHECK(std::abs(state.cells[i].U[k]-old.cells[i].U[k]-dt*evaluated.second[i][k]) < 1e-6);
    CHECK(state.face == old.face);
}

void config_and_adaptive_tests(const std::string& path)
{
    auto read = [&](const std::string& settings) {
        std::ofstream out(path);
        const char* axes[] = {"x", "y", "z"};
        for (int d = 0; d < hd_core::NDIM; ++d)
            out << "n" << axes[d] << "a=3\nn" << axes[d] << "g=2\n"
                << axes[d] << "_lower=1\n" << axes[d] << "_upper=2\n";
        out << settings;
        out.close();
        return config::ReadParameters(path);
    };
    CHECK(read("").time_integrator == "euler"); // Old inputs need no implicit settings.
    auto cfg = read("time_integrator=ImPlIcIt\ndt_implicit=0.1\nimplicit_adaptive=true\n"
                    "implicit_tolerance=1e-9\nimplicit_newton_max_iterations=10\n");
    CHECK(hd_time::parse_time_integrator(cfg.time_integrator) == hd_time::TimeIntegrator::Implicit);
    CHECK(cfg.implicit_adaptive && cfg.dt_implicit == 0.1 && cfg.implicit_tolerance == 1e-9);
    CHECK(cfg.implicit_newton_max_iterations == 10);
    CHECK(!read("time_integrator=implicit\ndt_implicit=0.1\nimplicit_adaptive=false\n").implicit_adaptive);
    const auto diagnostic = read("implicit_verbose=true\nimplicit_report_file=solver.csv\n");
    CHECK(diagnostic.implicit_verbose && diagnostic.implicit_report_file == "solver.csv");
    for (const auto& invalid : {"", "dt_implicit=0\n", "dt_implicit=-1\n", "dt_implicit=nan\n",
         "dt_implicit=0.1\nimplicit_adaptive=maybe\n", "dt_implicit=0.1\nimplicit_tolerance=0\n",
         "dt_implicit=0.1\nimplicit_tolerance=1\n", "dt_implicit=0.1\nimplicit_newton_max_iterations=0\n"})
        expect_throw<std::runtime_error>([&] { read(std::string("time_integrator=implicit\n")+invalid); });

    // E'=20E: dt=0.1 has no positive backward Euler solution, 0.05 is singular,
    // and 0.025 is admissible. Verify shrinking, history and final-time clipping.
    auto state = mesh();
    hd_source::SourceList sources{{"growth", hd_source::SourceKind::Thermal,
        [](const hd_core::MeshBlock& s, double, const config::SimulationConfig&,
           std::vector<hd_core::Conserved>& rhs, double&) {
            for (std::size_t i = 0; i < s.cells.size(); ++i)
                rhs[i][hd_core::ENERGY] += 20*s.cells[i].U[hd_core::ENERGY];
        }}};
    hd_implicit::Integrator workspace(state, cfg);
    double next = std::numeric_limits<double>::infinity();
    auto advance = [&](double time, double remaining) {
        return hd_time::advance_one_step(state, time, remaining, hd_time::TimeIntegrator::Implicit,
            hydro::SpatialScheme::Godunov, ZeroFlux{}, sources, cfg, {}, cfg.cfl, next,
            nullptr, nullptr, &workspace);
    };
    const auto first = state.flat_index(hd_core::Index{}, true);
    double energy = state.cells[first].U[hd_core::ENERGY];
    double time = 0;
    for (int step = 0; step < 2; ++step)
    {
        const double desired = next;
        const double dt = advance(time, 1.0);
        CHECK(dt > 0 && dt < 0.05);
        if (step == 0) CHECK(near(dt, 0.025));
        else CHECK(near(dt, desired));
        energy /= 1-20*dt;
        CHECK(near(state.cells[first].U[hd_core::ENERGY], energy, 1e-8));
        CHECK(next > dt && next <= cfg.dt_implicit);
        time += dt;
    }
    CHECK(advance(time, 0.001) == 0.001);
}
}

int main(int argc, char** argv)
{
    return run_test("standalone implicit analytical cooling, low density, rollback, nonlocal boundary and gravity", [&] {
        CHECK(argc > 1);
        source_tests(); hydro_test(false);
        if constexpr (hd_core::NDIM == 1) hydro_test(true);
        config_and_adaptive_tests(std::string(argv[1])+".txt");
    });
}
