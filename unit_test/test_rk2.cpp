#include "test_support.h"
#include "time/integrator.h"
#include "hydro/wave_timestep.h"
#include "hydro/riemann/riemann.h"
#include "user/user_api.h"
#include <fstream>
#include <limits>

namespace
{
// Zero flux isolates source evolution. The prescribed wave speed gives rate=rho
// on our mesh, so the stage CFL can be varied independently of the exact update.
class TestSolver : public hydro::RiemannSolver
{
public:
    bool zero_speed = false;
    hd_core::Conserved flux(const hd_core::Primitives& left,
                           const hd_core::Primitives& right, int, double, double) const override
    {
        hydro::wave_timestep.face_speed = zero_speed ? 0.0
            : std::max(left[hd_core::RHO], right[hd_core::RHO]) / (2.0 * hd_core::NDIM);
        return {};
    }
};

double density_source = 0.0;
double energy_source = 0.0;
bool time_dependent_energy = false;
double initial_source_limit = std::numeric_limits<double>::infinity();
double stage_source_limit = std::numeric_limits<double>::infinity();
int source_calls = 0;

void source(const hd_core::MeshBlock& mesh, double time,
            const config::SimulationConfig&, std::vector<hd_core::Conserved>& rhs,
            double& dt_limit)
{
    ++source_calls;
    // This also checks that stage ghosts were refilled before evaluating the RHS.
    const auto& reference = mesh.at_active(hd_core::Index{});
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
    {
        CHECK(near(mesh.cells[i].U[hd_core::RHO], reference.U[hd_core::RHO]));
        rhs[i][hd_core::RHO] += density_source;
        rhs[i][hd_core::ENERGY] += time_dependent_energy
            ? (time == 0.0 ? 0.0 : energy_source) : energy_source;
    }
    dt_limit = time == 0.0 ? initial_source_limit : stage_source_limit;
}

void loose_source(const hd_core::MeshBlock&, double,
                  const config::SimulationConfig&, std::vector<hd_core::Conserved>&,
                  double& dt_limit)
{
    dt_limit = 100.0; // Must not erase a previous source's tighter bound.
}

hd_core::MeshBlock make_mesh()
{
    source_calls = 0;
    density_source = energy_source = 0.0;
    time_dependent_energy = false;
    initial_source_limit = stage_source_limit = std::numeric_limits<double>::infinity();
    hd_core::MeshBlock mesh;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        mesh.nactive[d] = 2;
        mesh.nghost[d] = 1;
        mesh.set_uniform_faces(d, 0.0, 1.0);
        mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
    }
    mesh.allocate();
    mesh.build_geometry();
    for (auto& cell : mesh.cells)
    {
        cell.U[hd_core::RHO] = 1.0;
        cell.U[hd_core::ENERGY] = 10.0;
        cell.gamma = 1.4;
    }
    return mesh;
}

void same_cells(const hd_core::MeshBlock& mesh, const hd_core::MeshBlock& original)
{
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
    {
        CHECK(mesh.cells[i].U.data == original.cells[i].U.data);
        CHECK(mesh.cells[i].gamma == original.cells[i].gamma);
    }
}
}

int main(int argc, char** argv)
{
    return run_test("RK2 CFL targets, source limits and rejection", [&] {
        config::SimulationConfig cfg;
        hd_user::UserHooks hooks;
        TestSolver solver;
        hd_source::SourceList sources{{"test", hd_source::SourceKind::General, source}};
        const auto scheme = hydro::SpatialScheme::Godunov;
        auto advance = [&](hd_core::MeshBlock& mesh, double& dt_next,
                           double limit = 1.0) {
            return hd_time::advance_one_step(mesh, 0.0, limit, hd_time::TimeIntegrator::RK2,
                scheme, solver, sources, cfg, hooks, cfg.cfl, dt_next);
        };

        // Above desired but below/equal max is accepted, including tiny perturbations.
        for (double growth : {1e-12, 0.1, 0.625})
        {
            auto mesh = make_mesh();
            density_source = growth;
            double dt_next = std::numeric_limits<double>::infinity();
            CHECK(near(advance(mesh, dt_next), 0.4));
            CHECK(source_calls == 2); // No rejected attempt.
            CHECK(near(dt_next, 0.4 / (1.0 + 0.4 * growth)));
            for (const auto& cell : mesh.cells)
                CHECK(near(cell.U[hd_core::RHO], 1.0 + 0.4 * growth));
            const double suggested = dt_next;
            CHECK(near(advance(mesh, dt_next), suggested));
        }

        // Genuine CFL excess: scale by desired/observed, not by a fixed half.
        auto mesh = make_mesh();
        density_source = 2.0;
        double dt_next = std::numeric_limits<double>::infinity();
        const double dt = advance(mesh, dt_next);
        CHECK(near(dt, 0.4 / 1.8));
        CHECK(source_calls == 3); // One rejected CFL attempt.
        for (const auto& cell : mesh.cells)
            CHECK(near(cell.U[hd_core::RHO], 1.0 + 2.0 * dt));

        // Exhaustion never commits a stage or a new suggested timestep.
        mesh = make_mesh();
        density_source = 2.0;
        const auto original = mesh;
        cfg.integrator_retry_max = 0;
        dt_next = std::numeric_limits<double>::infinity();
        expect_throw<std::runtime_error>([&] { advance(mesh, dt_next); });
        same_cells(mesh, original);
        CHECK(std::isinf(dt_next));
        cfg.integrator_retry_max = 20;

        // A stage source cap is hard even when the CFL is acceptable.
        mesh = make_mesh();
        density_source = 1.0;
        initial_source_limit = 0.3;
        stage_source_limit = 0.1;
        dt_next = std::numeric_limits<double>::infinity();
        CHECK(near(advance(mesh, dt_next), 0.09));
        CHECK(source_calls == 3); // One rejected source-limit attempt.
        CHECK(dt_next <= 0.1);
        CHECK(near(mesh.at_active(hd_core::Index{}).U[hd_core::RHO], 1.09));
        sources.push_back({"loose", hd_source::SourceKind::General, loose_source});
        double rate, source_limit;
        const auto evaluated = hd_time::evaluate_RHS(mesh, 0.0, scheme, sources, solver, cfg,
                                                      cfg.cfl, rate, source_limit);
        CHECK(near(source_limit, 0.3));
        CHECK(near(evaluated.first, 0.3));
        CHECK(evaluated.second.size() == mesh.cells.size());
        sources.pop_back();

        // Invalid first stage and invalid final combination both retry from U^n.
        for (bool final_only : {false, true})
        {
            mesh = make_mesh();
            energy_source = final_only ? -60.0 : -30.0;
            time_dependent_energy = final_only;
            dt_next = std::numeric_limits<double>::infinity();
            CHECK(near(advance(mesh, dt_next), 0.2));
            CHECK(source_calls == (final_only ? 3 : 2));
            CHECK(near(mesh.at_active(hd_core::Index{}).U[hd_core::ENERGY], 4.0));
        }

        // External clipping must not poison the next suggested timestep.
        mesh = make_mesh();
        dt_next = std::numeric_limits<double>::infinity();
        CHECK(near(advance(mesh, dt_next, 0.01), 0.01));
        CHECK(near(dt_next, 0.4));
        CHECK(near(advance(mesh, dt_next), 0.4));
        solver.zero_speed = true;
        dt_next = std::numeric_limits<double>::infinity();
        CHECK(near(advance(mesh, dt_next, 0.125), 0.125));
        CHECK(std::isinf(dt_next));
        initial_source_limit = std::numeric_limits<double>::quiet_NaN();
        expect_throw<std::invalid_argument>([&] { advance(mesh, dt_next); });

        // Parse new names and old input alias, with validation independent of case.
        CHECK(argc > 1);
        const std::string path = std::string(argv[1]) + ".txt";
        auto read_config = [&](const std::string& settings) {
            std::ofstream out(path);
            const char* axes[] = {"x", "y", "z"};
            for (int d = 0; d < hd_core::NDIM; ++d)
                out << "n" << axes[d] << "a=2\nn" << axes[d] << "g=1\n"
                    << axes[d] << "_lower=0\n" << axes[d] << "_upper=1\n";
            out << "time_integrator=RK2\n" << settings;
            out.close();
            return config::ReadParameters(path);
        };
        CHECK(near(read_config("cfl=0.3\n").cfl, 0.3));
        const auto parsed = read_config("cfl_desired=0.35\ncfl_max=0.45\n");
        CHECK(near(parsed.cfl, 0.35) && near(parsed.cfl_max, 0.45));
        for (const std::string settings : {"cfl=0.3\ncfl_desired=0.3\n",
             "cfl_desired=0.45\ncfl_max=0.4\n", "cfl_max=0.6\n",
             "cfl_desired=0\n", "integrator_retry_max=-1\n"})
            expect_throw<std::runtime_error>([&] { read_config(settings); });
    });
}
