#include "test_support.h"
#include "user/user_api.h"
#include "source/source.h"
#include "time/integrator.h"
#include "hydro/riemann/riemann.h"
#include "hydro/wave_timestep.h"
#include <limits>

namespace
{
class ZeroFluxSolver : public hydro::RiemannSolver
{
public:
    hd_core::Conserved flux(const hd_core::Primitives&, const hd_core::Primitives&,
                           int, double, double) const override
    {
        hydro::wave_timestep.face_speed = 0.0;
        return {};
    }
};

hd_core::MeshBlock make_mesh()
{
    hd_core::MeshBlock mesh;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        mesh.nactive[d] = 3;
        mesh.nghost[d] = 1;
        mesh.set_linear_uniform_faces(d, 1.0, 2.0);
        mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
    }
    mesh.allocate();
    mesh.build_geometry();
    for (auto& cell : mesh.cells)
    {
        cell.gamma = 1.4;
        cell.U[hd_core::RHO] = 2.0;
        cell.U[hd_core::ENERGY] = 20.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
            cell.U[hd_core::MOMDEN(d)] = 0.5 * (d + 1);
    }
    return mesh;
}

void heating_tests()
{
    config::SimulationConfig cfg;
    auto mesh = make_mesh();
    cfg.heating_rate = 5.0;
    CHECK(hd_source::make_sources(cfg, mesh, {}).empty()); // Disabled by default.
    cfg.constant_heating_enabled = true;
    const auto sources = hd_source::make_sources(cfg, mesh, {});
    CHECK(sources.size() == 1 && sources[0].name == "constant_heating");
    CHECK(sources[0].kind == hd_source::SourceKind::Thermal);
    cfg.heating_rate = 500.0; // The source must retain its captured rate of 5.
    const auto cold = mesh.flat_index(hd_core::Index{}, true);
    std::vector<bool> active(mesh.cells.size(), true);
    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
    {
        const auto index = mesh.coordinates_active(flat);
        for (int d = 0; d < hd_core::NDIM; ++d)
            active[flat] = active[flat] && index[d] >= 0 && index[d] < mesh.nactive[d];
        if (!active[flat])
        {
            // Uninitialized ghosts must neither be heated nor restrict dt.
            mesh.cells[flat].U.data.fill(std::numeric_limits<double>::quiet_NaN());
            continue;
        }
        hd_core::Primitives state{};
        state[hd_core::RHO] = 2.0;
        state[hd_core::PRE] = flat == cold ? 0.4 : 4.0;
        state[hd_core::VEL(0)] = 3.0; // Kinetic energy must not inflate the timescale.
        mesh.cells[flat] = hd_core::Cell(state, 1.4);
    }
    std::vector<hd_core::Conserved> rhs(mesh.cells.size());
    for (auto& value : rhs) value.data.fill(7.0);
    double limit = std::numeric_limits<double>::infinity();
    sources[0].func(mesh, 0.0, cfg, rhs, limit);
    CHECK(near(limit, 0.01)); // 0.1 * (0.4 / 0.4) / (2 * 5).
    for (std::size_t flat = 0; flat < rhs.size(); ++flat)
        for (int k = 0; k < hd_core::NVAR; ++k)
            CHECK(near(rhs[flat][k], 7.0 + (active[flat] && k == hd_core::ENERGY ? 10.0 : 0.0)));
    limit = 0.001;
    sources[0].func(mesh, 0.0, cfg, rhs, limit);
    CHECK(limit == 0.001); // Never loosen a tighter existing bound.
    const auto before = rhs;
    hd_source::make_constant_heating(0.0).func(mesh, 0.0, cfg, rhs, limit);
    CHECK(limit == 0.001);
    for (std::size_t flat = 0; flat < rhs.size(); ++flat)
        CHECK(rhs[flat].data == before[flat].data);

    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        if (active[flat]) mesh.cells[flat].U = {}; // Vacuum: no 0/0 timestep.
    limit = std::numeric_limits<double>::infinity();
    sources[0].func(mesh, 0.0, cfg, rhs, limit);
    CHECK(std::isinf(limit));
    for (std::size_t flat = 0; flat < rhs.size(); ++flat)
        CHECK(rhs[flat].data == before[flat].data);
    for (double bad : {-1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()})
        expect_throw<std::invalid_argument>([&] { hd_source::make_constant_heating(bad); });

    // With zero hydrodynamic wave speed, the heating timescale selects the step.
    ZeroFluxSolver solver;
    const hd_source::SourceList heating{hd_source::make_constant_heating(50.0)};
    for (auto method : {hd_time::TimeIntegrator::Euler, hd_time::TimeIntegrator::RK2})
    {
        auto trial = make_mesh();
        const auto initial = trial.at_active(hd_core::Index{}).U;
        double kinetic = 0.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
            kinetic += 0.5 * initial[hd_core::MOMDEN(d)] * initial[hd_core::MOMDEN(d)]
                / initial[hd_core::RHO];
        const double expected_dt = 0.1 * (20.0 - kinetic) / 100.0;
        double next = std::numeric_limits<double>::infinity();
        const double dt = hd_time::advance_one_step(trial, 0.0, 1.0, method,
            hydro::SpatialScheme::Godunov, solver, heating, cfg, {}, cfg.cfl, next);
        CHECK(near(dt, expected_dt));
        for (const auto& cell : trial.cells)
            for (int k = 0; k < hd_core::NVAR; ++k)
                CHECK(near(cell.U[k], initial[k] + (k == hd_core::ENERGY ? 100.0 * dt : 0.0)));
    }
}
}

int main()
{
    return run_test("source parameters, built-ins and time integration", [] {
        heating_tests();
        auto mesh = make_mesh();
        config::SimulationConfig cfg;
        hd_user::UserHooks hooks;
        CHECK(hd_source::make_sources(cfg, mesh, hooks).empty());

        // The list owns its parameters after the input configuration dies.
        auto sources = [&] {
            config::SimulationConfig local;
            local.constant_gravity_enabled = true;
            for (int d = 0; d < hd_core::NDIM; ++d)
                local.gravity[d] = -0.25 * (d + 1);
            return hd_source::make_sources(local, mesh, hooks);
        }();
        CHECK(sources.size() == 1 && sources[0].name == "constant_gravity");
        CHECK(sources[0].kind == hd_source::SourceKind::Force);
        cfg.gravity.fill(100.0); // Per-call config cannot replace owned parameters.
        std::vector<hd_core::Conserved> rhs(mesh.cells.size());
        for (auto& value : rhs) value.data.fill(1.0);
        double limit = 0.125;
        sources[0].func(mesh, 0.0, cfg, rhs, limit);
        CHECK(limit == 0.125);
        double energy_rate = 0.0, g_squared = 0.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const double g = -0.25 * (d + 1);
            energy_rate += 0.5 * (d + 1) * g;
            g_squared += g * g;
        }
        for (const auto& value : rhs)
        {
            CHECK(value[hd_core::RHO] == 1.0);
            CHECK(near(value[hd_core::ENERGY], 1.0 + energy_rate));
            for (int d = 0; d < hd_core::NDIM; ++d)
                CHECK(near(value[hd_core::MOMDEN(d)], 1.0 - 0.5 * (d + 1)));
        }

        // Instances keep independent values, including after copying a list.
        auto copied = sources;
        hd_source::ConstantGravityParams params;
        params.acceleration.fill(0.75);
        copied.push_back(hd_source::make_constant_gravity(params));
        params.acceleration.fill(50.0);
        rhs.assign(mesh.cells.size(), {});
        for (const auto& term : copied) term.func(mesh, 0.0, cfg, rhs, limit);
        for (int d = 0; d < hd_core::NDIM; ++d)
            CHECK(near(rhs[0][hd_core::MOMDEN(d)], 1.5 - 0.5 * (d + 1)));

        // Check gravity evolution through both integrators against analytic rates.
        ZeroFluxSolver solver;
        for (auto integrator : {hd_time::TimeIntegrator::Euler, hd_time::TimeIntegrator::RK2})
        {
            auto trial = mesh;
            double next = std::numeric_limits<double>::infinity();
            const double dt = hd_time::advance_one_step(trial, 0.0, 0.01, integrator,
                hydro::SpatialScheme::Godunov, solver, sources, cfg, hooks, cfg.cfl, next);
            CHECK(near(dt, 0.01));
            const double quadratic = integrator == hd_time::TimeIntegrator::RK2
                ? g_squared * dt * dt : 0.0;
            for (const auto& cell : trial.cells)
            {
                CHECK(near(cell.U[hd_core::ENERGY], 20.0 + energy_rate * dt + quadratic));
                for (int d = 0; d < hd_core::NDIM; ++d)
                    CHECK(near(cell.U[hd_core::MOMDEN(d)], (0.5 - 0.5 * dt) * (d + 1)));
            }
        }

        // Custom capture: stage-time-dependent rates and an independent timestep cap.
        const auto custom = [] {
            const std::vector<double> coefficients{2.0, 3.0};
            return hd_source::SourceTerm{"custom", hd_source::SourceKind::Thermal,
                [coefficients](const auto&, double time, const auto&, auto& rates, double& cap) {
                    for (auto& value : rates)
                        value[hd_core::ENERGY] += coefficients[0] + coefficients[1] * time;
                    cap = std::min(cap, 0.05);
                }};
        }();
        double wave_rate, source_limit;
        const auto result = hd_time::evaluate_RHS(mesh, 2.0, hydro::SpatialScheme::Godunov,
            {custom}, solver, cfg, cfg.cfl, wave_rate, source_limit);
        CHECK(near(result.first, 0.05) && near(source_limit, 0.05));
        CHECK(near(result.second[0][hd_core::ENERGY], 8.0));
        expect_throw<std::invalid_argument>([&] {
            hd_time::evaluate_RHS(mesh, 0.0, hydro::SpatialScheme::Godunov,
                {{"empty", hd_source::SourceKind::General, {}}}, solver, cfg,
                cfg.cfl, wave_rate, source_limit);
        });

        // One user hook can append multiple parameter-owning sources after built-ins.
        auto user_cfg = cfg;
        user_cfg.constant_gravity_enabled = true;
        user_cfg.gravity.fill(0.0);
        hd_user::UserHooks user_hooks;
        user_hooks.register_sources = [](const config::SimulationConfig& input,
                                        const hd_core::MeshBlock& input_mesh,
                                        hd_source::SourceList& registered) {
            CHECK(input.constant_gravity_enabled);
            CHECK(input_mesh.geometry == hd_core::Geometry::Cartesian);
            CHECK(registered.size() == 1 && registered[0].name == "constant_gravity");
            for (double heating : {2.0, 3.0})
                registered.push_back({"user_heating", hd_source::SourceKind::Thermal,
                    [heating](const auto&, double, const auto&, auto& rates, double& cap) {
                        for (auto& value : rates) value[hd_core::ENERGY] += heating;
                        cap = std::min(cap, 0.05);
                    }});
        };
        const auto registered = hd_source::make_sources(user_cfg, mesh, user_hooks);
        CHECK(registered.size() == 3);
        const auto user_result = hd_time::evaluate_RHS(mesh, 0.0,
            hydro::SpatialScheme::Godunov, registered, solver, user_cfg,
            user_cfg.cfl, wave_rate, source_limit);
        CHECK(near(user_result.first, 0.05));
        for (const auto& value : user_result.second)
            CHECK(near(value[hd_core::ENERGY], 5.0));

        // Geometry supports spherical 1D. A different stage mesh must supply pressure.
        if constexpr (hd_core::NDIM == 1)
        {
            mesh.geometry = hd_core::Geometry::Spherical;
            mesh.build_geometry();
            auto geometric = hd_source::make_sources(cfg, mesh, hooks);
            CHECK(geometric.size() == 1 && geometric[0].name == "geometric");
            cfg.constant_gravity_enabled = true;
            auto both = hd_source::make_sources(cfg, mesh, hooks);
            CHECK(both.size() == 2 && both[0].name == "constant_gravity"
                  && both[1].name == "geometric");
            auto stage = mesh;
            for (auto& cell : stage.cells)
            {
                cell.U[hd_core::MOMDEN(0)] = 0.0;
                cell.U[hd_core::ENERGY] = 12.5; // p = 5.
            }
            rhs.assign(stage.cells.size(), {});
            geometric[0].func(stage, 0.0, cfg, rhs, limit);
            CHECK(limit == 0.125);
            for (std::size_t i = 0; i < stage.cells.size(); ++i)
            {
                const int x = stage.coordinates(i)[0];
                const bool active = x >= stage.nghost[0]
                    && x < stage.nghost[0] + stage.nactive[0];
                const double left = stage.face[0][x], right = stage.face[0][x + 1];
                const double expected = active ? 15.0 * (right * right - left * left)
                    / (right * right * right - left * left * left) : 0.0;
                CHECK(near(rhs[i][hd_core::MOMDEN(0)], expected));
                CHECK(rhs[i][hd_core::RHO] == 0.0 && rhs[i][hd_core::ENERGY] == 0.0);
            }
        }
    });
}
