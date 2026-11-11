#include "test_support.h"
#include "time/integrator.h"
#include "hydro/riemann/exact.h"
#include "user/user_api.h"
#include <limits>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <array>
#include <cstdint>
#include <cstring>
#include <random>
namespace {
using hd_core::MeshBlock;
using hd_core::Conserved;
constexpr double infinity = std::numeric_limits<double>::infinity();
MeshBlock make_mesh(bool nonuniform = true, int narrow = 7, bool spherical = false,
                    int ghosts = 2, int active_cells = 20, double narrow_width = 0.01)
{
    MeshBlock mesh;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        mesh.nactive[d] = d == 0 ? active_cells : 1;
        mesh.nghost[d] = ghosts;
        mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
        mesh.set_linear_uniform_faces(d, 1.0, 3.0);
    }
    if (nonuniform)
    {
        auto& faces = mesh.face[0];
        const int n = mesh.nactive[0];
        faces[ghosts] = 1.0;
        for (int i = 0; i < n; ++i)
            faces[ghosts + i + 1] = faces[ghosts + i] + (i == narrow ? narrow_width : 0.1);
        for (int g = 1; g <= ghosts; ++g)
        {
            const int low = spherical ? ghosts + g - 1 : ghosts + n - g;
            const int high = spherical ? ghosts + n - g : ghosts + g - 1;
            faces[ghosts - g] = faces[ghosts - g + 1] - (faces[low + 1] - faces[low]);
            faces[ghosts + n + g] = faces[ghosts + n + g - 1] + (faces[high + 1] - faces[high]);
        }
    }
    if (spherical)
    {
        mesh.geometry = hd_core::Geometry::Spherical;
        mesh.bc_lower[0] = mesh.bc_upper[0] = boundary::BoundaryCondition::Reflective;
    }
    mesh.allocate();
    mesh.build_geometry();
    for (auto& cell : mesh.cells)
    {
        hd_core::Primitives W{};
        W[hd_core::RHO] = 1.0;
        W[hd_core::PRE] = 1.0;
        W[hd_core::VEL(0)] = spherical ? 0.0 : 0.2;
        cell = hd_core::Cell(W, 1.4);
    }
    boundary::FillBoundary(mesh, 0.0, {});
    return mesh;
}

Conserved integral(const MeshBlock& mesh)
{
    Conserved result{};
    for (std::size_t k = 0; k < mesh.cells.size(); ++k)
        if (mesh.is_active(k)) result += mesh.cells[k].U * mesh.volume_cell[k];
    return result;
}

void same(const MeshBlock& a, const MeshBlock& b, double tolerance = 2e-13)
{
    CHECK(a.cells.size() == b.cells.size());
    for (std::size_t k = 0; k < a.cells.size(); ++k)
        for (int v = 0; v < hd_core::NVAR; ++v)
            if (!near(a.cells[k].U[v], b.cells[k].U[v], tolerance))
            {
                std::ostringstream message;
                message << std::setprecision(17) << "State mismatch: cell=" << k
                    << " component=" << v << " actual=" << a.cells[k].U[v]
                    << " reference=" << b.cells[k].U[v] << " tolerance=" << tolerance;
                throw std::runtime_error(message.str());
            }
}

void contact(MeshBlock& mesh)
{
    const int ng = mesh.nghost[0];
    const double lo = mesh.face[0][ng], length = mesh.active_length(0);
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        hd_core::Primitives W{};
        const double x = mesh.cell_center(0, i + ng);
        W[hd_core::RHO] = 1.0 + 0.2 * std::sin(6.283185307179586 * (x - lo) / length);
        W[hd_core::PRE] = 1.0;
        W[hd_core::VEL(0)] = 0.2;
        mesh.cells[i + ng] = hd_core::Cell(W, 1.4);
    }
    boundary::FillBoundary(mesh, 0.0, {});
}

// Constant-speed transport isolates temporal error on a fixed spatial operator.
class Transport : public hydro::RiemannSolver
{
public:
    mutable int calls = 0;
    bool zero_flux = false;
    int fail_at = -1;
    int jump_after = -1;
    Conserved flux(const hd_core::Primitives& L, const hd_core::Primitives&,
                   int, double gamma, double, double& face_speed, double& face_velocity) const override
    {
        ++calls;
        if (calls == fail_at) throw std::runtime_error("injected flux failure");
        face_velocity = 0.0;
        face_speed = jump_after >= 0 && calls > jump_after ? 4.0 : 1.0;
        return zero_flux ? Conserved{} : L.to_Conserved(gamma);
    }
    Conserved flux_with_pressure_hint(const hd_core::Primitives& L, const hd_core::Primitives& R,
                   int d, double gL, double gR, double& hint, double& face_speed, double& face_velocity) const override
    {
        hint += 1.0;
        return flux(L, R, d, gL, gR, face_speed, face_velocity);
    }
};

double step(MeshBlock& mesh, const hydro::RiemannSolver& solver, config::SimulationConfig cfg,
            double time, double limit, const hd_source::SourceList& sources = {},
            hd_time::TimeIntegrator method = hd_time::TimeIntegrator::LTSRK2,
            hydro::StarPressureCache* hints = nullptr,
            const hd_user::UserHooks& hooks = {})
{
    double next = infinity;
    return hd_time::advance_one_step(mesh, time, limit, method, hydro::SpatialScheme::Godunov,
        solver, sources, cfg, hooks, cfg.cfl, next, nullptr, hints);
}

void evolve(MeshBlock& mesh, const hydro::RiemannSolver& solver, config::SimulationConfig cfg,
            double end, const hd_source::SourceList& sources = {},
            hd_time::TimeIntegrator method = hd_time::TimeIntegrator::LTSRK2, double cap = 1.0)
{
    double t = 0.0;
    while (t < end) t += step(mesh, solver, cfg, t, std::min(cap, end - t), sources, method);
}

double error(const MeshBlock& a, const MeshBlock& b)
{
    double result = 0.0;
    for (int i = 0; i < a.nactive[0]; ++i)
    {
        const int k = i + a.nghost[0];
        result += a.volume_cell[k] * std::abs(a.cells[k].U[hd_core::RHO] - b.cells[k].U[hd_core::RHO]);
    }
    return result;
}

void equivalence(const std::string& reconstruction)
{
    hydro::ExactSolver solver(1e-10, 200);
    config::SimulationConfig cfg;
    cfg.cfl = 0.2; cfg.reconstructor = reconstruction; cfg.limiter = "MC";
    for (bool uniform : {false, true})
    {
        auto a = make_mesh(!uniform);
        contact(a);
        auto b = a;
        // All-level-zero truncation must reproduce global Heun, including sources.
        hd_source::ConstantGravityParams g{};
        g.acceleration[0] = 0.02;
        auto sources = hd_source::SourceList{hd_source::make_constant_gravity(g)};
        const double h = step(a, solver, cfg, 0.0, 1e-4, sources);
        CHECK(h == step(b, solver, cfg, 0.0, h, sources, hd_time::TimeIntegrator::RK2));
        same(a, b);
        cfg.lts_max_levels = 1;
        a = make_mesh(!uniform); contact(a); b = a;
        const double H = step(a, solver, cfg, 0.0, 1.0);
        CHECK(near(H, step(b, solver, cfg, 0.0, H, {}, hd_time::TimeIntegrator::RK2)));
        same(a, b);
        cfg.lts_max_levels = 4;
    }
    // Uniform widths and fixed speeds: every cell takes eight fine RK2 steps.
    Transport transport;
    auto a = make_mesh(false); contact(a); auto b = a;
    cfg.lts_max_levels = 4;
    const double H = step(a, transport, cfg, 0.0, 1.0);
    evolve(b, transport, cfg, H, {}, hd_time::TimeIntegrator::RK2, H / 8.0);
    same(a, b);
}

void convergence(const std::string& reconstruction, bool with_source)
{
    auto initial = make_mesh(true, 31, false, 2, 64, 0.001);
    contact(initial);
    Transport solver;
    config::SimulationConfig cfg;
    cfg.lts_max_levels = 5; cfg.cfl = 0.2;
    cfg.reconstructor = reconstruction; cfg.limiter = "MC";
    hd_source::SourceList sources;
    if (with_source)
        sources.push_back({"nonlocal time-dependent growth", hd_source::SourceKind::General,
            [](const MeshBlock& mesh, double t, const config::SimulationConfig&,
               std::vector<Conserved>& rhs, double&) {
                const auto mean = integral(mesh) / mesh.active_length(0);
                for (std::size_t k = 0; k < mesh.cells.size(); ++k)
                    if (mesh.is_active(k)) rhs[k] += (0.1 + t) * (0.5 * mean + 0.5 * mesh.cells[k].U);
            }});
    auto reference = initial;
    evolve(reference, solver, cfg, 0.064, sources, hd_time::TimeIntegrator::RK2, 0.0000125);
    double previous = 0.0;
    for (double cfl : {0.2, 0.1, 0.05})
    {
        auto current = initial; cfg.cfl = cfl;
        evolve(current, solver, cfg, 0.064, sources);
        const double e = error(current, reference);
        std::cout << reconstruction << " sources=" << with_source << " CFL=" << cfl
            << " error=" << std::setprecision(12) << e;
        if (previous > 0.0)
        {
            std::cout << " ratio=" << previous / e;
            CHECK(previous / e > 3.5);
            CHECK(previous / e < 4.6);
        }
        std::cout << '\n';
        previous = e;
    }
}

void conservation_and_geometry(const std::string& reconstruction)
{
    hydro::ExactSolver solver(1e-10, 200);
    config::SimulationConfig cfg;
    cfg.cfl = 0.2; cfg.lts_max_levels = 5; cfg.reconstructor = reconstruction; cfg.limiter = "MC";
    for (int narrow : {0, 31, 63})
    {
        auto mesh = make_mesh(true, narrow, false, 2, 64, 0.001);
        contact(mesh);
        const auto before = integral(mesh);
        evolve(mesh, solver, cfg, 0.03);
        const auto after = integral(mesh);
        for (int v = 0; v < hd_core::NVAR; ++v) CHECK(near(before[v], after[v], 2e-12));
    }
    auto spherical = make_mesh(true, 7, true);
    const double mass = integral(spherical)[hd_core::RHO];
    evolve(spherical, solver, cfg, 0.02, {hd_source::make_geometric_source(),
        hd_source::make_spherical_self_gravity({0.01, true})});
    CHECK(near(integral(spherical)[hd_core::RHO], mass, 2e-12));
    for (const auto& c : spherical.cells) c.U.validate(c.gamma);
    Transport fast, global;
    auto a = make_mesh(true, 31, false, 2, 64, 0.001); contact(a); auto b = a;
    evolve(a, fast, cfg, 0.032);
    evolve(b, global, cfg, 0.032, {}, hd_time::TimeIntegrator::RK2);
    std::cout << reconstruction << " face solves local/global=" << fast.calls << '/' << global.calls << '\n';
    CHECK(fast.calls < global.calls / 2);
}

void retries_and_sources()
{
    config::SimulationConfig cfg;
    cfg.cfl = 0.2; cfg.lts_max_levels = 4;
    auto mesh = make_mesh(true, 7); auto before = mesh;
    Transport solver; solver.zero_flux = true; solver.jump_after = mesh.nactive[0];
    // Increased stage waves reject the first trial; accepted zero-flux state is exact.
    const double accepted = step(mesh, solver, cfg, 0.0, 1.0);
    CHECK(accepted < 8 * cfg.cfl * 0.01);
    same(mesh, before);
    auto reference = before;
    Transport control; control.zero_flux = true; control.jump_after = reference.nactive[0];
    CHECK(near(accepted, step(reference, control, cfg, 0.0, accepted)));
    same(mesh, reference);
    // Failure after a partially advanced subtree must leave BOTH mesh and hints intact.
    hydro::StarPressureCache hints;
    hints[0].assign(mesh.nactive[0] + 1, 3.0);
    auto old_hints = hints;
    Transport failure; failure.fail_at = 2 * mesh.nactive[0];
    expect_throw<std::runtime_error>([&] { step(mesh, failure, cfg, 0.0, 1.0, {},
        hd_time::TimeIntegrator::LTSRK2, &hints); });
    same(mesh, before, 0.0); CHECK(hints == old_hints);
    cfg.integrator_retry_max = 0;
    Transport rejected; rejected.zero_flux = true; rejected.jump_after = mesh.nactive[0];
    expect_throw<std::runtime_error>([&] { step(mesh, rejected, cfg, 0.0, 1.0); });
    same(mesh, before, 0.0);
    cfg.integrator_retry_max = 20;
    // Stage-time hard source bound must trigger a retry, independent of CFL.
    Transport source_solver; source_solver.zero_flux = true;
    int stages = 0;
    hd_source::SourceList sources{{"cap", hd_source::SourceKind::General,
        [&](const MeshBlock&, double t, const config::SimulationConfig&, std::vector<Conserved>&, double& cap) {
            ++stages; cap = t > 0.0 ? 0.002 : 1.0;
        }}};
    CHECK(step(mesh, source_solver, cfg, 0.0, 0.01, sources) <= 0.002);
    CHECK(stages > 2); same(mesh, before);
    // Time-dependent source integral: trapezoidal integration is exact for S_E=t.
    sources = {{"time heating", hd_source::SourceKind::Thermal,
        [](const MeshBlock& m, double t, const config::SimulationConfig&, std::vector<Conserved>& rhs, double&) {
            for (std::size_t k = 0; k < m.cells.size(); ++k) if (m.is_active(k)) rhs[k][hd_core::ENERGY] += t;
        }}};
    evolve(mesh, source_solver, cfg, 0.02, sources);
    for (int i = 0; i < mesh.nactive[0]; ++i)
        CHECK(near(mesh.cells[i+2].U[hd_core::ENERGY] - before.cells[i+2].U[hd_core::ENERGY], 0.0002, 2e-13));
}

void source_cache_stage_tests()
{
    Transport solver; solver.zero_flux = true;
    config::SimulationConfig cfg; cfg.cfl = 0.2; cfg.lts_max_levels = 5;
    auto mesh = make_mesh(false);
    const auto before = mesh;
    auto energy_source = [](double rate) -> hd_source::SourceTerm {
        return {"constant energy rate", hd_source::SourceKind::Thermal,
            [rate](const MeshBlock& m, double, const config::SimulationConfig&,
                   std::vector<Conserved>& rhs, double&) {
                for (std::size_t k = 0; k < m.cells.size(); ++k)
                    if (m.is_active(k)) rhs[k][hd_core::ENERGY] += rate;
            }};
    };
    // The zero source fills the cache; the second source must still be checked.
    // Initially internal energy is 2.5 and the cooling cap is 0.25. The first
    // evolved stage has a smaller cap and must reject H=0.25. Reusing the old
    // stage's thermodynamics would incorrectly accept it. The retry must refill
    // the cache again and accept H=0.125.
    const double h = step(mesh, solver, cfg, 0.0, 0.25,
                          {energy_source(0.0), energy_source(-1.0)});
    CHECK(h == 0.125);
    for (int i = 0; i < mesh.nactive[0]; ++i)
        CHECK(near(mesh.cells[i+2].U[hd_core::ENERGY], before.cells[i+2].U[hd_core::ENERGY] - h, 2e-13));

    // Equal/opposite sources cannot hide their independent 0.25 limits.
    mesh = before;
    CHECK(near(step(mesh, solver, cfg, 0.0, 1.0,
        {energy_source(1.0), energy_source(-1.0)}), 0.25, 2e-13));
    same(mesh, before);

    // Cached thermodynamics do not bypass per-source validity/vacuum checks.
    mesh = before;
    expect_throw<std::runtime_error>([&] {
        step(mesh, solver, cfg, 0.0, 1.0,
             {energy_source(0.0), energy_source(std::numeric_limits<double>::quiet_NaN())});
    });
    same(mesh, before, 0.0);
    for (auto& cell : mesh.cells) cell.U = Conserved{};
    const auto vacuum = mesh;
    CHECK(step(mesh, solver, cfg, 0.0, 0.01, {energy_source(0.0)}) == 0.01);
    same(mesh, vacuum, 0.0);
    expect_throw<std::invalid_argument>([&] {
        step(mesh, solver, cfg, 0.0, 1.0, {energy_source(0.0), energy_source(1.0)});
    });
    same(mesh, vacuum, 0.0);
}

void source_stage_failure_checks()
{
    Transport solver; solver.zero_flux = true;
    config::SimulationConfig cfg; cfg.cfl = 0.2; cfg.lts_max_levels = 5;
    auto mesh = make_mesh(false);
    const auto initial = mesh;
    const hd_source::SourceTerm excess{"stage excess", hd_source::SourceKind::Thermal,
        [](const MeshBlock& m, double t, const config::SimulationConfig&,
           std::vector<Conserved>& rhs, double&) {
            if (t > 0.0) rhs[m.nghost[0]][hd_core::ENERGY] = 1e5;
        }};
    const hd_source::SourceTerm bad{"later invalid stage rate", hd_source::SourceKind::Thermal,
        [](const MeshBlock& m, double t, const config::SimulationConfig&,
           std::vector<Conserved>& rhs, double&) {
            if (t > 0.0) rhs[m.nghost[0]+m.nactive[0]-1][hd_core::ENERGY] =
                std::numeric_limits<double>::quiet_NaN();
        }};
    // An earlier limit failure must not hide exceptions in later callbacks.
    expect_throw<std::runtime_error>([&] { step(mesh, solver, cfg, 0.0, 0.01, {excess, bad}); });
    same(mesh, initial, 0.0);
    auto combined = excess;
    combined.func = [excess, bad](const MeshBlock& m, double t, const config::SimulationConfig& c,
                                  std::vector<Conserved>& rhs, double& cap) {
        excess.func(m, t, c, rhs, cap); bad.func(m, t, c, rhs, cap);
    };
    // Nor may it hide a bad later cell within the same source scan.
    expect_throw<std::runtime_error>([&] { step(mesh, solver, cfg, 0.0, 0.01, {combined}); });
    same(mesh, initial, 0.0);
}

void shock(const std::string& reconstruction)
{
    hydro::ExactSolver solver(1e-10, 200);
    auto mesh = make_mesh(true, 7, true);
    mesh.geometry = hd_core::Geometry::Cartesian; mesh.build_geometry();
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        hd_core::Primitives W{};
        W[hd_core::RHO] = i < 6 ? 1.0 : 0.125; W[hd_core::PRE] = i < 6 ? 1.0 : 0.1;
        mesh.cells[i+2] = hd_core::Cell(W, 1.4);
    }
    boundary::FillBoundary(mesh, 0.0, {});
    auto reference = mesh; const auto before = integral(mesh);
    config::SimulationConfig cfg;
    cfg.cfl = 0.2; cfg.lts_max_levels = 4; cfg.reconstructor = reconstruction; cfg.limiter = "MC";
    evolve(mesh, solver, cfg, 0.2);
    cfg.cfl = 0.05;
    evolve(reference, solver, cfg, 0.2, {}, hd_time::TimeIntegrator::RK2);
    CHECK(error(mesh, reference) / before[hd_core::RHO] < 0.01);
    CHECK(near(integral(mesh)[hd_core::RHO], before[hd_core::RHO], 2e-12));
    CHECK(near(integral(mesh)[hd_core::ENERGY], before[hd_core::ENERGY], 2e-12));
}

double boundary_energy0 = 0.0;
int boundary_calls = 0;
void same_time_boundary(MeshBlock& mesh, double t, int, boundary::BoundarySide side)
{
    CHECK(side == boundary::BoundarySide::Lower);
    ++boundary_calls;
    // Even remote slow cells must be at the callback's stage time.
    for (int i = 0; i < mesh.nactive[0]; ++i)
        CHECK(near(mesh.cells[i + mesh.nghost[0]].U[hd_core::ENERGY], boundary_energy0 + t, 2e-13));
    for (int g = 0; g < mesh.nghost[0]; ++g) mesh.cells[g] = mesh.cells[mesh.nghost[0]];
}

void boundaries_and_inputs()
{
    auto mesh = make_mesh(true, 0, false, 2, 64, 0.001);
    mesh.bc_lower[0] = boundary::BoundaryCondition::User;
    mesh.bc_upper[0] = boundary::BoundaryCondition::Reflective;
    boundary_energy0 = mesh.cells[2].U[hd_core::ENERGY]; boundary_calls = 0;
    Transport solver; solver.zero_flux = true;
    config::SimulationConfig cfg; cfg.cfl = 0.2; cfg.lts_max_levels = 5; cfg.reconstructor = "PLM";
    hd_user::UserHooks hooks; hooks.userbc = same_time_boundary;
    hd_source::SourceList sources{{"uniform heating", hd_source::SourceKind::Thermal,
        [](const MeshBlock& m, double, const config::SimulationConfig&, std::vector<Conserved>& rhs, double&) {
            for (std::size_t k = 0; k < m.cells.size(); ++k) if (m.is_active(k)) rhs[k][hd_core::ENERGY] += 1.0;
        }}};
    CHECK(step(mesh, solver, cfg, 0.0, 1.0, sources, hd_time::TimeIntegrator::LTSRK2, nullptr, hooks) > 0.0);
    CHECK(boundary_calls > 16);
    mesh = make_mesh();
    auto before = mesh;
    for (int invalid : {0, 17})
    {
        cfg.lts_max_levels = invalid;
        expect_throw<std::invalid_argument>([&] { step(mesh, solver, cfg, 0.0, 1.0); });
    }
    cfg.lts_max_levels = 5; cfg.reconstructor = "typo";
    expect_throw<std::invalid_argument>([&] { step(mesh, solver, cfg, 0.0, 1.0); });
    cfg.reconstructor = "PLM";
    expect_throw<std::invalid_argument>([&] { step(mesh, solver, cfg, 0.0, -1.0); });
    same(mesh, before, 0.0);
    mesh = make_mesh(true, 7, false, 1);
    expect_throw<std::invalid_argument>([&] { step(mesh, solver, cfg, 0.0, 1.0); });
}
}
int main(int argc, char** argv)
{
    return run_test("multilevel conservative RK2", [&] {
        CHECK(hd_time::parse_time_integrator("LTS_RK2") == hd_time::TimeIntegrator::LTSRK2);
        if (argc > 1)
        {
            const std::string path = std::string(argv[1]) + ".txt";
            std::ofstream out(path);
            const char* axes[] = {"x", "y", "z"};
            for (int d = 0; d < hd_core::NDIM; ++d)
                out << "n" << axes[d] << "a=4\nn" << axes[d] << "g=2\n"
                    << axes[d] << "_lower=0\n" << axes[d] << "_upper=1\n";
            out << "time_integrator=lts_rk2\nlts_max_levels=8\n"; out.close();
            const auto cfg = config::ReadParameters(path);
            CHECK(cfg.time_integrator == "lts_rk2" && cfg.lts_max_levels == 8);
        }
        if (hd_core::NDIM != 1)
        {
            auto mesh = make_mesh(false); Transport solver;
            expect_throw<std::invalid_argument>([&] { step(mesh, solver, {}, 0.0, 1.0); });
            return;
        }
        for (const std::string method : {"PCM", "PLM", "PLM_corrected"})
        {
            equivalence(method);
            convergence(method, false);
            convergence(method, true);
            conservation_and_geometry(method);
            shock(method);
        }
        retries_and_sources();
        source_cache_stage_tests();
        source_stage_failure_checks();
        boundaries_and_inputs();
    });
}
