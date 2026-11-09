#include "test_support.h"
#include "time/integrator.h"
#include "hydro/riemann/exact.h"
#include "hydro/wave_timestep.h"
#include "user/user_api.h"
#include <limits>
#include <iomanip>
#include <sstream>

namespace
{
using hd_core::MeshBlock;
using hd_core::Conserved;
constexpr double infinity = std::numeric_limits<double>::infinity();

class CountSolver : public hydro::RiemannSolver
{
public:
    mutable int calls = 0;
    hydro::ExactSolver exact{1e-10, 200};
    Conserved flux(const hd_core::Primitives& L, const hd_core::Primitives& R,
                   int d, double gL, double gR) const override
    {
        ++calls;
        return exact.flux(L, R, d, gL, gR);
    }
};

class JumpSpeedSolver : public hydro::RiemannSolver
{
public:
    mutable int calls = 0;
    int initial_faces = 0;
    Conserved flux(const hd_core::Primitives&, const hd_core::Primitives&,
                   int, double, double) const override
    {
        hydro::wave_timestep.face_speed = ++calls <= initial_faces ? 1.0 : 100.0;
        return {};
    }
};

class HintSolver : public hydro::RiemannSolver
{
public:
    hydro::ExactSolver exact{1e-10, 200};
    mutable std::vector<double> inputs;
    int throw_on_call = -1;
    Conserved flux(const hd_core::Primitives& L, const hd_core::Primitives& R,
                   int d, double gL, double gR) const override
    { return exact.flux(L, R, d, gL, gR); }
    Conserved flux_with_pressure_hint(const hd_core::Primitives& L, const hd_core::Primitives& R,
                   int d, double gL, double gR, double& pressure) const override
    {
        inputs.push_back(pressure);
        if (static_cast<int>(inputs.size()) == throw_on_call)
            throw std::runtime_error("Injected failure after partial pressure-history update");
        return exact.flux_with_pressure_hint(L, R, d, gL, gR, pressure);
    }
};

MeshBlock make_mesh(bool nonuniform = true, int narrow = 7, bool spherical = false,
                    int ghosts = 2, int active_cells = 20)
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
            faces[ghosts + i + 1] = faces[ghosts + i] + (i == narrow ? 0.01 : 0.1);
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

double advance(MeshBlock& mesh, const hydro::RiemannSolver& solver,
               double time = 0.0, double limit = 1.0,
               const hd_source::SourceList& sources = {},
               config::SimulationConfig cfg = {}, hd_time::LTSEulerWorkspace* workspace = nullptr)
{
    double next = 0.123;
    cfg.cfl = 0.2;
    const double dt = hd_time::advance_one_step(mesh, time, limit,
        hd_time::TimeIntegrator::LTSEuler, hydro::SpatialScheme::Godunov,
        solver, sources, cfg, {}, cfg.cfl, next, workspace);
    CHECK(std::isinf(next));
    return dt;
}

Conserved integral(const MeshBlock& mesh)
{
    Conserved result{};
    for (auto k : mesh.get_all_active_flats()) result += mesh.cells[k].U * mesh.volume_cell[k];
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

void equivalence_tests(const std::string& method)
{
    CountSolver solver;
    config::SimulationConfig cfg;
    cfg.cfl = 0.2;
    cfg.reconstructor = method;
    cfg.limiter = "MC";
    // Truncated macro step: every cell is slow, identical to one Euler update.
    auto a = make_mesh();
    contact(a);
    auto b = a;
    const double dt = advance(a, solver, 0.0, 1e-4, {}, cfg);
    hd_time::advance_euler(b, 0.0, dt, hydro::SpatialScheme::Godunov, solver, {}, cfg, {}, cfg.cfl);
    boundary::FillBoundary(b, dt, {});
    same(a, b);

    // Uniform mesh/state: every cell is fast, identical to two Euler steps.
    a = make_mesh(false);
    b = a;
    const double H = advance(a, solver, 0.0, 1.0, {}, cfg);
    const double h1 = hd_time::advance_euler(b, 0.0, H / 2, hydro::SpatialScheme::Godunov,
                                            solver, {}, cfg, {}, cfg.cfl);
    boundary::FillBoundary(b, h1, {});
    const double h2 = hd_time::advance_euler(b, h1, H / 2, hydro::SpatialScheme::Godunov,
                                            solver, {}, cfg, {}, cfg.cfl);
    boundary::FillBoundary(b, h1 + h2, {});
    CHECK(near(h1 + h2, H));
    same(a, b);
}

void conservation_and_cost_tests(const std::string& method)
{
    for (int narrow : {0, 7, 19})
    {
        CountSolver solver;
        auto mesh = make_mesh(true, narrow);
        contact(mesh);
        const auto initial = integral(mesh);
        auto reference = mesh;
        config::SimulationConfig cfg;
        cfg.cfl = 0.2;
        cfg.reconstructor = method;
        cfg.limiter = "MC";
        double time = 0.0;
        hd_time::LTSEulerWorkspace workspace;
        for (int step = 0; step < 40; ++step)
        {
            const int before = solver.calls;
            const double H = advance(mesh, solver, time, 1.0, {}, cfg, &workspace);
            // One global pass and only the fast-adjacent faces in pass two.
            CHECK(solver.calls - before < 2 * mesh.nactive[0]);
            double ref_time = time;
            while (ref_time < time + H)
            {
                const double h = hd_time::advance_euler(reference, ref_time,
                    std::min(H / 2, time + H - ref_time), hydro::SpatialScheme::Godunov,
                    solver, {}, cfg, {}, cfg.cfl);
                ref_time += h;
                boundary::FillBoundary(reference, ref_time, {});
            }
            time += H;
            for (auto k : mesh.get_all_active_flats()) mesh.cells[k].U.validate(mesh.cells[k].gamma);
        }
        const auto final = integral(mesh);
        for (int v = 0; v < hd_core::NVAR; ++v) CHECK(near(final[v], initial[v], 1e-12));
        same(mesh, reference, 3e-3);
    }
}

void source_and_geometry_tests(const std::string& method)
{
    CountSolver solver;
    config::SimulationConfig cfg;
    cfg.reconstructor = method;
    cfg.limiter = "MC";
    auto mesh = make_mesh(true, 7, true);
    const auto original = mesh;
    // Spherical geometric source must cancel uniform pressure flux, even at
    // a fast/slow interface. A flux-only conservation test cannot catch this.
    const hd_source::SourceList geometric{hd_source::make_geometric_source()};
    advance(mesh, solver, 0.0, 1.0, geometric, cfg);
    // Allow the exact solver's finite pressure-root tolerance (1e-10).
    same(mesh, original, 2e-11);

    mesh = make_mesh();
    const auto before = integral(mesh);
    int source_calls = 0;
    hd_source::SourceList heating{{"heating_test", hd_source::SourceKind::Thermal,
        [&](const MeshBlock& state, double time, const config::SimulationConfig&,
            std::vector<Conserved>& rhs, double& limit) {
            ++source_calls;
            CHECK(time == 0.0); // Never use asynchronous states to call a source.
            for (auto k : state.get_all_active_flats()) rhs[k][hd_core::ENERGY] = 3.0;
            limit = 1.0;
        }}};
    const double H = advance(mesh, solver, 0.0, 1.0, heating, cfg);
    CHECK(source_calls == 1);
    const auto after = integral(mesh);
    CHECK(near(after[hd_core::RHO], before[hd_core::RHO]));
    CHECK(near(after[hd_core::ENERGY], before[hd_core::ENERGY] + 3.0 * H * mesh.active_length(0)));

    // Source hard limits can force an all-slow step.
    heating.push_back({"limit", hd_source::SourceKind::General,
        [](const MeshBlock&, double, const config::SimulationConfig&,
           std::vector<Conserved>&, double& limit) { limit = 1e-5; }});
    mesh = make_mesh();
    CHECK(near(advance(mesh, solver, 0.0, 1.0, heating, cfg), 1e-5));

    // Gravity has no existing callback timestep limit: LTS must impose its own.
    mesh = make_mesh();
    const hd_source::SourceList gravity{hd_source::make_constant_gravity({{1e8}})};
    CHECK(advance(mesh, solver, 0.0, 1.0, gravity, cfg) < 1e-8);
}

void retry_and_input_tests(const std::string& method)
{
    auto mesh = make_mesh();
    const auto before = mesh;
    JumpSpeedSolver solver;
    solver.initial_faces = mesh.nactive[0]; // Periodic seam solved only once.
    config::SimulationConfig cfg;
    cfg.reconstructor = method;
    cfg.integrator_retry_max = 0;
    expect_throw<std::runtime_error>([&] { advance(mesh, solver, 0.0, 1.0, {}, cfg); });
    same(mesh, before, 0.0);

    solver.calls = 0;
    cfg.integrator_retry_max = 1;
    CHECK(near(advance(mesh, solver, 0.0, 1.0, {}, cfg), 0.002));
    same(mesh, before, 0.0);

    CountSolver exact;
    cfg.reconstructor = "not_a_reconstructor";
    expect_throw<std::invalid_argument>([&] { advance(mesh, exact, 0.0, 1.0, {}, cfg); });
    same(mesh, before, 0.0);
    cfg.reconstructor = "PLM";
    auto one_ghost = make_mesh(true, 7, false, 1);
    const auto original = one_ghost;
    expect_throw<std::invalid_argument>([&] { advance(one_ghost, exact, 0.0, 1.0, {}, cfg); });
    same(one_ghost, original, 0.0);
    expect_throw<std::invalid_argument>([&] { advance(mesh, exact, 0.0, -1.0); });
    cfg.reconstructor = "PCM";
    const hd_source::SourceList bad{{"bad", hd_source::SourceKind::General,
        [](const MeshBlock&, double, const config::SimulationConfig&,
           std::vector<Conserved>& rhs, double&) { rhs.clear(); }}};
    expect_throw<std::logic_error>([&] { advance(mesh, exact, 0.0, 1.0, bad); });
    same(mesh, before, 0.0);
}

void shock_test(const std::string& method)
{
    CountSolver solver;
    auto mesh = make_mesh(true, 7, true);
    mesh.geometry = hd_core::Geometry::Cartesian;
    mesh.build_geometry();
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        hd_core::Primitives W{};
        W[hd_core::RHO] = i < 6 ? 1.0 : 0.125;
        W[hd_core::PRE] = i < 6 ? 1.0 : 0.1;
        mesh.cells[i + mesh.nghost[0]] = hd_core::Cell(W, 1.4);
    }
    boundary::FillBoundary(mesh, 0.0, {});
    auto reference = mesh;
    const auto initial = integral(mesh);
    config::SimulationConfig cfg;
    cfg.cfl = 0.2;
    cfg.reconstructor = method;
    cfg.limiter = "MC";
    double t = 0.0;
    hd_time::LTSEulerWorkspace workspace;
    while (t < 0.2)
    {
        const double H = advance(mesh, solver, t, 0.2 - t, {}, cfg, &workspace);
        double ref = t;
        while (ref < t + H)
        {
            const double h = hd_time::advance_euler(reference, ref, std::min(H / 2, t + H - ref),
                hydro::SpatialScheme::Godunov, solver, {}, cfg, {}, cfg.cfl);
            ref += h;
            boundary::FillBoundary(reference, ref, {});
        }
        t += H;
    }
    const auto final = integral(mesh);
    CHECK(near(final[hd_core::RHO], initial[hd_core::RHO], 1e-12));
    CHECK(near(final[hd_core::ENERGY], initial[hd_core::ENERGY], 1e-12));
    double difference = 0.0, mass = 0.0;
    for (auto k : mesh.get_all_active_flats())
    {
        mesh.cells[k].U.validate(mesh.cells[k].gamma);
        difference += mesh.volume_cell[k] * std::abs(mesh.cells[k].U[hd_core::RHO]
                                                    - reference.cells[k].U[hd_core::RHO]);
        mass += mesh.volume_cell[k] * reference.cells[k].U[hd_core::RHO];
    }
    CHECK(difference / mass < 0.02);
}

void midpoint_stencil_test(const std::string& method)
{
    // Independent oracle: at the second substep, all stencil cells must agree
    // with a global Euler half-step, including slow neighbours and source terms.
    class RecordingSolver : public CountSolver
    {
    public:
        mutable std::vector<std::pair<hd_core::Primitives, hd_core::Primitives>> states;
        Conserved flux(const hd_core::Primitives& L, const hd_core::Primitives& R,
                       int d, double gL, double gR) const override
        {
            states.emplace_back(L, R);
            return CountSolver::flux(L, R, d, gL, gR);
        }
    } recorded;
    CountSolver reference_solver;
    auto mesh = make_mesh();
    contact(mesh);
    auto midpoint = mesh;
    config::SimulationConfig cfg;
    cfg.cfl = 0.2;
    cfg.reconstructor = method;
    cfg.limiter = "MC";
    const hd_source::SourceList sources{hd_source::make_constant_heating(0.3)};
    double rate, source_limit;
    auto [h, rhs] = hd_time::evaluate_RHS(midpoint, 0.0, hydro::SpatialScheme::Godunov,
        sources, reference_solver, cfg, cfg.cfl, rate, source_limit);
    hd_time::add_scaled_rhs(midpoint, rhs, h);
    boundary::FillBoundary(midpoint, h, {});

    // Deliberately wrong mixed-time stencil to ensure this test distinguishes
    // the predictor from simply enabling PLM on the original PCM state view.
    auto stale = mesh;
    for (int i = 5; i <= 9; ++i) stale.cells[i + stale.nghost[0]] = midpoint.cells[i + midpoint.nghost[0]];
    boundary::FillBoundary(stale, h, {});

    const double H = advance(mesh, recorded, 0.0, 1.0, sources, cfg);
    CHECK(near(H, 2 * h));
    CHECK(recorded.states.size() == 26); // 20 unique initial + six fine faces.
    double stale_difference = 0.0;
    for (int f = 5; f <= 10; ++f)
    {
        hd_core::Index index{};
        index[0] = midpoint.nghost[0] + f;
        const auto expected = hydro::Reconstruct(hydro::parse_reconstructor(method),
            midpoint, 0, index, hydro::Limiter::MC, cfg.limiter_beta);
        const auto wrong = hydro::Reconstruct(hydro::parse_reconstructor(method),
            stale, 0, index, hydro::Limiter::MC, cfg.limiter_beta);
        const auto& actual = recorded.states[20 + f - 5];
        for (int k = 0; k < hd_core::NHYDRO; ++k)
        {
            CHECK(near(actual.first[k], expected.first[k], 1e-12));
            CHECK(near(actual.second[k], expected.second[k], 1e-12));
            stale_difference += std::abs(expected.first[k] - wrong.first[k])
                              + std::abs(expected.second[k] - wrong.second[k]);
        }
    }
    CHECK(stale_difference > 1e-6);
}

void workspace_tests()
{
    CountSolver solver;
    hd_time::LTSEulerWorkspace workspace;
    config::SimulationConfig cfg;
    // Reuse across distinct meshes, geometry rebuilds, states and reconstructors.
    // Fresh-workspace and reused-workspace paths must give the same results.
    for (int pass = 0; pass < 6; ++pass)
    {
        const bool spherical = pass % 2 != 0;
        auto reused = make_mesh(true, pass == 2 ? 0 : 7, spherical, 2 + pass % 2, 16 + pass);
        if (!spherical) contact(reused);
        if (pass == 4)
        {
            reused.face[0][reused.nghost[0] + 4] += 0.001;
            reused.build_geometry();
        }
        auto fresh = reused;
        cfg.reconstructor = pass % 3 == 0 ? "PCM" : pass % 3 == 1 ? "PLM" : "PLM_corrected";
        const auto sources = spherical
            ? hd_source::SourceList{hd_source::make_geometric_source(), hd_source::make_constant_heating(0.2)}
            : hd_source::SourceList{hd_source::make_constant_heating(0.2)};
        double t = 0.0;
        for (int step = 0; step < 3; ++step)
        {
            const double a = advance(reused, solver, t, 1.0, sources, cfg, &workspace);
            const double b = advance(fresh, solver, t, 1.0, sources, cfg);
            CHECK(a == b);
            same(reused, fresh, 0.0);
            t += a;
        }
    }

    // A failed trial must not poison the workspace for a subsequent call.
    JumpSpeedSolver jump;
    auto trial = make_mesh();
    jump.initial_faces = trial.nactive[0];
    cfg.integrator_retry_max = 0;
    const auto before = trial;
    expect_throw<std::runtime_error>([&] { advance(trial, jump, 0.0, 1.0, {}, cfg, &workspace); });
    same(trial, before, 0.0);
    auto reference = before;
    const double H = advance(trial, solver, 0.0, 1.0, {}, cfg, &workspace);
    CHECK(H == advance(reference, solver, 0.0, 1.0, {}, cfg));
    same(trial, reference, 0.0);
}

void pressure_history_test()
{
    auto mesh = make_mesh();
    contact(mesh);
    auto reference = mesh;
    HintSolver solver;
    hd_time::LTSEulerWorkspace workspace, reference_workspace;
    config::SimulationConfig cfg;
    cfg.reconstructor = "PLM";
    const int faces = mesh.nactive[0]; // Periodic seam is solved only once.
    double t = advance(mesh, solver, 0.0, 1.0, {}, cfg, &workspace);
    CHECK(t == advance(reference, solver.exact, 0.0, 1.0, {}, cfg, &reference_workspace));
    same(mesh, reference, 0.0);
    CHECK(solver.inputs.size() > static_cast<std::size_t>(faces));
    for (int f = 0; f < faces; ++f) CHECK(solver.inputs[f] == 0.0);
    for (std::size_t f = faces; f < solver.inputs.size(); ++f) CHECK(solver.inputs[f] > 0.0);

    solver.inputs.clear();
    solver.throw_on_call = faces + 2; // Fail during the fine step, after one updated face.
    const auto before = mesh;
    expect_throw<std::runtime_error>([&] { advance(mesh, solver, t, 1.0, {}, cfg, &workspace); });
    same(mesh, before, 0.0);
    const auto failed_inputs = solver.inputs;
    for (int f = 0; f < faces; ++f) CHECK(failed_inputs[f] > 0.0);

    solver.throw_on_call = -1;
    solver.inputs.clear();
    CHECK(advance(mesh, solver, t, 1.0, {}, cfg, &workspace) ==
          advance(reference, solver.exact, t, 1.0, {}, cfg, &reference_workspace));
    same(mesh, reference, 0.0);
    CHECK(solver.inputs.size() >= failed_inputs.size());
    for (std::size_t f = 0; f < failed_inputs.size(); ++f)
        CHECK(solver.inputs[f] == failed_inputs[f]);

    // Even same-sized coordinate edits invalidate the old face identities.
    mesh.face[0][mesh.nghost[0] + 4] += 0.001;
    mesh.build_geometry();
    solver.inputs.clear();
    advance(mesh, solver, 0.0, 1.0, {}, cfg, &workspace);
    for (int f = 0; f < faces; ++f) CHECK(solver.inputs[f] == 0.0);
}

double boundary_half_time = 0.0, boundary_expected_density = 0.0;
int boundary_half_calls = 0;
void global_read_boundary(MeshBlock& mesh, double time, int, boundary::BoundarySide side)
{
    CHECK(side == boundary::BoundarySide::Lower);
    const int ng = mesh.nghost[0];
    if (time == boundary_half_time)
    {
        ++boundary_half_calls;
        // This remote slow cell is outside every fine PLM stencil.
        CHECK(near(mesh.cells[ng + 12].U[hd_core::RHO], boundary_expected_density, 1e-12));
    }
    for (int g = 0; g < ng; ++g) mesh.cells[g] = mesh.cells[ng];
}

void custom_boundary_prediction_test()
{
    auto mesh = make_mesh(true, 0);
    contact(mesh);
    mesh.bc_lower[0] = boundary::BoundaryCondition::User;
    mesh.bc_upper[0] = boundary::BoundaryCondition::Reflective;
    hd_user::UserHooks hooks;
    hooks.userbc = global_read_boundary;
    boundary_half_time = -1.0;
    boundary_half_calls = 0;
    boundary::FillBoundary(mesh, 0.0, hooks);
    config::SimulationConfig cfg;
    cfg.reconstructor = "PLM";
    cfg.cfl = 0.2;
    CountSolver solver;
    double rate, source_limit;
    auto [h, rhs] = hd_time::evaluate_RHS(mesh, 0.0, hydro::SpatialScheme::Godunov,
        {}, solver, cfg, cfg.cfl, rate, source_limit);
    const auto middle = mesh.cells[mesh.nghost[0] + 12].U + h * rhs[mesh.nghost[0] + 12];
    boundary_half_time = h;
    boundary_expected_density = middle[hd_core::RHO];
    CHECK(std::abs(boundary_expected_density - mesh.cells[mesh.nghost[0] + 12].U[hd_core::RHO]) > 1e-7);
    hd_time::LTSEulerWorkspace workspace;
    double next = infinity;
    const double H = hd_time::advance_one_step(mesh, 0.0, 1.0, hd_time::TimeIntegrator::LTSEuler,
        hydro::SpatialScheme::Godunov, solver, {}, cfg, hooks, cfg.cfl, next, &workspace);
    CHECK(H == 2 * h);
    CHECK(boundary_half_calls == 1);
}
}

int main()
{
    return run_test("two-level LTS Euler conservation, synchronization and rollback", [] {
        CHECK(hd_time::parse_time_integrator("LTS_EULER") == hd_time::TimeIntegrator::LTSEuler);
        if (hd_core::NDIM != 1)
        {
            auto mesh = make_mesh(false);
            CountSolver solver;
            expect_throw<std::invalid_argument>([&] { advance(mesh, solver); });
            return;
        }
        for (const std::string method : {"PCM", "PLM", "PLM_corrected"})
        {
            std::cout << method << ": Euler equivalence\n"; equivalence_tests(method);
            std::cout << method << ": Conservation and face cost\n"; conservation_and_cost_tests(method);
            std::cout << method << ": Sources and geometry\n"; source_and_geometry_tests(method);
            std::cout << method << ": Retry and inputs\n"; retry_and_input_tests(method);
            std::cout << method << ": Shock crossing\n"; shock_test(method);
            if (method != "PCM")
            {
                std::cout << method << ": Same-time reconstruction stencil\n";
                midpoint_stencil_test(method);
            }
        }
        workspace_tests();
        pressure_history_test();
        custom_boundary_prediction_test();
    });
}
