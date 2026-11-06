#include "core/mesh.h"
#include "config/config.h"
#include "hydro/flux.h"
#include "hydro/riemann/exact.h"
#include "hydro/spatial.h"
#include "hydro/timestep_CFL.h"
#include "hydro/wave_timestep.h"
#include "time/integrator.h"
#include "user/user_api.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>

namespace
{
constexpr double gamma_gas = 1.4;

bool near(double a, double b)
{
    return std::abs(a - b) <= 1e-9 * std::max({1.0, std::abs(a), std::abs(b)});
}

template<class Exception, class Function>
void expect_throw(Function function)
{
    bool caught = false;
    try { function(); }
    catch (const Exception&) { caught = true; }
    assert(caught);
}

hd_core::Primitives state(double rho, double pressure, double velocity = 0.0, int d = 0)
{
    hd_core::Primitives w{};
    w[hd_core::RHO] = rho;
    w[hd_core::PRE] = pressure;
    w[hd_core::VEL(d)] = velocity;
    return w;
}

bool active(const hd_core::MeshBlock& mesh, std::size_t flat)
{
    const auto i = mesh.coordinates_active(flat);
    for (int d = 0; d < hd_core::NDIM; ++d)
        if (i[d] < 0 || i[d] >= mesh.nactive[d]) return false;
    return true;
}

hd_core::MeshBlock make_mesh(bool nonuniform = false)
{
    hd_core::MeshBlock mesh;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        mesh.nactive[d] = 2 + d;
        mesh.nghost[d] = 1;
        mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
        mesh.face[d].push_back(0.0);
        for (int j = 0; j < mesh.ntotal(d); ++j)
            mesh.face[d].push_back(mesh.face[d].back()
                + 0.1 * (d + 1) * (nonuniform ? 1.0 + 0.2 * j : 1.0));
    }
    mesh.allocate();
    mesh.build_geometry();
    for (auto& cell : mesh.cells)
        cell = hd_core::Cell(state(1.0, 1.0), gamma_gas);
    return mesh;
}

void same_cells(const std::vector<hd_core::Cell>& a, const std::vector<hd_core::Cell>& b)
{
    assert(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        assert(a[i].gamma == b[i].gamma);
        assert(a[i].U.data == b[i].U.data);
    }
}

void exact_wave_tests()
{
    hydro::ExactSolver solver(1e-12, 200);
    const double c = std::sqrt(gamma_gas);
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        const auto w = state(1.0, 1.0, 20.0, d);
        const auto uniform = solver.flux(w, w, d, gamma_gas, gamma_gas);
        assert(near(hydro::wave_timestep.face_speed, 20.0 + c));
        const auto physical = hydro::physical_flux(w, d, gamma_gas);
        for (int k = 0; k < hd_core::NVAR; ++k)
            assert(near(uniform[k], physical[k]));

        solver.flux(state(1000.0, 1000.0), state(1.0, 0.01), d, gamma_gas, gamma_gas);
        // Independent reference: p*=9.832118791607545, right shock speed below.
        assert(near(hydro::wave_timestep.face_speed, 3.4351917777511423));
        solver.flux(state(1.0, 0.01), state(1000.0, 1000.0), d, gamma_gas, gamma_gas);
        assert(near(hydro::wave_timestep.face_speed, 3.4351917777511423));

        const auto gas = state(1.0, 1.0);
        const auto vacuum = state(0.0, 0.0);
        solver.flux(gas, vacuum, d, gamma_gas, gamma_gas);
        assert(near(hydro::wave_timestep.face_speed, 2.0 * c / (gamma_gas - 1.0)));
        solver.flux(vacuum, gas, d, gamma_gas, gamma_gas);
        assert(near(hydro::wave_timestep.face_speed, 2.0 * c / (gamma_gas - 1.0)));

        const auto generated = solver.flux(state(1.0, 1.0, -10.0, d),
            state(1.0, 1.0, 10.0, d), d, gamma_gas, gamma_gas);
        assert(near(hydro::wave_timestep.face_speed, 10.0 + c));
        for (double f : generated.data) assert(f == 0.0);
        const auto empty = solver.flux(vacuum, vacuum, d, gamma_gas, gamma_gas);
        assert(hydro::wave_timestep.face_speed == 0.0);
        for (double f : empty.data) assert(f == 0.0);
    }
    hydro::ExactSolver unconverged(1e-12, 1);
    expect_throw<std::runtime_error>([&] {
        unconverged.flux(state(1000., 1000.), state(1., 0.01), 0, gamma_gas, gamma_gas);
    });
}

// A controlled face speed lets the test independently check the mesh CFL aggregation.
struct CountingSolver final : hydro::RiemannSolver
{
    mutable std::size_t calls = 0;
    bool invalid_flux = false;
    bool publish_speed = true;
    hd_core::Conserved flux(const hd_core::Primitives&, const hd_core::Primitives&,
                              int d, double, double) const override
    {
        ++calls;
        hd_core::Conserved r{};
        if (publish_speed) hydro::wave_timestep.face_speed = 2.0 * (d + 1.0);
        if (invalid_flux) r[hd_core::RHO] = std::numeric_limits<double>::quiet_NaN();
        return r;
    }
};

void integration_tests()
{
    const config::SimulationConfig cfg{};
    const hd_user::UserHooks hooks{};
    const hd_source::SourceList sources{}; // Existing checks exercise source-free hydro.
    CountingSolver counting;
    auto mesh = make_mesh(true);
    const auto original = mesh.cells;
    std::vector<hd_core::Conserved> rhs;
    hydro::wave_timestep.max_inverse_dt = 1e100; // Must reset for this new RHS.
    hydro::compute_hydro_rhs(mesh, hydro::SpatialScheme::Godunov, counting, rhs);
    double expected_rate = 0.0;
    for (std::size_t f = 0; f < mesh.cells.size(); ++f)
    {
        for (double value : rhs[f].data) assert(value == 0.0);
        if (!active(mesh, f)) continue;
        const auto idx = mesh.coordinates(f);
        double rate = 0.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
            rate += 2.0 * (d + 1.0) / (mesh.face[d][idx[d] + 1] - mesh.face[d][idx[d]]);
        expected_rate = std::max(expected_rate, rate);
    }
    assert(near(hydro::wave_timestep.max_inverse_dt, expected_rate));
    std::size_t face_count = 0;
    for (int n = 0; n < hd_core::NDIM; ++n)
    {
        std::size_t count = 1;
        for (int d = 0; d < hd_core::NDIM; ++d) count *= mesh.nactive[d] + (d == n);
        face_count += count;
    }
    assert(counting.calls == face_count);
    counting.calls = 0;
    double dt_next = std::numeric_limits<double>::infinity();
    const double dt = hd_time::advance_one_step(mesh, 0.0, 1.0,
        hd_time::TimeIntegrator::Euler, hydro::SpatialScheme::Godunov, counting, sources, cfg, hooks, 0.4, dt_next);
    assert(near(dt, 0.4 / expected_rate));
    assert(counting.calls == face_count); // No duplicate Riemann solves for dt.
    same_cells(mesh.cells, original);
    assert(hd_time::advance_euler(mesh, 0.0, dt / 2,
        hydro::SpatialScheme::Godunov, counting, sources, cfg, hooks) == dt / 2);
    counting.invalid_flux = true;
    expect_throw<std::invalid_argument>([&] {
        hd_time::advance_euler(mesh, 0.0, 1.0, hydro::SpatialScheme::Godunov, counting, sources, cfg, hooks);
    });
    same_cells(mesh.cells, original); // Failed candidate update is not committed.
    counting.invalid_flux = false;
    counting.publish_speed = false;
    expect_throw<std::runtime_error>([&] {
        hd_time::advance_euler(mesh, 0.0, 1.0, hydro::SpatialScheme::Godunov, counting, sources, cfg, hooks);
    });
    same_cells(mesh.cells, original);

    hydro::ExactSolver exact(1e-12, 200);
    auto shock = make_mesh();
    for (std::size_t f = 0; f < shock.cells.size(); ++f)
    {
        const auto idx = shock.coordinates(f);
        shock.cells[f] = hd_core::Cell(idx[0] <= shock.nghost[0]
            ? state(1000., 1000.) : state(1., 0.01), gamma_gas);
    }
    const auto before = shock.cells;
    const double old_estimate = hydro::compute_cfl_dt(shock);
    hydro::compute_hydro_rhs(shock, hydro::SpatialScheme::Godunov, exact, rhs);
    const double expected_shock_dt = hydro::compute_cfl_dt();
    const double shock_dt = hd_time::advance_euler(shock, 0.0, 1.0,
        hydro::SpatialScheme::Godunov, exact, sources, cfg, hooks);
    assert(shock_dt < old_estimate);
    assert(near(shock_dt, expected_shock_dt));
    for (std::size_t f = 0; f < shock.cells.size(); ++f)
    {
        if (!active(shock, f)) assert(shock.cells[f].U.data == before[f].U.data);
        else for (int k = 0; k < hd_core::NVAR; ++k)
            assert(near(shock.cells[f].U[k], before[f].U[k] + shock_dt * rhs[f][k]));
    }

    for (auto& cell : mesh.cells) cell = hd_core::Cell{};
    const auto empty_before = mesh.cells;
    assert(hd_time::advance_euler(mesh, 0.0, 0.1,
        hydro::SpatialScheme::Godunov, exact, sources, cfg, hooks) == 0.1);
    same_cells(mesh.cells, empty_before);
    assert(std::isinf(hydro::compute_cfl_dt()));
    expect_throw<std::invalid_argument>([&] {
        hd_time::advance_euler(mesh, 0., 1., hydro::SpatialScheme::Godunov, exact, sources, cfg, hooks, 0.0);
    });
    expect_throw<std::invalid_argument>([&] {
        hd_time::advance_euler(mesh, 0., -1., hydro::SpatialScheme::Godunov, exact, sources, cfg, hooks);
    });
    assert(hd_time::parse_time_integrator("euler") == hd_time::TimeIntegrator::Euler);
    assert(hd_time::parse_time_integrator("rk2") == hd_time::TimeIntegrator::RK2);
    expect_throw<std::invalid_argument>([] { hd_time::parse_time_integrator("unknown"); });
}
}

int main()
{
    exact_wave_tests();
    integration_tests();
    std::cout << "NDIM=" << hd_core::NDIM << ": shared wave state, CFL and Euler checks passed\n";
}
