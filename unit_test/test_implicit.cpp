#include "test_support.h"
#include "time/implicit.h"
#include "time/implicit_variables.h"
#include "hydro/riemann/exact.h"
#include "microphysics/composition/species.h"
#include "microphysics/eos/eos.h"
#include "physcon/constants.h"
#include "user/user_api.h"
#include <array>
#include <filesystem>
#include <fstream>
#include <limits>

namespace
{
using namespace hd_core;
int rejected_states = 0;

std::array<double, 9> metadata(const Cell& cell)
{
    return {cell.gamma, cell.temperature, cell.mean_molecular_weight,
            cell.n_baryon, cell.n_Ion, cell.n_e_free, cell.Y_Ion,
            cell.Y_e_free, cell.dual_energy_correction};
}

void stale(Cell& cell)
{
    cell.temperature = -1.0;
    cell.n_baryon = cell.n_Ion = cell.n_e_free = -1.0;
}

void check_thermodynamics(const Cell& cell)
{
    const double rho = cell.U[RHO];
    double mu = cell.mean_molecular_weight;
    if constexpr (NCOMP > 0)
    {
        double nuclei = 0.0;
        for (int k = 0; k < NCOMP; ++k)
            nuclei += cell.U[COMP(k)] / (hd_microphysics::species_table.species_active[k].A * physcon::m_u);
        mu = rho / (nuclei * physcon::m_u);
        CHECK(near(cell.n_baryon, rho / physcon::m_u));
        CHECK(near(cell.n_Ion, nuclei));
        CHECK(cell.n_e_free == 0.0);
        CHECK(near(cell.mean_molecular_weight, mu));
    }
    const double expected = (cell.gamma - 1.0) * cell.U.internal_energy() * mu
                          / (rho * physcon::R_g);
    CHECK(std::isfinite(cell.temperature) && cell.temperature > 0.0);
    CHECK(std::abs(cell.temperature / expected - 1.0) < 1e-10);
}

MeshBlock make_mesh()
{
    MeshBlock mesh;
    for (int d = 0; d < NDIM; ++d)
    {
        mesh.nactive[d] = 2;
        mesh.nghost[d] = 1;
        mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::User;
        mesh.set_linear_uniform_faces(d, 0.0, 1.0);
    }
    mesh.allocate();
    mesh.build_geometry();
    Primitives w{};
    w[RHO] = 2.0;
    w[PRE] = 3.0;
    if constexpr (NCOMP == 2)
    {
        w[COMP(0)] = 0.7;
        w[COMP(1)] = 0.3;
    }
    for (auto& cell : mesh.cells)
    {
        cell = Cell(w, 1.4);
        stale(cell);
    }
    return mesh;
}

// Mimic a user boundary which writes U but leaves thermodynamic caches stale.
void fill_ghosts(MeshBlock& mesh, double, int, boundary::BoundarySide)
{
    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
    {
        if (mesh.is_active(flat)) continue;
        auto index = mesh.coordinates(flat);
        for (int d = 0; d < NDIM; ++d)
            index[d] = std::clamp(index[d], mesh.nghost[d], mesh.nghost[d] + mesh.nactive[d] - 1);
        mesh.cells[flat] = mesh.at(index);
        stale(mesh.cells[flat]);
    }
}

config::SimulationConfig configuration(bool adaptive)
{
    config::SimulationConfig cfg;
    cfg.reconstructor = "PCM";
    cfg.dt_implicit = 0.1;
    cfg.implicit_adaptive = adaptive;
    cfg.implicit_tolerance = 1e-10;
    return cfg;
}

void temperature_source_and_line_search()
{
    auto mesh = make_mesh();
    const auto initial = mesh.cells;
    const auto cfg = configuration(false);
    hydro::ExactSolver solver(1e-10, 200);
    int evaluations = 0;
    rejected_states = 0;
    hd_user::UserHooks hooks;
    hooks.userbc = [](MeshBlock& state, double time, int direction, boundary::BoundarySide side) {
        fill_ghosts(state, time, direction, side);
        // The first full log-density Newton correction overshoots this bound;
        // the true backward-Euler solution rho=2/0.9 lies inside it.
        if (state.at_active(Index{}).U[RHO] > 2.228)
        {
            state.cells.front().U[RHO] = std::numeric_limits<double>::quiet_NaN();
            ++rejected_states;
        }
    };
    hd_source::SourceList sources{{"temperature-dependent heating", hd_source::SourceKind::General,
        [&](const MeshBlock& state, double, const config::SimulationConfig&,
            std::vector<Conserved>& rhs, double&) {
            ++evaluations;
            for (std::size_t flat = 0; flat < state.cells.size(); ++flat)
            {
                const auto& cell = state.cells[flat];
                check_thermodynamics(cell); // Includes newly populated user ghosts.
                if (!state.is_active(flat)) continue;
                const double heat = 0.5 * cell.U[RHO] * physcon::R_g * cell.temperature
                                  / (cell.mean_molecular_weight * (cell.gamma - 1.0));
                rhs[flat][RHO] += cell.U[RHO];
                rhs[flat][ENERGY] += heat;
                rhs[flat][INTERNAL_ENERGY] += heat;
                for (int k = 0; k < NCOMP; ++k)
                    rhs[flat][COMP(k)] += cell.U[COMP(k)];
            }
        }, true}};
    hd_implicit::Integrator integrator(mesh, cfg);
    const double dt = integrator.advance(mesh, 0.0, 0.1, hydro::SpatialScheme::Godunov,
                                         solver, sources, cfg, hooks);
    CHECK(dt == 0.1); // Recover through line search, without reducing the time step.
    CHECK(evaluations > 3 && rejected_states > 0);
    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
    {
        check_thermodynamics(mesh.cells[flat]);
        CHECK(near(mesh.cells[flat].U[RHO], 2.0 / (1.0 - dt), 1e-8));
        CHECK(near(mesh.cells[flat].U[ENERGY], initial[flat].U[ENERGY] / (1.0 - 0.5*dt), 1e-8));
        if constexpr (NCOMP > 0)
            hd_microphysics::validate_cell_composition(mesh.cells[flat]);
    }
}

void retry_and_rollback()
{
    hydro::ExactSolver solver(1e-10, 200);
    hd_user::UserHooks hooks;
    hooks.userbc = [](MeshBlock& state, double time, int direction, boundary::BoundarySide side) {
        fill_ghosts(state, time, direction, side);
        if (time > 0.075)
            state.cells.front().U[RHO] = std::numeric_limits<double>::quiet_NaN();
    };
    for (bool adaptive : {false, true})
    {
        auto mesh = make_mesh();
        const auto before = mesh.cells;
        const auto cfg = configuration(adaptive);
        hd_implicit::Integrator integrator(mesh, cfg);
        auto advance = [&] {
            return integrator.advance(mesh, 0.0, 0.1, hydro::SpatialScheme::Godunov,
                                      solver, {}, cfg, hooks);
        };
        if (adaptive)
        {
            CHECK(advance() == 0.05);
            for (const auto& cell : mesh.cells) check_thermodynamics(cell);
        }
        else
        {
            expect_throw<std::runtime_error>(advance);
            for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
            {
                CHECK(mesh.cells[flat].U.data == before[flat].U.data);
                CHECK(metadata(mesh.cells[flat]) == metadata(before[flat]));
            }
        }
    }
}

void recover_each_incoming_state()
{
    auto mesh = make_mesh();
    const auto cfg = configuration(false);
    hydro::ExactSolver solver(1e-10, 200);
    hd_user::UserHooks hooks;
    hooks.userbc = fill_ghosts;
    hd_source::SourceList sources{{"linear heating", hd_source::SourceKind::General,
        [](const MeshBlock& state, double, const config::SimulationConfig&,
           std::vector<Conserved>& rhs, double&) {
            for (std::size_t flat = 0; flat < state.cells.size(); ++flat)
            {
                check_thermodynamics(state.cells[flat]);
                if (!state.is_active(flat)) continue;
                const double heat = 0.5 * state.cells[flat].U.internal_energy();
                rhs[flat][ENERGY] += heat;
                rhs[flat][INTERNAL_ENERGY] += heat;
            }
        }, true}};
    hd_implicit::Integrator integrator(mesh, cfg);
    for (int step = 0; step < 2; ++step)
    {
        // Change U after workspace construction and between advances. Finite
        // but incorrect caches must also be recovered, not merely validated.
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            auto& cell = mesh.cells[flat];
            if (!mesh.is_active(flat))
            {
                cell.U[RHO] = std::numeric_limits<double>::quiet_NaN();
                continue; // The boundary hook must replace these before EoS.
            }
            cell.U[ENERGY] *= 1.5;
            cell.U[INTERNAL_ENERGY] *= 1.5;
            cell.temperature = 1e200;
            cell.n_baryon = cell.n_Ion = cell.n_e_free = 1e200;
            if constexpr (NCOMP > 0) cell.mean_molecular_weight = 1e-100;
        }
        const auto before = mesh.at_active(Index{}).U;
        const double dt = integrator.advance(mesh, step * 0.1, 0.1,
            hydro::SpatialScheme::Godunov, solver, sources, cfg, hooks);
        CHECK(dt == 0.1);
        for (const auto& cell : mesh.cells)
        {
            check_thermodynamics(cell);
            CHECK(near(cell.U[RHO], before[RHO], 1e-8));
            CHECK(near(cell.U[ENERGY], before[ENERGY] / (1.0 - 0.5*dt), 1e-8));
            for (int k = 0; k < NCOMP; ++k)
                CHECK(near(cell.U[COMP(k)], before[COMP(k)], 1e-8));
        }
    }
}

void eos_energy_coordinates()
{
    using hd_microphysics::Equation_of_State;
    for (auto model : {Equation_of_State::IdealGas_Default, Equation_of_State::IdealGas_Mixed,
                       Equation_of_State::Ionized, Equation_of_State::Ionized_Radiative})
    {
        if (NCOMP == 0 && model != Equation_of_State::IdealGas_Default) continue;
        hd_microphysics::eos_type = model;
        auto mesh = make_mesh();
        // These states use cgs stellar temperatures/sound speeds; use a
        // macroscopic box so this checks EoS coupling, not extreme acoustic CFL.
        for (int d = 0; d < NDIM; ++d) mesh.set_linear_uniform_faces(d, 0.0, 1e8);
        mesh.build_geometry();
        for (auto& cell : mesh.cells)
        {
            cell.U = {};
            cell.gamma = 5.0/3.0;
            cell.U[RHO] = 1e-10;
            if constexpr (NCOMP == 2)
            {
                cell.U[COMP(0)] = 0.65*cell.U[RHO];
                cell.U[COMP(1)] = 0.35*cell.U[RHO];
            }
            cell.temperature = 1e4;
            double pressure = 0;
            hd_microphysics::evaluate_eos(cell, hd_microphysics::EoSInput::RhoT, model, pressure);
        }
        const auto flat = mesh.active_to_storage.front();
        const auto& base = mesh.cells[flat];
        hd_implicit::detail::LogVariables variables(mesh, mesh, {flat});
        std::array<double, NVAR> dq{};
        dq[ENERGY] = 0.2;
        auto changed = variables.update(0, base, dq.data(), 1.0);
        CHECK(near(changed[ENERGY]/base.U[ENERGY], std::exp(0.2)));
        // The analytic direction must be the derivative of the actual update,
        // including density, velocity and independent species coordinates.
        dq[RHO] = 0.1;
        dq[VEL(0)] = 0.3;
        dq[INTERNAL_ENERGY] = 0.12;
        if constexpr (NCOMP > 0) dq[COMP(0)] = 0.05;
        const auto direction = variables.direction(0, base, dq.data());
        const double h = 1e-5;
        const auto plus = variables.update(0, base, dq.data(), h);
        const auto minus = variables.update(0, base, dq.data(), -h);
        for (int k = 0; k < NVAR; ++k)
            CHECK(std::abs((plus[k]-minus[k])/(2*h)-direction[k])
                  < 2e-9 * variables.row_scale(k));

        // A complete implicit step must add exactly rho*q*dt to both energies.
        const auto initial = base.U;
        const auto cfg = configuration(false);
        hydro::ExactSolver solver(1e-10, 256);
        hd_user::UserHooks hooks; hooks.userbc = fill_ghosts;
        constexpr double heating = 1e12;
        hd_implicit::Integrator integrator(mesh, cfg);
        const double dt = integrator.advance(mesh, 0.0, 0.01,
            hydro::SpatialScheme::Godunov, solver,
            {hd_source::make_constant_heating(heating)}, cfg, hooks);
        CHECK(dt == 0.01);
        for (auto index : mesh.active_to_storage)
        {
            const auto& cell = mesh.cells[index];
            CHECK(near(cell.U[ENERGY], initial[ENERGY] + initial[RHO]*heating*dt, 2e-9));
            CHECK(near(cell.U[INTERNAL_ENERGY], initial[INTERNAL_ENERGY] + initial[RHO]*heating*dt, 2e-9));
            CHECK(cell.temperature > 1e4);
            for (int k = 0; k < NCOMP; ++k) CHECK(near(cell.U[COMP(k)]/initial[COMP(k)], 1, 2e-9));
        }
    }
}
}

int main(int argc, char** argv)
{
    return run_test("implicit thermodynamics, line search and rollback", [&] {
        CHECK(NCOMP == 0 || NCOMP == 2);
        config::SimulationConfig cfg;
        if constexpr (NCOMP == 2)
        {
            CHECK(argc == 2);
            const auto directory = std::filesystem::absolute(argv[1]);
            std::filesystem::create_directories(directory);
            cfg.composition_activate_file = (directory / "active.txt").string();
            std::ofstream(cfg.composition_activate_file) << "h1 0 0\nhe4 0 0\n";
            hd_microphysics::eos_type = hd_microphysics::Equation_of_State::IdealGas_Mixed;
        }
        hd_microphysics::initialize_species_table(cfg);
        temperature_source_and_line_search();
        retry_and_rollback();
        recover_each_incoming_state();
        if constexpr (NCOMP == 2)
        {
            std::ofstream(cfg.composition_activate_file) << "h1 0 1\n";
            hd_microphysics::initialize_species_table(cfg);
        }
        eos_energy_coordinates();
    });
}
