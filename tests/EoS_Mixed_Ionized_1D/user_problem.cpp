#include "user_config.h"
#include "user/user_api.h"
#include "config/config.h"
#include "microphysics/composition/species.h"
#include "microphysics/eos/eos.h"
#include "physcon/constants.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

static_assert(hd_core::NDIM == 1 && hd_core::NCOMP == 2);
namespace
{
void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    using namespace hd_core;
    const bool heating = cfg.constant_heating_enabled;
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        auto& cell = mesh.at_active({i});
        cell.gamma = cfg.gamma;
        const double x = mesh.cell_center(0, i + mesh.nghost[0]);
        const double phase = 2*physcon::pi*(x-cfg.lower[0])/mesh.active_length(0);
        Primitives w{};
        w[RHO] = heating ? 1e-10 : 1e-10*(1+0.2*std::sin(phase));
        w[VEL(0)] = heating ? 0 : 1e5;
        w[COMP(0)] = heating ? 0.65 : 0.5+0.3*std::sin(phase);
        w[COMP(1)] = 1-w[COMP(0)];
        w[PRE] = 100;
        cell.Set_Conserved_From_Primitives(w);
        if (heating)
        {
            cell.temperature = 1e4;
            double pressure = 0;
            hd_microphysics::evaluate_eos(cell, hd_microphysics::EoSInput::RhoT,
                                         hd_microphysics::eos_type, pressure);
        }
    }
}

// Check the refreshed cache independently, and expose it without changing the
// framework CSV format. This file always corresponds to the latest snapshot.
void thermodynamics(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    using namespace hd_core;
    std::ofstream out(cfg.output_directory + "/thermodynamics_last.csv");
    if (!out) throw std::runtime_error("Cannot write thermodynamic diagnostics");
    out << std::setprecision(17) << "x,temperature,mu,n_e\n";
    for (auto flat : mesh.active_to_storage)
    {
        const auto& cell = mesh.cells[flat];
        double particles = 0, electrons = 0, ion = 0;
        for (int k = 0; k < NCOMP; ++k)
        {
            const auto& sp = hd_microphysics::species_table.species_active[k];
            const double n = cell.U[COMP(k)]/(sp.A*physcon::m_u);
            particles += (1+sp.charge)*n;
            electrons += sp.charge*n;
            if (hd_microphysics::eos_type == hd_microphysics::Equation_of_State::Ionized)
                ion += n*sp.energy_eV*physcon::eV;
        }
        const double expected = (cell.U.internal_energy()-ion)*(cell.gamma-1)/(particles*physcon::k_B);
        if (std::abs(cell.temperature/expected-1)>1e-9 ||
            std::abs(cell.n_e_free-electrons)>1e-9*particles ||
            std::abs(cell.mean_molecular_weight*particles*physcon::m_u/cell.U[RHO]-1)>1e-9)
            throw std::runtime_error("EoS cache disagrees with independent particle/energy calculation");
        out << mesh.cell_center(0, mesh.coordinates(flat)) << ',' << cell.temperature
            << ',' << cell.mean_molecular_weight << ',' << cell.n_e_free << '\n';
    }
}
}
hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    hooks.useroutput = thermodynamics;
    return hooks;
}
