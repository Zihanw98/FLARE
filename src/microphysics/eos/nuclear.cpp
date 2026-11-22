#include "nuclear.h"
#include "core/cell.h"
#include "physcon/constants.h"
#include "physcon/unit.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hd_microphysics
{
void nuclear_cold_state(double rho, double gamma_low,
                        double& pressure, double& energy_density)
{
    if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(gamma_low)
        || gamma_low <= 1.0 || gamma_low > nuclear_gamma_stiff)
        throw std::invalid_argument("Nuclear requires rho > 0 and 1 < gamma <= gamma_stiff");
    const double rho_n = nuclear_density_cgs / physcon::units.unit_density;
    const double p_n = nuclear_K_cgs * std::pow(nuclear_density_cgs, gamma_low)
                     / physcon::units.unit_pressure;
    const bool stiff = rho > rho_n;
    const double g = stiff ? nuclear_gamma_stiff : gamma_low;
    pressure = p_n * std::pow(rho / rho_n, g);
    energy_density = pressure / (g - 1.0);
    if (stiff)
        // Integrating de_cold/drho=P_cold/rho^2 gives this specific-energy offset.
        energy_density += rho * (p_n / rho_n)
            * (1.0 / (gamma_low - 1.0) - 1.0 / (nuclear_gamma_stiff - 1.0));
    if (!std::isfinite(pressure) || !std::isfinite(energy_density))
        throw std::invalid_argument("Nuclear cold state overflow");
}

void evaluate_nuclear_eos(hd_core::Cell& cell, EoSInput input, double& pressure)
{
    const double rho = cell.U[hd_core::RHO];
    double pc, uc;
    nuclear_cold_state(rho, cell.gamma, pc, uc);
    const double b = rho * physcon::R_g / cell.mean_molecular_weight;
    const double gm1 = nuclear_gamma_thermal - 1.0;
    double thermal;
    if (input == EoSInput::RhoT)
        thermal = b * cell.temperature / gm1;
    else if (input == EoSInput::RhoP)
        thermal = (pressure - pc) / gm1;
    else
        thermal = cell.U.internal_energy() - uc;

    // Only absorb subtraction roundoff at T=0; no physical energy floor or U repair.
    if (!std::isfinite(thermal) || thermal < -1e-12 * uc)
        throw std::invalid_argument("Nuclear requires nonnegative thermal energy (u >= u_cold, P >= P_cold)");
    thermal = std::max(0.0, thermal);
    cell.temperature = thermal * gm1 / b;
    if (input != EoSInput::RhoP)
        pressure = pc + gm1 * thermal;
    if (input != EoSInput::RhoE)
    {
        cell.U[hd_core::INTERNAL_ENERGY] = uc + thermal;
        cell.U[hd_core::ENERGY] = uc + thermal + cell.U.kinetic_energy();
    }
}
}
