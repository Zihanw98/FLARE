#pragma once
#include "eos.h"

namespace hd_microphysics
{
    // Toy hybrid EoS, not a tabulated nuclear-matter EoS. Constants are in CGS.
    // cell.gamma is the low-density cold exponent; the demo uses 4/3.
    inline constexpr double nuclear_density_cgs = 2.0e14;
    inline constexpr double nuclear_K_cgs = 4.0e14;
    inline constexpr double nuclear_gamma_stiff = 2.5;
    inline constexpr double nuclear_gamma_thermal = 1.5;

    // Returns cold pressure and cold INTERNAL ENERGY DENSITY in code units.
    // Both are continuous at nuclear_density_cgs; the derivative changes there.
    void nuclear_cold_state(double rho, double gamma_low,
                            double& pressure, double& energy_density);

    // Called by evaluate_eos after common input validation. RhoE never changes U.
    // T is a toy thermal temperature, with cV = R_g/[mu*(gamma_th-1)].
    // The exactly cold state T=0 is admissible. Negative thermal energy is not.
    void evaluate_nuclear_eos(hd_core::Cell& cell, EoSInput input, double& pressure);
}
