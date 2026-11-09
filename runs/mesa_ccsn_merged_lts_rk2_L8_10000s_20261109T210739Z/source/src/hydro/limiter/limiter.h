/**
 * @file limiter.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include <string>
namespace hydro
{
    /** @brief Available one-dimensional slope-limiter formulas. */
    enum struct Limiter
    {
        none, ///< Unlimited centered derivative.
        minmod, ///< Minimum-magnitude same-sign derivative.
        MC, ///< Monotonized-central slope limiter.
        van_Leer, ///< Harmonic-mean same-sign slope limiter.
        MCbeta, ///< Generalized MC limiter with beta in [1, 2].
        superbee ///< Superbee compressive slope limiter.
    };

    // Parse once during setup; names are case-insensitive.
    /**
     * @brief Parse a case-insensitive slope-limiter name.
     * @param[in] name none, minmod, mc, vanleer/van_leer, mcbeta/mc_beta, or superbee.
     * @return The selected Limiter enum.
     * @throws std::invalid_argument If the name is unknown.
     */
    Limiter parse_limiter(const std::string& name);

    // Left/right derivatives dq/dx on a uniform grid. Inputs must be finite.
    // none returns the unlimited centered slope, not zero (PCM).
    // beta is used only by MCbeta and must be in [1, 2].
    /**
     * @brief Combine one-sided derivatives using the selected slope limiter.
     * @param[in] left Finite derivative on the left side of the cell.
     * @param[in] right Finite derivative on the right side of the cell.
     * @param[in] limiter Limiter formula; none returns the centered average.
     * @param[in] beta MCbeta parameter; finite and in [1, 2] when used.
     * @return Limited derivative, with zero at sign changes for all limiters except none.
     * @throws std::invalid_argument If derivatives, limiter enum, or the required beta are invalid.
     * @note none means an unlimited centered slope, not piecewise-constant reconstruction.
     */
    double limited_slope(double left, double right, Limiter limiter,
                         double beta = 1.5);
}
