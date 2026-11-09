/**
 * @file flux.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"

namespace hydro
{
    /**
     * @brief Compute the ideal-gas Euler flux along one coordinate axis.
     * @param[in] W Primitive state; valid gas or zero-density vacuum.
     * @param[in] direction Flux-normal axis in [0, hd_core::NDIM).
     * @param[in] gamma Finite ideal-gas adiabatic index greater than one.
     * @return Mass, momentum, and energy fluxes; zero for vacuum. Extra-component fluxes remain zero.
     * @throws std::out_of_range If direction is invalid.
     * @throws std::invalid_argument If the state or gas conversion fails validation.
     */
    hd_core::Conserved physical_flux(
        const hd_core::Primitives& W,
        int direction,
        double gamma
    );
}