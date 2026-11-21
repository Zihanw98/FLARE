/**
 * @file flux.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"

namespace hydro
{
    /**
     * @brief Convert W through the global EoS, then compute the Euler flux.
     * @param[in] W Primitive state; valid gas or zero-density vacuum.
     * @param[in] direction Flux-normal axis in [0, hd_core::NDIM).
     * @param[in] gamma Finite ideal-gas adiabatic index greater than one.
     * @return Mass, momentum, energy and composition fluxes; zero for vacuum.
     * User-defined field fluxes remain zero.
     * @throws std::out_of_range If direction is invalid.
     * @throws std::invalid_argument If the state or gas conversion fails validation.
     */
    hd_core::Conserved physical_flux(
        const hd_core::Primitives& W,
        int direction,
        double gamma
    );

    /**
     * @brief Compute the Euler flux from an already converted state, without EoS evaluation.
     * @pre W and U describe the same validated density, velocity, pressure and composition.
     * U supplies both energy densities; W's specific internal energy is not read.
     * @param[in] direction Flux-normal axis in [0, hd_core::NDIM).
     * @return Mass, momentum, energy and composition fluxes; zero for vacuum.
     * User-defined field fluxes remain zero. Neither input is modified.
     * @throws std::out_of_range If direction is invalid.
     * @throws std::invalid_argument If U fails generic conserved-state validation.
     */
    hd_core::Conserved physical_flux(
        const hd_core::Primitives& W,
        const hd_core::Conserved& U,
        int direction
    );
}
