/**
 * @file cell.h
 * @brief Physical states and mesh storage.
 */
#pragma once
#include "core.h"
namespace hd_core
{
    /** @brief A conserved cell state and its ideal-gas adiabatic index. */
    struct Cell
    {
        Conserved U;
        double gamma = 5.0/3;
        // Accumulated density of energy added/removed by cold-cell synchronization.
        double dual_energy_correction = 0.0;
        void synchronize_dual_energy()
        {
            U.validate(gamma);
            if (U.total_energy_reliable())
                U[INTERNAL_ENERGY] = U[ENERGY] - U.kinetic_energy();
            else
            {
                const double energy = U.kinetic_energy() + U[INTERNAL_ENERGY];
                dual_energy_correction += energy - U[ENERGY];
                U[ENERGY] = energy;
            }
        }
        /**
         * @brief Construct a cell with zero conserved components and gamma = 5/3.
         * @note Default construction does not validate the state.
         */
        Cell() = default;
        /**
         * @brief Construct and validate a cell from conserved variables.
         * @param[in] state Conserved state to copy.
         * @param[in] g Finite ideal-gas adiabatic index greater than one.
         * @throws std::invalid_argument If the conserved state or adiabatic index is invalid.
         */
        explicit Cell(const Conserved& state, double g = 5.0/3)
            : U(state), gamma(g)
        {
            U.validate(gamma);
        }
        /**
         * @brief Construct a cell by converting primitive variables.
         * @param[in] state Primitive state to convert.
         * @param[in] g Ideal-gas adiabatic index greater than one.
         * @pre The input is a valid primitive state for g.
         * @note Uses Primitives::to_Conserved(); a zero-density input produces a zero conserved state and extra components are not transferred.
         * @throws std::invalid_argument If the non-vacuum conversion fails validation.
         */
        explicit Cell(const Primitives& state, double g = 5.0/3)
            : U(state.to_Conserved(g)), gamma(g)
        {
        }
    };
}
