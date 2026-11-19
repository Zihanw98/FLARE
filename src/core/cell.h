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
        double temperature = 1.0;
        double n_baryon = 0.0;              //baryon number density  
        double Y_e_free = 0.0;             //mu_e ---> how much molecular weight does the gas has when it has one electron, or baryons per electron    Y_e = 1/ mu_e : free electron per mean molecular weight(per baryon)
        //mu_e can be blow up if it's purely neutral gas so instead we store Y_e which is its inverse
        double Y_Ion = 0.0;                 // number of ions per baryon(molecular weight)
        double n_e_free = 0.0;
        double n_Ion = 0.0;
        double mean_molecular_weight = 1.0;   //mu, mean molecular weight
        // Accumulated density of energy added/removed by cold-cell synchronization.
        double dual_energy_correction = 0.0;
        void synchronize_dual_energy()
        {
            U.validate();
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
         * @brief Construct a cell and check its conserved-state arithmetic.
         * @param[in] state Conserved state to copy.
         * @param[in] g Adiabatic index to store; model-specific validation is separate.
         * @throws std::invalid_argument If the conserved state fails the generic checks.
         */
        explicit Cell(const Conserved& state, double g = 5.0/3)
            : U(state), gamma(g)
        {
            U.validate();
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
