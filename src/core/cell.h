/**
 * @file cell.h
 * @brief Physical states and mesh storage.
 */
#pragma once
#include "core.h"

namespace hd_microphysics
{
    enum struct Equation_of_State;
    extern Equation_of_State eos_type;
}

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
         * @note For positive density, initializes composition-derived members when
         * HD_NCOMP > 0; the shared active species table must be initialized first.
         * @throws std::invalid_argument If the conserved state fails the generic checks.
         */
        explicit Cell(const Conserved& state, double g = 5.0/3);
        /**
         * @brief Construct a cell by converting primitive variables.
         * @param[in] state Primitive state to convert.
         * @param[in] g Ideal-gas adiabatic index greater than one.
         * @pre The input is a valid primitive state for g.
         * @note Sets gamma before building U through the Setter with the global EoS model.
         * For positive density, then initializes composition-derived members when
         * HD_NCOMP > 0; the shared active species table must be initialized first. Temperature
         * retains its default value; user-defined conserved fields are copied unchanged.
         * @throws std::invalid_argument If the non-vacuum conversion fails validation.
         */
        explicit Cell(const Primitives& state, double g = 5.0/3);
        /**
         * @brief Recover primitives through RhoE on a private Cell copy.
         * @note All members of this Cell remain unchanged. The returned specific
         * internal energy uses the dual-energy selection, including EoS energy
         * contributions. User-defined conserved fields are copied unchanged.
         */
        Primitives Get_Primitives(
            hd_microphysics::Equation_of_State model = hd_microphysics::eos_type) const;

        /**
         * @brief Recover primitives using the stored temperature and molecular weight.
         * @pre For non-vacuum cells, temperature and mean_molecular_weight must
         * already match the current U and model (e.g. after RhoE on this Cell).
         * @note Does not evaluate EoS, recompute number densities, or modify any
         * members. Pressure includes radiation only for Ionized_Radiative;
         * specific internal energy still uses U's dual-energy selection.
         * User-defined conserved fields are copied unchanged. Cache freshness
         * is the caller's responsibility and is not checked here.
         */
        Primitives Get_Primitives_From_Updated_State(
            hd_microphysics::Equation_of_State model = hd_microphysics::eos_type) const;

        /**
         * @brief Return pressure from rho and the stored temperature and molecular weight.
         * @pre For non-vacuum cells, temperature and mean_molecular_weight must
         * already match the current U and model; freshness is not checked here.
         * @note Adds radiation pressure for Ionized_Radiative. Returns zero for
         * vacuum. Does not evaluate EoS, scan composition, or modify any members.
         */
        double Get_Pressure_From_Updated_State(
            hd_microphysics::Equation_of_State model = hd_microphysics::eos_type) const;

        /**
         * @brief Build conserved variables through RhoP on a private Cell copy.
         * @param state Density, velocity, pressure and mass fractions to convert.
         * Its stored specific internal energy is ignored and recomputed by EoS.
         * @note Uses this Cell's gamma and default-model molecular weight without
         * modifying any members. User-defined conserved fields are copied unchanged.
         */
        Conserved to_Conserved(
            const Primitives& state,
            hd_microphysics::Equation_of_State model = hd_microphysics::eos_type) const;

        /**
         * @brief Replace only U after successful primitive-to-conserved conversion.
         * @note Temperature, molecular weight, number densities, gamma and the
         * dual-energy correction remain unchanged, even when conversion throws.
         * Thermodynamic caches may therefore need an explicit EoS refresh later.
         */
        void Set_Conserved_From_Primitives(
            const Primitives& state,
            hd_microphysics::Equation_of_State model = hd_microphysics::eos_type);
    };
}
