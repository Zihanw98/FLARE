#include "cell.h"
#include "microphysics/composition/species.h"
#include "microphysics/eos/eos.h"
#include "physcon/constants.h"
#include <cmath>
#include <stdexcept>

namespace hd_core
{
    Cell::Cell(const Conserved& state, double g)
        : U(state), gamma(g)
    {
        U.validate();
        if constexpr (NCOMP > 0)
            if (U[RHO] > 0.0)
                hd_microphysics::evaluate_cell_density(*this);
    }

    Cell::Cell(const Primitives& state, double g)
        : gamma(g)
    {
        Set_Conserved_From_Primitives(state);
        if constexpr (NCOMP > 0)
            if (U[RHO] > 0.0)
                hd_microphysics::evaluate_cell_density(*this);
    }

    Primitives Cell::Get_Primitives(hd_microphysics::Equation_of_State model) const
    {
        U.validate();
        Primitives result{};
        for (int k = 0; k < NUSER; ++k)
            result[CONUSR(k)] = U[CONUSR(k)];

        const double rho = U[RHO];
        if (rho == 0.0) return result; // Vacuum does not enter the positive-density EoS.

        Cell work = *this;
        double pressure = 0.0;
        hd_microphysics::evaluate_eos(work, hd_microphysics::EoSInput::RhoE, model, pressure);

        result[RHO] = rho;
        for (int d = 0; d < NDIM; ++d)
            result[VEL(d)] = work.U[MOMDEN(d)] / rho;
        result[PRE] = pressure;
        result[INTERNAL_ENERGY] = work.U.internal_energy() / rho;
        for (int k = 0; k < NCOMP; ++k)
            result[COMP(k)] = work.U[COMP(k)] / rho;
        result.validate(gamma);
        return result;
    }

    double Cell::Get_Pressure_From_Updated_State(
        hd_microphysics::Equation_of_State model) const
    {
        const double rho = U[RHO];
        if (!std::isfinite(rho) || rho < 0.0)
            throw std::invalid_argument(
                "Get_Pressure_From_Updated_State requires finite nonnegative rho");
        if (rho == 0.0) return 0.0;

        using hd_microphysics::Equation_of_State;
        switch (model)
        {
            case Equation_of_State::IdealGas_Default:
            case Equation_of_State::IdealGas_Mixed:
            case Equation_of_State::Ionized:
            case Equation_of_State::Ionized_Radiative:
                break;
            default:
                throw std::invalid_argument("Invalid equation of state enum");
        }
        if (!std::isfinite(temperature) || temperature <= 0.0
            || !std::isfinite(mean_molecular_weight) || mean_molecular_weight <= 0.0)
            throw std::invalid_argument(
                "Get_Pressure_From_Updated_State requires finite positive T and mu");

        double pressure = rho / mean_molecular_weight * physcon::R_g * temperature;
        if (model == Equation_of_State::Ionized_Radiative)
        {
            const double T2 = temperature * temperature;
            pressure += physcon::a_Rad * T2 * T2 / 3.0;
        }
        if (!std::isfinite(pressure) || pressure <= 0.0)
            throw std::invalid_argument(
                "Get_Pressure_From_Updated_State requires finite positive pressure");
        return pressure;
    }

    Primitives Cell::Get_Primitives_From_Updated_State(
        hd_microphysics::Equation_of_State model) const
    {
        U.validate();
        Primitives result{};
        for (int k = 0; k < NUSER; ++k)
            result[CONUSR(k)] = U[CONUSR(k)];

        const double rho = U[RHO];
        if (rho == 0.0) return result;

        result[RHO] = rho;
        for (int d = 0; d < NDIM; ++d)
            result[VEL(d)] = U[MOMDEN(d)] / rho;
        result[PRE] = Get_Pressure_From_Updated_State(model);
        result[INTERNAL_ENERGY] = U.internal_energy() / rho;
        for (int k = 0; k < NCOMP; ++k)
            result[COMP(k)] = U[COMP(k)] / rho;
        result.validate(gamma);
        return result;
    }

    Conserved Cell::to_Conserved(const Primitives& state,
                                hd_microphysics::Equation_of_State model) const
    {
        // Only rho, v, P, X and user fields are inputs; stale specific energy is ignored.
        Primitives input = state;
        input[INTERNAL_ENERGY] = 0.0;
        input.validate(gamma);

        Cell work = *this;
        work.U = Conserved{};
        for (int k = 0; k < NUSER; ++k)
            work.U[CONUSR(k)] = input[CONUSR(k)];

        const double rho = input[RHO];
        if (rho == 0.0) return work.U;

        work.U[RHO] = rho;
        for (int d = 0; d < NDIM; ++d)
            work.U[MOMDEN(d)] = rho * input[VEL(d)];
        for (int k = 0; k < NCOMP; ++k)
            work.U[COMP(k)] = rho * input[COMP(k)];

        double pressure = input[PRE];
        hd_microphysics::evaluate_eos(work, hd_microphysics::EoSInput::RhoP, model, pressure);
        work.U.validate();
        return work.U;
    }

    void Cell::Set_Conserved_From_Primitives(const Primitives& state,
                                            hd_microphysics::Equation_of_State model)
    {
        U = to_Conserved(state, model);
    }

    double Cell::get_sound_speed(hd_microphysics::Equation_of_State model) const
    {
        if(model==hd_microphysics::Equation_of_State::IdealGas_Default || model == hd_microphysics::Equation_of_State::IdealGas_Mixed || model==hd_microphysics::Equation_of_State::Ionized)
        {
            double pressure = Get_Pressure_From_Updated_State();
            return std::sqrt(gamma * pressure / U[hd_core::RHO]);
        }
        if(model==hd_microphysics::Equation_of_State::Ionized_Radiative)
        {
            double pressure_gas = U[hd_core::RHO] / mean_molecular_weight * physcon::R_g * temperature;
            double pressure_rad = physcon::a_Rad / 3 * std::pow(temperature,4);
            const double A = pressure_gas + 4.0 * pressure_rad;
            const double B = pressure_gas / (gamma - 1) + 12.0 * pressure_rad;
            
            return std::sqrt((pressure_gas + A* A /B )/ U[hd_core::RHO]);
        }
        throw std::runtime_error("Haven't implemented sound speed for this EoS");
    }
}
