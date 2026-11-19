#include "eos.h"
#include "core/cell.h"
#include "species.h"
#include "microphysics/composition/species.h"
#include "utils/utils.h"
#include "physcon/constants.h"
#include "hdmath.h"
#include <stdexcept>

namespace hd_microphysics
{
Equation_of_State eos_type = Equation_of_State::IdealGas_Default;

namespace
{
    void evaluate_eos_rho_temperature(hd_core::Cell& cell, Equation_of_State model, double& pressure)
    {
        auto rho = cell.U[hd_core::RHO];
        auto temperature = cell.temperature;
        const double kinetic = cell.U.kinetic_energy();
        if(model==Equation_of_State::IdealGas_Default)
        {
            const auto mu = cell.mean_molecular_weight;
            const auto gamma = cell.gamma;
            const auto pre = rho / mu * physcon::R_g * temperature;
            const auto internal_energy = pre / (gamma - 1);
            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            pressure = pre;
            return;
        }
        hd_microphysics::evaluate_cell_density(cell);

        if(model==Equation_of_State::IdealGas_Mixed)
        {
            const auto gamma = cell.gamma;
            auto pre = rho / cell.mean_molecular_weight * physcon::R_g * temperature;
            const auto internal_energy = pre / (gamma - 1);
            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            pressure = pre;
            return;
        }
        if(model==Equation_of_State::Ionized)
        {
            const auto gamma = cell.gamma;
            auto pre = rho / cell.mean_molecular_weight * physcon::R_g * temperature;
            const auto thermal_energy = pre / (gamma - 1);
            double ionization_energy = 0.0;
            for(int k=0; k<hd_core::NCOMP; k++)
            {
                const auto& species = species_table.species_active[k];
                const double n_k = cell.U[hd_core::COMP(k)] / (species.A * physcon::m_u);
                ionization_energy += n_k * species.energy_eV;
            }
            const auto internal_energy = thermal_energy + ionization_energy;
            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            pressure = pre;
            return;
        }
        if(model==Equation_of_State::Ionized_Radiative)
        {
            const auto gamma = cell.gamma;
            auto pre = rho / cell.mean_molecular_weight * physcon::R_g * temperature;
            const auto thermal_energy = pre / (gamma - 1);
            double ionization_energy = 0.0;
            for(int k=0; k<hd_core::NCOMP; k++)
            {
                const auto& species = species_table.species_active[k];
                const double n_k = cell.U[hd_core::COMP(k)] / (species.A * physcon::m_u);
                ionization_energy += n_k * species.energy_eV;
            }
            const double radiation_energy = physcon::a_Rad * std::pow(temperature,4);
            const auto internal_energy = thermal_energy + ionization_energy;
            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            pressure = pre + radiation_energy / 3.0;
            return;
        }
        
    }
    void evaluate_eos_rho_pressure(hd_core::Cell& cell, Equation_of_State model, double pressure)
    {
        const auto rho = cell.U[hd_core::RHO];
        const double kinetic = cell.U.kinetic_energy();
        if(model==Equation_of_State::IdealGas_Default)
        {
            const auto mu = cell.mean_molecular_weight;
            const auto gamma = cell.gamma;
            const auto temperature = pressure * mu / (rho * physcon::R_g);  //P = rho / mu * R_g(k/m_u) T
            const auto internal_energy = pressure / (gamma - 1);
            cell.temperature = temperature;
            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            return;
        }
        hd_microphysics::evaluate_cell_density(cell);
        if(model==Equation_of_State::IdealGas_Mixed)
        {
            const auto gamma = cell.gamma;
            const auto mu = cell.mean_molecular_weight;
            const auto temperature = pressure * mu / (rho * physcon::R_g); 
            const auto internal_energy = pressure / (gamma - 1);
            cell.temperature = temperature;
            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            return;
        }
        if(model==Equation_of_State::Ionized)
        {
            const auto gamma = cell.gamma;
            const auto mu = cell.mean_molecular_weight;
            const auto temperature = pressure * mu / (rho * physcon::R_g); 

            cell.temperature = temperature;

            const auto thermal_energy = pressure / (gamma - 1);
            double ionization_energy = 0.0;
            for(int k=0; k<hd_core::NCOMP; k++)
            {
                const auto& species = species_table.species_active[k];
                const double n_k = cell.U[hd_core::COMP(k)] / (species.A * physcon::m_u);
                ionization_energy += n_k * species.energy_eV;
            }
            const auto internal_energy = thermal_energy + ionization_energy;

            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            return;
        }
        if(model==Equation_of_State::Ionized_Radiative)
        {
            const auto gamma = cell.gamma;
            const auto mu = cell.mean_molecular_weight;
            const auto b = rho / mu * physcon::R_g;

            const double T0 = std::min(pressure/b,std::pow(3.0*pressure/physcon::a_Rad,0.25));
            // P(T) = rho/mu * R_g * T + 1/3 a T^4
            //f(T) = rho/mu *R_g  * T + 1/3 a T^4 - P
            //df(T)/dT = rho/mu * R_g + 4/3 a T^3
            auto f = [=](double T){
                return (b * T + physcon::a_Rad / 3.0 * std::pow(T,4) - pressure)/pressure;
            };
            auto df = [=](double T){
                return (b + 4.0/3.0 * physcon::a_Rad * std::pow(T,3))/pressure;
            };
            double ionization_energy = 0.0;
            for(int k=0; k<hd_core::NCOMP; k++)
            {
                const auto& species = species_table.species_active[k];
                const double n_k = cell.U[hd_core::COMP(k)] / (species.A * physcon::m_u);
                ionization_energy += n_k * species.energy_eV;
            }

            const auto temperature = hdmath::newton_raphson(f,df,T0,1e-10,100);
            cell.temperature = temperature;

        }
        throw std::logic_error("EoS RhoP is not implemented");
    }
    void evaluate_eos_rho_internal_energy(hd_core::Cell& cell, Equation_of_State model, double& pressure)
    {
        throw std::logic_error("EoS RhoE is not implemented");
    }
}
void evaluate_eos(hd_core::Cell& cell, EoSInput input,
                  Equation_of_State model, double& pressure)
{
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

    switch (input)
    {
        case EoSInput::RhoT:
            evaluate_eos_rho_temperature(cell, model, pressure);
            return;
        case EoSInput::RhoP:
            evaluate_eos_rho_pressure(cell, model, pressure);
            return;
        case EoSInput::RhoE:
            evaluate_eos_rho_internal_energy(cell, model, pressure);
            return;
    }
    throw std::invalid_argument("Invalid EoS input enum");
}

Equation_of_State parse_equation_of_state(const std::string& name)
{
    if (utils::iequals(name, "IdealGas_Default"))
        return Equation_of_State::IdealGas_Default;
    if (utils::iequals(name, "IdealGas_Mixed"))
        return Equation_of_State::IdealGas_Mixed;
    if (utils::iequals(name, "Ionized"))
        return Equation_of_State::Ionized;
    if (utils::iequals(name, "Ionized_Radiative"))
        return Equation_of_State::Ionized_Radiative;

    throw std::invalid_argument("Unknown equation of state: " + name);
}
}
