#include "eos.h"
#include "core/cell.h"
#include "microphysics/composition/species.h"
#include "utils/utils.h"
#include "physcon/constants.h"
#include <stdexcept>

namespace hd_microphysics
{
    
namespace
{
    // 在下面对应函数中补充物理计算，保持 EoSFunction 的公共签名不变。
    void ideal_gas_default(hd_core::Cell&)
    {
        throw std::logic_error("EoS IdealGas_Default is not implemented");
    }

    void ideal_gas_mixed(hd_core::Cell&)
    {
        throw std::logic_error("EoS IdealGas_Mixed is not implemented");
    }

    void ionized(hd_core::Cell&)
    {
        throw std::logic_error("EoS Ionized is not implemented");
    }

    void ionized_radiative(hd_core::Cell&)
    {
        throw std::logic_error("EoS Ionized_Radiative is not implemented");
    }
}
    void evaluate_eos(hd_core::Cell& cell,EoSInput input,Equation_of_State model)
    {
        
    }
    void evaluate_eos_rho_T(hd_core::Cell& cell, Equation_of_State model)
    {
        auto rho = cell.U[hd_core::RHO];
        auto temperature = cell.temperature;
        
        if(model==Equation_of_State::IdealGas_Default)
        {
            const auto mu = cell.mean_molecular_weight;
            const double kinetic = cell.U.kinetic_energy();
            const auto gamma = cell.gamma;
            const auto pre = rho / mu * physcon::R_g * temperature;
            const auto internal_energy = pre / (gamma - 1);
            cell.U[hd_core::INTERNAL_ENERGY] = internal_energy;
            cell.U[hd_core::ENERGY] = internal_energy + kinetic;
            return;
        }
        double mu_inv=0.0;
        //Calculate mu first for every three EoS left
        double n_e = 0.0;
        double n_I = 0.0;
        double mu_e = 0.0;
        double mu_I = 0.0;
        for(int k=0; k<hd_core::NCOMP; k++)
        {
            const auto A_k = hd_microphysics::species_table.species_active[k].A;
            const auto Z_k = hd_microphysics::species_table.species_active[k].charge;
            const auto X_k = cell.U[hd_core::COMP(k)]/rho;
            n_e += Z_k*X_k/A_k;
            mu_inv += (1+Z_k)*X_k/A_k;
            //1/mu = \Sigma (1+Z_k)*X_k/A_k
        }
        cell.mean_molecular_weight = 1.0/mu_inv;
        
        
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

EoSFunction get_eos_function(Equation_of_State type)
{
    switch (type)
    {
        case Equation_of_State::IdealGas_Default: return &ideal_gas_default;
        case Equation_of_State::IdealGas_Mixed: return &ideal_gas_mixed;
        case Equation_of_State::Ionized: return &ionized;
        case Equation_of_State::Ionized_Radiative: return &ionized_radiative;
    }

    throw std::invalid_argument("Invalid equation of state enum");
}
}
