#include "eos.h"
#include "core/cell.h"
#include "microphysics/composition/species.h"
#include "utils/utils.h"
#include <stdexcept>

namespace hd_microphysics
{
namespace
{
    // 在下面对应函数中补充物理计算，保持 EoSFunction 的公共签名不变。
    void ideal_gas_default(hd_core::Cell&, const SpeciesTable&)
    {
        throw std::logic_error("EoS IdealGas_Default is not implemented");
    }

    void ideal_gas_single(hd_core::Cell&, const SpeciesTable&)
    {
        throw std::logic_error("EoS IdealGas_Single is not implemented");
    }

    void ideal_gas_mixed(hd_core::Cell&, const SpeciesTable&)
    {
        throw std::logic_error("EoS IdealGas_Mixed is not implemented");
    }

    void ionized(hd_core::Cell&, const SpeciesTable&)
    {
        throw std::logic_error("EoS Ionized is not implemented");
    }

    void ionized_radiative(hd_core::Cell&, const SpeciesTable&)
    {
        throw std::logic_error("EoS Ionized_Radiative is not implemented");
    }
}

Equation_of_State parse_equation_of_state(const std::string& name)
{
    if (utils::iequals(name, "IdealGas_Default"))
        return Equation_of_State::IdealGas_Default;
    if (utils::iequals(name, "IdealGas_Single"))
        return Equation_of_State::IdealGas_Single;
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
        case Equation_of_State::IdealGas_Single: return &ideal_gas_single;
        case Equation_of_State::IdealGas_Mixed: return &ideal_gas_mixed;
        case Equation_of_State::Ionized: return &ionized;
        case Equation_of_State::Ionized_Radiative: return &ionized_radiative;
    }

    throw std::invalid_argument("Invalid equation of state enum");
}
}
