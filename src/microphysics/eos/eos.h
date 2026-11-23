#pragma once
#include <string>
namespace hd_core
{
    struct Cell;
    struct MeshBlock;
}

namespace hd_microphysics
{
    enum struct Equation_of_State
    {
        IdealGas_Default,
        IdealGas_Mixed,
        Ionized,
        Ionized_Radiative,
        Nuclear         // Toy hybrid: cold piecewise polytrope + thermal gamma law.
    };
    enum class EoSInput
    {
        RhoT,   // density and temperature are known
        RhoP,   // density and pressure are known
        RhoE    // density and internal energy are known (internal energy calculated from dual energy method). Also the way to get primitives from convserved quantities.
    };

    // global eos
    extern Equation_of_State eos_type;

    // 启动时解析配置名称（忽略大小写），未知名称抛出 invalid_argument。
    Equation_of_State parse_equation_of_state(const std::string& name);

    //calculate 
    void evaluate_eos(
        hd_core::Cell& cell,
        EoSInput input,
        Equation_of_State model,
        double& pressure);

    //update thermodynamics states for certain cell by using rhoE, only temperature & number densities are updated
    void refresh_thermodynamics(
        hd_core::Cell& cell,
        Equation_of_State model = eos_type);

    //update thermo states of entire mesh of cells
    void refresh_thermodynamics(
        hd_core::MeshBlock& mesh,
        Equation_of_State model = eos_type);
}
