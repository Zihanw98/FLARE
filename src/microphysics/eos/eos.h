#pragma once
#include <string>

namespace hd_core
{
    struct Cell;
}

namespace hd_microphysics
{
    enum struct Equation_of_State
    {
        IdealGas_Default,
        IdealGas_Mixed,
        Ionized,
        Ionized_Radiative
    };
    enum class EoSInput
    {
        RhoT,   // 密度、温度是已知量
        RhoP,   // 密度、压力是已知量
        RhoE    // 密度、比内能是已知量
    };

    // EOS implementations that need composition read the shared species_table.
    using EoSFunction = void (*)(hd_core::Cell& cell);

    // 启动时解析配置名称（忽略大小写），未知名称抛出 invalid_argument。
    Equation_of_State parse_equation_of_state(const std::string& name);

    // 按枚举绑定具体实现；返回值非空，非法枚举抛出 invalid_argument。
    // 当前五种物理实现均为占位函数，调用时抛出 logic_error，等待后续实现。
    EoSFunction get_eos_function(Equation_of_State type);
    void evaluate_eos(
        hd_core::Cell& cell,
        EoSInput input,
        Equation_of_State model);
    void evaluate_eos_rho_T(hd_core::Cell& cell, Equation_of_State model);
}
