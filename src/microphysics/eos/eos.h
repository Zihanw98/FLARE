#pragma once
#include <string>

namespace hd_core
{
    struct Cell;
}

namespace hd_microphysics
{
    struct SpeciesTable;

    enum struct Equation_of_State
    {
        IdealGas_Default,
        IdealGas_Single,
        IdealGas_Mixed,
        Ionized,
        Ionized_Radiative
    };

    // 所有 EoS 使用相同签名：更新当前 cell，物种表只读。
    // 调用者保存此指针，之后直接 eos(cell, species)，无需再判断 EoS 类型。
    using EoSFunction = void (*)(hd_core::Cell& cell, const SpeciesTable& species);

    // 启动时解析配置名称（忽略大小写），未知名称抛出 invalid_argument。
    Equation_of_State parse_equation_of_state(const std::string& name);

    // 按枚举绑定具体实现；返回值非空，非法枚举抛出 invalid_argument。
    // 当前五种物理实现均为占位函数，调用时抛出 logic_error，等待后续实现。
    EoSFunction get_eos_function(Equation_of_State type);
}
