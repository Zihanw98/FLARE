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
        RhoT,   // 密度、温度是已知量
        RhoP,   // 密度、压力是已知量
        RhoE    // 密度、内能密度是已知量（从 U 按双能量规则选择）
    };

    // 启动时从配置设置；演化期间保持不变，需要时显式传给 evaluate_eos。
    extern Equation_of_State eos_type;

    // 启动时解析配置名称（忽略大小写），未知名称抛出 invalid_argument。
    Equation_of_State parse_equation_of_state(const std::string& name);

    // 唯一计算入口。RhoT/RhoE 输出 pressure；RhoP 读取 pressure，不修改它。
    // 密度和组分来自 cell；RhoT 的温度来自 cell.temperature。
    // 非真空状态要求 rho、mu、T、P 为有限正数，gamma > 1；RhoE 要求
    // 扣除电离能后仍有正的能量。非法状态抛出 invalid_argument。
    // Nuclear 不依赖组分；允许 T=0，但要求 u>=u_cold、P>=P_cold。
    // 尚未实现的输入模式抛出 logic_error。
    void evaluate_eos(
        hd_core::Cell& cell,
        EoSInput input,
        Equation_of_State model,
        double& pressure);

    // 校验 U，并通过 RhoE 更新温度和模型需要的数密度、mu；不修改 U。
    // 真空跳过 EoS；不做组分归一性检查，也不进行双能量同步。
    void refresh_thermodynamics(
        hd_core::Cell& cell,
        Equation_of_State model = eos_type);

    // 更新所有 cell（包括 ghost）；调用前应先填充当前阶段的边界。
    // 状态校验或数值求解异常附带 flat 下标；已更新的热力学量不回滚。
    void refresh_thermodynamics(
        hd_core::MeshBlock& mesh,
        Equation_of_State model = eos_type);
}
