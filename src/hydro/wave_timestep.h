/**
 * @file wave_timestep.h
 * @brief Shared mesh-wide CFL rate.
 */
#pragma once

namespace hydro
{
    // 当前网格的 max_cells sum_d(max(a_left, a_right) / dx_d)。
    // 单位为 1 / 时间；最终 dt = cfl / max_inverse_dt。
    // 在 timestep_CFL.cpp 中唯一定义；每次空间 RHS 计算时重置。
    // 仅用于串行网格汇总，界面波速和速度由求解器引用参数返回。
    extern double max_inverse_dt;
}
