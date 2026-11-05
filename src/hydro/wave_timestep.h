/**
 * @file wave_timestep.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once

namespace hydro
{
    /** @brief Shared wave-speed measurements for the most recent serial RHS evaluation. */
    struct WaveTimestep
    {
        // 最近一次 flux() 求解的 max(|S_left|, |S_right|)，每个面覆盖一次。
        double face_speed = 0.0;

        // 当前网格的 max_cells sum_d(max(a_left, a_right) / dx_d)。
        // 单位为 1 / 时间；最终 dt = cfl / max_inverse_dt。
        double max_inverse_dt = 0.0;
    };

    // 在 timestep_CFL.cpp 中唯一定义；需要访问时包含本头文件。
    // 每次 compute_godunov_rhs() 开始时重置，仅适用于串行网格计算。
    extern WaveTimestep wave_timestep;
}
