#pragma once

namespace hd_core
{
    struct MeshBlock;
}

namespace hydro
{
    // 读取 wave_timestep；必须先完成本次 RHS 计算。
    // 零传播速率返回正无穷，需由外层 dt_limit 限制。
    double compute_cfl_dt(double cfl = 0.4);

    // 原来的单元状态估计仍可用于初始猜测；积分器使用上面的波速版本。
    double compute_cfl_dt(
        const hd_core::MeshBlock& mesh,
        double cfl = 0.4
    );
}
