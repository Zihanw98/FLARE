/**
 * @file timestep_CFL.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once

namespace hd_core
{
    struct MeshBlock;
}

namespace hydro
{
    // 读取 max_inverse_dt；必须先完成本次 RHS 计算。
    // 零传播速率返回正无穷，需由外层 dt_limit 限制。
    /**
     * @brief Compute the CFL timestep from the most recently recorded face-wave rate.
     * @param[in] cfl Finite positive dimensionless CFL coefficient.
     * @return cfl / max_inverse_dt, or positive infinity when the rate is zero.
     * @pre The current RHS evaluation has completed and updated max_inverse_dt.
     * @throws std::invalid_argument If cfl is invalid or the recorded rate is non-finite or negative.
     * @note Does not inspect the mesh. The caller must cap an infinite result with an external timestep limit.
     */
    double compute_cfl_dt(double cfl = 0.4);

    // 单元状态估计用于初始化诊断；实际积分器从界面波速计算步长限制。
    /**
     * @brief Estimate a CFL timestep from cell-centered velocities and sound speeds.
     * @param[in] mesh Allocated mesh with valid active states and positive cell widths.
     * @param[in] cfl Finite positive dimensionless CFL coefficient.
     * @return cfl divided by the largest active-cell sum of (absolute velocity + sound speed) / cell width.
     * @pre Geometry and active cell states are initialized.
     * @throws std::invalid_argument If cfl or a converted gas state is invalid.
     * @throws std::runtime_error If a rate is non-finite or no positive active-cell rate is found.
     * @note Skips vacuum cells. Unlike the wave-record overload, an all-vacuum mesh raises an error; time integrators use the face-wave rate instead.
     */
    double compute_cfl_dt(
        const hd_core::MeshBlock& mesh,
        double cfl = 0.4
    );
}
