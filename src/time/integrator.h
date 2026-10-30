#pragma once
#include "core.h"
#include "hydro/spatial.h"
#include <string>

namespace hd_core
{
    struct MeshBlock;
}

namespace hd_time
{
    enum class TimeIntegrator
    {
        Euler
        // 实现 RK2 后再增加 RK2。
    };

    TimeIntegrator parse_time_integrator(const std::string& name);

    // 统一入口：将活跃单元从 time 推进到 time + dt。
    void advance_one_step(
        hd_core::MeshBlock& mesh,
        double time,
        double dt,
        TimeIntegrator integrator,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver
    );

    // 显式 Euler。
    void advance_euler(
        hd_core::MeshBlock& mesh,
        double time,
        double dt,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver
    );

    // 以后在这里声明 advance_rk2()，使用相同参数。
}