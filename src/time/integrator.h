#pragma once
#include "hydro/spatial.h"
#include "source/source.h"
#include "config/config.h"
#include <string>
#include <array>
#include <vector>
#include <utility>
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

    // 调用前填好当前 time 的 ghost。dt_limit 是外层给出的有限正上限。
    // 返回实际使用的 dt；主循环必须 time += 返回值。
    // 当前 Euler 使用保守的 0 < cfl <= 0.5；默认 0.4。
    double advance_one_step(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        TimeIntegrator integrator,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        hd_source::SourceList source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl = 0.4
    );

    // 同一次计算获取 RHS 和波速，先限制 dt，再更新活跃单元。
    // 更新失败时保持原 mesh；ghost 的 U 和 gamma 不由积分器修改。
    double advance_euler(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        hd_source::SourceList source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl = 0.4
    );

    double advance_rk2(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        hd_source::SourceList source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl=0.4
    );
    std::pair<double,std::vector<hd_core::Conserved>> evaluate_RHS(
        const hd_core::MeshBlock& mesh,
        hydro::SpatialScheme scheme,
        const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver,
        const config::SimulationConfig& cfg,
        double cfl=0.4
    );
}
