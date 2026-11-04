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
        Euler,
        RK2
    };

    TimeIntegrator parse_time_integrator(const std::string& name);

    // 对活跃单元执行 U += dt * rhs，并逐个校验候选状态。
    // rhs 与 mesh.cells 等长；dt 必须有限，零表示不更新。
    // 全部校验成功后才提交；失败时保持原 mesh，ghost、gamma 和几何不变。
    void add_scaled_rhs(
        hd_core::MeshBlock& mesh,
        const std::vector<hd_core::Conserved>& rhs,
        double dt
    );

    // 调用前填好当前 time 的 ghost。dt_limit 是外层给出的有限正上限。
    // 返回实际使用的 dt；主循环必须 time += 返回值。
    // cfl 为目标值；dt_next 初始设为 +infinity，RK 方法可返回下一步建议。
    // Euler 每步独立估计，dt_next 返回 +infinity 表示没有额外建议。
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
        double cfl,
        double& dt_next
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

    // cfl 用于选步长；cfg.cfl_max 用于接受/拒绝。两个阶段共用同一 dt。
    // dt_next：输入上一步建议（首次设为 +infinity），成功后写入下一步建议。
    double advance_rk2(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        hd_source::SourceList source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl,
        double& dt_next
    );
    // 保留 {目标步长, RHS} 返回形式；引用输出用于分别检查 CFL 和源项硬限制。
    std::pair<double, std::vector<hd_core::Conserved>> evaluate_RHS(
        const hd_core::MeshBlock& mesh,
        double time,
        hydro::SpatialScheme scheme,
        const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver,
        const config::SimulationConfig& cfg,
        double cfl,
        double& hydro_inverse_dt,
        double& source_dt_limit
    );
}
