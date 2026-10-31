#pragma once
#include "core.h"
#include <vector>

namespace hd_core
{
    struct MeshBlock;
}
namespace hydro
{
    class RiemannSolver;
    // ghost 状态必须预先填好，rhs 必须已按 mesh.cells.size() 分配。
    // 重置并更新 wave_timestep 中的本次 CFL 限制。
    void compute_godunov_rhs(
        const hd_core::MeshBlock& mesh,
        const RiemannSolver& solver,
        std::vector<hd_core::Conserved>& rhs
    );
}
