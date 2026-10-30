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
    void compute_godunov_rhs(
        const hd_core::MeshBlock& mesh,
        const RiemannSolver& solver,
        std::vector<hd_core::Conserved>& rhs
    );
}