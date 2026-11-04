#pragma once
#include "core.h"
#include "reconstruction.h"
#include "limiter.h"
#include <string>
#include <vector>

namespace hd_core
{
    struct MeshBlock;
}

namespace hydro
{
class RiemannSolver;

enum class SpatialScheme
{
    Godunov
};
SpatialScheme parse_spatial_scheme(const std::string& name);
// 填充 rhs，同时通过 wave_timestep 更新本次波速限制。
void compute_hydro_rhs(
    const hd_core::MeshBlock& mesh,
    SpatialScheme scheme,
    const RiemannSolver& solver,
    std::vector<hd_core::Conserved>& rhs,
    const Reconstructor& reconstructor = Reconstructor::PCM,
    const Limiter& limiter = Limiter::none
);

}
