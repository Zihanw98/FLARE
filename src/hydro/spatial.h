#pragma once
#include "core.h"
#include <string>
#include <vector>

namespace hd_core
{
    struct MeshBlock;
}

namespace hydro
{
class RiemannSolver;    //Link to spatial.cpp->godunov.h->godunov.cpp->riemman.h->RiemannSolver

enum class SpatialScheme
{
    Godunov
};
SpatialScheme parse_spatial_scheme(const std::string& name);
void compute_hydro_rhs(
    const hd_core::MeshBlock& mesh,
    SpatialScheme scheme,
    const RiemannSolver& solver,
    std::vector<hd_core::Conserved>& rhs
);

}
