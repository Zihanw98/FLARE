#include "spatial.h"
#include "godunov.h"
#include "mesh.h"

#include <stdexcept>

namespace hydro
{

SpatialScheme parse_spatial_scheme(const std::string& name)
{
    if (name == "godunov")
    {
        return SpatialScheme::Godunov;
    }

    throw std::invalid_argument(
        "Unknown or unavailable spatial scheme: " + name
    );
}

void compute_hydro_rhs(
    const hd_core::MeshBlock& mesh,
    SpatialScheme scheme,
    const RiemannSolver& solver,
    std::vector<hd_core::Conserved>& rhs
)
{
    //Clear the results and re-calculate
    rhs.assign(
        mesh.cells.size(),
        hd_core::Conserved{}
    );

    switch (scheme)
    {
        case SpatialScheme::Godunov:
            compute_godunov_rhs(mesh, solver, rhs);
            return;
    }

    throw std::invalid_argument("Invalid spatial scheme");
}

}