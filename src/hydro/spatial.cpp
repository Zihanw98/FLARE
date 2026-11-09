/**
 * @file spatial.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "spatial.h"
#include "godunov.h"
#include "mesh.h"
#include "utils.h"
#include "utils/performance.h"
#include <stdexcept>

namespace hydro
{

SpatialScheme parse_spatial_scheme(const std::string& name)
{
    if (utils::iequals(name,"godunov"))
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
    std::vector<hd_core::Conserved>& rhs,
    const Reconstructor& reconstructor,
    const Limiter& limiter,
    const double limiter_beta,
    hydro::StarPressureCache* star_pressure_cache
)
{
    performance::Scope timer(performance::Section::Hydro);
    //Clear the results and re-calculate
    rhs.assign(
        mesh.cells.size(),
        hd_core::Conserved{}
    );

    switch (scheme)
    {
        case SpatialScheme::Godunov:
            compute_godunov_rhs(mesh, solver, rhs,reconstructor,limiter,limiter_beta,star_pressure_cache);
            return;
    }

    throw std::invalid_argument("Invalid spatial scheme");
}

}
