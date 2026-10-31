#include "integrator.h"
#include "core.h"
#include "mesh.h"
#include "facefield.h"
#include "riemann.h"
#include "exact.h"
#include "spatial.h"
#include <vector>
namespace time
{
    void advance_euler(
        hd_core::MeshBlock& mesh,
        double time,
        double dt,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver
    )
    {
        std::vector<hd_core::Conserved> rhs;
        hydro::compute_hydro_rhs(mesh,spatial_scheme,solver,rhs);
    }
}