#include "core.h"
#include "mesh.h"
#include "integrator.h"
namespace time
{
    void advance_one_step(
        hd_core::MeshBlock& mesh,
        double time,
        double dt,
        TimeIntegrator integrator,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver
    )
    {
        if(integrator == TimeIntegrator::Euler)
        {
            advance_euler(mesh,time,dt,spatial_scheme,solver);
        }
    }
}