#include "integrator.h"
#include <stdexcept>

namespace hd_time
{
    TimeIntegrator parse_time_integrator(const std::string& name)
    {
        if (name == "euler")
            return TimeIntegrator::Euler;
        throw std::invalid_argument("Unknown or unavailable time integrator: " + name);
    }

    double advance_one_step(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        TimeIntegrator integrator,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        double cfl
    )
    {
        switch (integrator)
        {
            case TimeIntegrator::Euler:
                return advance_euler(mesh, time, dt_limit, spatial_scheme, solver, cfl);
        }
        throw std::invalid_argument("Invalid time integrator");
    }
}
