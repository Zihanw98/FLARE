#include "integrator.h"
namespace hd_time
{
    double advance_rk2(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        hd_source::SourceList source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl=0.4
    )
    {
        if (!std::isfinite(time) || !std::isfinite(dt_limit) || dt_limit <= 0.0)
            throw std::invalid_argument("Euler requires finite time and a finite positive dt_limit");

        if (!std::isfinite(cfl) || cfl <= 0.0 || cfl > 0.5)
            throw std::invalid_argument("Euler requires 0 < CFL <= 0.5");

        std::vector<hd_core::Conserved> rhs;

        
    }
}