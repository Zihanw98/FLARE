#include "integrator.h"
#include "utils.h"
#include "timestep_CFL.h"
#include <stdexcept>

namespace hd_time
{
    std::pair<double,std::vector<hd_core::Conserved>> evaluate_RHS(
        const hd_core::MeshBlock& mesh,
        double time,
        hydro::SpatialScheme scheme,
        const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver,
        const config::SimulationConfig& cfg,
        double cfl=0.4
    )
    {
        std::vector<hd_core::Conserved> rhs;
        hydro::compute_hydro_rhs(mesh, scheme, solver,rhs); 
        double dt_rhs;
        dt_rhs = hydro::compute_cfl_dt(cfg.cfl);
        for(auto i=0;i<source_list.size();i++)
        {
            const auto& source = source_list[i];
            source.func(mesh,time,cfg,rhs,dt_rhs);
        }
        return {dt_rhs,rhs}; 
    }
    TimeIntegrator parse_time_integrator(const std::string& name)
    {
        if (utils::iequals(name,"euler"))
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
        hd_source::SourceList source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl = 0.4
    )
    {
        double dt;
        switch (integrator)
        {
            case TimeIntegrator::Euler:
                dt =  advance_euler(mesh, time, dt_limit, spatial_scheme, solver, source_list, cfg, hooks);
                break;
            default:
                throw std::invalid_argument("Invalid time integrator");
        }
        boundary::FillBoundary(mesh,time + dt,hooks);
        return dt;
    }
}
