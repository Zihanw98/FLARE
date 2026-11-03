#include "integrator.h"
#include "core/mesh.h"
#include "hydro/timestep_CFL.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace hd_time
{
    double advance_euler(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        hd_source::SourceList source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl = 0.4
    )
    {
        if (!std::isfinite(time) || !std::isfinite(dt_limit) || dt_limit <= 0.0)
            throw std::invalid_argument("Euler requires finite time and a finite positive dt_limit");

        if (!std::isfinite(cfl) || cfl <= 0.0 || cfl > 0.5)
            throw std::invalid_argument("Euler requires 0 < CFL <= 0.5");

        std::vector<hd_core::Conserved> rhs;
        hydro::compute_hydro_rhs(mesh, spatial_scheme, solver, rhs);
        double dt_h = std::min(dt_limit, hydro::compute_cfl_dt(cfl));
        double dt_s = dt_h;
        for(int i=0;i < source_list.size(); i++)
        {
            auto source_function = source_list.at(i).func;
            source_function(mesh,time,cfg,rhs,dt_s);
        }
        const double dt = std::min(dt_h,dt_s);
        if (!std::isfinite(dt) || dt <= 0.0 ||
            !std::isfinite(time + dt) || time + dt <= time)
            throw std::runtime_error("Euler time step cannot advance time");

        
        auto next_cells = mesh.cells;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            const auto index = mesh.coordinates_active(flat);
            bool is_active = true;
            for (int d = 0; d < hd_core::NDIM; ++d)
            {
                if (index[d] < 0 || index[d] >= mesh.nactive[d])
                {
                    is_active = false;
                    break;
                }
            }
            if (!is_active)
                continue;

            next_cells[flat].U += dt * rhs[flat];
            next_cells[flat].U.validate(next_cells[flat].gamma);
        }
        mesh.cells.swap(next_cells);
        return dt;
    }
}
