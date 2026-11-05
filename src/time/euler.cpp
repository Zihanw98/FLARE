/**
 * @file euler.cpp
 * @brief Time integration.
 */
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
        double cfl
    )
    {
        if (!std::isfinite(time) || !std::isfinite(dt_limit) || dt_limit <= 0.0)
            throw std::invalid_argument("Euler requires finite time and a finite positive dt_limit");

        if (!std::isfinite(cfl) || cfl <= 0.0 || cfl > 0.5)
            throw std::invalid_argument("Euler requires 0 < CFL <= 0.5");
        double hydro_inverse_dt;
        double source_dt_limit;
        auto [dt_rhs, rhs] = evaluate_RHS(
            mesh,
            time,
            spatial_scheme,
            source_list,
            solver,
            cfg,
            cfl,
            hydro_inverse_dt,
            source_dt_limit
        );

        const double dt = std::min(dt_limit,dt_rhs);
        if (!std::isfinite(dt) || dt <= 0.0 ||
            !std::isfinite(time + dt) || time + dt <= time)
            throw std::runtime_error("Euler time step cannot advance time");

        add_scaled_rhs(mesh, rhs, dt);
        return dt;
    }
}
