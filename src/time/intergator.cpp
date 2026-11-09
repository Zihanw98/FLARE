/**
 * @file intergator.cpp
 * @brief Time integration.
 */
#include "integrator.h"
#include "core/mesh.h"
#include "utils.h"
#include "timestep_CFL.h"
#include "limiter.h"
#include "reconstruction.h"
#include "wave_timestep.h"
#include "utils/performance.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hd_time
{
    void add_scaled_rhs(
        hd_core::MeshBlock& mesh,
        const std::vector<hd_core::Conserved>& rhs,
        double dt
    )
    {
        performance::Scope timer(performance::Section::StateUpdate);
        if (!std::isfinite(dt))
            throw std::invalid_argument("RHS scale must be finite");
        if (rhs.size() != mesh.cells.size())
            throw std::invalid_argument("RHS size must match mesh cells");
        (void)mesh.flat_index(hd_core::Index{}, true);
        if (dt == 0.0)
            return;

        auto next_cells = performance::measure(performance::Section::MeshCopy, [&] { return mesh.cells; });
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
    }

    std::pair<double, std::vector<hd_core::Conserved>> evaluate_RHS(
        const hd_core::MeshBlock& mesh,
        double time,
        hydro::SpatialScheme scheme,
        const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver,
        const config::SimulationConfig& cfg,
        double cfl,
        double& hydro_inverse_dt,
        double& source_dt_limit,
        bool use_source_limits,
        std::array<std::vector<double>, hd_core::NDIM>* pressure_hints
    )
    {
        performance::Scope timer(performance::Section::RHS);
        if (!std::isfinite(cfl) || cfl <= 0.0)
            throw std::invalid_argument("RHS requires a finite positive CFL target");
        std::vector<hd_core::Conserved> rhs;
        hydro::Limiter limiter = hydro::parse_limiter(cfg.limiter);
        hydro::Reconstructor reconstructor = hydro::parse_reconstructor(cfg.reconstructor);
        hydro::compute_hydro_rhs(mesh, scheme, solver, rhs,reconstructor,limiter,cfg.limiter_beta,pressure_hints);
        const double inverse_dt = hydro::wave_timestep.max_inverse_dt;
        if (!std::isfinite(inverse_dt) || inverse_dt < 0.0)
            throw std::invalid_argument("RHS requires a finite nonnegative CFL rate");
        const double hydro_limit = inverse_dt == 0.0
            ? std::numeric_limits<double>::infinity() : cfl / inverse_dt;
        double source_limit_all = std::numeric_limits<double>::infinity();
        performance::Scope source_timer(performance::Section::Sources);
        for (const auto& source : source_list)
        {
            if (!source.func)
                throw std::invalid_argument("Missing source function");
            double source_limit = std::numeric_limits<double>::infinity();
            source.func(mesh, time, cfg, rhs, source_limit);
            if (!use_source_limits)
                continue;
            if (std::isnan(source_limit) || source_limit <= 0.0)
                throw std::invalid_argument("Source timestep must be positive (infinity allowed)");
            source_limit_all = std::min(source_limit_all, source_limit);
        }
        source_timer.stop();
        if (rhs.size() != mesh.cells.size())
            throw std::logic_error("Source changed RHS storage size");
        hydro_inverse_dt = inverse_dt;
        source_dt_limit = source_limit_all;
        return {std::min(hydro_limit, source_limit_all), std::move(rhs)};
    }
    TimeIntegrator parse_time_integrator(const std::string& name)
    {
        if (utils::iequals(name,"euler"))
            return TimeIntegrator::Euler;
        if (utils::iequals(name,"lts_euler"))
            return TimeIntegrator::LTSEuler;
        if(utils::iequals(name,"rk2"))
            return TimeIntegrator::RK2;
        if (utils::iequals(name, "backward_euler"))
            return TimeIntegrator::BackwardEuler;
        throw std::invalid_argument("Unknown or unavailable time integrator: " + name);
    }

    double advance_one_step(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        TimeIntegrator integrator,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl,
        double& dt_next,
        LTSEulerWorkspace* lts_workspace,
        std::array<std::vector<double>, hd_core::NDIM>* pressure_history
    )
    {
        performance::Scope timer(performance::Section::Advance);
        double dt;
        switch (integrator)
        {
            case TimeIntegrator::Euler:
                dt =  advance_euler(mesh, time, dt_limit, spatial_scheme, solver, source_list, cfg, hooks,cfl,pressure_history);
                dt_next = std::numeric_limits<double>::infinity(); // Euler 每步重新估计，无跨步建议。
                break;
            case TimeIntegrator::RK2:
                dt = advance_rk2(mesh,time,dt_limit,spatial_scheme,solver,source_list,cfg,hooks,cfl,dt_next,pressure_history);
                break;
            case TimeIntegrator::LTSEuler:
                dt = advance_lts_euler(mesh, time, dt_limit, spatial_scheme, solver,
                                       source_list, cfg, hooks, cfl, lts_workspace);
                dt_next = std::numeric_limits<double>::infinity();
                return dt; // LTS commits only after the synchronized boundary fill.
            case TimeIntegrator::BackwardEuler:
                dt = advance_backward_euler(mesh, time, dt_limit, spatial_scheme, solver,
                                            source_list, cfg, hooks, cfl);
                dt_next = std::numeric_limits<double>::infinity();
                return dt; // The accepted residual already filled the new-time ghosts.
            default:
                throw std::invalid_argument("Invalid time integrator");
        }
        boundary::FillBoundary(mesh,time + dt,hooks);
        return dt;
    }
}
