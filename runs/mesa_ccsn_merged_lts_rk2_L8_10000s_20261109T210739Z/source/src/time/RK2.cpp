/**
 * @file RK2.cpp
 * @brief Time integration.
 */
#include "integrator.h"
#include "core/mesh.h"
#include "utils/performance.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
namespace hd_time
{
    namespace
    {
        /**
         * @brief Convert a nonnegative inverse timescale to a CFL timestep.
         * @param[in] cfl Finite positive CFL coefficient.
         * @param[in] inverse_dt Finite nonnegative inverse timescale.
         * @return cfl / inverse_dt, or positive infinity for a zero rate.
         * @pre Inputs have already been validated by the caller.
         */
        double cfl_step(double cfl, double inverse_dt)
        {
            return inverse_dt == 0.0
                ? std::numeric_limits<double>::infinity() : cfl / inverse_dt;
        }

        /**
         * @brief Combine the original state and Euler stage into an SSP-RK2 candidate.
         * @param[in] mesh Original mesh at the beginning of the timestep.
         * @param[in] stage First-stage mesh with matching dimensions and storage layout.
         * @param[in] rhs Second-stage RHS with one entry per stored cell.
         * @param[in] dt Timestep shared by both stages.
         * @return Validated mesh copy with active states equal to 0.5 * U0 + 0.5 * (U1 + dt * rhs).
         * @pre dt is finite and both stages refer to the same physical grid.
         * @throws std::logic_error If dimensions or storage layouts differ.
         * @throws std::invalid_argument If a candidate active state fails validation.
         * @note Retains the original mesh geometry, ghost states, and cell adiabatic indices; inputs are unchanged.
         */
        hd_core::MeshBlock combine_rk2(
            const hd_core::MeshBlock& mesh,
            const hd_core::MeshBlock& stage,
            const std::vector<hd_core::Conserved>& rhs,
            double dt
        )
        {
            performance::Scope timer(performance::Section::StateUpdate);
            if (stage.cells.size() != mesh.cells.size()
                || stage.nactive != mesh.nactive || stage.nghost != mesh.nghost
                || rhs.size() != mesh.cells.size())
                throw std::logic_error("RK2 combination requires matching mesh and RHS layouts");

            auto candidate = performance::measure(performance::Section::MeshCopy, [&] { return mesh; });
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

                candidate.cells[flat].U =
                    0.5 * mesh.cells[flat].U
                    + 0.5 * (stage.cells[flat].U + dt * rhs[flat]);
            }
            candidate.validate_active();
            return candidate;
        }
    }

    double advance_rk2(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl,
        double& dt_next,
        hydro::StarPressureCache* star_pressure_cache
    )
    {
        if (!std::isfinite(time) || !std::isfinite(dt_limit) || dt_limit <= 0.0)
            throw std::invalid_argument("RK2 requires finite time and a finite positive dt_limit");

        if (!std::isfinite(cfl) || !std::isfinite(cfg.cfl_max)
            || cfl <= 0.0 || cfl > cfg.cfl_max || cfg.cfl_max > 0.5)
            throw std::invalid_argument("RK2 requires 0 < cfl_desired <= cfl_max <= 0.5");
        if (cfg.integrator_retry_max < 0)
            throw std::invalid_argument("RK2 retry count must be nonnegative");
        if (std::isnan(dt_next) || dt_next <= 0.0)
            throw std::invalid_argument("RK2 next timestep must be positive (infinity allowed)");
        mesh.validate_active();
        hydro::StarPressureCache initial_pressure, trial_pressure;
        if (star_pressure_cache) initial_pressure = *star_pressure_cache;
        double rate0, source_limit0;
        auto [dt_rhs_0, rhs_0] = evaluate_RHS(mesh,time,spatial_scheme,source_list,solver,cfg,
                                             cfl,rate0,source_limit0,
                                             star_pressure_cache ? &initial_pressure : nullptr);
        // 使用目标 CFL 选步长，同时遵守当前源项上限、上一步建议和外部时间上限。
        double dt = std::min({dt_limit, dt_next, dt_rhs_0});

        for(long long attempt = 0; attempt <= cfg.integrator_retry_max; attempt++)
        {
            performance::count(performance::Counter::RK2Attempts);
            if (!std::isfinite(dt) || dt <= 0.0
                || !std::isfinite(time + dt) || time + dt <= time)
                throw std::runtime_error("RK2 timestep cannot advance time");
            auto mesh_trial = performance::measure(performance::Section::MeshCopy, [&] { return mesh; });
            try {
                add_scaled_rhs(mesh_trial, rhs_0, dt);
            }
            catch (const std::invalid_argument&) {
                performance::count(performance::Counter::RK2Rejected);
                dt *= 0.5;
                continue;
            }
            boundary::FillBoundary(mesh_trial,time+dt,hooks);

            double rate1, source_limit1;
            // Every retry starts from stage 0, never the rejected stage's roots.
            if (star_pressure_cache) trial_pressure = initial_pressure;
            auto [dt_rhs_1, rhs_1] = evaluate_RHS(mesh_trial,time+dt,spatial_scheme,source_list,solver,cfg,
                                                 cfl,rate1,source_limit1,
                                                 star_pressure_cache ? &trial_pressure : nullptr);
            const double inverse_dt = std::max(rate0, rate1);
            const double actual_cfl = dt * inverse_dt;
            const double source_limit = std::min(source_limit0, source_limit1);
            const bool cfl_exceeded = actual_cfl > cfg.cfl_max;
            const bool source_exceeded = dt > source_limit;

            if (cfl_exceeded || source_exceeded)
            {
                performance::count(performance::Counter::RK2Rejected);
                // CFL 部分等价于 dt * cfl_desired / actual_cfl，直接除速率避免溢出。
                // 源项是独立硬限制；只有源项超限时为重试预留余量。
                const double source_retry_limit = source_exceeded ? 0.9 * source_limit : source_limit;
                const double retry_dt = std::min({dt_limit, cfl_step(cfl, inverse_dt), source_retry_limit});
                // 即使除法舍入到旧 dt，拒绝后的步长也必须严格变小。
                dt = std::min(retry_dt, std::nextafter(dt, 0.0));
                continue;
            }

            try {
                auto candidate = combine_rk2(mesh, mesh_trial, rhs_1, dt);
                mesh.cells.swap(candidate.cells);
                if (star_pressure_cache) star_pressure_cache->swap(trial_pressure);
                // 建议值不受本步 dt_limit 截断，避免碰到输出/终止时刻后一直使用小步长。
                dt_next = std::min(cfl_step(cfl, inverse_dt), source_limit);
                return dt;
            }
            catch (const std::invalid_argument&) {
                performance::count(performance::Counter::RK2Rejected);
                dt *= 0.5;
                continue;
            }
        }
        throw std::runtime_error("RK2 retry limit exceeded");
    }
}
