#include "integrator.h"
#include "core/mesh.h"
#include "core/core.h"
#include "hydro/timestep_CFL.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <optional>
namespace hd_time
{
    namespace
    {
        int max_newton_iteration = 15;
        std::vector<hd_core::Conserved> Residual(hd_core::MeshBlock& mesh_trial, const hd_core::MeshBlock& mesh_old,double t_new, double dt,hydro::SpatialScheme scheme, const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver, const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks, double cfl);

        std::vector<hd_core::Conserved> Residual(hd_core::MeshBlock& mesh_trial, const hd_core::MeshBlock& mesh_old,double t_new, double dt,hydro::SpatialScheme scheme, const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver, const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks, double cfl)
        {
            std::vector<hd_core::Conserved> residual;
            residual.resize(mesh_trial.cells.size());
            boundary::FillBoundary(mesh_trial,t_new,hooks);
            double dt1=0,dt2=0;
            auto [_,rhs]  = evaluate_RHS(mesh_trial,t_new,scheme,source_list,solver,cfg,cfl,dt1,dt2);
            for(size_t flat = 0; flat < mesh_trial.cells.size(); flat++)
            {
                if(!mesh_trial.is_active(flat))
                {
                    continue;
                }
                auto R = mesh_trial.cells[flat].U - mesh_old.cells[flat].U - dt * rhs[flat];
                residual[flat] = R;
            }
            return residual;
        }
        // active contains the flat indices of the active cells, in Newton-vector order.
        // R is the unperturbed residual computed with the same step parameters.
        std::vector<std::vector<double>> finite_difference_jacobian()
        {
            
        }

        double scaled_norm(
            const std::vector<hd_core::Conserved>& R,
            const hd_core::MeshBlock& mesh_old,
            const hd_core::MeshBlock& mesh_trial,
            const std::array<double, hd_core::NVAR>& atol,
            double rtol)
        {
            double worst = 0.0;

            for (std::size_t flat = 0; flat < R.size(); ++flat)
            {
                if (!mesh_trial.is_active(flat))
                    continue;

                for (int k = 0; k < hd_core::NVAR; ++k)
                {
                    const double old_value = mesh_old.cells[flat].U[k];
                    const double trial_value = mesh_trial.cells[flat].U[k];
                    const double scale = std::max(std::abs(old_value), std::abs(trial_value));
                    const double allowed = atol[k] + rtol * scale;

                    if (!std::isfinite(allowed) || allowed <= 0.0)
                        throw std::invalid_argument("Invalid Newton residual tolerance");

                    const double ratio = std::abs(R[flat][k]) / allowed;
                    if (!std::isfinite(ratio))
                        return std::numeric_limits<double>::infinity();

                    worst = std::max(worst, ratio);
                }
            }
            return worst;
        }
        std::optional<hd_core::MeshBlock> Solve_Newton(const hd_core::MeshBlock& mesh_old, double t_new, double dt,hydro::SpatialScheme scheme, const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver, const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks, double cfl)
        {
            constexpr double rtol = 1e-6;
            std::array<double, hd_core::NVAR> atol;
            atol.fill(1e-10);
            hd_core::MeshBlock mesh_trial = mesh_old;
            for(int k=0; k<max_newton_iteration; k++)
            {
                auto R = Residual(mesh_trial,mesh_old,t_new,dt,scheme,source_list,solver,cfg,hooks,cfl);
                if(scaled_norm(R,mesh_old,mesh_trial,atol,rtol) <= 1.0)
                {
                    //If we have some sources that lasts for finite time, then here we will check the time...(hard upper limit)---> Add a t_stop, and cut off before call this function is the best way
                    return mesh_trial;
                }
            }
            return std::nullopt;
        }
    }
    double advance_backward_euler(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        TimeIntegrator integrator,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl
    )
    {
        throw std::logic_error("Backward Euler time stepping is not implemented");
    }
}
