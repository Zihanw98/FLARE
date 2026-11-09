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
        // Dense reference Jacobian of R(U) = U - U_old - dt * RHS(U).
        // active lists every active cell once, in Newton-vector order; R is the
        // baseline Residual evaluated at mesh_trial with these same parameters.
        // Row/column = active-cell position * NVAR + conserved-component index.
        // Cost: N residual evaluations and O(N^2) storage, N = active.size()*NVAR.
        // The unit perturbation scale assumes normalized conserved variables.
        // Exact vacuum has no admissible coordinate-wise derivative and throws.
        [[maybe_unused]] std::vector<std::vector<double>> finite_difference_jacobian(
            const hd_core::MeshBlock& mesh_trial,
            const hd_core::MeshBlock& mesh_old,
            const std::vector<std::size_t>& active,
            const std::vector<hd_core::Conserved>& R,
            double t_new,
            double dt,
            hydro::SpatialScheme scheme,
            const hd_source::SourceList& source_list,
            const hydro::RiemannSolver& solver,
            const config::SimulationConfig& cfg,
            const hd_user::UserHooks& hooks,
            double cfl)
        {
            mesh_trial.validate_active();
            if (R.size() != mesh_trial.cells.size()
                || mesh_old.cells.size() != mesh_trial.cells.size())
                throw std::invalid_argument("Jacobian mesh/residual size mismatch");

            std::vector<bool> seen(mesh_trial.cells.size(), false);
            for (const auto flat : active)
            {
                if (flat >= seen.size() || !mesh_trial.is_active(flat) || seen[flat])
                    throw std::invalid_argument("Invalid Jacobian active-cell ordering");
                seen[flat] = true;
                for (int k = 0; k < hd_core::NVAR; ++k)
                    if (!std::isfinite(R[flat][k]))
                        throw std::invalid_argument("Non-finite Jacobian baseline residual");
            }
            for (std::size_t flat = 0; flat < seen.size(); ++flat)
                if (mesh_trial.is_active(flat) != seen[flat])
                    throw std::invalid_argument("Incomplete Jacobian active-cell ordering");

            const std::size_t n = active.size() * hd_core::NVAR;
            std::vector<std::vector<double>> J(n, std::vector<double>(n));
            hd_core::MeshBlock perturbed = mesh_trial;
            const double relative_step = std::sqrt(std::numeric_limits<double>::epsilon());

            for (std::size_t column = 0; column < n; ++column)
            {
                const auto flat = active[column / hd_core::NVAR];
                const int k = static_cast<int>(column % hd_core::NVAR);
                const double value = mesh_trial.cells[flat].U[k];
                double step = relative_step * std::max(1.0, std::abs(value));

                // Restore all cells, including ghosts modified by the last probe.
                perturbed.cells = mesh_trial.cells;
                auto& cell = perturbed.cells[flat];
                double h = 0.0;
                // Near the pressure floor, use the opposite direction or reduce
                // the step. Only state-validation errors are retried; errors in
                // boundary/source callbacks and the Riemann solver propagate.
                for (int attempt = 0; attempt < 20 && h == 0.0; ++attempt, step *= 0.5)
                    for (const double sign : {1.0, -1.0})
                    {
                        cell.U[k] = value + sign * step;
                        if (!std::isfinite(cell.U[k]) || cell.U[k] == value)
                            continue;
                        try { cell.U.validate(cell.gamma); }
                        catch (const std::invalid_argument&) { continue; }
                        h = cell.U[k] - value; // Actual representable increment.
                        break;
                    }
                if (h == 0.0)
                    throw std::runtime_error("No admissible finite-difference Jacobian perturbation");

                // Residual refills boundaries at t_new, differentiating their
                // dependence on active states (including periodic wraparound).
                const auto R_perturbed = Residual(perturbed, mesh_old, t_new, dt,
                    scheme, source_list, solver, cfg, hooks, cfl);
                for (std::size_t row = 0; row < n; ++row)
                {
                    const auto row_flat = active[row / hd_core::NVAR];
                    const int component = static_cast<int>(row % hd_core::NVAR);
                    J[row][column] = (R_perturbed[row_flat][component]
                                     - R[row_flat][component]) / h;
                    if (!std::isfinite(J[row][column]))
                        throw std::runtime_error("Non-finite finite-difference Jacobian entry");
                }
            }
            return J;
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
