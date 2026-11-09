#include "integrator.h"
#include "core/mesh.h"
#include "core/core.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <optional>
#include <utility>
namespace hd_time
{
    namespace
    {
        constexpr int max_newton_iteration = 15;

        std::vector<hd_core::Conserved> Residual(hd_core::MeshBlock& mesh_trial, const hd_core::MeshBlock& mesh_old,double t_new, double dt,hydro::SpatialScheme scheme, const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver, const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks, double cfl)
        {
            std::vector<hd_core::Conserved> residual;
            residual.resize(mesh_trial.cells.size());
            boundary::FillBoundary(mesh_trial,t_new,hooks);
            double hydro_inverse_dt, source_limit;
            auto [_,rhs] = evaluate_RHS(mesh_trial,t_new,scheme,source_list,solver,cfg,cfl,
                                        hydro_inverse_dt,source_limit, /*use_source_limits=*/false);
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
        std::vector<std::vector<double>> finite_difference_jacobian(
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
                double step = relative_step * std::max(1.0, std::abs(value));       //At least it's sqrt(episilon), at order of 1e-8
                //Need further modification when working with composition...

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
        // Dense Gaussian elimination with partial pivoting. Row swaps preserve unknown order.
        std::optional<std::vector<double>> linear_solve(
            std::vector<std::vector<double>> A, std::vector<double> b)
        {
            const std::size_t n = b.size();
            if (A.size() != n)
                throw std::invalid_argument("linear_solve: matrix and RHS sizes differ");
            for (const auto& row : A)
                if (row.size() != n)
                    throw std::invalid_argument("linear_solve: expected a square matrix");

            for (std::size_t k = 0; k < n; ++k)
            {
                std::size_t pivot = k;
                for (std::size_t i = k + 1; i < n; ++i)
                    if (std::abs(A[i][k]) > std::abs(A[pivot][k]))
                        pivot = i;

                if (A[pivot][k] == 0.0 || !std::isfinite(A[pivot][k]))
                    return std::nullopt;
                if (pivot != k)
                {
                    std::swap(A[k], A[pivot]);
                    std::swap(b[k], b[pivot]);
                }

                for (std::size_t i = k + 1; i < n; ++i)
                {
                    const double factor = A[i][k] / A[k][k];
                    A[i][k] = 0.0;
                    for (std::size_t j = k + 1; j < n; ++j)
                        A[i][j] -= factor * A[k][j];
                    b[i] -= factor * b[k];
                }
            }

            std::vector<double> dX(n);
            for (std::size_t end = n; end > 0; --end)
            {
                const std::size_t i = end - 1;
                double rhs = b[i];
                for (std::size_t j = i + 1; j < n; ++j)
                    rhs -= A[i][j] * dX[j];
                dX[i] = rhs / A[i][i];
                if (!std::isfinite(dX[i]))
                    return std::nullopt;
            }
            return dX;
        }

        std::optional<hd_core::MeshBlock> Solve_Newton(const hd_core::MeshBlock& mesh_old, double t_new, double dt,hydro::SpatialScheme scheme, const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver, const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks, double cfl)
        {
            constexpr double rtol = 1e-6;
            std::array<double, hd_core::NVAR> atol;
            atol.fill(1e-10);
            hd_core::MeshBlock mesh_trial = mesh_old;
            const auto active = mesh_old.get_all_active_flats();

            auto R = Residual(mesh_trial, mesh_old, t_new, dt, scheme, source_list,
                              solver, cfg, hooks, cfl);
            for (int iter = 0; ; ++iter)
            {
                const double norm = scaled_norm(R, mesh_old, mesh_trial, atol, rtol);
                if (norm <= 1.0)
                    return mesh_trial;
                if (!std::isfinite(norm) || iter == max_newton_iteration)
                    return std::nullopt;

                auto J = finite_difference_jacobian(
                    mesh_trial, mesh_old, active, R, t_new, dt,
                    scheme, source_list, solver, cfg, hooks, cfl);
                std::vector<double> b;
                b.reserve(J.size());
                for (const auto flat : active)
                    for (int k = 0; k < hd_core::NVAR; ++k)
                        b.push_back(-R[flat][k]);
                const auto result = linear_solve(std::move(J), std::move(b));
                if (!result)
                    return std::nullopt;
                const auto& dX = *result;

                bool accepted = false;
                for (double alpha = 1.0; alpha >= 1.0 / 4096; alpha *= 0.5)
                {
                    auto candidate = mesh_trial;
                    std::size_t j = 0;
                    for (const auto flat : active)
                        for (int k = 0; k < hd_core::NVAR; ++k)
                            candidate.cells[flat].U[k] += alpha * dX[j++];
                    try {
                        candidate.validate_active();
                    } catch (const std::invalid_argument&) {
                        continue;
                    }

                    auto next_R = Residual(candidate, mesh_old, t_new, dt, scheme,
                        source_list, solver, cfg, hooks, cfl);
                    // Compare both residuals using the current iterate's scales.
                    if (scaled_norm(next_R, mesh_old, mesh_trial, atol, rtol) < norm)
                    {
                        mesh_trial.cells.swap(candidate.cells);
                        R.swap(next_R);
                        accepted = true;
                        break;
                    }
                }
                if (!accepted)
                    return std::nullopt;
            }
        }
    }
    double advance_backward_euler(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl
    )
    {
        if (!std::isfinite(time) || !std::isfinite(dt_limit) || dt_limit <= 0.0
            || cfg.integrator_retry_max < 0)
            throw std::invalid_argument("Invalid backward Euler timestep or retry settings");
        mesh.validate_active();
        double dt = dt_limit;
        for (long long attempt = 0; attempt <= cfg.integrator_retry_max; ++attempt)
        {
            if (!std::isfinite(time + dt) || time + dt <= time)
                throw std::runtime_error("Backward Euler timestep cannot advance time");
            auto solution = Solve_Newton(mesh, time + dt, dt, spatial_scheme, source_list,
                                         solver, cfg, hooks, cfl);
            if (!solution)
            {
                dt *= 0.5;
                continue;
            }
            mesh.cells.swap(solution->cells); // Includes the converged new-time ghosts.
            return dt;
        }
        throw std::runtime_error("Backward Euler retry limit exceeded");
    }
}
