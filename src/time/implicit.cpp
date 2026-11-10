#include "implicit.h"
#include "implicit_linear.h"
#include "implicit_variables.h"
#include "user/user_api.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace hd_implicit
{
using detail::Vector;
using detail::LogVariables;
namespace
{
constexpr int nv = hd_core::NVAR;
constexpr double difference_step = 2e-6;
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}
double merit(const Vector& r, const LogVariables& coordinates)
{
    double sum = 0;
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        const double value = r[i] / coordinates.row_scale(i);
        sum += value * value;
    }
    return std::sqrt(sum);
}
}

struct Integrator::Impl
{
    const config::SimulationConfig settings;
    hydro::Reconstructor reconstruction;
    hydro::Limiter limiter;
    hd_core::MeshBlock trial, probe, minus, candidate;
    std::vector<std::size_t> active, color;
    std::vector<std::vector<std::size_t>> columns, groups;
    detail::SparseMatrix preconditioner;
    std::ofstream report;
    double next_dt, residual_norm = 0;
    int newton = 0, linear = 0, rhs_calls = 0, preconditioner_calls = 0;
    double matrix_seconds = 0, linear_seconds = 0;
    std::string reason;

    Impl(const hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
        : settings(cfg), reconstruction(hydro::parse_reconstructor(cfg.reconstructor)),
          limiter(hydro::parse_limiter(cfg.limiter)), trial(mesh), probe(mesh), minus(mesh), candidate(mesh),
          next_dt(settings.dt_implicit)
    {
        if (!std::isfinite(settings.dt_implicit) || settings.dt_implicit <= 0 || !std::isfinite(settings.implicit_tolerance)
            || settings.implicit_tolerance <= 0 || settings.implicit_tolerance >= 1 || settings.implicit_newton_max_iterations < 1)
            throw std::invalid_argument("Implicit requires dt_implicit > 0, 0 < implicit_tolerance < 1, and implicit_newton_max_iterations > 0");
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
            if (mesh.is_active(flat)) active.push_back(flat);
        if (active.empty()) throw std::invalid_argument("Implicit requires active cells");
        build_pattern(mesh);
        if (!settings.implicit_report_file.empty())
        {
            report.open(settings.implicit_report_file);
            if (!report) throw std::runtime_error("Cannot open implicit report: " + settings.implicit_report_file);
            report << "time,dt,accepted,newton,gmres,rhs_calls,preconditioner_calls,residual,wall_s,matrix_s,linear_s,reason\n";
        }
        const auto& a = preconditioner;
        const auto bytes = (a.value.size()*2 + (30+6)*a.size())*sizeof(double)
                         + (a.column.size()+a.row.size()+a.diagonal.size())*sizeof(std::size_t);
        std::cerr << "[implicit] unknowns=" << a.size() << " colors=" << groups.size()
                  << " preconditioner_nnz=" << a.value.size()
                  << " matrix_and_krylov_bytes_estimate=" << bytes << '\n';
    }

    void build_pattern(const hd_core::MeshBlock& mesh)
    {
        // Only the preconditioner freezes ghosts and omits physical sources.
        // Thus its graph is local regardless of the real boundary/source graph.
        const int radius = reconstruction == hydro::Reconstructor::PCM ? 1 : 2;
        const auto missing = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> mapping(mesh.cells.size(), missing);
        for (std::size_t i = 0; i < active.size(); ++i) mapping[active[i]] = i;
        columns.resize(active.size());
        std::vector<std::vector<std::size_t>> rows(active.size());
        for (std::size_t i = 0; i < active.size(); ++i)
        {
            auto& neighbors = columns[i];
            neighbors.push_back(i);
            const auto center = mesh.coordinates_active(active[i]);
            for (int d = 0; d < hd_core::NDIM; ++d)
                for (int offset = -radius; offset <= radius; ++offset)
                {
                    auto index = center;
                    index[d] += offset;
                    if (offset && index[d] >= 0 && index[d] < mesh.nactive[d])
                        neighbors.push_back(mapping[mesh.flat_index(index, true)]);
                }
            std::sort(neighbors.begin(), neighbors.end());
            for (auto j : neighbors) rows[j].push_back(i);
        }
        color.assign(active.size(), missing);
        std::vector<std::size_t> marked(active.size(), missing);
        for (std::size_t j = 0; j < active.size(); ++j)
        {
            for (auto row : rows[j])
                for (auto k : columns[row])
                    if (color[k] != missing) marked[color[k]] = j;
            std::size_t c = 0;
            while (marked[c] == j) ++c;
            color[j] = c;
            if (c == groups.size()) groups.emplace_back();
            groups[c].push_back(j);
        }
        auto& a = preconditioner;
        a.row.push_back(0);
        for (std::size_t i = 0; i < active.size(); ++i)
            for (int k = 0; k < nv; ++k)
            {
                for (auto j : columns[i])
                    for (int c = 0; c < nv; ++c)
                    {
                        if (i == j && k == c) a.diagonal.push_back(a.column.size());
                        a.column.push_back(j*nv+c);
                    }
                a.row.push_back(a.column.size());
            }
        a.value.resize(a.column.size());
    }

    void local_rhs(const hd_core::MeshBlock& state, const hydro::RiemannSolver& solver,
                   hydro::SpatialScheme scheme, const config::SimulationConfig& cfg,
                   std::vector<hd_core::Conserved>& rhs)
    {
        hydro::compute_hydro_rhs(state, scheme, solver, rhs, reconstruction, limiter, cfg.limiter_beta);
        double limit = std::numeric_limits<double>::infinity();
        if (state.geometry == hd_core::Geometry::Spherical)
            hd_source::geometric_source(state, 0.0, cfg, rhs, limit);
        ++preconditioner_calls;
    }

    void build_preconditioner(const LogVariables& coordinates, double dt,
                              const hydro::RiemannSolver& solver, hydro::SpatialScheme scheme,
                              const config::SimulationConfig& cfg)
    {
        const auto start = Clock::now();
        std::vector<hd_core::Conserved> baseline, perturbed;
        local_rhs(trial, solver, scheme, cfg, baseline);
        auto& a = preconditioner;
        for (std::size_t c = 0; c < groups.size(); ++c)
            for (int component = 0; component < nv; ++component)
            {
                hd_core::Conserved dq{};
                dq[component] = 1;
                probe.cells = trial.cells;
                for (auto j : groups[c])
                    probe.cells[active[j]].U = coordinates.update(j, trial.cells[active[j]], dq.data.data(), difference_step);
                local_rhs(probe, solver, scheme, cfg, perturbed);
                for (std::size_t i = 0; i < active.size(); ++i)
                    for (std::size_t block = 0; block < columns[i].size(); ++block)
                    {
                        const auto j = columns[i][block];
                        if (color[j] != c) continue;
                        hd_core::Conserved du{};
                        if (i == j) du = coordinates.direction(i, trial.cells[active[i]], dq.data.data());
                        for (int k = 0; k < nv; ++k)
                        {
                            const auto row = i*nv+k;
                            a.value[a.row[row]+block*nv+component] = (du[k] - dt*
                                (perturbed[active[i]][k]-baseline[active[i]][k])/difference_step)
                                / coordinates.row_scale(row);
                        }
                    }
            }
        matrix_seconds += seconds(start);
    }

    bool solve(const hd_core::MeshBlock& old, double time, double dt, hydro::SpatialScheme scheme,
               const hydro::RiemannSolver& solver, const hd_source::SourceList& sources,
               const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks)
    {
        trial.cells = old.cells;
        newton = linear = rhs_calls = preconditioner_calls = 0;
        matrix_seconds = linear_seconds = 0;
        reason = "newton_limit";
        auto rhs = [&](hd_core::MeshBlock& state) {
            boundary::FillBoundary(state, time+dt, hooks);
            double rate, source_limit;
            auto value = hd_time::evaluate_RHS(state, time+dt, scheme, sources, solver, cfg,
                                              cfg.cfl, rate, source_limit);
            ++rhs_calls;
            // Keep user-reported source limits meaningful without changing the API.
            if (dt > source_limit) throw std::runtime_error("Implicit step exceeds source timestep limit");
            return std::move(value.second);
        };
        auto residual = [&](const hd_core::MeshBlock& state, const std::vector<hd_core::Conserved>& rates) {
            Vector r(active.size()*nv);
            for (std::size_t i = 0; i < active.size(); ++i)
                for (int k = 0; k < nv; ++k)
                    r[i*nv+k] = state.cells[active[i]].U[k] - old.cells[active[i]].U[k] - dt*rates[active[i]][k];
            return r;
        };
        auto rates = rhs(trial);
        auto r = residual(trial, rates);
        for (; newton <= settings.implicit_newton_max_iterations; ++newton)
        {
            LogVariables coordinates(old, trial, active);
            residual_norm = coordinates.norm(r);
            if (settings.implicit_verbose)
                std::cerr << "[implicit] Newton=" << newton << " t=" << time << " dt=" << dt
                          << " residual=" << residual_norm << '\n';
            if (residual_norm <= settings.implicit_tolerance) { reason = "converged"; return true; }
            if (!std::isfinite(residual_norm) || newton == settings.implicit_newton_max_iterations) return false;
            build_preconditioner(coordinates, dt, solver, scheme, cfg);
            const detail::LinearAction action = [&](const Vector& v, Vector& result) {
                double largest = 0;
                for (double value : v) largest = std::max(largest, std::abs(value));
                result.assign(v.size(), 0.0);
                if (largest == 0) return;
                const double h = difference_step/largest;
                probe.cells = minus.cells = trial.cells;
                for (std::size_t i = 0; i < active.size(); ++i)
                {
                    probe.cells[active[i]].U = coordinates.update(i, trial.cells[active[i]], &v[i*nv], h);
                    minus.cells[active[i]].U = coordinates.update(i, trial.cells[active[i]], &v[i*nv], -h);
                }
                const auto plus_rhs = rhs(probe), minus_rhs = rhs(minus);
                for (std::size_t i = 0; i < active.size(); ++i)
                {
                    const auto du = coordinates.direction(i, trial.cells[active[i]], &v[i*nv]);
                    for (int k = 0; k < nv; ++k)
                        result[i*nv+k] = (du[k]-dt*(plus_rhs[active[i]][k]-minus_rhs[active[i]][k])/(2*h))
                                        / coordinates.row_scale(i*nv+k);
                }
            };
            Vector b(r.size()), correction;
            for (std::size_t i = 0; i < r.size(); ++i) b[i] = -r[i]/coordinates.row_scale(i);
            int iterations = 0;
            const auto start = Clock::now();
            // Inexact Newton: avoid solving each early linear system to 1e-9.
            const bool solved = detail::solve_gmres(preconditioner, action, b, correction, 0.01, iterations);
            linear_seconds += seconds(start);
            linear += iterations;
            if (!solved) { reason = "linear_solve"; return false; }
            double largest = 0;
            for (double value : correction) largest = std::max(largest, std::abs(value));
            const double start_alpha = largest > 2 ? 2/largest : 1;
            const double old_merit = merit(r, coordinates);
            bool accepted = false;
            for (double alpha = start_alpha; alpha >= start_alpha/1048576; alpha *= 0.5)
            {
                candidate.cells = trial.cells;
                try
                {
                    for (std::size_t i = 0; i < active.size(); ++i)
                        candidate.cells[active[i]].U = coordinates.update(i, trial.cells[active[i]], &correction[i*nv], alpha);
                }
                catch (const std::invalid_argument&) { continue; }
                auto next_rates = rhs(candidate);
                auto next_r = residual(candidate, next_rates);
                if (merit(next_r, coordinates) <= (1-1e-4*alpha)*old_merit)
                {
                    trial.cells.swap(candidate.cells);
                    rates.swap(next_rates);
                    r.swap(next_r);
                    accepted = true;
                    break;
                }
            }
            if (!accepted) { reason = "line_search"; return false; }
        }
        return false;
    }
};

Integrator::Integrator(const hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    : impl_(std::make_unique<Impl>(mesh, cfg)) {}
Integrator::~Integrator() = default;
double Integrator::suggested_dt() const { return impl_->next_dt; }

double Integrator::advance(hd_core::MeshBlock& mesh, double time, double remaining,
    hydro::SpatialScheme scheme, const hydro::RiemannSolver& solver,
    const hd_source::SourceList& sources, const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks)
{
    auto& p = *impl_;
    if (!std::isfinite(time) || !std::isfinite(remaining) || remaining <= 0
        || mesh.face != p.trial.face || mesh.nactive != p.trial.nactive || mesh.nghost != p.trial.nghost
        || mesh.geometry != p.trial.geometry)
        throw std::invalid_argument("Implicit requires finite time and unchanged mesh geometry");
    double dt = std::min(remaining, p.settings.implicit_adaptive ? p.next_dt : p.settings.dt_implicit);
    for (int attempt = 0; attempt < (p.settings.implicit_adaptive ? 16 : 1); ++attempt)
    {
        if (!std::isfinite(time+dt) || time+dt <= time)
            throw std::runtime_error("Implicit timestep cannot advance time");
        const auto start = Clock::now();
        const bool accepted = p.solve(mesh, time, dt, scheme, solver, sources, cfg, hooks);
        const double elapsed = seconds(start);
        std::cerr << "[implicit] " << (accepted ? "accepted" : "rejected") << " t=" << time
                  << " dt=" << dt << " Newton=" << p.newton << " GMRES=" << p.linear
                  << " residual=" << p.residual_norm << " wall_s=" << elapsed << " " << p.reason << '\n';
        if (p.report)
        {
            p.report << std::setprecision(17) << time << ',' << dt << ',' << accepted << ','
                     << p.newton << ',' << p.linear << ',' << p.rhs_calls << ',' << p.preconditioner_calls
                     << ',' << p.residual_norm << ',' << elapsed << ',' << p.matrix_seconds << ','
                     << p.linear_seconds << ',' << p.reason << '\n';
            p.report.flush();
        }
        if (accepted)
        {
            p.trial.synchronize_dual_energy();
            boundary::FillBoundary(p.trial, time + dt, hooks);
            mesh.cells.swap(p.trial.cells);
            p.next_dt = std::min(p.settings.dt_implicit, dt*(p.newton <= 8 ? 1.5 : 1.0));
            return dt;
        }
        dt *= 0.5;
    }
    throw std::runtime_error("Standalone implicit failed ("+p.reason+"); accepted state unchanged");
}
}
