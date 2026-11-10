/**
 * @file lts_euler.cpp
 * @brief Experimental multilevel, conservative, first-order explicit subcycling.
 *
 * Fixed spatial grid, 1D Godunov with PCM/PLM/PLM_corrected. Each macro step starts
 * synchronized. Time levels take H/2^level, with neighbour buffers. PCM holds
 * coarse boundary data at their local step's start state. PLM predicts coarse
 * states to the fine-stage time before
 * reconstructing the entire stencil. Time integration remains first order.
 * Flux registers include face area and time, and synchronize both sides
 * of every interface, including the periodic seam. Failed trials never commit.
 *
 * All source rates are evaluated on the synchronized start state and held for
 * H, including geometric pressure and self-gravity. This deliberate first-order
 * approximation avoids mixed-time enclosed masses and repeated global source
 * scans. Source callbacks' hard limits and independent source-change estimates
 * cap H. Refinement in H is still required to assess source accuracy.
 *
 * Dynamic wave checks include the slow side's FULL step; they are necessary
 * safeguards, not a proof of nonlinear LTS stability for arbitrary shocks.
 */
#include "integrator.h"
#include "hydro/riemann/riemann.h"
#include "hydro/wave_timestep.h"
#include "utils.h"
#include "utils/performance.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace hd_time
{
namespace
{
using hd_core::Conserved;
using hd_core::MeshBlock;
constexpr double infinity = std::numeric_limits<double>::infinity();

struct FaceResult
{
    Conserved flux{};
    double speed = 0.0;
    double velocity = 0.0;
};

// Both interface outputs are stored directly in the caller-owned result.
FaceResult solve_face(const MeshBlock& mesh, hydro::ReconstructionCache& reconstruction, int face,
                      const hydro::RiemannSolver& solver, double& pressure_hint)
{
    const int right = mesh.nghost[0] + face;
    const auto& L = mesh.cells[right - 1];
    const auto& R = mesh.cells[right];
    hd_core::Index index{};
    index[0] = right;
    const auto [WL, WR] = performance::measure(performance::Section::Reconstruction,
        [&] { return reconstruction.face(0, index); }, true);
    FaceResult result;
    result.speed = -1.0;
    result.flux = performance::measure(performance::Section::Riemann, [&] {
        return solver.flux_with_pressure_hint(WL, WR, 0, L.gamma, R.gamma, pressure_hint, result.speed, result.velocity);
    }, true);
    if (!std::isfinite(result.speed) || result.speed < 0.0)
        throw std::runtime_error("LTS received an invalid Riemann wave speed");
    for (double x : result.flux.data)
        if (!std::isfinite(x))
            throw std::runtime_error("LTS received a non-finite face flux");
    return result;
}

bool physical(const Conserved& U, double gamma)
{
    try { U.validate(gamma); }
    catch (const std::invalid_argument&) { return false; }
    return true;
}

// Independent estimate per source: do not hide large terms by cancellation.
// Bound density change, momentum relative to rho*(|v|+cs), and internal-energy
// change. The acoustic momentum scale also bounds acceleration starting at v=0.
double source_change_limit(const MeshBlock& mesh, const std::vector<Conserved>& S,
                          std::vector<std::array<double, 3>>& scales, bool& cache_ready)
{
    double limit = infinity;
    auto bound = [&](double scale, double rate) {
        if (rate != 0.0)
            limit = std::min(limit, 0.1 * scale / std::abs(rate));
    };
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        const int k = mesh.nghost[0] + i;
        const auto& s = S[k];
        for (double x : s.data)
            if (!std::isfinite(x))
                throw std::runtime_error("LTS received a non-finite source rate");
        // rho, momentum scale, internal energy.
        // Without sources the existing single zero-rate scan needs no storage.
        std::array<double, 3> uncached{};
        auto& cached = scales.empty() ? uncached : scales[i];
        if (!cache_ready || scales.empty())
        {
            const auto& cell = mesh.cells[k];
            const auto W = cell.U.to_Primtivies(cell.gamma);
            const double rho = W[hd_core::RHO];
            cached = {rho, 0.0, 0.0};
            if (rho != 0.0)
            {
                const double v = W[hd_core::VEL(0)];
                const double cs = std::sqrt(cell.gamma * W[hd_core::PRE] / rho);
                cached = {rho, rho * (std::abs(v) + cs),
                          W[hd_core::PRE] / (cell.gamma - 1.0)};
            }
        }
        const auto [rho, momentum, internal] = cached;
        if (rho == 0.0)
        {
            for (double x : s.data)
                if (x != 0.0)
                    throw std::invalid_argument("LTS does not support sources creating gas in vacuum");
            continue;
        }
        bound(rho, s[hd_core::RHO]);
        bound(momentum, s[hd_core::MOMDEN(0)]);
        bound(internal, s[hd_core::INTERNAL_ENERGY]);
    }
    cache_ready = true;
    return limit;
}
}

struct LTSEulerWorkspace::Impl
{
    MeshBlock start, work, candidate;
    std::unique_ptr<hydro::ReconstructionCache> initial_reconstruction, sub_reconstruction;
    std::vector<double> local_dt, step_fraction;
    // One scalar per face, plus transactional copies for macro steps/retries.
    hydro::StarPressureCache star_pressure_cache, initial_pressure, trial_pressure;
    std::vector<FaceResult> initial, stage;
    std::vector<Conserved> source, contribution, pending_flux, step_start;
    // Reuse storage across macro steps; values belong only to their start state.
    std::vector<std::array<double, 3>> source_scales;
    std::vector<int> level, face_level;
    std::vector<unsigned char> predicted, lower_ghosts, upper_ghosts;
    std::vector<std::vector<int>> level_cells, level_faces, prediction_cells;

    void prepare(const MeshBlock& mesh)
    {
        if (start.nactive != mesh.nactive || start.nghost != mesh.nghost ||
            start.geometry != mesh.geometry || start.face != mesh.face ||
            start.volume_cell != mesh.volume_cell || start.face_area[0].faces != mesh.face_area[0].faces)
        {
            mesh.validate_geometry();
            start = mesh;
            work = mesh;
            candidate = mesh;
            star_pressure_cache.reset(0, mesh.nactive[0] + 1);
        }
        start.bc_lower = work.bc_lower = candidate.bc_lower = mesh.bc_lower;
        start.bc_upper = work.bc_upper = candidate.bc_upper = mesh.bc_upper;
        start.cells = mesh.cells; // Transaction snapshot: geometry is not copied.
        const auto n = static_cast<std::size_t>(mesh.nactive[0]);
        local_dt.resize(n); step_fraction.resize(n);
        initial.resize(n + 1); stage.resize(n + 1);
        source.resize(mesh.cells.size()); contribution.resize(mesh.cells.size());
        pending_flux.resize(n); step_start.resize(n);
        level.resize(n); face_level.resize(n + 1); predicted.resize(n);
    }
};

LTSEulerWorkspace::LTSEulerWorkspace() : impl(std::make_unique<Impl>()) {}
LTSEulerWorkspace::~LTSEulerWorkspace() = default;

double advance_lts_euler(
    MeshBlock& mesh, double time, double dt_limit,
    hydro::SpatialScheme spatial_scheme, const hydro::RiemannSolver& solver,
    const hd_source::SourceList& source_list, const config::SimulationConfig& cfg,
    const hd_user::UserHooks& hooks, double cfl, LTSEulerWorkspace* workspace)
{
    if (hd_core::NDIM != 1 || spatial_scheme != hydro::SpatialScheme::Godunov)
        throw std::invalid_argument("lts_euler requires 1D Godunov");
    if (!utils::iequals(cfg.reconstructor, "PCM") &&
        !utils::iequals(cfg.reconstructor, "PLM") &&
        !utils::iequals(cfg.reconstructor, "PLM_corrected"))
        throw std::invalid_argument("lts_euler requires PCM, PLM or PLM_corrected reconstruction");
    const auto method = hydro::parse_reconstructor(cfg.reconstructor);
    const auto limiter = hydro::parse_limiter(cfg.limiter);
    const bool linear = method != hydro::Reconstructor::PCM;
    if (!std::isfinite(time) || !std::isfinite(dt_limit) || dt_limit <= 0.0 ||
        !std::isfinite(cfl) || !std::isfinite(cfg.cfl_max) ||
        cfl <= 0.0 || cfl > cfg.cfl_max || cfg.cfl_max > 0.5 ||
        cfg.integrator_retry_max < 0 || cfg.lts_max_levels < 1 || cfg.lts_max_levels > 16)
        throw std::invalid_argument("Invalid LTS time, CFL or retry settings");
    (void)mesh.flat_index(hd_core::Index{}, true);
    if (mesh.nghost[0] < (linear ? 2 : 1))
        throw std::invalid_argument("LTS requires one ghost cell for PCM, two for PLM variants");

    const int n = mesh.nactive[0], ng = mesh.nghost[0];
    const bool periodic = mesh.bc_lower[0] == boundary::BoundaryCondition::Periodic;
    if (periodic != (mesh.bc_upper[0] == boundary::BoundaryCondition::Periodic))
        throw std::invalid_argument("LTS requires paired periodic boundaries");
    if (periodic && mesh.geometry != hd_core::Geometry::Cartesian)
        throw std::invalid_argument("LTS periodic boundaries require Cartesian geometry");

    std::optional<LTSEulerWorkspace> temporary_workspace;
    if (!workspace)
    {
        temporary_workspace.emplace();
        workspace = &*temporary_workspace;
    }
    auto& scratch = *workspace->impl;
    scratch.prepare(mesh);
    auto& start = scratch.start;
    auto& work = scratch.work;
    auto& candidate = scratch.candidate;
    boundary::FillBoundary(start, time, hooks);
    const double* area = start.face_area[0].faces.data() + ng;
    const double* volume = start.volume_cell.data() + ng;
    auto& local_dt = scratch.local_dt;
    auto& initial = scratch.initial;
    if (!scratch.initial_reconstruction)
        scratch.initial_reconstruction = std::make_unique<hydro::ReconstructionCache>(
            start, method, limiter, cfg.limiter_beta);
    else
        scratch.initial_reconstruction->reset(start, method, limiter, cfg.limiter_beta);
    scratch.initial_pressure = scratch.star_pressure_cache;
    for (int f = 0; f <= n; ++f)
    {
        if (periodic && f == n)
        {
            initial[f] = initial[0];
            scratch.initial_pressure[0][f] = scratch.initial_pressure[0][0];
        }
        else
            initial[f] = solve_face(start, *scratch.initial_reconstruction, f, solver,
                                    scratch.initial_pressure[0][f]);
    }
    double min_dt = infinity;
    hydro::max_inverse_dt = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double rate = std::max(initial[i].speed, initial[i + 1].speed) / start.cell_width(0, ng + i);
        if (!std::isfinite(rate))
            throw std::runtime_error("Non-finite LTS CFL rate");
        hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
        local_dt[i] = rate > 0.0 ? cfl / rate : infinity;
        min_dt = std::min(min_dt, local_dt[i]);
    }

    auto& source = scratch.source;
    std::fill(source.begin(), source.end(), Conserved{});
    double source_limit = infinity;
    scratch.source_scales.resize(source_list.empty() ? 0 : n);
    // Refresh on every macro step, including after failures or workspace reuse.
    // Retries keep the same start state and already computed frozen source rates.
    bool source_cache_ready = false;
    for (const auto& term : source_list)
    {
        if (!term.func) throw std::invalid_argument("Missing LTS source function");
        auto& contribution = scratch.contribution;
        std::fill(contribution.begin(), contribution.end(), Conserved{});
        double limit = infinity;
        term.func(start, time, cfg, contribution, limit);
        hd_source::complete_internal_source(start, term, contribution);
        if (contribution.size() != source.size())
            throw std::logic_error("Source changed LTS RHS storage size");
        if (std::isnan(limit) || limit <= 0.0)
            throw std::invalid_argument("Invalid LTS source timestep limit");
        source_limit = std::min({source_limit, limit,
            source_change_limit(start, contribution, scratch.source_scales, source_cache_ready)});
        for (int i = 0; i < n; ++i) source[ng + i] += contribution[ng + i];
    }
    source_limit = std::min(source_limit,
        source_change_limit(start, source, scratch.source_scales, source_cache_ready));

    const int max_depth = cfg.lts_max_levels - 1;
    double H = std::min({dt_limit, std::ldexp(min_dt, max_depth), source_limit});
    auto& level = scratch.level;
    auto& face_level = scratch.face_level;
    auto& stage = scratch.stage;
    auto& pending = scratch.pending_flux;
    auto& step_start = scratch.step_start;
    auto& step_fraction = scratch.step_fraction;
    std::vector<double> pressure_work(n);
    auto left = [&](int f) { return f == 0 ? (periodic ? n - 1 : -1) : f - 1; };
    auto right = [&](int f) { return f == n ? (periodic ? 0 : -1) : f; };
    const int unique_faces = periodic ? n : n + 1;
    const int buffer = linear ? 2 : 1;
    for (long long attempt = 0; attempt <= cfg.integrator_retry_max; ++attempt, H *= 0.5)
    {
        if (!std::isfinite(H) || H <= 0.0 || !std::isfinite(time + H) || time + H <= time)
            throw std::runtime_error("LTS macro step cannot advance time");
        std::fill(level.begin(), level.end(), 0);
        for (int i = 0; i < n; ++i)
        {
            int depth = 0;
            while (depth < max_depth && std::ldexp(H, -depth) > local_dt[i]) ++depth;
            // Preserve the original fast-region buffer at every required depth.
            for (int offset = -buffer; offset <= buffer; ++offset)
            {
                int j = i + offset;
                if (periodic) j = (j % n + n) % n;
                if (j >= 0 && j < n) level[j] = std::max(level[j], depth);
            }
        }
        const int deepest = *std::max_element(level.begin(), level.end());
        // Grade transitions: every reconstruction neighbour is at most one
        // level coarser. Descending passes propagate buffers without substep scans.
        for (int depth = deepest; depth > 1; --depth)
            for (int i = 0; i < n; ++i)
                if (level[i] == depth)
                    for (int offset = -buffer; offset <= buffer; ++offset)
                    {
                        int j = i + offset;
                        if (periodic) j = (j % n + n) % n;
                        if (j >= 0 && j < n) level[j] = std::max(level[j], depth - 1);
                    }

        auto& cells = scratch.level_cells;
        auto& faces = scratch.level_faces;
        auto& predictions = scratch.prediction_cells;
        cells.resize(deepest + 1); faces.resize(deepest + 1); predictions.resize(deepest + 1);
        scratch.lower_ghosts.assign(deepest + 1, false);
        scratch.upper_ghosts.assign(deepest + 1, false);
        for (int depth = 0; depth <= deepest; ++depth)
        {
            cells[depth].clear(); faces[depth].clear(); predictions[depth].clear();
        }
        for (int i = 0; i < n; ++i) cells[level[i]].push_back(i);
        for (int f = 0; f < unique_faces; ++f)
        {
            face_level[f] = std::max(left(f) >= 0 ? level[left(f)] : 0,
                                     right(f) >= 0 ? level[right(f)] : 0);
            for (int depth = 0; depth <= face_level[f]; ++depth) faces[depth].push_back(f);
        }
        if (periodic) face_level[n] = face_level[0];

        for (int depth = 1; depth <= deepest; ++depth)
        {
            std::fill(scratch.predicted.begin(), scratch.predicted.end(), false);
            bool global_boundary = false;
            auto mark = [&](int i) {
                if (level[i] < depth && !scratch.predicted[i])
                {
                    scratch.predicted[i] = true;
                    predictions[depth].push_back(i);
                }
            };
            for (int f : faces[depth])
                for (int i = f - (linear ? 2 : 1); i <= f + (linear ? 1 : 0); ++i)
                {
                    if (i >= 0 && i < n) { mark(i); continue; }
                    const bool lower = i < 0;
                    scratch.lower_ghosts[depth] |= lower;
                    scratch.upper_ghosts[depth] |= !lower;
                    const auto bc = lower ? start.bc_lower[0] : start.bc_upper[0];
                    if (bc == boundary::BoundaryCondition::Periodic) mark((i % n + n) % n);
                    else if (bc == boundary::BoundaryCondition::Reflective)
                        mark(lower ? -i - 1 : 2 * n - i - 1);
                    else if (bc == boundary::BoundaryCondition::User) global_boundary = true;
                }
            // A custom boundary can inspect any active cell, not just its stencil.
            if (global_boundary) for (int i = 0; i < n; ++i) mark(i);
        }

        work.cells = start.cells;
        candidate.cells = start.cells;
        scratch.trial_pressure = scratch.initial_pressure;
        stage = initial;
        std::fill(pending.begin(), pending.end(), Conserved{});

        // Recursively take two child steps. No extra integrator/level class:
        // lists, local predictors and flux registers live in the existing workspace.
        auto advance_level = [&](auto&& self, int depth, double fraction, bool reuse_flux) -> bool {
            const double h = std::ldexp(H, -depth);
            const double span = std::ldexp(1.0, -depth);
            const double stage_time = time + fraction * H;
            if (time + (fraction + span) * H <= stage_time)
                throw std::runtime_error("LTS substep cannot advance time");
            if (!reuse_flux)
            {
                for (int i : predictions[depth])
                {
                    auto& cell = work.cells[ng + i];
                    const double theta = linear ? std::ldexp(fraction - step_fraction[i], level[i]) : 0.0;
                    cell.U = (1.0 - theta) * step_start[i] + theta * candidate.cells[ng + i].U;
                    if (!physical(cell.U, cell.gamma)) return false;
                }
                if (scratch.lower_ghosts[depth])
                    boundary::ApplyBoundaryCondition(work, 0, boundary::BoundarySide::Lower, stage_time, hooks);
                if (scratch.upper_ghosts[depth])
                    boundary::ApplyBoundaryCondition(work, 0, boundary::BoundarySide::Upper, stage_time, hooks);
                if (!scratch.sub_reconstruction)
                    scratch.sub_reconstruction = std::make_unique<hydro::ReconstructionCache>(
                        work, method, limiter, cfg.limiter_beta);
                else
                    scratch.sub_reconstruction->reset(work, method, limiter, cfg.limiter_beta);
                for (int f : faces[depth])
                    stage[f] = solve_face(work, *scratch.sub_reconstruction, f, solver,
                                          scratch.trial_pressure[0][f]);
                if (periodic)
                {
                    stage[n] = stage[0];
                    scratch.trial_pressure[0][n] = scratch.trial_pressure[0][0];
                }
            }
            for (int f : faces[depth])
            {
                for (int i : {left(f), right(f)})
                    if (i >= 0)
                    {
                        const double rate = stage[f].speed / start.cell_width(0, ng + i);
                        hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
                        // Check EACH side with its own full local step, not h.
                        if (!std::isfinite(rate) || std::ldexp(H, -level[i]) * rate > cfg.cfl_max)
                            return false;
                    }
                if (face_level[f] == depth)
                {
                    // Integrate a physical face only at the finer adjacent level.
                    // Both sides receive exactly the same area*time-weighted flux.
                    const auto integral = stage[f].flux * (h * area[f]);
                    if (left(f) >= 0) pending[left(f)] -= integral;
                    if (right(f) >= 0) pending[right(f)] += integral;
                }
            }
            for (int i : cells[depth])
            {
                step_start[i] = work.cells[ng + i].U;
                step_fraction[i] = fraction;
                pressure_work[i] = -(work.cells[ng + i].gamma - 1.0) * step_start[i].internal_energy()
                    * (area[i + 1] * stage[i + 1].velocity - area[i] * stage[i].velocity) / volume[i];
                auto& cell = candidate.cells[ng + i];
                cell.U = step_start[i]
                    + (stage[i].flux * area[i] - stage[i + 1].flux * area[i + 1]) * (h / volume[i])
                    + source[ng + i] * h;
                cell.U[hd_core::INTERNAL_ENERGY] += h * pressure_work[i];
                if (!physical(cell.U, cell.gamma)) return false;
            }
            if (depth < deepest)
            {
                if (!self(self, depth + 1, fraction, true)) return false;
                if (!self(self, depth + 1, fraction + span / 2.0, false)) return false;
            }
            for (int i : cells[depth])
            {
                auto& cell = work.cells[ng + i];
                cell.U = step_start[i] + pending[i] / volume[i] + source[ng + i] * h;
                cell.U[hd_core::INTERNAL_ENERGY] += h * pressure_work[i];
                if (!physical(cell.U, cell.gamma)) return false;
                cell.synchronize_dual_energy();
                pending[i] = Conserved{};
            }
            return true;
        };
        if (!advance_level(advance_level, 0, 0.0, true)) continue;

        // All cells are synchronized; exceptions cannot commit partial states/history.
        boundary::FillBoundary(work, time + H, hooks);
        mesh.cells.swap(work.cells);
        scratch.star_pressure_cache.swap(scratch.trial_pressure);
        return H;
    }
    throw std::runtime_error("LTS Euler retry limit exceeded");
}
}
