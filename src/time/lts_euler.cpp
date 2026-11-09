/**
 * @file lts_euler.cpp
 * @brief Experimental two-level, conservative, first-order explicit subcycling.
 *
 * Fixed spatial grid, 1D Godunov with PCM/PLM/PLM_corrected. Each macro step starts
 * synchronized. Fast cells and a reconstruction-dependent neighbour buffer take
 * two half steps; other cells take one full step. PCM holds slow boundary data
 * at their start state. PLM predicts slow states to the fine-stage time before
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
};

// Keep the existing solver ABI; capture its wave speed immediately, locally.
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
    hydro::wave_timestep.face_speed = -1.0;
    FaceResult result;
    result.flux = performance::measure(performance::Section::Riemann, [&] {
        return solver.flux_with_pressure_hint(WL, WR, 0, L.gamma, R.gamma, pressure_hint);
    }, true);
    result.speed = hydro::wave_timestep.face_speed;
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
double source_change_limit(const MeshBlock& mesh, const std::vector<Conserved>& S)
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
        const auto& cell = mesh.cells[k];
        const auto W = cell.U.to_Primtivies(cell.gamma);
        const double rho = W[hd_core::RHO];
        if (rho == 0.0)
        {
            for (double x : s.data)
                if (x != 0.0)
                    throw std::invalid_argument("LTS does not support sources creating gas in vacuum");
            continue;
        }
        const double v = W[hd_core::VEL(0)];
        const double cs = std::sqrt(cell.gamma * W[hd_core::PRE] / rho);
        const double internal = W[hd_core::PRE] / (cell.gamma - 1.0);
        bound(rho, s[hd_core::RHO]);
        bound(rho * (std::abs(v) + cs), s[hd_core::MOMDEN(0)]);
        bound(internal, s[hd_core::ENERGY] - v * s[hd_core::MOMDEN(0)]
                       + 0.5 * v * v * s[hd_core::RHO]);
    }
    return limit;
}
}

struct LTSEulerWorkspace::Impl
{
    MeshBlock start, work, candidate;
    std::unique_ptr<hydro::ReconstructionCache> initial_reconstruction, sub_reconstruction;
    std::vector<double> local_dt;
    // One scalar per face, plus transactional copies for macro steps/retries.
    std::array<std::vector<double>, hd_core::NDIM> pressure_history, initial_pressure, trial_pressure;
    std::vector<FaceResult> initial;
    std::vector<Conserved> source, contribution, fine_integral, coarse_integral, delta, next_fast;
    std::vector<unsigned char> fast, predicted;
    std::vector<int> fast_cells, fast_faces, prediction_cells;

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
            pressure_history[0].assign(mesh.nactive[0] + 1, 0.0);
        }
        start.bc_lower = work.bc_lower = candidate.bc_lower = mesh.bc_lower;
        start.bc_upper = work.bc_upper = candidate.bc_upper = mesh.bc_upper;
        start.cells = mesh.cells; // Transaction snapshot: geometry is not copied.
        const auto n = static_cast<std::size_t>(mesh.nactive[0]);
        local_dt.resize(n); initial.resize(n + 1);
        initial_pressure[0].resize(n + 1); trial_pressure[0].resize(n + 1);
        source.resize(mesh.cells.size()); contribution.resize(mesh.cells.size());
        fine_integral.resize(n + 1); coarse_integral.resize(n + 1); delta.resize(n + 1);
        next_fast.resize(n); fast.resize(n); predicted.resize(n);
        fast_cells.reserve(n); fast_faces.reserve(n + 1); prediction_cells.reserve(n);
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
        cfg.integrator_retry_max < 0)
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
    scratch.initial_pressure = scratch.pressure_history;
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
    hydro::wave_timestep.max_inverse_dt = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double rate = std::max(initial[i].speed, initial[i + 1].speed) / start.cell_width(0, ng + i);
        if (!std::isfinite(rate))
            throw std::runtime_error("Non-finite LTS CFL rate");
        hydro::wave_timestep.max_inverse_dt = std::max(hydro::wave_timestep.max_inverse_dt, rate);
        local_dt[i] = rate > 0.0 ? cfl / rate : infinity;
        min_dt = std::min(min_dt, local_dt[i]);
    }

    auto& source = scratch.source;
    std::fill(source.begin(), source.end(), Conserved{});
    double source_limit = infinity;
    for (const auto& term : source_list)
    {
        if (!term.func) throw std::invalid_argument("Missing LTS source function");
        auto& contribution = scratch.contribution;
        std::fill(contribution.begin(), contribution.end(), Conserved{});
        double limit = infinity;
        term.func(start, time, cfg, contribution, limit);
        if (contribution.size() != source.size())
            throw std::logic_error("Source changed LTS RHS storage size");
        if (std::isnan(limit) || limit <= 0.0)
            throw std::invalid_argument("Invalid LTS source timestep limit");
        source_limit = std::min({source_limit, limit, source_change_limit(start, contribution)});
        for (int i = 0; i < n; ++i) source[ng + i] += contribution[ng + i];
    }
    source_limit = std::min(source_limit, source_change_limit(start, source));

    double H = std::min({dt_limit, 2.0 * min_dt, source_limit});
    auto& fine_integral = scratch.fine_integral;
    auto& coarse_integral = scratch.coarse_integral;
    auto& delta = scratch.delta;
    auto& next_fast = scratch.next_fast;
    for (long long attempt = 0; attempt <= cfg.integrator_retry_max; ++attempt, H *= 0.5)
    {
        if (!std::isfinite(H) || H <= 0.0 || !std::isfinite(time + H) || time + H <= time)
            throw std::runtime_error("LTS macro step cannot advance time");
        auto& fast = scratch.fast;
        std::fill(fast.begin(), fast.end(), false);
        for (int i = 0; i < n; ++i)
            if (H > local_dt[i])
            {
                // PLM faces depend on one additional cell on each side.
                const int buffer = linear ? 2 : 1;
                for (int offset = -buffer; offset <= buffer; ++offset)
                {
                    int j = i + offset;
                    if (periodic) j = (j % n + n) % n;
                    if (j >= 0 && j < n) fast[j] = true;
                }
            }
        auto& fast_cells = scratch.fast_cells;
        auto& fast_faces = scratch.fast_faces;
        fast_cells.clear(); fast_faces.clear();
        for (int i = 0; i < n; ++i) if (fast[i]) fast_cells.push_back(i);
        if (!fast_cells.empty() && time + H / 2.0 <= time)
            throw std::runtime_error("LTS half step cannot advance time");
        auto left = [&](int f) { return f == 0 ? (periodic ? n - 1 : -1) : f - 1; };
        auto right = [&](int f) { return f == n ? (periodic ? 0 : -1) : f; };
        for (int f = 0; f < (periodic ? n : n + 1); ++f)
            if ((left(f) >= 0 && fast[left(f)]) || (right(f) >= 0 && fast[right(f)]))
                fast_faces.push_back(f);

        // Union of active slow cells actually read by fine-face stencils.
        // Include the donors of any required periodic/reflecting ghost cells.
        auto& predicted = scratch.predicted;
        auto& prediction_cells = scratch.prediction_cells;
        std::fill(predicted.begin(), predicted.end(), false);
        prediction_cells.clear();
        bool lower_ghosts = false, upper_ghosts = false, global_boundary = false;
        auto mark_prediction = [&](int i) {
            if (!fast[i] && !predicted[i])
            {
                predicted[i] = true;
                prediction_cells.push_back(i);
            }
        };
        for (int f : fast_faces)
            for (int i = f - (linear ? 2 : 1); i <= f + (linear ? 1 : 0); ++i)
            {
                if (i >= 0 && i < n)
                {
                    if (linear) mark_prediction(i);
                    continue;
                }
                const bool lower = i < 0;
                lower_ghosts = lower_ghosts || lower;
                upper_ghosts = upper_ghosts || !lower;
                const auto bc = lower ? start.bc_lower[0] : start.bc_upper[0];
                if (!linear) continue;
                if (bc == boundary::BoundaryCondition::Periodic)
                    mark_prediction((i % n + n) % n);
                else if (bc == boundary::BoundaryCondition::Reflective)
                    mark_prediction(lower ? -i - 1 : 2 * n - i - 1);
                else if (bc == boundary::BoundaryCondition::User)
                    global_boundary = true; // Arbitrary userbc may read any active cell.
            }
        if (global_boundary)
            for (int i = 0; i < n; ++i) mark_prediction(i);

        work.cells = start.cells;
        candidate.cells = start.cells;
        scratch.trial_pressure = scratch.initial_pressure;
        std::fill(fine_integral.begin(), fine_integral.end(), Conserved{});
        for (int f = 0; f <= n; ++f) coarse_integral[f] = initial[f].flux * (H * area[f]);
        bool valid = true;
        for (int i = 0; i < n; ++i)
            if (!fast[i])
            {
                auto& cell = candidate.cells[ng + i];
                cell.U += (coarse_integral[i] - coarse_integral[i + 1]) / volume[i]
                          + source[ng + i] * H;
                valid = valid && physical(cell.U, cell.gamma);
            }
        if (!valid) continue;

        const double h = H / 2.0;
        for (int sub = 0; sub < 2 && valid && !fast_cells.empty(); ++sub)
        {
            if (sub != 0)
            {
                if (linear)
                {
                    // Only cells read by PLM/boundaries need to represent t+h:
                    // fast cells completed Euler(h); slow cells interpolate
                    // between their start and provisional Euler(H) states.
                    // This is prediction only, never a second committed update.
                    for (int i : prediction_cells)
                    {
                        auto& cell = work.cells[ng + i];
                        cell.U = 0.5 * start.cells[ng + i].U
                               + 0.5 * candidate.cells[ng + i].U;
                        if (!physical(cell.U, cell.gamma)) valid = false;
                    }
                    if (!valid) break;
                }
                // Unused ghosts do not need a half-time fill. Required custom
                // boundaries above retain the full synchronized prediction view.
                if (lower_ghosts)
                    boundary::ApplyBoundaryCondition(work, 0, boundary::BoundarySide::Lower, time + h, hooks);
                if (upper_ghosts)
                    boundary::ApplyBoundaryCondition(work, 0, boundary::BoundarySide::Upper, time + h, hooks);
                // No slope/primitive cache may survive a state or ghost update.
                if (!scratch.sub_reconstruction)
                    scratch.sub_reconstruction = std::make_unique<hydro::ReconstructionCache>(
                        work, method, limiter, cfg.limiter_beta);
                else
                    scratch.sub_reconstruction->reset(work, method, limiter, cfg.limiter_beta);
            }
            for (int f : fast_faces)
            {
                const auto result = sub == 0 ? initial[f]
                    : solve_face(work, *scratch.sub_reconstruction, f, solver, scratch.trial_pressure[0][f]);
                for (int i : {left(f), right(f)})
                    if (i >= 0)
                    {
                        const double rate = result.speed / start.cell_width(0, ng + i);
                        hydro::wave_timestep.max_inverse_dt =
                            std::max(hydro::wave_timestep.max_inverse_dt, rate);
                        // Crucially check a slow neighbour with H, not h.
                        if (!std::isfinite(rate) || (fast[i] ? h : H) * rate > cfg.cfl_max)
                            valid = false;
                    }
                if (!valid) break;
                delta[f] = result.flux * (h * area[f]);
                fine_integral[f] += delta[f];
                if (periodic && f == 0)
                {
                    delta[n] = delta[0];
                    fine_integral[n] = fine_integral[0];
                    scratch.trial_pressure[0][n] = scratch.trial_pressure[0][0];
                }
            }
            if (!valid) break;
            for (int i : fast_cells)
            {
                next_fast[i] = work.cells[ng + i].U
                    + (delta[i] - delta[i + 1]) / volume[i] + source[ng + i] * h;
                if (!physical(next_fast[i], work.cells[ng + i].gamma)) valid = false;
            }
            if (!valid) break;
            for (int i : fast_cells) work.cells[ng + i].U = next_fast[i];
        }
        if (!valid) continue;

        // Replace the provisional slow-side integral at each mixed interface.
        for (int f : fast_faces)
        {
            const auto correction = fine_integral[f] - coarse_integral[f];
            const int L = left(f), R = right(f);
            if (L >= 0 && !fast[L]) candidate.cells[ng + L].U -= correction / volume[L];
            if (R >= 0 && !fast[R]) candidate.cells[ng + R].U += correction / volume[R];
        }
        for (int i : fast_cells) candidate.cells[ng + i].U = work.cells[ng + i].U;
        for (int i = 0; i < n; ++i)
            if (!physical(candidate.cells[ng + i].U, candidate.cells[ng + i].gamma)) valid = false;
        if (!valid) continue;

        // Even boundary exceptions must not leave a partially committed step.
        boundary::FillBoundary(candidate, time + H, hooks);
        mesh.cells.swap(candidate.cells);
        scratch.pressure_history.swap(scratch.trial_pressure);
        return H;
    }
    throw std::runtime_error("LTS Euler retry limit exceeded");
}
}
