/**
 * Conservative recursive Heun subcycling. See docs/lts_rk2.md.
 * No additional workspace/level class: local predictors and face registers are
 * arrays. Source callbacks see a complete approximation at their stage time.
 */
#include "integrator.h"
#include "hydro/godunov.h"
#include "hydro/wave_timestep.h"
#include "utils.h"
#include "utils/performance.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hd_time::detail
{
// Same decision as H > allowance/abs(rate), for a finite positive trial H.
// allowance is the cached, rounded 0.1*scale; twice_H is computed once as 2*H.
// A generous factor-two margin makes the multiplication-only acceptance path
// insensitive to rounding at the 10% boundary. Near that boundary (and for
// overflow/underflow or NaNs not resolved by the fast path), use the original
// division/comparison, preserving its floating-point decision. A zero rate
// imposes no bound, even when the scale is zero.
bool source_change_exceeds(double H, double twice_H, double allowance, double rate)
{
    const double magnitude = std::abs(rate);
    if (twice_H * magnitude < allowance) return false;
    return rate != 0.0 && H > allowance / magnitude;
}
}

namespace hd_time
{
double advance_lts_rk2(
    hd_core::MeshBlock& mesh, double time, double dt_limit,
    hydro::SpatialScheme spatial_scheme, const hydro::RiemannSolver& solver,
    const hd_source::SourceList& sources, const config::SimulationConfig& cfg,
    const hd_user::UserHooks& hooks, double cfl,
    hydro::StarPressureCache* star_pressure_cache)
{
    using hd_core::Conserved;
    constexpr double inf = std::numeric_limits<double>::infinity();
    if (hd_core::NDIM != 1 || spatial_scheme != hydro::SpatialScheme::Godunov)
        throw std::invalid_argument("lts_rk2 requires 1D Godunov");
    if (!utils::iequals(cfg.reconstructor, "PCM") && !utils::iequals(cfg.reconstructor, "PLM") &&
        !utils::iequals(cfg.reconstructor, "PLM_corrected"))
        throw std::invalid_argument("lts_rk2 requires PCM, PLM or PLM_corrected reconstruction");
    const auto method = hydro::parse_reconstructor(cfg.reconstructor);
    const auto limiter = hydro::parse_limiter(cfg.limiter);
    const bool linear = method != hydro::Reconstructor::PCM;
    if (!std::isfinite(time) || !std::isfinite(dt_limit) || dt_limit <= 0.0 ||
        !std::isfinite(cfl) || !std::isfinite(cfg.cfl_max) ||
        cfl <= 0.0 || cfl > cfg.cfl_max || cfg.cfl_max > 0.5 ||
        cfg.integrator_retry_max < 0 || cfg.lts_max_levels < 1 || cfg.lts_max_levels > 16)
        throw std::invalid_argument("Invalid LTS RK2 time, CFL, level or retry settings");
    mesh.validate_geometry();
    (void)mesh.flat_index(hd_core::Index{}, true);
    if (mesh.nghost[0] < (linear ? 2 : 1))
        throw std::invalid_argument("LTS RK2 needs one ghost cell for PCM, two for PLM");
    const int n = mesh.nactive[0], ng = mesh.nghost[0];
    const bool periodic = mesh.bc_lower[0] == boundary::BoundaryCondition::Periodic;
    if (periodic != (mesh.bc_upper[0] == boundary::BoundaryCondition::Periodic) ||
        (periodic && mesh.geometry != hd_core::Geometry::Cartesian))
        throw std::invalid_argument("LTS RK2 needs paired Cartesian periodic boundaries");
    const int nf = periodic ? n : n + 1;
    auto left = [&](int f) { return f == 0 ? (periodic ? n - 1 : -1) : f - 1; };
    auto right = [&](int f) { return f == n ? -1 : f; };
    auto physical = [](const Conserved& U, double gamma) {
        try { U.validate(gamma); }
        catch (const std::invalid_argument&) { return false; }
        return true;
    };
    auto start = mesh, state = mesh, view = mesh;
    boundary::FillBoundary(start, time, hooks);
    const double* area = mesh.face_area[0].faces.data() + ng;
    const double* volume = mesh.volume_cell.data() + ng;
    hydro::StarPressureCache hints;
    if (star_pressure_cache) hints = *star_pressure_cache;
    hints.ensure_size(0, n + 1);
    std::vector<Conserved> flux(n + 1), initial_flux, first_flux(n + 1);
    std::vector<double> face_velocity(n + 1), initial_velocity, work0(n);
    std::vector<double> speed(n + 1), initial_speed, local_dt(n), start_fraction(n);
    std::vector<Conserved> source(mesh.cells.size()), initial_source, contribution(mesh.cells.size());
    // Source-limit stage cache: allocate once per macro step, reuse across its
    // substeps. Each stage's first change_limit call fills it; later independent
    // sources and their sum reuse the same thermodynamic scales.
    const std::size_t source_cache_size = sources.empty() ? 0 : n;
    std::vector<double> source_rho(source_cache_size), source_density_allowance(source_cache_size),
        source_momentum_allowance(source_cache_size), source_internal_allowance(source_cache_size);
    bool source_cache_ready = false;
    std::vector<Conserved> base(n), slope(n), source0(n), pending(n);
    std::vector<int> level(n), all_faces(nf);
    for (int f = 0; f < nf; ++f) all_faces[f] = f;
    hydro::ReconstructionCache reconstruction(start, method, limiter, cfg.limiter_beta);
    auto solve = [&](const hd_core::MeshBlock& stage, const std::vector<int>& faces) {
        reconstruction.reset(stage, method, limiter, cfg.limiter_beta);
        for (int f : faces)
        {
            hd_core::Index index{};
            index[0] = ng + f;
            const auto [WL, WR] = performance::measure(performance::Section::Reconstruction,
                [&] { return reconstruction.face(0, index); }, true);
            speed[f] = -1.0;
            face_velocity[f] = 0.0;
            flux[f] = performance::measure(performance::Section::Riemann, [&] {
                return solver.flux_with_pressure_hint(WL, WR, 0,
                    stage.cells[ng + f - 1].gamma, stage.cells[ng + f].gamma, hints[0][f], speed[f], face_velocity[f]);
            }, true);
            if (!std::isfinite(speed[f]) || speed[f] < 0.0)
                throw std::runtime_error("LTS RK2 received an invalid wave speed");
            for (double x : flux[f].data)
                if (!std::isfinite(x)) throw std::runtime_error("LTS RK2 received a non-finite flux");
        }
        if (periodic)
        {
            flux[n] = flux[0]; speed[n] = speed[0]; face_velocity[n] = face_velocity[0]; hints[0][n] = hints[0][0];
        }
    };
    // Independent source-change caps prevent cancellation from hiding stiff terms.
    // Specialize the scan for its bound callback, so the choice between initial
    // step selection and stage acceptance is outside the per-cell loop.
    auto scan_source_changes = [&](const hd_core::MeshBlock& stage,
                                   const std::vector<Conserved>& rates, auto bound) {
        for (int i = 0; i < n; ++i)
        {
            const auto& s = rates[ng + i];
            for (double x : s.data)
                if (!std::isfinite(x)) throw std::runtime_error("Non-finite LTS RK2 source");
            if (!source_cache_ready)
            {
                // The stage is immutable across source callbacks. Convert and
                // validate its state, and compute sqrt(gamma*p/rho), only once.
                const auto& cell = stage.cells[ng + i];
                const auto W = cell.U.to_Primtivies(cell.gamma);
                const double rho = W[hd_core::RHO], p = W[hd_core::PRE];
                source_rho[i] = rho;
                source_density_allowance[i] = 0.1 * rho;
                source_momentum_allowance[i] = source_internal_allowance[i] = 0.0;
                if (rho != 0.0)
                {
                    const double v = W[hd_core::VEL(0)], cs = std::sqrt(cell.gamma * p / rho);
                    // Cache the original rounded 10% numerators once per stage.
                    source_momentum_allowance[i] = 0.1 * (rho * (std::abs(v) + cs));
                    source_internal_allowance[i] = 0.1 * (p / (cell.gamma - 1.0));
                }
            }
            const double rho = source_rho[i];
            if (rho == 0.0)
            {
                for (double x : s.data)
                    if (x != 0.0) throw std::invalid_argument("LTS RK2 cannot create gas in vacuum");
                continue;
            }
            bound(source_density_allowance[i], s[hd_core::RHO]);
            bound(source_momentum_allowance[i], s[hd_core::MOMDEN(0)]);
            bound(source_internal_allowance[i], s[hd_core::INTERNAL_ENERGY]);
        }
        source_cache_ready = true;
    };
    // trial_H == 0 selects the exact initial limit. Stage checks return only
    // infinity/zero for pass/fail, but still scan every cell for invalid rates.
    auto change_limit = [&](const hd_core::MeshBlock& stage,
                            const std::vector<Conserved>& rates, double trial_H) {
        if (trial_H > 0.0)
        {
            const double twice_H = 2.0 * trial_H;
            bool exceeded = false;
            scan_source_changes(stage, rates, [&](double allowance, double rate) {
                exceeded |= detail::source_change_exceeds(trial_H, twice_H, allowance, rate);
            });
            return exceeded ? 0.0 : inf;
        }
        double limit = inf;
        scan_source_changes(stage, rates, [&](double allowance, double rate) {
            if (rate != 0.0) limit = std::min(limit, allowance / std::abs(rate));
        });
        return limit;
    };
    auto evaluate_sources = [&](const hd_core::MeshBlock& stage, double t, double trial_H) {
        performance::Scope timer(performance::Section::Sources);
        // Refill even when the MeshBlock address is unchanged: stage_view
        // replaces its states at each RK stage/substep, including retry stages.
        source_cache_ready = false;
        std::fill(source.begin(), source.end(), Conserved{});
        double limit = inf;
        for (const auto& term : sources)
        {
            if (!term.func) throw std::invalid_argument("Missing LTS RK2 source function");
            std::fill(contribution.begin(), contribution.end(), Conserved{});
            double cap = inf;
            term.func(stage, t, cfg, contribution, cap);
            hd_source::complete_internal_source(stage, term, contribution);
            if (contribution.size() != source.size()) throw std::logic_error("Source resized LTS RK2 storage");
            if (std::isnan(cap) || cap <= 0.0) throw std::invalid_argument("Invalid LTS RK2 source limit");
            limit = std::min({limit, cap, change_limit(stage, contribution, trial_H)});
            for (int i = 0; i < n; ++i) source[ng + i] += contribution[ng + i];
        }
        return sources.empty() ? inf : std::min(limit, change_limit(stage, source, trial_H));
    };
    solve(start, all_faces);
    initial_flux = flux; initial_speed = speed; initial_velocity = face_velocity;
    const auto initial_hints = hints;
    const double source_limit = evaluate_sources(start, time, 0.0);
    initial_source = source;
    double min_dt = inf;
    hydro::max_inverse_dt = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double rate = std::max(speed[i], speed[i + 1]) / mesh.cell_width(0, ng + i);
        if (!std::isfinite(rate)) throw std::runtime_error("Non-finite LTS RK2 CFL rate");
        hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
        local_dt[i] = rate > 0.0 ? cfl / rate : inf;
        min_dt = std::min(min_dt, local_dt[i]);
    }
    const int max_depth = cfg.lts_max_levels - 1, buffer = linear ? 2 : 1;
    double H = std::min({dt_limit, std::ldexp(min_dt, max_depth), source_limit});
    for (long long attempt = 0; attempt <= cfg.integrator_retry_max; ++attempt, H *= 0.5)
    {
        performance::count(performance::Counter::RK2Attempts);
        if (!(H > 0.0) || !std::isfinite(time + H) || time + H <= time)
            throw std::runtime_error("LTS RK2 macro step cannot advance time");
        std::fill(level.begin(), level.end(), 0);
        auto promote = [&](int i, int depth) {
            for (int offset = -buffer; offset <= buffer; ++offset)
            {
                int j = i + offset;
                if (periodic) j = (j % n + n) % n;
                if (j >= 0 && j < n) level[j] = std::max(level[j], depth);
            }
        };
        for (int i = 0; i < n; ++i)
        {
            int depth = 0;
            while (depth < max_depth && std::ldexp(H, -depth) > local_dt[i]) ++depth;
            promote(i, depth);
        }
        const int deepest = *std::max_element(level.begin(), level.end());
        for (int d = deepest; d > 1; --d)
            for (int i = 0; i < n; ++i) if (level[i] == d) promote(i, d - 1);
        std::vector<std::vector<int>> cells(deepest + 1), faces(deepest + 1), owned(deepest + 1);
        for (int i = 0; i < n; ++i) cells[level[i]].push_back(i);
        for (int f = 0; f < nf; ++f)
        {
            const int d = std::max(left(f) >= 0 ? level[left(f)] : 0,
                                   right(f) >= 0 ? level[right(f)] : 0);
            owned[d].push_back(f);
            for (int l = 0; l <= d; ++l) faces[l].push_back(f);
        }
        state.cells = start.cells;
        face_velocity = initial_velocity;
        hints = initial_hints; flux = initial_flux; speed = initial_speed; source = initial_source;
        std::fill(pending.begin(), pending.end(), Conserved{});
        auto check_waves = [&](const std::vector<int>& list) {
            for (int f : list)
                for (int i : {left(f), right(f)}) if (i >= 0)
                {
                    const double rate = speed[f] / mesh.cell_width(0, ng + i);
                    hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
                    if (!std::isfinite(rate) || std::ldexp(H, -level[i]) * rate > cfg.cfl_max) return false;
                }
            return true;
        };
        // Build a complete same-time stage view, including all PLM stencils,
        // user boundary donors and nonlocal source inputs. Accepted state stays intact.
        auto stage_view = [&](int depth, double fraction, bool endpoint) {
            view.cells = state.cells;
            for (int i = 0; i < n; ++i)
                if (level[i] < depth || (endpoint && level[i] == depth))
                {
                    auto& cell = view.cells[ng + i];
                    cell.U = base[i] + slope[i] * ((fraction - start_fraction[i]) * H);
                    if (!physical(cell.U, cell.gamma)) return false;
                }
            boundary::FillBoundary(view, time + fraction * H, hooks);
            return true;
        };
        auto advance = [&](auto&& self, int depth, double fraction, bool reuse) -> bool {
            const double span = std::ldexp(1.0, -depth), h = H * span;
            const double end = fraction + span;
            if (time + end * H <= time + fraction * H)
                throw std::runtime_error("LTS RK2 substep cannot advance time");
            if (!reuse)
            {
                if (!stage_view(depth, fraction, false)) return false;
                solve(view, faces[depth]);
                if (H > evaluate_sources(view, time + fraction * H, H)) return false;
            }
            if (!check_waves(faces[depth])) return false;
            for (int f : owned[depth]) first_flux[f] = flux[f];
            for (int i : cells[depth])
            {
                base[i] = state.cells[ng + i].U;
                start_fraction[i] = fraction;
                source0[i] = source[ng + i];
                // Local pressure work uses the same sampled face velocities and areas as Godunov.
                work0[i] = -(state.cells[ng + i].gamma - 1.0) * base[i].internal_energy()
                    * (area[i + 1] * face_velocity[i + 1] - area[i] * face_velocity[i]) / volume[i];
                slope[i] = (flux[i] * area[i] - flux[i + 1] * area[i + 1]) / volume[i] + source0[i];
                slope[i][hd_core::INTERNAL_ENERGY] += work0[i];
                if (!physical(base[i] + h * slope[i], state.cells[ng + i].gamma)) return false;
            }
            if (depth < deepest)
            {
                if (!self(self, depth + 1, fraction, true)) return false;
                if (!self(self, depth + 1, fraction + span / 2.0, false)) return false;
            }
            if (cells[depth].empty()) return true;
            if (!stage_view(depth, end, true)) return false;
            solve(view, owned[depth]);
            if (!check_waves(owned[depth]) || H > evaluate_sources(view, time + end * H, H)) return false;
            for (int f : owned[depth])
            {
                // Exactly one trapezoidal time integral per face and local step.
                // Fine/coarse interfaces deposit the same integral on both sides.
                const auto integral = (first_flux[f] + flux[f]) * (0.5 * h * area[f]);
                if (left(f) >= 0) pending[left(f)] -= integral;
                if (right(f) >= 0) pending[right(f)] += integral;
            }
            for (int i : cells[depth])
            {
                auto& cell = state.cells[ng + i];
                cell.U = base[i] + pending[i] / volume[i] + (source0[i] + source[ng + i]) * (0.5 * h);
                const double work1 = -(cell.gamma - 1.0) * view.cells[ng + i].U.internal_energy()
                    * (area[i + 1] * face_velocity[i + 1] - area[i] * face_velocity[i]) / volume[i];
                cell.U[hd_core::INTERNAL_ENERGY] += 0.5 * h * (work0[i] + work1);
                if (!physical(cell.U, cell.gamma)) return false;
                cell.synchronize_dual_energy();
                pending[i] = Conserved{};
            }
            return true;
        };
        if (!advance(advance, 0, 0.0, true))
        {
            performance::count(performance::Counter::RK2Rejected);
            continue;
        }
        boundary::FillBoundary(state, time + H, hooks);
        mesh.cells.swap(state.cells);
        if (star_pressure_cache) star_pressure_cache->swap(hints);
        return H;
    }
    throw std::runtime_error("LTS RK2 retry limit exceeded");
}
}
