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
    std::vector<double> speed(n + 1), initial_speed, local_dt(n), start_fraction(n);
    std::vector<Conserved> source(mesh.cells.size()), initial_source, contribution(mesh.cells.size());
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
            hydro::wave_timestep.face_speed = -1.0;
            flux[f] = performance::measure(performance::Section::Riemann, [&] {
                return solver.flux_with_pressure_hint(WL, WR, 0,
                    stage.cells[ng + f - 1].gamma, stage.cells[ng + f].gamma, hints[0][f]);
            }, true);
            speed[f] = hydro::wave_timestep.face_speed;
            if (!std::isfinite(speed[f]) || speed[f] < 0.0)
                throw std::runtime_error("LTS RK2 received an invalid wave speed");
            for (double x : flux[f].data)
                if (!std::isfinite(x)) throw std::runtime_error("LTS RK2 received a non-finite flux");
        }
        if (periodic)
        {
            flux[n] = flux[0]; speed[n] = speed[0]; hints[0][n] = hints[0][0];
        }
    };
    // Independent source-change caps prevent cancellation from hiding stiff terms.
    auto change_limit = [&](const hd_core::MeshBlock& stage, const std::vector<Conserved>& rates) {
        double limit = inf;
        auto bound = [&](double scale, double rate) {
            if (rate != 0.0) limit = std::min(limit, 0.1 * scale / std::abs(rate));
        };
        for (int i = 0; i < n; ++i)
        {
            const auto& s = rates[ng + i];
            for (double x : s.data)
                if (!std::isfinite(x)) throw std::runtime_error("Non-finite LTS RK2 source");
            const auto& cell = stage.cells[ng + i];
            const auto W = cell.U.to_Primtivies(cell.gamma);
            const double rho = W[hd_core::RHO], p = W[hd_core::PRE];
            if (rho == 0.0)
            {
                for (double x : s.data)
                    if (x != 0.0) throw std::invalid_argument("LTS RK2 cannot create gas in vacuum");
                continue;
            }
            const double v = W[hd_core::VEL(0)], cs = std::sqrt(cell.gamma * p / rho);
            bound(rho, s[hd_core::RHO]);
            bound(rho * (std::abs(v) + cs), s[hd_core::MOMDEN(0)]);
            bound(p / (cell.gamma - 1.0), s[hd_core::ENERGY] - v * s[hd_core::MOMDEN(0)]
                + 0.5 * v * v * s[hd_core::RHO]);
        }
        return limit;
    };
    auto evaluate_sources = [&](const hd_core::MeshBlock& stage, double t) {
        performance::Scope timer(performance::Section::Sources);
        std::fill(source.begin(), source.end(), Conserved{});
        double limit = inf;
        for (const auto& term : sources)
        {
            if (!term.func) throw std::invalid_argument("Missing LTS RK2 source function");
            std::fill(contribution.begin(), contribution.end(), Conserved{});
            double cap = inf;
            term.func(stage, t, cfg, contribution, cap);
            if (contribution.size() != source.size()) throw std::logic_error("Source resized LTS RK2 storage");
            if (std::isnan(cap) || cap <= 0.0) throw std::invalid_argument("Invalid LTS RK2 source limit");
            limit = std::min({limit, cap, change_limit(stage, contribution)});
            for (int i = 0; i < n; ++i) source[ng + i] += contribution[ng + i];
        }
        return sources.empty() ? inf : std::min(limit, change_limit(stage, source));
    };
    solve(start, all_faces);
    initial_flux = flux; initial_speed = speed;
    const auto initial_hints = hints;
    const double source_limit = evaluate_sources(start, time);
    initial_source = source;
    double min_dt = inf;
    hydro::wave_timestep.max_inverse_dt = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double rate = std::max(speed[i], speed[i + 1]) / mesh.cell_width(0, ng + i);
        if (!std::isfinite(rate)) throw std::runtime_error("Non-finite LTS RK2 CFL rate");
        hydro::wave_timestep.max_inverse_dt = std::max(hydro::wave_timestep.max_inverse_dt, rate);
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
        hints = initial_hints; flux = initial_flux; speed = initial_speed; source = initial_source;
        std::fill(pending.begin(), pending.end(), Conserved{});
        auto check_waves = [&](const std::vector<int>& list) {
            for (int f : list)
                for (int i : {left(f), right(f)}) if (i >= 0)
                {
                    const double rate = speed[f] / mesh.cell_width(0, ng + i);
                    hydro::wave_timestep.max_inverse_dt = std::max(hydro::wave_timestep.max_inverse_dt, rate);
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
                if (H > evaluate_sources(view, time + fraction * H)) return false;
            }
            if (!check_waves(faces[depth])) return false;
            for (int f : owned[depth]) first_flux[f] = flux[f];
            for (int i : cells[depth])
            {
                base[i] = state.cells[ng + i].U;
                start_fraction[i] = fraction;
                source0[i] = source[ng + i];
                slope[i] = (flux[i] * area[i] - flux[i + 1] * area[i + 1]) / volume[i] + source0[i];
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
            if (!check_waves(owned[depth]) || H > evaluate_sources(view, time + end * H)) return false;
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
                if (!physical(cell.U, cell.gamma)) return false;
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
