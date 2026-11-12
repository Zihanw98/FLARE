/**
 * @file lts_euler.cpp
 * @brief 守恒的多层一阶 Euler 局部时间步。
 *
 * 总体 pipeline：
 * 1. 校验配置，把正式网格复制到 workspace，并填充同步起点边界；
 * 2. 求全部独立面的初始通量和波速，计算每个单元的 CFL 局部步；
 * 3. 在同步起点计算一次源项，并在整个同步大步 H 内冻结；
 * 4. 由 CFL、源项和外部上限共同选择 H，再按 H/2^level 给单元分层；
 * 5. 每个面交给相邻单元中更细的层，并预先建立细层需要的慢单元预测表；
 * 6. advance_level() 递归执行两个半步，累计两侧共享的守恒面通量积分；
 * 7. 每层结束时结算本层单元并同步双能量；失败则 H 减半、恢复起点并重试；
 * 8. 所有层成功后才提交正式网格和星区压力缓存。
 *
 * 时间离散仍是一阶。PLM 只提高空间重构阶数，并用线性预测把慢层模板
 * 移到细层阶段时刻；源项（包括几何项和自引力）在整个 H 内保持起点值。
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

struct ValidatedOptions
{
    hydro::Reconstructor method;
    hydro::Limiter limiter;
};

ValidatedOptions validate_lts_euler_inputs(
    const MeshBlock& mesh, double time, double dt_limit,
    hydro::SpatialScheme spatial_scheme, const config::SimulationConfig& cfg,
    double cfl)
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

    const bool periodic = mesh.bc_lower[0] == boundary::BoundaryCondition::Periodic;
    if (periodic != (mesh.bc_upper[0] == boundary::BoundaryCondition::Periodic))
        throw std::invalid_argument("LTS requires paired periodic boundaries");
    if (periodic && mesh.geometry != hd_core::Geometry::Cartesian)
        throw std::invalid_argument("LTS periodic boundaries require Cartesian geometry");

    return {method, limiter};
}
}

// 跨同步大步复用的存储。这里不保存一次调用的控制流程。
struct LTSEulerWorkspace::Impl
{
    MeshBlock start, work, candidate;
    std::unique_ptr<hydro::ReconstructionCache> initial_reconstruction, sub_reconstruction;
    std::vector<double> local_dt, step_fraction, pressure_work;

    // 星区压力的持久缓存、同步起点快照和当前试算副本。
    hydro::StarPressureCache star_pressure_cache, initial_pressure, trial_pressure;
    std::vector<FaceResult> initial_faces, stage_faces;

    // 冻结源项、跨层面通量寄存器和局部 Euler 预测器。
    std::vector<Conserved> source, source_contribution, pending_flux_integral, step_start;
    std::vector<std::array<double, 3>> source_scales;

    // 每次重试重新生成的层级表、面表和预测模板。
    std::vector<int> level, face_level;
    std::vector<unsigned char> predicted, lower_ghosts, upper_ghosts;
    std::vector<std::vector<int>> cells_by_level, faces_needed_by_level, prediction_cells;

    void prepare(const MeshBlock& mesh)
    {
        const bool geometry_changed =
            start.nactive != mesh.nactive || start.nghost != mesh.nghost ||
            start.geometry != mesh.geometry || start.face != mesh.face ||
            start.volume_cell != mesh.volume_cell ||
            start.face_area[0].faces != mesh.face_area[0].faces;
        if (geometry_changed)
        {
            mesh.validate_geometry();
            start = mesh;
            work = mesh;
            candidate = mesh;
            star_pressure_cache.reset(0, mesh.nactive[0] + 1);
        }

        start.bc_lower = work.bc_lower = candidate.bc_lower = mesh.bc_lower;
        start.bc_upper = work.bc_upper = candidate.bc_upper = mesh.bc_upper;
        start.cells = mesh.cells;

        const auto cell_count = static_cast<std::size_t>(mesh.nactive[0]);
        local_dt.resize(cell_count);
        step_fraction.resize(cell_count);
        pressure_work.resize(cell_count);
        initial_faces.resize(cell_count + 1);
        stage_faces.resize(cell_count + 1);
        source.resize(mesh.cells.size());
        source_contribution.resize(mesh.cells.size());
        pending_flux_integral.resize(cell_count);
        step_start.resize(cell_count);
        level.resize(cell_count);
        face_level.resize(cell_count + 1);
        predicted.resize(cell_count);
    }
};

LTSEulerWorkspace::LTSEulerWorkspace() : impl(std::make_unique<Impl>()) {}
LTSEulerWorkspace::~LTSEulerWorkspace() = default;

namespace
{
/**
 * 一次同步 LTS Euler 大步的控制对象。
 *
 * 算法通过下面的具名成员直接访问状态；这些成员引用 workspace 中的持久存储，
 * 因而既保持与 LTS RK2 一致的阅读方式，也避免每个同步大步重复分配数组。
 */
class LTSEulerStep
{
public:
    LTSEulerStep(
        MeshBlock& mesh, double time, double dt_limit,
        const hydro::RiemannSolver& solver, const hd_source::SourceList& sources,
        const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks,
        double cfl, LTSEulerWorkspace::Impl& storage,
        hydro::Reconstructor method, hydro::Limiter limiter);

    double run();

private:
    // ---------- 网格拓扑与状态有效性 ----------
    int left_cell(int face) const;
    int right_cell(int face) const;
    static bool is_physical(const Conserved& state, double gamma);

    // ---------- 面重构、Riemann 求解与初始 CFL ----------
    hydro::ReconstructionCache& reset_reconstruction(
        std::unique_ptr<hydro::ReconstructionCache>& cache, const MeshBlock& stage);
    FaceResult solve_face(
        const MeshBlock& stage, hydro::ReconstructionCache& reconstruction,
        int face, double& pressure_hint) const;
    double solve_initial_faces_and_cfl();

    // ---------- 同步起点的冻结源项 ----------
    double source_change_limit(const std::vector<Conserved>& rates, bool& cache_ready);
    double evaluate_frozen_sources();

    // ---------- CFL 分层、面所有权与预测模板 ----------
    void promote_neighbors(int cell, int depth);
    void assign_time_levels();
    void assign_faces_to_levels();
    void mark_prediction_cell(int depth, int cell);
    void assign_prediction_cells();
    void reset_attempt();

    // ---------- 阶段预测与递归 Euler 子循环 ----------
    bool refresh_stage_faces(int depth, double fraction, double stage_time);
    bool check_waves_and_accumulate_flux(int depth, double local_step);
    bool build_level_predictors(int depth, double fraction, double local_step);
    bool commit_level(int depth, double local_step);
    bool advance_level(int depth, double fraction, bool reuse_start_flux);

    MeshBlock& mesh_;
    const hydro::RiemannSolver& solver_;
    const hd_source::SourceList& sources_;
    const config::SimulationConfig& cfg_;
    const hd_user::UserHooks& hooks_;

    // 同步起点、当前已接受状态和 Euler 预测终点。
    MeshBlock& start_;
    MeshBlock& state_;
    MeshBlock& candidate_;
    std::unique_ptr<hydro::ReconstructionCache>& initial_reconstruction_;
    std::unique_ptr<hydro::ReconstructionCache>& stage_reconstruction_;

    // 面量与星区压力：initial_* 是重试快照，其余是当前试算值。
    hydro::StarPressureCache& star_pressure_cache_;
    hydro::StarPressureCache& initial_pressure_;
    hydro::StarPressureCache& pressure_hints_;
    std::vector<FaceResult>& initial_faces_;
    std::vector<FaceResult>& stage_faces_;

    // CFL、源项、局部预测器和跨层守恒通量寄存器。
    std::vector<double>& local_dt_;
    std::vector<double>& step_fraction_;
    std::vector<double>& pressure_work_;
    std::vector<Conserved>& source_;
    std::vector<Conserved>& source_contribution_;
    std::vector<Conserved>& pending_flux_integral_;
    std::vector<Conserved>& step_start_;
    std::vector<std::array<double, 3>>& source_scales_;

    // LTS 层级、面所有权与慢层预测模板。
    std::vector<int>& level_;
    std::vector<int>& face_level_;
    std::vector<unsigned char>& predicted_;
    std::vector<unsigned char>& lower_ghosts_;
    std::vector<unsigned char>& upper_ghosts_;
    std::vector<std::vector<int>>& cells_by_level_;
    std::vector<std::vector<int>>& faces_needed_by_level_;
    std::vector<std::vector<int>>& prediction_cells_;

    double time_;
    double dt_limit_;
    double cfl_;
    hydro::Reconstructor method_;
    hydro::Limiter limiter_;
    bool linear_;
    int n_;
    int ng_;
    bool periodic_;
    int independent_face_count_;
    int max_depth_;
    int level_buffer_;

    // mesh 几何不随局部子步改变；指针从第一个 active 单元/面开始。
    const double* area_ = nullptr;
    const double* volume_ = nullptr;
    int deepest_ = 0;
    double macro_step_ = 0.0;
};

LTSEulerStep::LTSEulerStep(
    MeshBlock& mesh, double time, double dt_limit,
    const hydro::RiemannSolver& solver, const hd_source::SourceList& sources,
    const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks,
    double cfl, LTSEulerWorkspace::Impl& storage,
    hydro::Reconstructor method, hydro::Limiter limiter)
    : mesh_(mesh), solver_(solver), sources_(sources), cfg_(cfg), hooks_(hooks),
      start_(storage.start), state_(storage.work), candidate_(storage.candidate),
      initial_reconstruction_(storage.initial_reconstruction),
      stage_reconstruction_(storage.sub_reconstruction),
      star_pressure_cache_(storage.star_pressure_cache),
      initial_pressure_(storage.initial_pressure), pressure_hints_(storage.trial_pressure),
      initial_faces_(storage.initial_faces), stage_faces_(storage.stage_faces),
      local_dt_(storage.local_dt), step_fraction_(storage.step_fraction),
      pressure_work_(storage.pressure_work), source_(storage.source),
      source_contribution_(storage.source_contribution),
      pending_flux_integral_(storage.pending_flux_integral),
      step_start_(storage.step_start), source_scales_(storage.source_scales),
      level_(storage.level), face_level_(storage.face_level),
      predicted_(storage.predicted), lower_ghosts_(storage.lower_ghosts),
      upper_ghosts_(storage.upper_ghosts), cells_by_level_(storage.cells_by_level),
      faces_needed_by_level_(storage.faces_needed_by_level),
      prediction_cells_(storage.prediction_cells),
      time_(time), dt_limit_(dt_limit), cfl_(cfl),
      method_(method), limiter_(limiter), linear_(method != hydro::Reconstructor::PCM),
      n_(mesh.nactive[0]), ng_(mesh.nghost[0]),
      periodic_(mesh.bc_lower[0] == boundary::BoundaryCondition::Periodic),
      independent_face_count_(periodic_ ? n_ : n_ + 1),
      max_depth_(cfg.lts_max_levels - 1), level_buffer_(linear_ ? 2 : 1)
{
    storage.prepare(mesh_);
    boundary::FillBoundary(start_, time_, hooks_);
    area_ = start_.face_area[0].faces.data() + ng_;
    volume_ = start_.volume_cell.data() + ng_;
}

int LTSEulerStep::left_cell(int face) const
{
    return face == 0 ? (periodic_ ? n_ - 1 : -1) : face - 1;
}

int LTSEulerStep::right_cell(int face) const
{
    return face == n_ ? (periodic_ ? 0 : -1) : face;
}

bool LTSEulerStep::is_physical(const Conserved& state, double gamma)
{
    try { state.validate(gamma); }
    catch (const std::invalid_argument&) { return false; }
    return true;
}

hydro::ReconstructionCache& LTSEulerStep::reset_reconstruction(
    std::unique_ptr<hydro::ReconstructionCache>& cache, const MeshBlock& stage)
{
    if (!cache)
        cache = std::make_unique<hydro::ReconstructionCache>(
            stage, method_, limiter_, cfg_.limiter_beta);
    else
        cache->reset(stage, method_, limiter_, cfg_.limiter_beta);
    return *cache;
}

FaceResult LTSEulerStep::solve_face(
    const MeshBlock& stage, hydro::ReconstructionCache& reconstruction,
    int face, double& pressure_hint) const
{
    const int right = ng_ + face;
    hd_core::Index right_index{};
    right_index[0] = right;
    const auto [left_state, right_state] = performance::measure(
        performance::Section::Reconstruction,
        [&] { return reconstruction.face(0, right_index); }, true);

    FaceResult result;
    result.speed = -1.0;
    result.flux = performance::measure(performance::Section::Riemann, [&] {
        return solver_.flux_with_pressure_hint(
            left_state, right_state, 0,
            stage.cells[right - 1].gamma, stage.cells[right].gamma,
            pressure_hint, result.speed, result.velocity);
    }, true);

    if (!std::isfinite(result.speed) || result.speed < 0.0)
        throw std::runtime_error("LTS received an invalid Riemann wave speed");
    for (double value : result.flux.data)
        if (!std::isfinite(value))
            throw std::runtime_error("LTS received a non-finite face flux");
    return result;
}

double LTSEulerStep::solve_initial_faces_and_cfl()
{
    auto& reconstruction = reset_reconstruction(
        initial_reconstruction_, start_);
    initial_pressure_ = star_pressure_cache_;
    for (int face = 0; face <= n_; ++face)
    {
        if (periodic_ && face == n_)
        {
            initial_faces_[face] = initial_faces_[0];
            initial_pressure_[0][face] = initial_pressure_[0][0];
        }
        else
        {
            initial_faces_[face] = solve_face(
                start_, reconstruction, face, initial_pressure_[0][face]);
        }
    }

    double minimum_local_dt = infinity;
    hydro::max_inverse_dt = 0.0;
    for (int cell = 0; cell < n_; ++cell)
    {
        const double rate =
            std::max(initial_faces_[cell].speed,
                     initial_faces_[cell + 1].speed)
            / start_.cell_width(0, ng_ + cell);
        if (!std::isfinite(rate))
            throw std::runtime_error("Non-finite LTS CFL rate");
        hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
        local_dt_[cell] = rate > 0.0 ? cfl_ / rate : infinity;
        minimum_local_dt = std::min(minimum_local_dt, local_dt_[cell]);
    }
    return minimum_local_dt;
}

double LTSEulerStep::source_change_limit(
    const std::vector<Conserved>& rates, bool& cache_ready)
{
    double limit = infinity;
    auto bound = [&](double scale, double rate) {
        if (rate != 0.0)
            limit = std::min(limit, 0.1 * scale / std::abs(rate));
    };

    for (int cell_index = 0; cell_index < n_; ++cell_index)
    {
        const int storage_index = ng_ + cell_index;
        const auto& rate = rates[storage_index];
        for (double value : rate.data)
            if (!std::isfinite(value))
                throw std::runtime_error("LTS received a non-finite source rate");

        std::array<double, 3> uncached{};
        auto& scale = source_scales_.empty()
            ? uncached
            : source_scales_[cell_index];
        if (!cache_ready || source_scales_.empty())
        {
            const auto& cell = start_.cells[storage_index];
            const auto primitive = cell.U.to_Primtivies(cell.gamma);
            const double rho = primitive[hd_core::RHO];
            scale = {rho, 0.0, 0.0};
            if (rho != 0.0)
            {
                const double velocity = primitive[hd_core::VEL(0)];
                const double sound_speed =
                    std::sqrt(cell.gamma * primitive[hd_core::PRE] / rho);
                scale = {rho, rho * (std::abs(velocity) + sound_speed),
                         primitive[hd_core::PRE] / (cell.gamma - 1.0)};
            }
        }

        const auto [rho, momentum, internal_energy] = scale;
        if (rho == 0.0)
        {
            for (double value : rate.data)
                if (value != 0.0)
                    throw std::invalid_argument(
                        "LTS does not support sources creating gas in vacuum");
            continue;
        }
        bound(rho, rate[hd_core::RHO]);
        bound(momentum, rate[hd_core::MOMDEN(0)]);
        bound(internal_energy, rate[hd_core::INTERNAL_ENERGY]);
    }
    cache_ready = true;
    return limit;
}

double LTSEulerStep::evaluate_frozen_sources()
{
    std::fill(source_.begin(), source_.end(), Conserved{});
    source_scales_.resize(sources_.empty() ? 0 : n_);
    bool cache_ready = false;
    double source_limit = infinity;

    for (const auto& term : sources_)
    {
        if (!term.func) throw std::invalid_argument("Missing LTS source function");
        std::fill(source_contribution_.begin(),
                  source_contribution_.end(), Conserved{});
        double callback_limit = infinity;
        term.func(start_, time_, cfg_, source_contribution_, callback_limit);
        hd_source::complete_internal_source(
            start_, term, source_contribution_);
        if (source_contribution_.size() != source_.size())
            throw std::logic_error("Source changed LTS RHS storage size");
        if (std::isnan(callback_limit) || callback_limit <= 0.0)
            throw std::invalid_argument("Invalid LTS source timestep limit");

        source_limit = std::min({
            source_limit, callback_limit,
            source_change_limit(source_contribution_, cache_ready)});
        for (int cell = 0; cell < n_; ++cell)
            source_[ng_ + cell] += source_contribution_[ng_ + cell];
    }
    return std::min(
        source_limit, source_change_limit(source_, cache_ready));
}

void LTSEulerStep::promote_neighbors(int cell, int depth)
{
    for (int offset = -level_buffer_; offset <= level_buffer_; ++offset)
    {
        int neighbor = cell + offset;
        if (periodic_) neighbor = (neighbor % n_ + n_) % n_;
        if (neighbor >= 0 && neighbor < n_)
            level_[neighbor] = std::max(level_[neighbor], depth);
    }
}

void LTSEulerStep::assign_time_levels()
{
    std::fill(level_.begin(), level_.end(), 0);
    for (int cell = 0; cell < n_; ++cell)
    {
        int depth = 0;
        while (depth < max_depth_ &&
               std::ldexp(macro_step_, -depth) > local_dt_[cell])
            ++depth;
        promote_neighbors(cell, depth);
    }

    deepest_ = *std::max_element(level_.begin(), level_.end());
    for (int depth = deepest_; depth > 1; --depth)
        for (int cell = 0; cell < n_; ++cell)
            if (level_[cell] == depth)
                promote_neighbors(cell, depth - 1);
}

void LTSEulerStep::assign_faces_to_levels()
{
    cells_by_level_.resize(deepest_ + 1);
    faces_needed_by_level_.resize(deepest_ + 1);
    prediction_cells_.resize(deepest_ + 1);
    lower_ghosts_.assign(deepest_ + 1, false);
    upper_ghosts_.assign(deepest_ + 1, false);
    for (int depth = 0; depth <= deepest_; ++depth)
    {
        cells_by_level_[depth].clear();
        faces_needed_by_level_[depth].clear();
        prediction_cells_[depth].clear();
    }

    for (int cell = 0; cell < n_; ++cell)
        cells_by_level_[level_[cell]].push_back(cell);

    for (int face = 0; face < independent_face_count_; ++face)
    {
        const int left = left_cell(face);
        const int right = right_cell(face);
        face_level_[face] = std::max(
            left >= 0 ? level_[left] : 0,
            right >= 0 ? level_[right] : 0);
        for (int depth = 0; depth <= face_level_[face]; ++depth)
            faces_needed_by_level_[depth].push_back(face);
    }
    if (periodic_) face_level_[n_] = face_level_[0];
}

void LTSEulerStep::mark_prediction_cell(int depth, int cell)
{
    if (level_[cell] < depth && !predicted_[cell])
    {
        predicted_[cell] = true;
        prediction_cells_[depth].push_back(cell);
    }
}

void LTSEulerStep::assign_prediction_cells()
{
    for (int depth = 1; depth <= deepest_; ++depth)
    {
        std::fill(predicted_.begin(), predicted_.end(), false);
        bool global_boundary = false;
        for (int face : faces_needed_by_level_[depth])
        {
            const int stencil_left = face - (linear_ ? 2 : 1);
            const int stencil_right = face + (linear_ ? 1 : 0);
            for (int cell = stencil_left; cell <= stencil_right; ++cell)
            {
                if (cell >= 0 && cell < n_)
                {
                    mark_prediction_cell(depth, cell);
                    continue;
                }

                const bool lower = cell < 0;
                lower_ghosts_[depth] |= lower;
                upper_ghosts_[depth] |= !lower;
                const auto condition = lower
                    ? start_.bc_lower[0]
                    : start_.bc_upper[0];
                if (condition == boundary::BoundaryCondition::Periodic)
                    mark_prediction_cell(depth, (cell % n_ + n_) % n_);
                else if (condition == boundary::BoundaryCondition::Reflective)
                    mark_prediction_cell(depth, lower ? -cell - 1 : 2 * n_ - cell - 1);
                else if (condition == boundary::BoundaryCondition::User)
                    global_boundary = true;
            }
        }
        if (global_boundary)
            for (int cell = 0; cell < n_; ++cell)
                mark_prediction_cell(depth, cell);
    }
}

void LTSEulerStep::reset_attempt()
{
    assign_time_levels();
    assign_faces_to_levels();
    assign_prediction_cells();
    state_.cells = start_.cells;
    candidate_.cells = start_.cells;
    pressure_hints_ = initial_pressure_;
    stage_faces_ = initial_faces_;
    std::fill(pending_flux_integral_.begin(),
              pending_flux_integral_.end(), Conserved{});
}

bool LTSEulerStep::refresh_stage_faces(
    int depth, double fraction, double stage_time)
{
    for (int cell_index : prediction_cells_[depth])
    {
        auto& cell = state_.cells[ng_ + cell_index];
        const double theta = linear_
            ? std::ldexp(
                fraction - step_fraction_[cell_index],
                level_[cell_index])
            : 0.0;
        cell.U = (1.0 - theta) * step_start_[cell_index]
            + theta * candidate_.cells[ng_ + cell_index].U;
        if (!is_physical(cell.U, cell.gamma)) return false;
    }

    if (lower_ghosts_[depth])
        boundary::ApplyBoundaryCondition(
            state_, 0, boundary::BoundarySide::Lower, stage_time, hooks_);
    if (upper_ghosts_[depth])
        boundary::ApplyBoundaryCondition(
            state_, 0, boundary::BoundarySide::Upper, stage_time, hooks_);

    auto& reconstruction = reset_reconstruction(
        stage_reconstruction_, state_);
    for (int face : faces_needed_by_level_[depth])
        stage_faces_[face] = solve_face(
            state_, reconstruction, face, pressure_hints_[0][face]);
    if (periodic_)
    {
        stage_faces_[n_] = stage_faces_[0];
        pressure_hints_[0][n_] = pressure_hints_[0][0];
    }
    return true;
}

bool LTSEulerStep::check_waves_and_accumulate_flux(int depth, double local_step)
{
    for (int face : faces_needed_by_level_[depth])
    {
        for (int cell : {left_cell(face), right_cell(face)})
        {
            if (cell < 0) continue;
            const double rate =
                stage_faces_[face].speed
                / start_.cell_width(0, ng_ + cell);
            hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
            const double cell_step = std::ldexp(macro_step_, -level_[cell]);
            if (!std::isfinite(rate) || cell_step * rate > cfg_.cfl_max) return false;
        }

        if (face_level_[face] == depth)
        {
            const auto integral =
                stage_faces_[face].flux * (local_step * area_[face]);
            const int left = left_cell(face);
            const int right = right_cell(face);
            if (left >= 0) pending_flux_integral_[left] -= integral;
            if (right >= 0) pending_flux_integral_[right] += integral;
        }
    }
    return true;
}

bool LTSEulerStep::build_level_predictors(
    int depth, double fraction, double local_step)
{
    for (int cell_index : cells_by_level_[depth])
    {
        const auto& current_cell = state_.cells[ng_ + cell_index];
        step_start_[cell_index] = current_cell.U;
        step_fraction_[cell_index] = fraction;
        pressure_work_[cell_index] =
            -(current_cell.gamma - 1.0)
            * step_start_[cell_index].internal_energy()
            * (area_[cell_index + 1] * stage_faces_[cell_index + 1].velocity
               - area_[cell_index] * stage_faces_[cell_index].velocity)
            / volume_[cell_index];

        auto& candidate = candidate_.cells[ng_ + cell_index];
        candidate.U = step_start_[cell_index]
            + (stage_faces_[cell_index].flux * area_[cell_index]
               - stage_faces_[cell_index + 1].flux * area_[cell_index + 1])
                * (local_step / volume_[cell_index])
            + source_[ng_ + cell_index] * local_step;
        candidate.U[hd_core::INTERNAL_ENERGY] +=
            local_step * pressure_work_[cell_index];
        if (!is_physical(candidate.U, candidate.gamma)) return false;
    }
    return true;
}

bool LTSEulerStep::commit_level(int depth, double local_step)
{
    for (int cell_index : cells_by_level_[depth])
    {
        auto& cell = state_.cells[ng_ + cell_index];
        cell.U = step_start_[cell_index]
            + pending_flux_integral_[cell_index] / volume_[cell_index]
            + source_[ng_ + cell_index] * local_step;
        cell.U[hd_core::INTERNAL_ENERGY] +=
            local_step * pressure_work_[cell_index];
        if (!is_physical(cell.U, cell.gamma)) return false;
        cell.synchronize_dual_energy();
        pending_flux_integral_[cell_index] = Conserved{};
    }
    return true;
}

bool LTSEulerStep::advance_level(int depth, double fraction, bool reuse_start_flux)
{
    const double local_step = std::ldexp(macro_step_, -depth);
    const double span = std::ldexp(1.0, -depth);
    const double stage_time = time_ + fraction * macro_step_;
    if (time_ + (fraction + span) * macro_step_ <= stage_time)
        throw std::runtime_error("LTS substep cannot advance time");

    if (!reuse_start_flux && !refresh_stage_faces(depth, fraction, stage_time))
        return false;
    if (!check_waves_and_accumulate_flux(depth, local_step)) return false;
    if (!build_level_predictors(depth, fraction, local_step)) return false;

    if (depth < deepest_)
    {
        if (!advance_level(depth + 1, fraction, true)) return false;
        if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
    }
    return commit_level(depth, local_step);
}

double LTSEulerStep::run()
{
    // Pipeline 1：同步起点的面量、CFL 局部步和冻结源项。
    const double minimum_local_dt = solve_initial_faces_and_cfl();
    const double source_limit = evaluate_frozen_sources();

    // Pipeline 2：选择所有限制共同允许的最大同步大步 H。
    macro_step_ = std::min({
        dt_limit_, std::ldexp(minimum_local_dt, max_depth_), source_limit});
    for (long long attempt = 0; attempt <= cfg_.integrator_retry_max;
         ++attempt, macro_step_ *= 0.5)
    {
        if (!std::isfinite(macro_step_) || macro_step_ <= 0.0 ||
            !std::isfinite(time_ + macro_step_) || time_ + macro_step_ <= time_)
            throw std::runtime_error("LTS macro step cannot advance time");

        // Pipeline 3：每次重试重新分层，并恢复同步起点的全部试算状态。
        reset_attempt();
        if (!advance_level(0, 0.0, true)) continue;

        // Pipeline 4：递归全部成功后才提交正式网格和压力历史。
        boundary::FillBoundary(state_, time_ + macro_step_, hooks_);
        mesh_.cells.swap(state_.cells);
        star_pressure_cache_.swap(pressure_hints_);
        return macro_step_;
    }
    throw std::runtime_error("LTS Euler retry limit exceeded");
}
}

double advance_lts_euler(
    MeshBlock& mesh, double time, double dt_limit,
    hydro::SpatialScheme spatial_scheme, const hydro::RiemannSolver& solver,
    const hd_source::SourceList& source_list, const config::SimulationConfig& cfg,
    const hd_user::UserHooks& hooks, double cfl, LTSEulerWorkspace* workspace)
{
    const auto options = validate_lts_euler_inputs(
        mesh, time, dt_limit, spatial_scheme, cfg, cfl);

    std::optional<LTSEulerWorkspace> temporary_workspace;
    if (!workspace)
    {
        temporary_workspace.emplace();
        workspace = &*temporary_workspace;
    }

    LTSEulerStep step(
        mesh, time, dt_limit, solver, source_list, cfg, hooks, cfl,
        *workspace->impl, options.method, options.limiter);
    return step.run();
}
}
