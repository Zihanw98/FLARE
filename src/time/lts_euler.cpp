/**
 * @file lts_euler.cpp
 * @brief 守恒的多层一阶 Euler 局部时间步。
 *
 * 总体 pipeline：
 * 1. 校验配置，建立同步起点和本次推进状态，并填充同步起点边界；
 * 2. 求全部独立面的初始通量和波速，计算每个单元的 CFL 局部步；
 * 3. 在同步起点计算一次源项，并在整个同步大步 H 内冻结；
 * 4. 由 CFL、源项和外部上限共同选择 H，再按 H/2^level 给单元分层；
 * 5. 每个面交给相邻单元中更细的层，通过同一时刻的阶段视图重构；
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
    if (spatial_scheme != hydro::SpatialScheme::Godunov)
        throw std::invalid_argument("lts_euler requires Godunov");
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

    mesh.validate_geometry();
    (void)mesh.flat_index(hd_core::Index{}, true);
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        if (mesh.nghost[d] < (linear ? 2 : 1))
            throw std::invalid_argument(
                "LTS requires one ghost cell for PCM, two for PLM variants on every axis");
        const bool periodic = mesh.bc_lower[d] == boundary::BoundaryCondition::Periodic;
        if (periodic != (mesh.bc_upper[d] == boundary::BoundaryCondition::Periodic))
            throw std::invalid_argument("LTS requires paired periodic boundaries on every axis");
        if (periodic && mesh.geometry != hd_core::Geometry::Cartesian)
            throw std::invalid_argument("LTS periodic boundaries require Cartesian geometry");
    }
    const auto max_count = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (mesh.active_count > max_count || mesh.face_info.size() > max_count)
        throw std::length_error("LTS cell or face count exceeds int range");

    return {method, limiter};
}

MeshBlock make_start_stage(
    const MeshBlock& mesh, double time, const hd_user::UserHooks& hooks)
{
    auto start = mesh;
    boundary::FillBoundary(start, time, hooks);
    return start;
}

// 一次同步 Euler 大步的状态、缓存和递归推进全部由此类持有。
class LTSEulerStep
{
public:
    LTSEulerStep(
        MeshBlock& mesh, double time, double dt_limit,
        const hydro::RiemannSolver& solver, const hd_source::SourceList& sources,
        const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks,
        double cfl, hydro::StarPressureCache* star_pressure_cache,
        hydro::Reconstructor method, hydro::Limiter limiter);

    double run();

private:
    // ---------- 网格拓扑与状态有效性 ----------
    std::size_t storage_index(int cell) const { return mesh_.active_to_storage[cell]; }
    double face_area(int face) const;
    int left_cell(int face) const;
    int right_cell(int face) const;
    static bool is_physical(const Conserved& state, double gamma);

    // ---------- 面重构、Riemann 求解与多维 CFL ----------
    FaceResult solve_face(const MeshBlock& stage, int face);
    double cell_wave_rate(int cell) const;
    double solve_initial_faces_and_cfl();

    // ---------- 同步起点的冻结源项 ----------
    double source_change_limit(const std::vector<Conserved>& rates, bool& cache_ready);
    double evaluate_frozen_sources();

    // ---------- CFL 分层与面所有权 ----------
    void promote_neighbors(int cell, int depth);
    void assign_time_levels();
    void assign_faces_to_levels();
    void reset_attempt();

    // ---------- 同时刻阶段视图与递归 Euler 子循环 ----------
    bool build_stage_view(int depth, double fraction, double stage_time);
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
    double time_, dt_limit_, cfl_;
    hydro::Reconstructor method_;
    hydro::Limiter limiter_;
    bool linear_;
    int n_, independent_face_count_, max_depth_, level_buffer_;
    int deepest_ = 0;
    double macro_step_ = 0.0;

    // start_：同步起点；state_：各 cell 已接受的局部状态；view_：同一时刻的重构视图。
    MeshBlock start_, state_, view_;
    hydro::ReconstructionCache reconstruction_;

    // 只有可选的星区压力猜测跨大步传递；所有试算存储都属于本类。
    hydro::StarPressureCache* star_pressure_cache_;
    hydro::StarPressureCache initial_pressure_, pressure_hints_;
    std::vector<FaceResult> initial_faces_, stage_faces_;
    std::vector<double> max_seen_speed_;
    std::vector<double> local_dt_, step_fraction_, pressure_work_;
    std::vector<Conserved> source_, source_contribution_, pending_flux_integral_;
    std::vector<Conserved> step_start_, candidate_;
    // rho、各方向 momentum、internal energy 的源项变化尺度。
    std::vector<std::array<double, hd_core::NDIM + 2>> source_scales_;
    std::vector<int> level_, face_level_;
    std::vector<std::vector<int>> cells_by_level_, faces_needed_by_level_;
};

LTSEulerStep::LTSEulerStep(
    MeshBlock& mesh, double time, double dt_limit,
    const hydro::RiemannSolver& solver, const hd_source::SourceList& sources,
    const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks,
    double cfl, hydro::StarPressureCache* star_pressure_cache,
    hydro::Reconstructor method, hydro::Limiter limiter)
    : mesh_(mesh), solver_(solver), sources_(sources), cfg_(cfg), hooks_(hooks),
      time_(time), dt_limit_(dt_limit), cfl_(cfl), method_(method), limiter_(limiter),
      linear_(method != hydro::Reconstructor::PCM),
      n_(static_cast<int>(mesh.active_count)),
      independent_face_count_(static_cast<int>(mesh.face_info.size())),
      max_depth_(cfg.lts_max_levels - 1), level_buffer_(linear_ ? 2 : 1),
      start_(make_start_stage(mesh, time, hooks)), state_(start_), view_(start_),
      reconstruction_(start_, method, limiter, cfg.limiter_beta),
      star_pressure_cache_(star_pressure_cache),
      initial_faces_(independent_face_count_), stage_faces_(independent_face_count_),
      max_seen_speed_(independent_face_count_), local_dt_(n_), step_fraction_(n_),
      pressure_work_(n_), source_(mesh.cells.size()), source_contribution_(mesh.cells.size()),
      pending_flux_integral_(n_), step_start_(n_), candidate_(n_),
      source_scales_(sources.empty() ? 0 : n_), level_(n_), face_level_(independent_face_count_)
{
    if (star_pressure_cache_) pressure_hints_ = *star_pressure_cache_;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        std::size_t count = 1;
        for (int axis = 0; axis < hd_core::NDIM; ++axis)
            count *= static_cast<std::size_t>(mesh.nactive[axis]) + (axis == d ? 1 : 0);
        pressure_hints_.ensure_size(d, count);
    }
}

double LTSEulerStep::face_area(int face) const
{
    const auto& info = mesh_.face_info[face];
    return mesh_.face_area[info.direction].at(info.storage_coordinates);
}

int LTSEulerStep::left_cell(int face) const
{
    const auto cell = mesh_.face_info[face].left_active;
    return cell == MeshBlock::invalid_active_index ? -1 : static_cast<int>(cell);
}

int LTSEulerStep::right_cell(int face) const
{
    const auto cell = mesh_.face_info[face].right_active;
    return cell == MeshBlock::invalid_active_index ? -1 : static_cast<int>(cell);
}

bool LTSEulerStep::is_physical(const Conserved& state, double gamma)
{
    try { state.validate(gamma); }
    catch (const std::invalid_argument&) { return false; }
    return true;
}

FaceResult LTSEulerStep::solve_face(const MeshBlock& stage, int face)
{
    const auto& info = mesh_.face_info[face];
    const int d = info.direction;
    const auto& right_index = info.storage_coordinates;
    auto left_index = right_index;
    --left_index[d];
    const auto [left_state, right_state] = performance::measure(
        performance::Section::Reconstruction,
        [&] { return reconstruction_.face(d, right_index); }, true);

    // 压力猜测按各方向的 active 面场展平；通量则按全局独立 face_id 存储。
    std::size_t hint_index = 0, stride = 1, normal_stride = 1;
    for (int axis = 0; axis < hd_core::NDIM; ++axis)
    {
        if (axis == d) normal_stride = stride;
        hint_index += static_cast<std::size_t>(right_index[axis] - stage.nghost[axis]) * stride;
        stride *= static_cast<std::size_t>(stage.nactive[axis]) + (axis == d ? 1 : 0);
    }
    auto& pressure_hint = pressure_hints_[d][hint_index];
    FaceResult result;
    result.speed = -1.0;
    result.flux = performance::measure(performance::Section::Riemann, [&] {
        return solver_.flux_with_pressure_hint(
            left_state, right_state, d, stage.at(left_index).gamma, stage.at(right_index).gamma,
            pressure_hint, result.speed, result.velocity);
    }, true);

    if (!std::isfinite(result.speed) || result.speed < 0.0)
        throw std::runtime_error("LTS received an invalid Riemann wave speed");
    for (double value : result.flux.data)
        if (!std::isfinite(value))
            throw std::runtime_error("LTS received a non-finite face flux");
    if (stage.bc_lower[d] == boundary::BoundaryCondition::Periodic &&
        right_index[d] == stage.nghost[d])
        pressure_hints_[d][hint_index + stage.nactive[d] * normal_stride] = pressure_hint;
    return result;
}

double LTSEulerStep::cell_wave_rate(int cell) const
{
    const auto index = mesh_.coordinates(storage_index(cell));
    double rate = 0.0;
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        const auto& pair = mesh_.cell_faces[cell][d];
        rate += std::max(max_seen_speed_[pair.lower], max_seen_speed_[pair.upper])
              / mesh_.cell_width(d, index[d]);
    }
    return rate;
}

double LTSEulerStep::solve_initial_faces_and_cfl()
{
    for (int face = 0; face < independent_face_count_; ++face)
    {
        initial_faces_[face] = solve_face(start_, face);
        max_seen_speed_[face] = initial_faces_[face].speed;
    }
    initial_pressure_ = pressure_hints_;

    double minimum_local_dt = infinity;
    hydro::max_inverse_dt = 0.0;
    for (int cell = 0; cell < n_; ++cell)
    {
        const double rate = cell_wave_rate(cell);
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
        const auto flat = storage_index(cell_index);
        const auto& rate = rates[flat];
        for (double value : rate.data)
            if (!std::isfinite(value))
                throw std::runtime_error("LTS received a non-finite source rate");

        auto& scale = source_scales_[cell_index];
        if (!cache_ready)
        {
            const auto& cell = start_.cells[flat];
            const auto primitive = cell.U.to_Primtivies(cell.gamma);
            const double rho = primitive[hd_core::RHO];
            scale.fill(0.0);
            scale[0] = rho;
            if (rho != 0.0)
            {
                const double sound_speed =
                    std::sqrt(cell.gamma * primitive[hd_core::PRE] / rho);
                for (int d = 0; d < hd_core::NDIM; ++d)
                    scale[1 + d] = rho * (std::abs(primitive[hd_core::VEL(d)]) + sound_speed);
                scale[1 + hd_core::NDIM] = primitive[hd_core::PRE] / (cell.gamma - 1.0);
            }
        }

        if (scale[0] == 0.0)
        {
            for (double value : rate.data)
                if (value != 0.0)
                    throw std::invalid_argument("LTS does not support sources creating gas in vacuum");
            continue;
        }
        bound(scale[0], rate[hd_core::RHO]);
        for (int d = 0; d < hd_core::NDIM; ++d)
            bound(scale[1 + d], rate[hd_core::MOMDEN(d)]);
        bound(scale[1 + hd_core::NDIM], rate[hd_core::INTERNAL_ENERGY]);
    }
    cache_ready = true;
    return limit;
}

double LTSEulerStep::evaluate_frozen_sources()
{
    std::fill(source_.begin(), source_.end(), Conserved{});
    if (sources_.empty()) return infinity;
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
            source_[storage_index(cell)] += source_contribution_[storage_index(cell)];
    }
    return std::min(
        source_limit, source_change_limit(source_, cache_ready));
}

void LTSEulerStep::promote_neighbors(int cell, int depth)
{
    level_[cell] = std::max(level_[cell], depth);
    for (int d = 0; d < hd_core::NDIM; ++d)
        for (int side : {-1, 1})
        {
            int neighbor = cell;
            for (int distance = 0; distance < level_buffer_; ++distance)
            {
                const auto& pair = mesh_.cell_faces[neighbor][d];
                neighbor = side < 0 ? left_cell(static_cast<int>(pair.lower))
                                    : right_cell(static_cast<int>(pair.upper));
                if (neighbor < 0) break;
                level_[neighbor] = std::max(level_[neighbor], depth);
            }
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
    for (int depth = 0; depth <= deepest_; ++depth)
    {
        cells_by_level_[depth].clear();
        faces_needed_by_level_[depth].clear();
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
}

void LTSEulerStep::reset_attempt()
{
    assign_time_levels();
    assign_faces_to_levels();
    state_.cells = start_.cells;
    pressure_hints_ = initial_pressure_;
    stage_faces_ = initial_faces_;
    for (int face = 0; face < independent_face_count_; ++face)
        max_seen_speed_[face] = initial_faces_[face].speed;
    std::fill(pending_flux_integral_.begin(), pending_flux_integral_.end(), Conserved{});
}

bool LTSEulerStep::build_stage_view(int depth, double fraction, double stage_time)
{
    view_.cells = state_.cells;
    for (int cell_index = 0; cell_index < n_; ++cell_index)
    {
        if (level_[cell_index] >= depth) continue;
        auto& cell = view_.cells[storage_index(cell_index)];
        // PLM 使用慢层 Euler 预测；PCM 保持慢层局部步起点值。
        const double theta = linear_
            ? std::ldexp(fraction - step_fraction_[cell_index], level_[cell_index]) : 0.0;
        cell.U = (1.0 - theta) * step_start_[cell_index] + theta * candidate_[cell_index];
        if (!is_physical(cell.U, cell.gamma)) return false;
    }
    boundary::FillBoundary(view_, stage_time, hooks_);
    return true;
}

bool LTSEulerStep::refresh_stage_faces(int depth, double fraction, double stage_time)
{
    if (!build_stage_view(depth, fraction, stage_time)) return false;
    reconstruction_.reset(view_, method_, limiter_, cfg_.limiter_beta);
    for (int face : faces_needed_by_level_[depth])
        stage_faces_[face] = solve_face(view_, face);
    return true;
}

bool LTSEulerStep::check_waves_and_accumulate_flux(int depth, double local_step)
{
    // 先更新整批面，再用各方向历史峰值的总和检查 CFL；重试时清除本次峰值。
    for (int face : faces_needed_by_level_[depth])
        max_seen_speed_[face] = std::max(max_seen_speed_[face], stage_faces_[face].speed);
    for (int face : faces_needed_by_level_[depth])
    {
        for (int cell : {left_cell(face), right_cell(face)})
        {
            if (cell < 0) continue;
            const double rate = cell_wave_rate(cell);
            hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
            const double cell_step = std::ldexp(macro_step_, -level_[cell]);
            if (!std::isfinite(rate) || cell_step * rate > cfg_.cfl_max) return false;
        }
        if (face_level_[face] == depth)
        {
            const auto integral = stage_faces_[face].flux * (local_step * face_area(face));
            const int left = left_cell(face), right = right_cell(face);
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
        const auto flat = storage_index(cell_index);
        const auto& current_cell = state_.cells[flat];
        const double volume = mesh_.volume_cell[flat];
        step_start_[cell_index] = current_cell.U;
        step_fraction_[cell_index] = fraction;
        Conserved flux_balance{};
        double div_velocity_balance = 0.0;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const auto& pair = mesh_.cell_faces[cell_index][d];
            const int lower = static_cast<int>(pair.lower), upper = static_cast<int>(pair.upper);
            const double lower_area = face_area(lower), upper_area = face_area(upper);
            flux_balance += stage_faces_[lower].flux * lower_area - stage_faces_[upper].flux * upper_area;
            div_velocity_balance += upper_area * stage_faces_[upper].velocity
                                  - lower_area * stage_faces_[lower].velocity;
        }
        pressure_work_[cell_index] = -(current_cell.gamma - 1.0)
            * step_start_[cell_index].internal_energy() * div_velocity_balance / volume;
        auto& candidate = candidate_[cell_index];
        candidate = step_start_[cell_index] + flux_balance * (local_step / volume)
                  + source_[flat] * local_step;
        candidate[hd_core::INTERNAL_ENERGY] += local_step * pressure_work_[cell_index];
        if (!is_physical(candidate, current_cell.gamma)) return false;
    }
    return true;
}

bool LTSEulerStep::commit_level(int depth, double local_step)
{
    for (int cell_index : cells_by_level_[depth])
    {
        const auto flat = storage_index(cell_index);
        auto& cell = state_.cells[flat];
        cell.U = step_start_[cell_index]
            + pending_flux_integral_[cell_index] / mesh_.volume_cell[flat]
            + source_[flat] * local_step;
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
        if (star_pressure_cache_) star_pressure_cache_->swap(pressure_hints_);
        return macro_step_;
    }
    throw std::runtime_error("LTS Euler retry limit exceeded");
}
}

double advance_lts_euler(
    MeshBlock& mesh, double time, double dt_limit,
    hydro::SpatialScheme spatial_scheme, const hydro::RiemannSolver& solver,
    const hd_source::SourceList& source_list, const config::SimulationConfig& cfg,
    const hd_user::UserHooks& hooks, double cfl, hydro::StarPressureCache* star_pressure_cache)
{
    const auto options = validate_lts_euler_inputs(
        mesh, time, dt_limit, spatial_scheme, cfg, cfl);

    LTSEulerStep step(
        mesh, time, dt_limit, solver, source_list, cfg, hooks, cfl,
        star_pressure_cache, options.method, options.limiter);
    return step.run();
}
}
