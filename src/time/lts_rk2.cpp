/**
 * 守恒的多层 Heun/RK2 局部时间步。完整算法说明见 docs/lts_rk2.md。
 *
 * 总体 pipeline：
 * 1. 校验配置，复制大步起点状态并填充边界；
 * 2. 求所有独立面的初始通量、波速和源项，确定同步大步 H；
 * 3. 根据各单元 CFL 时间步分配 level，并给重构模板扩展缓冲层；
 * 4. 每个面交给左右单元中更细的 level，保证面通量只积分一次；
 * 5. advance_level() 递归执行两个半步，慢层用线性预测提供同一时刻状态；
 * 6. 用 Heun 梯形公式累计面通量、源项和压力功，再同步双能量；
 * 7. 任一阶段失败则 H 减半并完全重试，成功后才提交 mesh 和压力缓存。
 *
 * LTSRK2Step 只在本文件可见，集中管理预测量、面寄存器和重试状态。
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
namespace
{
using hd_core::Conserved;
constexpr double inf = std::numeric_limits<double>::infinity();

// 与 H > allowance/abs(rate) 相同。远离阈值时先用乘法快速接受，
// 接近阈值或出现极端浮点数时回到原除法，保留原来的舍入判定。
bool source_change_exceeds(double H, double twice_H, double allowance, double rate)
{
    const double magnitude = std::abs(rate);
    if (twice_H * magnitude < allowance) return false;
    return rate != 0.0 && H > allowance / magnitude;
}

struct ValidatedOptions
{
    hydro::Reconstructor method;
    hydro::Limiter limiter;
};

ValidatedOptions validate_lts_rk2_inputs(
    const hd_core::MeshBlock& mesh, double time, double dt_limit,
    hydro::SpatialScheme spatial_scheme, const config::SimulationConfig& cfg,
    double cfl)
{
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

    const bool periodic = mesh.bc_lower[0] == boundary::BoundaryCondition::Periodic;
    if (periodic != (mesh.bc_upper[0] == boundary::BoundaryCondition::Periodic) ||
        (periodic && mesh.geometry != hd_core::Geometry::Cartesian))
        throw std::invalid_argument("LTS RK2 needs paired Cartesian periodic boundaries");

    return {method, limiter};
}

hd_core::MeshBlock make_start_stage(
    const hd_core::MeshBlock& mesh, double time, const hd_user::UserHooks& hooks)
{
    auto start = mesh;
    boundary::FillBoundary(start, time, hooks);
    return start;
}

class LTSRK2Step
{
public:
    LTSRK2Step(
        hd_core::MeshBlock& mesh, double time, double dt_limit,
        const hydro::RiemannSolver& solver, const hd_source::SourceList& sources,
        const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks,
        double cfl, hydro::StarPressureCache* star_pressure_cache,
        hydro::Reconstructor method, hydro::Limiter limiter);

    double run();

private:
    // ---------- 网格拓扑与状态有效性 ----------
    int left_cell(int face) const;
    int right_cell(int face) const;
    static bool is_physical(const Conserved& state, double gamma);

    // ---------- 面重构与 Riemann 求解 ----------
    void solve_faces(const hd_core::MeshBlock& stage, const std::vector<int>& faces);

    // ---------- 源项计算与 10% 变化限制 ----------
    template<class Bound>
    void scan_source_changes(
        const hd_core::MeshBlock& stage, const std::vector<Conserved>& rates, Bound bound);
    double source_change_limit(
        const hd_core::MeshBlock& stage, const std::vector<Conserved>& rates, double trial_H);
    double evaluate_sources(const hd_core::MeshBlock& stage, double stage_time, double trial_H);

    // ---------- CFL 分层与面所有权 ----------
    void promote_neighbors(int cell, int depth);
    void assign_time_levels();
    void assign_faces_to_levels();
    void reset_attempt();

    // ---------- 同时刻阶段视图与递归 Heun 子循环 ----------
    bool check_wave_speeds(const std::vector<int>& faces);
    bool build_stage_view(int depth, double fraction, bool endpoint);
    bool advance_level(int depth, double fraction, bool reuse_start_stage);

    // 外部输入：引用的生命周期覆盖本次 run()。
    hd_core::MeshBlock& mesh_;
    const hydro::RiemannSolver& solver_;
    const hd_source::SourceList& sources_;
    const config::SimulationConfig& cfg_;
    const hd_user::UserHooks& hooks_;
    hydro::StarPressureCache* star_pressure_cache_;

    // 本次同步大步的固定配置和网格尺寸。
    double time_;
    double dt_limit_;
    double cfl_;
    hydro::Reconstructor method_;
    hydro::Limiter limiter_;
    int n_;
    int ng_;
    bool periodic_;
    int independent_face_count_;
    int max_depth_;
    int level_buffer_;

    // start_ 永不修改；state_ 保存已接受的局部步；view_ 是阶段预测快照。
    hd_core::MeshBlock start_;
    hd_core::MeshBlock state_;
    hd_core::MeshBlock view_;
    const double* area_;
    const double* volume_;

    // 面量：当前值、同步起点快照，以及 Heun 第一阶段通量。
    hydro::StarPressureCache hints_;
    hydro::StarPressureCache initial_hints_;
    std::vector<Conserved> flux_;
    std::vector<Conserved> initial_flux_;
    std::vector<Conserved> first_flux_;
    std::vector<double> face_velocity_;
    std::vector<double> initial_velocity_;
    std::vector<double> pressure_work_start_;
    std::vector<double> speed_;
    std::vector<double> initial_speed_;
    std::vector<double> local_dt_;
    std::vector<double> local_step_start_;

    // 源项工作区及同一阶段内复用的热力学尺度缓存。
    std::vector<Conserved> source_;
    std::vector<Conserved> initial_source_;
    std::vector<Conserved> source_contribution_;
    std::vector<double> source_rho_;
    std::vector<double> source_density_allowance_;
    std::vector<double> source_momentum_allowance_;
    std::vector<double> source_internal_allowance_;
    bool source_cache_ready_ = false;

    // 局部步预测器、跨层面通量寄存器和层级索引。
    std::vector<Conserved> base_;
    std::vector<Conserved> slope_;
    std::vector<Conserved> source_start_;
    std::vector<Conserved> pending_flux_integral_;
    std::vector<int> level_;
    std::vector<int> all_faces_;
    std::vector<std::vector<int>> cells_by_level_;
    std::vector<std::vector<int>> faces_needed_by_level_;
    std::vector<std::vector<int>> faces_owned_by_level_;
    int deepest_ = 0;
    double macro_step_ = 0.0;

    hydro::ReconstructionCache reconstruction_;
};

LTSRK2Step::LTSRK2Step(
    hd_core::MeshBlock& mesh, double time, double dt_limit,
    const hydro::RiemannSolver& solver, const hd_source::SourceList& sources,
    const config::SimulationConfig& cfg, const hd_user::UserHooks& hooks,
    double cfl, hydro::StarPressureCache* star_pressure_cache,
    hydro::Reconstructor method, hydro::Limiter limiter)
    : mesh_(mesh), solver_(solver), sources_(sources), cfg_(cfg), hooks_(hooks),
      star_pressure_cache_(star_pressure_cache), time_(time), dt_limit_(dt_limit), cfl_(cfl),
      method_(method), limiter_(limiter),
      n_(mesh.nactive[0]), ng_(mesh.nghost[0]),
      periodic_(mesh.bc_lower[0] == boundary::BoundaryCondition::Periodic),
      independent_face_count_(periodic_ ? n_ : n_ + 1),
      max_depth_(cfg.lts_max_levels - 1),
      level_buffer_(method_ == hydro::Reconstructor::PCM ? 1 : 2),
      start_(make_start_stage(mesh, time, hooks)), state_(mesh), view_(mesh),
      area_(mesh.face_area[0].faces.data() + ng_), volume_(mesh.volume_cell.data() + ng_),
      flux_(n_ + 1), first_flux_(n_ + 1),
      face_velocity_(n_ + 1), pressure_work_start_(n_),
      speed_(n_ + 1), local_dt_(n_), local_step_start_(n_),
      source_(mesh.cells.size()), source_contribution_(mesh.cells.size()),
      source_rho_(sources.empty() ? 0 : n_),
      source_density_allowance_(sources.empty() ? 0 : n_),
      source_momentum_allowance_(sources.empty() ? 0 : n_),
      source_internal_allowance_(sources.empty() ? 0 : n_),
      base_(n_), slope_(n_), source_start_(n_), pending_flux_integral_(n_),
      level_(n_), all_faces_(independent_face_count_),
      reconstruction_(start_, method_, limiter_, cfg.limiter_beta)
{
    if (star_pressure_cache_) hints_ = *star_pressure_cache_;
    hints_.ensure_size(0, n_ + 1);
    for (int face = 0; face < independent_face_count_; ++face) all_faces_[face] = face;
}

int LTSRK2Step::left_cell(int face) const
{
    return face == 0 ? (periodic_ ? n_ - 1 : -1) : face - 1;
}

int LTSRK2Step::right_cell(int face) const
{
    return face == n_ ? -1 : face;
}

bool LTSRK2Step::is_physical(const Conserved& state, double gamma)
{
    try { state.validate(gamma); }
    catch (const std::invalid_argument&) { return false; }
    return true;
}

void LTSRK2Step::solve_faces(
    const hd_core::MeshBlock& stage, const std::vector<int>& faces)
{
    // 将指定面重构到同一阶段时刻，再求数值通量和 CFL 波速。
    reconstruction_.reset(stage, method_, limiter_, cfg_.limiter_beta);
    for (int face : faces)
    {
        hd_core::Index right_index{};
        right_index[0] = ng_ + face;
        const auto [left_state, right_state] = performance::measure(
            performance::Section::Reconstruction,
            [&] { return reconstruction_.face(0, right_index); }, true);

        speed_[face] = -1.0;
        face_velocity_[face] = 0.0;
        flux_[face] = performance::measure(performance::Section::Riemann, [&] {
            return solver_.flux_with_pressure_hint(left_state, right_state, 0,
                stage.cells[ng_ + face - 1].gamma, stage.cells[ng_ + face].gamma,
                hints_[0][face], speed_[face], face_velocity_[face]);
        }, true);

        if (!std::isfinite(speed_[face]) || speed_[face] < 0.0)
            throw std::runtime_error("LTS RK2 received an invalid wave speed");
        for (double value : flux_[face].data)
            if (!std::isfinite(value))
                throw std::runtime_error("LTS RK2 received a non-finite flux");
    }

    if (periodic_)
    {
        flux_[n_] = flux_[0];
        speed_[n_] = speed_[0];
        face_velocity_[n_] = face_velocity_[0];
        hints_[0][n_] = hints_[0][0];
    }
}

template<class Bound>
void LTSRK2Step::scan_source_changes(
    const hd_core::MeshBlock& stage, const std::vector<Conserved>& rates, Bound bound)
{
    // 第一个源项填充尺度缓存；同阶段的后续源项及合计直接复用。
    for (int cell_index = 0; cell_index < n_; ++cell_index)
    {
        const auto& rate = rates[ng_ + cell_index];
        for (double value : rate.data)
            if (!std::isfinite(value)) throw std::runtime_error("Non-finite LTS RK2 source");

        if (!source_cache_ready_)
        {
            const auto& cell = stage.cells[ng_ + cell_index];
            const auto primitive = cell.U.to_Primtivies(cell.gamma);
            const double rho = primitive[hd_core::RHO];
            const double pressure = primitive[hd_core::PRE];
            source_rho_[cell_index] = rho;
            source_density_allowance_[cell_index] = 0.1 * rho;
            source_momentum_allowance_[cell_index] = 0.0;
            source_internal_allowance_[cell_index] = 0.0;
            if (rho != 0.0)
            {
                const double velocity = primitive[hd_core::VEL(0)];
                const double sound_speed = std::sqrt(cell.gamma * pressure / rho);
                source_momentum_allowance_[cell_index] =
                    0.1 * (rho * (std::abs(velocity) + sound_speed));
                source_internal_allowance_[cell_index] =
                    0.1 * (pressure / (cell.gamma - 1.0));
            }
        }

        if (source_rho_[cell_index] == 0.0)
        {
            for (double value : rate.data)
                if (value != 0.0)
                    throw std::invalid_argument("LTS RK2 cannot create gas in vacuum");
            continue;
        }

        bound(source_density_allowance_[cell_index], rate[hd_core::RHO]);
        bound(source_momentum_allowance_[cell_index], rate[hd_core::MOMDEN(0)]);
        bound(source_internal_allowance_[cell_index], rate[hd_core::INTERNAL_ENERGY]);
    }
    source_cache_ready_ = true;
}

double LTSRK2Step::source_change_limit(
    const hd_core::MeshBlock& stage, const std::vector<Conserved>& rates, double trial_H)
{
    if (trial_H > 0.0)
    {
        const double twice_H = 2.0 * trial_H;
        bool exceeded = false;
        scan_source_changes(stage, rates, [&](double allowance, double rate) {
            exceeded |= source_change_exceeds(trial_H, twice_H, allowance, rate);
        });
        return exceeded ? 0.0 : inf;
    }

    double limit = inf;
    scan_source_changes(stage, rates, [&](double allowance, double rate) {
        if (rate != 0.0) limit = std::min(limit, allowance / std::abs(rate));
    });
    return limit;
}

double LTSRK2Step::evaluate_sources(
    const hd_core::MeshBlock& stage, double stage_time, double trial_H)
{
    performance::Scope timer(performance::Section::Sources);
    // A stage view reuses its MeshBlock object while replacing the states.
    source_cache_ready_ = false;
    std::fill(source_.begin(), source_.end(), Conserved{});
    double limit = inf;

    for (const auto& term : sources_)
    {
        if (!term.func) throw std::invalid_argument("Missing LTS RK2 source function");
        std::fill(source_contribution_.begin(), source_contribution_.end(), Conserved{});
        double callback_limit = inf;
        term.func(stage, stage_time, cfg_, source_contribution_, callback_limit);
        hd_source::complete_internal_source(stage, term, source_contribution_);
        if (source_contribution_.size() != source_.size())
            throw std::logic_error("Source resized LTS RK2 storage");
        if (std::isnan(callback_limit) || callback_limit <= 0.0)
            throw std::invalid_argument("Invalid LTS RK2 source limit");

        limit = std::min({limit, callback_limit,
                          source_change_limit(stage, source_contribution_, trial_H)});
        for (int cell = 0; cell < n_; ++cell)
            source_[ng_ + cell] += source_contribution_[ng_ + cell];
    }

    return sources_.empty()
        ? inf
        : std::min(limit, source_change_limit(stage, source_, trial_H));
}

void LTSRK2Step::promote_neighbors(int cell, int depth)
{
    for (int offset = -level_buffer_; offset <= level_buffer_; ++offset)
    {
        int neighbor = cell + offset;
        if (periodic_) neighbor = (neighbor % n_ + n_) % n_;
        if (neighbor >= 0 && neighbor < n_)
            level_[neighbor] = std::max(level_[neighbor], depth);
    }
}

void LTSRK2Step::assign_time_levels()
{
    // level=0 使用 H，level=d 使用 H/2^d；缓冲层保护 PCM/PLM 模板。
    std::fill(level_.begin(), level_.end(), 0);
    for (int cell = 0; cell < n_; ++cell)
    {
        int depth = 0;
        while (depth < max_depth_ && std::ldexp(macro_step_, -depth) > local_dt_[cell])
            ++depth;
        promote_neighbors(cell, depth);
    }

    deepest_ = *std::max_element(level_.begin(), level_.end());
    for (int depth = deepest_; depth > 1; --depth)
        for (int cell = 0; cell < n_; ++cell)
            if (level_[cell] == depth) promote_neighbors(cell, depth - 1);
}

void LTSRK2Step::assign_faces_to_levels()
{
    // 一个面由相邻单元中更细的层负责积分，但粗层构造预测时仍可能需要它。
    cells_by_level_.assign(deepest_ + 1, {});
    faces_needed_by_level_.assign(deepest_ + 1, {});
    faces_owned_by_level_.assign(deepest_ + 1, {});

    for (int cell = 0; cell < n_; ++cell)
        cells_by_level_[level_[cell]].push_back(cell);

    for (int face = 0; face < independent_face_count_; ++face)
    {
        const int left = left_cell(face);
        const int right = right_cell(face);
        const int owner_depth = std::max(
            left >= 0 ? level_[left] : 0,
            right >= 0 ? level_[right] : 0);
        faces_owned_by_level_[owner_depth].push_back(face);
        for (int depth = 0; depth <= owner_depth; ++depth)
            faces_needed_by_level_[depth].push_back(face);
    }
}

void LTSRK2Step::reset_attempt()
{
    assign_time_levels();
    assign_faces_to_levels();
    state_.cells = start_.cells;
    face_velocity_ = initial_velocity_;
    hints_ = initial_hints_;
    flux_ = initial_flux_;
    speed_ = initial_speed_;
    source_ = initial_source_;
    std::fill(pending_flux_integral_.begin(), pending_flux_integral_.end(), Conserved{});
}

bool LTSRK2Step::check_wave_speeds(const std::vector<int>& faces)
{
    for (int face : faces)
    {
        for (int cell : {left_cell(face), right_cell(face)})
        {
            if (cell < 0) continue;
            const double rate = speed_[face] / mesh_.cell_width(0, ng_ + cell);
            hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
            const double cell_step = std::ldexp(macro_step_, -level_[cell]);
            if (!std::isfinite(rate) || cell_step * rate > cfg_.cfl_max) return false;
        }
    }
    return true;
}

bool LTSRK2Step::build_stage_view(int depth, double fraction, bool endpoint)
{
    // 把慢层线性预测到目标时刻；细层使用已经推进到该时刻的 state_。
    view_.cells = state_.cells;
    for (int cell_index = 0; cell_index < n_; ++cell_index)
    {
        const bool needs_prediction =
            level_[cell_index] < depth || (endpoint && level_[cell_index] == depth);
        if (!needs_prediction) continue;

        auto& cell = view_.cells[ng_ + cell_index];
        const double elapsed = (fraction - local_step_start_[cell_index]) * macro_step_;
        cell.U = base_[cell_index] + slope_[cell_index] * elapsed;
        if (!is_physical(cell.U, cell.gamma)) return false;
    }
    boundary::FillBoundary(view_, time_ + fraction * macro_step_, hooks_);
    return true;
}

bool LTSRK2Step::advance_level(int depth, double fraction, bool reuse_start_stage)
{
    // 当前层走一个 h=H/2^depth；更细层在中间递归走两个半步。
    const double span = std::ldexp(1.0, -depth);
    const double local_step = macro_step_ * span;
    const double end_fraction = fraction + span;
    if (time_ + end_fraction * macro_step_ <= time_ + fraction * macro_step_)
        throw std::runtime_error("LTS RK2 substep cannot advance time");

    if (!reuse_start_stage)
    {
        if (!build_stage_view(depth, fraction, false)) return false;
        solve_faces(view_, faces_needed_by_level_[depth]);
        if (macro_step_ > evaluate_sources(
                view_, time_ + fraction * macro_step_, macro_step_)) return false;
    }
    if (!check_wave_speeds(faces_needed_by_level_[depth])) return false;

    for (int face : faces_owned_by_level_[depth]) first_flux_[face] = flux_[face];
    for (int cell_index : cells_by_level_[depth])
    {
        const auto& current_cell = state_.cells[ng_ + cell_index];
        base_[cell_index] = current_cell.U;
        local_step_start_[cell_index] = fraction;
        source_start_[cell_index] = source_[ng_ + cell_index];

        pressure_work_start_[cell_index] =
            -(current_cell.gamma - 1.0) * base_[cell_index].internal_energy()
            * (area_[cell_index + 1] * face_velocity_[cell_index + 1]
               - area_[cell_index] * face_velocity_[cell_index])
            / volume_[cell_index];
        slope_[cell_index] =
            (flux_[cell_index] * area_[cell_index]
             - flux_[cell_index + 1] * area_[cell_index + 1])
            / volume_[cell_index] + source_start_[cell_index];
        slope_[cell_index][hd_core::INTERNAL_ENERGY] += pressure_work_start_[cell_index];

        if (!is_physical(base_[cell_index] + local_step * slope_[cell_index],
                         current_cell.gamma)) return false;
    }

    if (depth < deepest_)
    {
        if (!advance_level(depth + 1, fraction, true)) return false;
        if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
    }

    if (cells_by_level_[depth].empty()) return true;
    if (!build_stage_view(depth, end_fraction, true)) return false;
    solve_faces(view_, faces_owned_by_level_[depth]);
    if (!check_wave_speeds(faces_owned_by_level_[depth]) ||
        macro_step_ > evaluate_sources(
            view_, time_ + end_fraction * macro_step_, macro_step_)) return false;

    for (int face : faces_owned_by_level_[depth])
    {
        // 每个面每个局部步只产生一个梯形积分，左右单元共享同一份积分。
        const auto integral =
            (first_flux_[face] + flux_[face]) * (0.5 * local_step * area_[face]);
        const int left = left_cell(face);
        const int right = right_cell(face);
        if (left >= 0) pending_flux_integral_[left] -= integral;
        if (right >= 0) pending_flux_integral_[right] += integral;
    }

    for (int cell_index : cells_by_level_[depth])
    {
        auto& cell = state_.cells[ng_ + cell_index];
        cell.U = base_[cell_index]
            + pending_flux_integral_[cell_index] / volume_[cell_index]
            + (source_start_[cell_index] + source_[ng_ + cell_index])
                * (0.5 * local_step);
        const double pressure_work_end =
            -(cell.gamma - 1.0) * view_.cells[ng_ + cell_index].U.internal_energy()
            * (area_[cell_index + 1] * face_velocity_[cell_index + 1]
               - area_[cell_index] * face_velocity_[cell_index])
            / volume_[cell_index];
        cell.U[hd_core::INTERNAL_ENERGY] +=
            0.5 * local_step * (pressure_work_start_[cell_index] + pressure_work_end);
        if (!is_physical(cell.U, cell.gamma)) return false;
        cell.synchronize_dual_energy();
        pending_flux_integral_[cell_index] = Conserved{};
    }
    return true;
}

double LTSRK2Step::run()
{
    // Pipeline 1：同步起点的全体面量和源项，给选步与所有重试提供快照。
    solve_faces(start_, all_faces_);
    initial_flux_ = flux_;
    initial_speed_ = speed_;
    initial_velocity_ = face_velocity_;
    initial_hints_ = hints_;
    const double source_limit = evaluate_sources(start_, time_, 0.0);
    initial_source_ = source_;

    double minimum_local_dt = inf;
    hydro::max_inverse_dt = 0.0;
    for (int cell = 0; cell < n_; ++cell)
    {
        const double rate =
            std::max(speed_[cell], speed_[cell + 1])
            / mesh_.cell_width(0, ng_ + cell);
        if (!std::isfinite(rate)) throw std::runtime_error("Non-finite LTS RK2 CFL rate");
        hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
        local_dt_[cell] = rate > 0.0 ? cfl_ / rate : inf;
        minimum_local_dt = std::min(minimum_local_dt, local_dt_[cell]);
    }

    // Pipeline 2：最大层仍能满足最小 CFL 步，同时服从终止时间和源项上限。
    macro_step_ = std::min({
        dt_limit_, std::ldexp(minimum_local_dt, max_depth_), source_limit});
    for (long long attempt = 0; attempt <= cfg_.integrator_retry_max;
         ++attempt, macro_step_ *= 0.5)
    {
        // Pipeline 3：每次重试都重新分层，并从同步起点恢复全部可回滚状态。
        performance::count(performance::Counter::RK2Attempts);
        if (!(macro_step_ > 0.0) || !std::isfinite(time_ + macro_step_) ||
            time_ + macro_step_ <= time_)
            throw std::runtime_error("LTS RK2 macro step cannot advance time");

        reset_attempt();
        if (!advance_level(0, 0.0, true))
        {
            performance::count(performance::Counter::RK2Rejected);
            continue;
        }

        // Pipeline 4：递归全部成功后才提交正式网格和可选压力缓存。
        boundary::FillBoundary(state_, time_ + macro_step_, hooks_);
        mesh_.cells.swap(state_.cells);
        if (star_pressure_cache_) star_pressure_cache_->swap(hints_);
        return macro_step_;
    }
    throw std::runtime_error("LTS RK2 retry limit exceeded");
}
}

double advance_lts_rk2(
    hd_core::MeshBlock& mesh, double time, double dt_limit,
    hydro::SpatialScheme spatial_scheme, const hydro::RiemannSolver& solver,
    const hd_source::SourceList& sources, const config::SimulationConfig& cfg,
    const hd_user::UserHooks& hooks, double cfl,
    hydro::StarPressureCache* star_pressure_cache)
{
    const auto options = validate_lts_rk2_inputs(
        mesh, time, dt_limit, spatial_scheme, cfg, cfl);

    LTSRK2Step step(mesh, time, dt_limit, solver, sources, cfg, hooks, cfl,
                    star_pressure_cache, options.method, options.limiter);
    return step.run();
}
}
