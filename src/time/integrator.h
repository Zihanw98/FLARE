/**
 * @file integrator.h
 * @brief Time integration.
 */
#pragma once
#include "hydro/spatial.h"
#include "source/source.h"
#include "config/config.h"
#include <string>
#include <array>
#include <vector>
#include <utility>
#include <memory>
namespace hd_core
{
    struct MeshBlock;
}

namespace hd_time
{
    /** Reusable scratch storage for a serial LTS evolution. Own one per caller;
     * never share concurrently. Geometry changes invalidate its mesh views. */
    struct LTSEulerWorkspace
    {
        struct Impl;
        std::unique_ptr<Impl> impl;
        LTSEulerWorkspace();
        ~LTSEulerWorkspace();
        LTSEulerWorkspace(const LTSEulerWorkspace&) = delete;
        LTSEulerWorkspace& operator=(const LTSEulerWorkspace&) = delete;
    };
    /** @brief Available time integration methods. */
    enum class TimeIntegrator
    {
        Euler, ///< One-stage forward Euler update.
        LTSEuler, ///< Experimental 1D Euler with two-level conservative subcycling.
        RK2 ///< Two-stage SSP-RK2 update with timestep retries.
    };

    /**
     * @brief Parse a case-insensitive time-integrator name.
     * @param[in] name euler, lts_euler or rk2.
     * @return The corresponding TimeIntegrator enum.
     * @throws std::invalid_argument If the method is unknown or unavailable.
     */
    TimeIntegrator parse_time_integrator(const std::string& name);

    /**
     * @brief Advance one synchronized macro step using two-level explicit subcycling.
     * @details Restricted to 1D Godunov with PCM, PLM or PLM_corrected.
     * Slow cells take H and fast cells two
     * H/2 steps. Fine integrated face fluxes replace coarse fluxes at interfaces.
     * PLM variants require two ghost cells and predict slow stencil states to
     * the half-step time before reconstruction; time accuracy remains first order.
     * Source rates (including self-gravity) are frozen at the synchronized start;
     * source limits and a 10-percent source-change estimate cap H. This is a
     * first-order experimental method, not a general asynchronous RK scheme.
     * @return Accepted macro step, at most twice the initial minimum hydro step.
     * @post All cells and ghosts are synchronized. Failures leave mesh unchanged.
     * @throws std::invalid_argument For unsupported dimensions/reconstruction or inputs.
     * @throws std::runtime_error If retries are exhausted or time cannot advance.
     */
    double advance_lts_euler(
        hd_core::MeshBlock& mesh, double time, double dt_limit,
        hydro::SpatialScheme spatial_scheme, const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list, const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks, double cfl, LTSEulerWorkspace* workspace = nullptr);

    // 对活跃单元执行 U += dt * rhs，并逐个校验候选状态。
    // rhs 与 mesh.cells 等长；dt 必须有限，零表示不更新。
    // 全部校验成功后才提交；失败时保持原 mesh，ghost、gamma 和几何不变。
    /**
     * @brief Add a scaled RHS to active states and commit only after validation.
     * @param[in,out] mesh Mesh whose active conserved states receive dt * rhs.
     * @param[in] rhs One RHS value per stored cell; ghost entries are ignored.
     * @param[in] dt Finite scale; zero performs no state update and negative values are allowed.
     * @pre Cell storage matches the mesh dimensions.
     * @throws std::invalid_argument If dt, RHS size, or a candidate active state is invalid.
     * @throws std::logic_error If cell storage does not match mesh dimensions.
     * @post On success, active states are updated; ghost states, gamma, and geometry are unchanged.
     * @note If validation fails, no cell states are committed.
     */
    void add_scaled_rhs(
        hd_core::MeshBlock& mesh,
        const std::vector<hd_core::Conserved>& rhs,
        double dt
    );

    // 调用前填好当前 time 的 ghost。dt_limit 是外层给出的有限正上限。
    // 返回实际使用的 dt；主循环必须 time += 返回值。
    // cfl 为目标值；dt_next 初始设为 +infinity，RK 方法可返回下一步建议。
    // Euler 每步独立估计，dt_next 返回 +infinity 表示没有额外建议。
    /**
     * @brief Advance one timestep and fill boundaries at the accepted new time.
     * @param[in,out] mesh Allocated mesh whose active states are advanced.
     * @param[in] time Finite physical time at the start of the step.
     * @param[in] dt_limit Finite positive upper bound, such as time remaining to the end.
     * @param[in] integrator Time integration method to dispatch.
     * @param[in] spatial_scheme Hydrodynamic spatial operator.
     * @param[in] solver Riemann solver used for fluxes and CFL wave speeds.
     * @param[in] source_list Source callbacks that add unscaled rates and supply independent timestep limits.
     * @param[in] cfg Numerical settings, including reconstruction, limiter, and RK2 acceptance/retry controls.
     * @param[in] hooks Boundary callbacks used by RK2 stages and final boundary filling.
     * @param[in] cfl Finite positive target CFL; Euler requires at most 0.5, RK2 requires cfl <= cfg.cfl_max <= 0.5.
     * @param[in,out] dt_next RK2 timestep suggestion (initially infinity). Both Euler methods
     * ignore the incoming value and reset it to infinity after acceptance.
     * @param[in,out] pressure_history Optional ordinary face arrays for Euler/RK2
     * initial guesses; commit with accepted states. LTS uses its own workspace.
     * @return The finite positive timestep actually used; the caller advances time by this value.
     * @pre Active states and ghost states at time are initialized and consistent.
     * @throws std::invalid_argument If timestep/CFL settings or a required state are invalid.
     * @throws std::runtime_error If the timestep cannot advance time or RK2 exhausts its retries.
     * @note Euler/RK2 commit before the final boundary fill. LTS commits a
     * boundary-filled candidate.
     */
    double advance_one_step(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        TimeIntegrator integrator,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl,
        double& dt_next,
        LTSEulerWorkspace* lts_workspace = nullptr,
        std::array<std::vector<double>, hd_core::NDIM>* pressure_history = nullptr
    );

    // 同一次计算获取 RHS 和波速，先限制 dt，再更新活跃单元。
    // 更新失败时保持原 mesh；ghost 的 U 和 gamma 不由积分器修改。
    /**
     * @brief Advance active states by one forward-Euler step.
     * @param[in,out] mesh Allocated mesh whose active states are advanced.
     * @param[in] time Finite physical time at the start of the step.
     * @param[in] dt_limit Finite positive upper bound, such as time remaining to the end.
     * @param[in] spatial_scheme Hydrodynamic spatial operator.
     * @param[in] solver Riemann solver used for fluxes and CFL wave speeds.
     * @param[in] source_list Source callbacks that add unscaled rates and supply independent timestep limits.
     * @param[in] cfg Numerical settings, including reconstruction, limiter, and RK2 acceptance/retry controls.
     * @param[in] hooks Reserved boundary hooks; not used by this direct Euler routine.
     * @param[in] cfl Finite positive target CFL; Euler requires at most 0.5, RK2 requires cfl <= cfg.cfl_max <= 0.5.
     * @return The finite positive timestep actually used; the caller advances time by this value.
     * @pre Active states and ghost states at time are initialized and consistent.
     * @throws std::invalid_argument If timestep/CFL settings or a required state are invalid.
     * @throws std::runtime_error If the timestep cannot advance time.
     * @note Commits only after every candidate active state validates; direct Euler leaves ghost states unchanged.
     */
    double advance_euler(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl = 0.4,
        std::array<std::vector<double>, hd_core::NDIM>* pressure_history = nullptr
    );

    // cfl 用于选步长；cfg.cfl_max 用于接受/拒绝。两个阶段共用同一 dt。
    // dt_next：输入上一步建议（首次设为 +infinity），成功后写入下一步建议。
    /**
     * @brief Advance one accepted SSP-RK2 step with CFL checks and retries.
     * @param[in,out] mesh Allocated mesh whose active states are advanced.
     * @param[in] time Finite physical time at the start of the step.
     * @param[in] dt_limit Finite positive upper bound, such as time remaining to the end.
     * @param[in] spatial_scheme Hydrodynamic spatial operator.
     * @param[in] solver Riemann solver used for fluxes and CFL wave speeds.
     * @param[in] source_list Source callbacks that add unscaled rates and supply independent timestep limits.
     * @param[in] cfg Numerical settings, including reconstruction, limiter, and RK2 acceptance/retry controls.
     * @param[in] hooks Boundary callbacks used by RK2 stages and final boundary filling.
     * @param[in] cfl Finite positive target CFL; Euler requires at most 0.5, RK2 requires cfl <= cfg.cfl_max <= 0.5.
     * @param[in,out] dt_next Previous positive timestep suggestion (initially infinity); replaced by the next RK2 suggestion or infinity for Euler.
     * @return The finite positive timestep actually used; the caller advances time by this value.
     * @pre Active states and ghost states at time are initialized and consistent.
     * @throws std::invalid_argument If timestep/CFL settings or a required state are invalid.
     * @throws std::runtime_error If the timestep cannot advance time or RK2 exhausts its retries.
     * @note Rejected RK2 trials do not commit cell states; direct RK2 leaves final ghost states unchanged. The shared wave-speed record may still change.
     */
    double advance_rk2(
        hd_core::MeshBlock& mesh,
        double time,
        double dt_limit,
        hydro::SpatialScheme spatial_scheme,
        const hydro::RiemannSolver& solver,
        const hd_source::SourceList& source_list,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        double cfl,
        double& dt_next,
        std::array<std::vector<double>, hd_core::NDIM>* pressure_history = nullptr
    );
    // 保留 {目标步长, RHS} 返回形式；引用输出用于分别检查 CFL 和源项硬限制。
    /**
     * @brief Evaluate hydrodynamic and source rates together with timestep limits.
     * @param[in] mesh Allocated mesh with valid states and stage-time ghost values.
     * @param[in] time Stage time forwarded to source callbacks.
     * @param[in] scheme Hydrodynamic spatial operator.
     * @param[in] source_list Sources that add unscaled RHS rates without resizing the array.
     * @param[in] solver Riemann solver supplying fluxes and interface wave speeds.
     * @param[in] cfg Reconstruction, limiter, and source settings.
     * @param[in] cfl Finite positive target CFL coefficient.
     * @param[out] hydro_inverse_dt Largest cell inverse timescale measured from interface waves.
     * @param[out] source_dt_limit Minimum independent source limit, or infinity if unrestricted.
     * @param[in,out] pressure_hints Optional stage-local pressure guesses, updated
     * while solving faces. Integrators pass trial arrays to preserve rollback.
     * @return Pair {target timestep, RHS}, where the timestep is the minimum hydrodynamic and source limit and may be infinity.
     * @pre Each source callback is valid, returns a positive limit or infinity, and does not make its rate depend on the incoming limit.
     * @throws std::invalid_argument If CFL, recorded wave rate, callback, or source limit is invalid.
     * @throws std::logic_error If a source changes the RHS array size.
     * @note Each source receives a fresh infinity limit. Updates the shared hydro::wave_timestep record; source exceptions propagate.
     */
    std::pair<double, std::vector<hd_core::Conserved>> evaluate_RHS(
        const hd_core::MeshBlock& mesh,
        double time,
        hydro::SpatialScheme scheme,
        const hd_source::SourceList& source_list,
        const hydro::RiemannSolver& solver,
        const config::SimulationConfig& cfg,
        double cfl,
        double& hydro_inverse_dt,
        double& source_dt_limit,
        std::array<std::vector<double>, hd_core::NDIM>* pressure_hints = nullptr
    );
}
