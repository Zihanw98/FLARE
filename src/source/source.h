/**
 * @file source.h
 * @brief Physical source terms.
 */
#pragma once
#include <array>
#include <functional>
#include <string>
#include <vector>
#include "mesh.h"
#include "config.h"

namespace hd_user { struct UserHooks; }

namespace hd_source
{
    // 将未乘 dt 的源项加到 rhs；dt_limit 只能收紧，不能增大。
    // RK2 以 +infinity 初始化每个源项的独立限制；rhs 不应依赖 dt_limit。
    /**
     * @brief Callable adding unscaled source rates and only tightening a timestep limit.
     * @note Capture source-specific parameters by value (or own shared immutable data).
     *       RK2 may evaluate rejected stages; callbacks must not commit physical state
     *       outside rhs or make their rates depend on the number of calls.
     */
    using SourceFunction = std::function<void(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    )>;
    // 输运源项：独立面沿坐标轴正向的通量（不乘面积、dt）。输出初始为零。
    // inverse_dt[0/1] 是该面对左/右 cell 的非负逆时间尺度贡献，包含面积/体积，
    // 不含 CFL。F=-D*grad(q) 时为 A*D/(V_left/right*distance)。多过程贡献相加。
    // 从同一阶段的相邻 cell/ghost 计算梯度，不使用重构跳跃；不修改 mesh。
    // LTS 和隐式预条件器目前支持相邻单元模板。重试/扰动可重复调用，不能提交外部状态。
    // 必须显式填写受影响的总能量与辅助内能通量，不执行普通源项的内能自动补全。
    using TransportFunction = std::function<void(
        const hd_core::MeshBlock& stage, std::size_t face, double time,
        const config::SimulationConfig& cfg, hd_core::Conserved& flux,
        std::array<double, 2>& inverse_dt)>;
    /** @brief Physical categories used to label registered sources. */
    enum class SourceKind
    {
        Thermal, ///< Source categorized as thermal.
        Force, ///< Source categorized as a force.
        General, ///< Source without a specialized physical category.
        Transport ///< Conservative face transport; uses transport_flux instead of func.
    };

    /** @brief Named source callback and its physical category. */
    struct SourceTerm
    {
        std::string name;
        SourceKind kind;
        SourceFunction func;
        bool supplies_internal_energy = false; // legacy callbacks: derive thermal rate from E,m,rho
        TransportFunction transport_flux{}; // Only used when kind == Transport.
    };

    /** @brief Ordered collection of source terms evaluated for each RHS. */
    using SourceList = std::vector<SourceTerm>;
    // 输运接入：仅计算标为 Transport 的源项，叠加到已有面通量；覆盖 inverse_dt。
    // 左右 cell 的更新由积分器统一执行，回调中不能先计算成逐 cell 的源项。
    void add_transport_flux(
        const hd_core::MeshBlock& stage, std::size_t face, double time,
        const config::SimulationConfig& cfg, const SourceList& sources,
        hd_core::Conserved& flux, std::array<double, 2>& inverse_dt);
    // 全局显式/隐式 RHS 共用的守恒散度；返回各 cell 输运逆时间尺度的最大值。
    // 每个独立面只计算一次，物理边界只更新域内一侧。
    double add_transport_rhs(
        const hd_core::MeshBlock& stage, double time,
        const config::SimulationConfig& cfg, const SourceList& sources,
        std::vector<hd_core::Conserved>& rhs);
    void complete_internal_source(const hd_core::MeshBlock& mesh, const SourceTerm& term,
                                  std::vector<hd_core::Conserved>& rates);

    /**
     * @brief Assemble enabled built-in sources, then invoke the optional user registration hook.
     * @param[in] cfg Configuration controlling constant gravity and heating registration.
     * @param[in] mesh Mesh geometry controlling geometric-source registration.
     * @param[in] hooks Optional register_sources callback, called once while building this list.
     * @return Built-in sources followed by any sources appended by the user hook.
     * @note Source parameters are copied at construction; rebuild the list to change them.
     */
    SourceList make_sources(
        const config::SimulationConfig& cfg,
        const hd_core::MeshBlock& mesh,
        const hd_user::UserHooks& hooks
    );

    /** @brief Parameters owned by a constant-gravity source. */
    struct ConstantGravityParams
    {
        std::array<double, hd_core::NDIM> acceleration{};
    };
    struct SphericalGravityParams
    {
        double central_mass;
        bool self_gravity_enabled;
    };
    /** @brief Create a gravity source owning a copy of params. */
    SourceTerm make_constant_gravity(ConstantGravityParams params);

    /** @brief Create a heating source owning one finite nonnegative specific heating rate. */
    SourceTerm make_constant_heating(double heating_rate);

    /** @brief Create a geometric source using the mesh passed at each evaluation. */
    SourceTerm make_geometric_source();

    SourceTerm make_spherical_self_gravity(SphericalGravityParams params);
    // Built-in source implementations; factories adapt them to SourceFunction.
    /**
     * @brief Add constant gravitational acceleration to momentum and energy RHS terms.
     * @param[in] mesh Mesh whose conserved density and momentum supply the source terms.
     * @param[in] time Physical stage time; unused for constant gravity.
     * @param[in] cfg Shared simulation configuration; unused for constant gravity.
     * @param[in] params Constant acceleration vector for this source instance.
     * @param[in,out] rhs Existing RHS; gravity contributions are added without multiplying by dt.
     * @param[in,out] dt_limit Source timestep bound; left unchanged by this source.
     * @pre rhs has mesh.cells.size() entries; states and gravity are valid.
     * @note Updates all stored entries, including ghosts. Time integrators apply the resulting RHS only to active cells.
     */
    void constant_gravity(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        const ConstantGravityParams& params,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
    void geometric_source(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
    void spherical_self_gravity(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& /*cfg*/,
        const SphericalGravityParams& params,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
    /**
     * @brief Add rho * heating_rate to active-cell energy rates.
     * @param[in] mesh Current stage mesh; positive-density cells must have positive internal energy.
     * @param[in] time Stage time; unused for constant heating.
     * @param[in] cfg Shared configuration; unused by this source.
     * @param[in] heating_rate Energy added per unit mass per unit time, finite and nonnegative.
     * @param[in,out] rhs Existing RHS; only active-cell energy entries are changed, without multiplying by dt.
     * @param[in,out] dt_limit Tightened to min(0.1 * internal_energy / (rho * heating_rate)) over active gas cells.
     * @note Zero heating rate and vacuum cells add nothing and impose no timestep limit.
     * @throws std::invalid_argument If heating_rate is negative or nonfinite.
     */
    void constant_heating(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        double heating_rate,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
}
