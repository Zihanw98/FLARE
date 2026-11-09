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
    /** @brief Physical categories used to label registered sources. */
    enum class SourceKind
    {
        Thermal, ///< Source categorized as thermal.
        Force, ///< Source categorized as a force.
        General ///< Source without a specialized physical category.
    };

    /** @brief Named source callback and its physical category. */
    struct SourceTerm
    {
        std::string name;
        SourceKind kind;
        SourceFunction func;
    };

    /** @brief Ordered collection of source terms evaluated for each RHS. */
    using SourceList = std::vector<SourceTerm>;

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
