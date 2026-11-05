/**
 * @file source.h
 * @brief Physical source terms.
 */
#pragma once
#include <string>
#include "mesh.h"
#include "config.h"
namespace hd_source
{
    // 将未乘 dt 的源项加到 rhs；dt_limit 只能收紧，不能增大。
    // RK2 以 +infinity 初始化每个源项的独立限制；rhs 不应依赖 dt_limit。
    /** @brief Callback adding unscaled source rates to the RHS and only tightening a timestep limit. */
    using SourceFunction = void (*)(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
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
     * @brief Assemble the list of enabled built-in source terms.
     * @param[in] cfg Configuration controlling constant-gravity registration.
     * @param[in] hooks Reserved user hooks; currently unused.
     * @return A list containing constant gravity when enabled, otherwise an empty list.
     */
    SourceList make_sources(
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks
    );

    //Built-in Source Function
    /**
     * @brief Add constant gravitational acceleration to momentum and energy RHS terms.
     * @param[in] mesh Mesh whose conserved density and momentum supply the source terms.
     * @param[in] time Physical stage time; unused for constant gravity.
     * @param[in] cfg Configuration containing the acceleration vector gravity.
     * @param[in,out] rhs Existing RHS; gravity contributions are added without multiplying by dt.
     * @param[in,out] dt_limit Source timestep bound; left unchanged by this source.
     * @pre rhs has mesh.cells.size() entries; states and gravity are valid.
     * @note Updates all stored entries, including ghosts. Time integrators apply the resulting RHS only to active cells.
     */
    void constant_gravity(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
}
