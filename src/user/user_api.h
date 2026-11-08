/**
 * @file user_api.h
 * @brief Problem configuration and callback interfaces.
 */
#pragma once
#include "boundary.h"
#include "physcon/unit.h"
#include <vector>

namespace hd_source
{
    struct SourceTerm;
}
namespace hd_core
{
    struct MeshBlock;
}

namespace config
{
    struct SimulationConfig;
}

namespace hd_user
{
    /**
     * @brief Optional callback run once after parameter parsing and before any mesh or profile is read.
     * @param[in,out] units Simulation unit system. Select scales with init_unitsystem_MLT or init_unitsystem_RHOLV.
     * @param[in] cfg Parsed input parameters. Profile-derived mesh bounds and counts are populated later.
     */
    using UnitInitializer = void (*)(physcon::Unit& units, const config::SimulationConfig& cfg);
    /** @brief Problem callback that initializes cell states and their adiabatic indices. */
    using Initializer = void (*)(
        hd_core::MeshBlock&,
        const config::SimulationConfig&
    );
    /** @brief Problem callback that fills one boundary side at a specified physical time. */
    using UserBoundaryCondition = void (*)(
        hd_core::MeshBlock&,
        double time,
        int direction,
        boundary::BoundarySide side
    );
    /** @brief Optional problem output callback invoked after successful CSV output. */
    using UserOutput = void (*)(
        hd_core::MeshBlock&,
        const config::SimulationConfig&
    );
    /** @brief Append problem-specific sources after built-in source registration. */
    using RegisterSources = void (*)(
        const config::SimulationConfig& cfg,
        const hd_core::MeshBlock& mesh,
        std::vector<hd_source::SourceTerm>& sources
    );
    /** @brief Problem callbacks; initialization is required and other hooks are optional. */
    struct UserHooks
    {
        Initializer initialize = nullptr;
        UserBoundaryCondition userbc = nullptr;
        UserOutput useroutput = nullptr;
        RegisterSources register_sources = nullptr;
        UnitInitializer initialize_units = nullptr;
    };
    // void register_sources(/* 配置、网格、sources */)
    // {
    //     sources.push_back(make_heating(...));
    //     sources.push_back(make_cooling(...));
    //     sources.push_back(make_external_force(...));
    // }
    /**
     * @brief Provide the problem-specific callbacks linked into a simulation.
     * @return UserHooks containing the required initializer and optional boundary/output/source callbacks.
     * @note Implemented by the selected problem in user_problem.cpp, not by the framework. A user boundary requires userbc to be non-null.
     */
    UserHooks make_user_hooks();

}
