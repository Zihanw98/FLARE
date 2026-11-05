/**
 * @file user_api.h
 * @brief Problem configuration and callback interfaces.
 */
#pragma once
#include "boundary.h"
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
    /** @brief Problem callbacks; initialization is required and other hooks are optional. */
    struct UserHooks
    {
        Initializer initialize = nullptr;
        UserBoundaryCondition userbc = nullptr;
        UserOutput useroutput = nullptr;
        // 将来可以增加可选源项等接口。
    };

    /**
     * @brief Provide the problem-specific callbacks linked into a simulation.
     * @return UserHooks containing the required initialize callback and any optional boundary/output callbacks.
     * @note Implemented by the selected problem in user_problem.cpp, not by the framework. A user boundary requires userbc to be non-null.
     */
    UserHooks make_user_hooks();

}