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

    using Initializer = void (*)(
        hd_core::MeshBlock&,
        const config::SimulationConfig&
    );
    using UserBoundaryCondition = void (*)(
        hd_core::MeshBlock&,
        double time,
        int direction,
        boundary::BoundarySide side
    );
    using UserOutput = void (*)(
        hd_core::MeshBlock&,
        const config::SimulationConfig&
    );
    struct UserHooks
    {
        Initializer initialize = nullptr;
        UserBoundaryCondition userbc = nullptr;
        UserOutput useroutput = nullptr;
        // 将来可以增加可选源项等接口。
    };

    UserHooks make_user_hooks();

}