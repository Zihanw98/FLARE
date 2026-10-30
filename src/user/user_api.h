#pragma once

namespace hd_core
{
    struct MeshBlock;
}

namespace hd_config
{
    struct SimulationConfig;
}

namespace hd_user
{

using Initializer = void (*)(
    hd_core::MeshBlock&,
    const hd_config::SimulationConfig&
);

struct UserHooks
{
    Initializer initialize = nullptr;

    // 将来可以增加可选源项等接口。
};

UserHooks make_user_hooks();

}