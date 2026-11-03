#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <stdexcept>

static_assert(hd_core::NDIM == 2, "Gravity adjustment requires 2D");

namespace
{
void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (!cfg.constant_gravity_enabled || cfg.gravity[0] != 0.0 || cfg.gravity[1] >= 0.0)
        throw std::invalid_argument("This test requires enabled, downward vertical gravity");
    hd_core::Primitives w{};
    w[hd_core::RHO] = 1.0;
    w[hd_core::PRE] = 1.0;
    for (int j = 0; j < mesh.nactive[1]; ++j)
        for (int i = 0; i < mesh.nactive[0]; ++i)
            mesh.at_active({i, j}) = hd_core::Cell(w, cfg.gamma);
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    return hooks;
}
