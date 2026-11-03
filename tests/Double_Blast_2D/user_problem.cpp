#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <algorithm>
#include <stdexcept>

static_assert(hd_core::NDIM == 2, "Double blast requires HD_NDIM=2");

namespace
{
void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (cfg.constant_gravity_enabled)
        throw std::invalid_argument("Double blast requires gravity disabled");
    const double lx = mesh.active_length(0), ly = mesh.active_length(1);
    const double radius = 0.08 * std::min(lx, ly);
    const double cy = cfg.lower[1] + 0.5 * ly;
    const double cx_left = cfg.lower[0] + 0.375 * lx;
    const double cx_right = cfg.lower[0] + 0.625 * lx;
    for (int j = 0; j < mesh.nactive[1]; ++j)
        for (int i = 0; i < mesh.nactive[0]; ++i)
        {
            const hd_core::Index index{i + mesh.nghost[0], j + mesh.nghost[1]};
            const double x = mesh.cell_center(0, index);
            const double y = mesh.cell_center(1, index) - cy;
            const double dl = x - cx_left, dr = x - cx_right;
            const bool hot = std::min(dl * dl, dr * dr) + y * y < radius * radius;
            hd_core::Primitives w{};
            w[hd_core::RHO] = 1.0;
            w[hd_core::PRE] = hot ? 10.0 : 1.0;
            mesh.at(index) = hd_core::Cell(w, cfg.gamma);
        }
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    return hooks;
}
