#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <algorithm>
#include <stdexcept>

static_assert(hd_core::NDIM == 3, "Spherical blast requires HD_NDIM=3");

namespace
{
void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (cfg.constant_gravity_enabled)
        throw std::invalid_argument("Spherical blast requires gravity disabled");

    const double length = std::min({mesh.active_length(0), mesh.active_length(1),
                                    mesh.active_length(2)});
    const double radius = 0.12 * length;
    for (int k = 0; k < mesh.nactive[2]; ++k)
        for (int j = 0; j < mesh.nactive[1]; ++j)
            for (int i = 0; i < mesh.nactive[0]; ++i)
            {
                const hd_core::Index index{i + mesh.nghost[0], j + mesh.nghost[1],
                                           k + mesh.nghost[2]};
                double r2 = 0.0;
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    const double offset = mesh.cell_center(d, index)
                                        - 0.5 * (cfg.lower[d] + cfg.upper[d]);
                    r2 += offset * offset;
                }
                hd_core::Primitives w{};
                w[hd_core::RHO] = 1.0;
                // Finite-radius overpressure sphere, sampled at cell centers.
                // This is not a point-energy Sedov similarity benchmark.
                w[hd_core::PRE] = r2 < radius * radius ? 5.0 : 1.0;
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
