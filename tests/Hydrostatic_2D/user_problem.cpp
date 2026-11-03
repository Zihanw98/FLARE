#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <cmath>
#include <stdexcept>

static_assert(hd_core::NDIM == 2, "Hydrostatic equilibrium requires 2D");

namespace
{
void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (!cfg.constant_gravity_enabled || cfg.gravity[0] != 0.0 || cfg.gravity[1] >= 0.0)
        throw std::invalid_argument("This test requires enabled, downward vertical gravity");
    // Isothermal equilibrium: p/rho = 1, rho(y_lower)=p(y_lower)=1.
    // Thus dp/dy=rho*g_y. This stratification is convectively stable for gamma>1.
    // Isothermal describes the initial equilibrium, not an imposed cooling law.
    const double a = cfg.gravity[1];
    for (int j = 0; j < mesh.nactive[1]; ++j)
    {
        const int row = j + mesh.nghost[1];
        const double y = mesh.cell_center(1, row) - cfg.lower[1];
        const double z = 0.5 * a * mesh.cell_width(1, row);
        // Exact cell average of the exponential, avoiding a point-value bias.
        const double average_factor = std::abs(z) < 1e-6
            ? 1.0 + z*z/6.0 + z*z*z*z/120.0 : std::sinh(z)/z;
        hd_core::Primitives w{};
        w[hd_core::RHO] = std::exp(a*y) * average_factor;
        w[hd_core::PRE] = w[hd_core::RHO];
        for (int i = 0; i < mesh.nactive[0]; ++i)
            mesh.at_active({i, j}) = hd_core::Cell(w, cfg.gamma);
    }
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    return hooks;
}
