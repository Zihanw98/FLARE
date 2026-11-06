#include "user_config.h"  // Load this problem's settings before framework headers.
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <cmath>
#include "physcon/constants.h"
static_assert(hd_core::NDIM == 2, "KH requires a 2D setting");

namespace
{
    constexpr double A = 0.01;
    constexpr double sigma = 0.05; // Fraction of the domain height.
    constexpr int perturbation_modes = 4; // Four wavelengths along each shear layer.

    // Periodic double shear layer at 1/4 and 3/4 of the domain height.
    void initialize_KH2D(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        const double length_x = mesh.active_length(0);
        const double length_y = mesh.active_length(1);
        for(int y=0; y < mesh.nactive[1]; y++)
        {
            for (int x=0; x < mesh.nactive[0]; x++)
            {
                const hd_core::Index index{x+mesh.nghost[0], y+mesh.nghost[1]};
                const double x_fraction = (mesh.cell_center(0,index) - cfg.lower[0]) / length_x;
                const double y_fraction = (mesh.cell_center(1,index) - cfg.lower[1]) / length_y;
                const double d_1 = y_fraction - 0.25;
                const double d_2 = y_fraction - 0.75;
                const bool middle = y_fraction > 0.25 && y_fraction < 0.75;
                hd_core::Primitives W{};
                W[hd_core::RHO] = middle? 2 : 1;
                W[hd_core::VEL(0)] = middle? 0.5 : -0.5;
                W[hd_core::VEL(1)] = A * std::sin(2*physcon::pi * perturbation_modes * x_fraction)
                    * (std::exp(-(d_1*d_1)/(2*sigma*sigma))
                       + std::exp(-(d_2*d_2)/(2*sigma*sigma)));
                W[hd_core::PRE] = 2.5;
                mesh.at(index) = hd_core::Cell(W, cfg.gamma);

            }
        }
    }

}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize_KH2D;
    return hooks;
}
