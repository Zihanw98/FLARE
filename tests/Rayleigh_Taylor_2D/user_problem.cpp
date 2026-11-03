#include "user_config.h"  // Load this problem's settings before framework headers.
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <cmath>
static_assert(hd_core::NDIM == 2, "RT initialization requires 2D");

namespace
{
    /**
     * @brief Initialize hydrostatic RT layers with four horizontal perturbation periods.
     * @param[in,out] mesh Allocated 2D mesh; built-in boundaries fill the ghosts.
     * @param[in] cfg Supplies the domain, gas gamma, and vertical gravity.
     */
    void initialize_RTI2D(
        hd_core::MeshBlock& mesh,
        const config::SimulationConfig& cfg
    )
    {
        static_assert(hd_core::NDIM == 2, "RT initialization requires 2D");

        const double pi = std::acos(-1.0);
        const double rho_bottom = 1.0;
        const double rho_top = 2.0;
        const double p_interface = 2.5;
        const double amplitude = 0.01;
        constexpr int mode = 4; // Four wavelengths across Lx; 32 cells each at nxa=128.

        const double x_lower = cfg.lower[0];
        const double y_lower = cfg.lower[1];
        const double Lx = cfg.upper[0] - x_lower;
        const double Ly = cfg.upper[1] - y_lower;
        const double y_interface = y_lower + 0.5 * Ly;

        for (int j = 0; j < mesh.nactive[1]; ++j)
        {
            for (int i = 0; i < mesh.nactive[0]; ++i)
            {
                const hd_core::Index index{
                    i + mesh.nghost[0],
                    j + mesh.nghost[1]
                };

                const double x = mesh.cell_center(0, index);
                const double y = mesh.cell_center(1, index);

                const double rho =
                    y < y_interface ? rho_bottom : rho_top;

                const double envelope =
                    std::sin(pi * (y - y_lower) / Ly);

                hd_core::Primitives W{};
                W[hd_core::RHO] = rho;
                W[hd_core::PRE] =
                    p_interface + rho * cfg.gravity[1] * (y - y_interface);

                W[hd_core::VEL(0)] = 0.0;
                W[hd_core::VEL(1)] =
                    amplitude
                    * std::cos(2.0 * pi * mode * (x - x_lower) / Lx)
                    * envelope * envelope;

                mesh.at(index) = hd_core::Cell(W, cfg.gamma);
            }
        }
    }
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize_RTI2D;
    return hooks;
}
