#include "user_config.h"  // Load this problem's settings before framework headers.
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <stdexcept>
#include <math.h>
#include "constants/constants.h"
static_assert(hd_core::NDIM == 2, "KH rquires 2D setting");

namespace
{
    /**
     * @brief Initialize the active cells of the standard Sod shock tube.
     * @param[in,out] mesh Allocated one-dimensional mesh; ghosts are filled by userbc.
     * @param[in] cfg Supplies the gas gamma; the discontinuity is at x=0.5.
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
                    * std::cos(2.0 * pi * (x - x_lower) / Lx)
                    * envelope * envelope;

                mesh.at(index) = hd_core::Cell(W, cfg.gamma);
            }
        }
    }

    /**
     * @brief Fill every ghost layer by copying the closest active cell (zero gradient).
     * @param[in,out] mesh Copies both U and gamma; active cells remain unchanged.
     * @param[in] time Current time; this boundary is time independent.
     * @param[in] direction Must be 0 for this one-dimensional problem.
     * @param[in] side Lower or upper end of the tube.
     */
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize_RTI2D;
    return hooks;
}
