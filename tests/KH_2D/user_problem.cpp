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
    const double A = 0.01;
    const double sigma = 0.05;
    void initialize_KH2D(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        for(int y=0; y < mesh.nactive[1]; y++)
        {
            for (int x=0; x < mesh.nactive[0]; x++)
            {
                const hd_core::Index index{x+mesh.nghost[0], y+mesh.nghost[1]};
                const double x_c = mesh.cell_center(0,index);
                const double y_c = mesh.cell_center(1,index);
                const double d_1 = y_c - mesh.active_length(1) * 1/4.0 ;
                const double d_2 = y_c - mesh.active_length(1) * 3/4.0;
                const bool middle = ( y_c > mesh.active_length(1) * 1/4.0 && y_c < mesh.active_length(1) * 3/4.0);
                hd_core::Primitives W{};
                W[hd_core::RHO] = middle? 2 : 1;
                W[hd_core::VEL(0)] = middle? 0.5 : -0.5;
                W[hd_core::VEL(1)] = A * std::sin(4*hd_const::pi * x_c) * (std::exp(-(d_1*d_1)/(2*sigma*sigma))+std::exp(-(d_2*d_2)/(2*sigma*sigma)));
                W[hd_core::PRE] = 2.5;
                mesh.at(index) = hd_core::Cell(W, cfg.gamma);

            }
        }
        // for (int i = 0; i < mesh.nactive[0]; i++)
        // {
            
        //     for (int j =0; j < mesh.nactive[1]; j++)
        //     const hd_core::Index index{i + mesh.nghost[0]};
        //     const bool left = mesh.cell_center(0, index) < 0.5;
        //     hd_core::Primitives W{};
        //     W[hd_core::RHO] = left ? 1.0 : 0.125;
        //     W[hd_core::VEL(0)] = 0.0;
        //     W[hd_core::PRE] = left ? 1.0 : 0.1;
        //     mesh.at(index) = hd_core::Cell(W, cfg.gamma);
        // }
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
    hooks.initialize = initialize_KH2D;
    return hooks;
}
