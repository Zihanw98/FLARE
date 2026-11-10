#include "user_config.h"  // Load this problem's settings before framework headers.
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <algorithm>
#include <stdexcept>
static_assert(hd_core::NDIM == 2, "This test-suite requires 2D setting");

namespace
{
    /**
     * @brief Initialize a denser, higher-pressure patch in the central third of each axis.
     * @param[in,out] mesh Allocated two-dimensional mesh; ghosts are filled by userbc.
     * @param[in] cfg Supplies the gas gamma.
     */
    void initialize_expansion_2D(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        const double x_lower = mesh.face[0][mesh.nghost[0]];
        const double y_lower = mesh.face[1][mesh.nghost[1]];
        const double x_length = mesh.active_length(0);
        const double y_length = mesh.active_length(1);
        for(int y=0; y < mesh.nactive[1]; y++)
        {
            for (int x=0; x < mesh.nactive[0]; x++)
            {
                const hd_core::Index index{x+mesh.nghost[0], y+mesh.nghost[1]};
                const double x_c = mesh.cell_center(0,index);
                const double y_c = mesh.cell_center(1,index);
                const bool middle = y_c > y_lower + y_length / 3.0 && y_c < y_lower + 2.0 * y_length / 3.0
                                 && x_c > x_lower + x_length / 3.0 && x_c < x_lower + 2.0 * x_length / 3.0;
                hd_core::Primitives W{};
                W[hd_core::RHO] = middle? 2.0 : 1.0;
                W[hd_core::VEL(0)] = 0;
                W[hd_core::VEL(1)] = 0;
                W[hd_core::PRE] = middle? 1.0 : 0.2;
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
     * @param[in] direction 0 for x, 1 for y.
     * @param[in] side Lower or upper side of the selected direction.
     */
    void expansion_outflow_boundary(hd_core::MeshBlock& mesh, double /*time*/,
                                    int direction, boundary::BoundarySide side)
    {
        if (direction < 0 || direction >= hd_core::NDIM)
            throw std::invalid_argument("Invalid expansion boundary direction");
        if (side != boundary::BoundarySide::Lower && side != boundary::BoundarySide::Upper)
            throw std::invalid_argument("Invalid expansion boundary side");

        const int first = mesh.nghost[direction];
        const int end = first + mesh.nactive[direction];
        const bool lower = side == boundary::BoundarySide::Lower;
        const int begin_ghost = lower ? 0 : end;
        const int end_ghost = lower ? first : mesh.ntotal(direction);
        const int transverse = 1 - direction;
        const int transverse_first = mesh.nghost[transverse];
        const int transverse_last = transverse_first + mesh.nactive[transverse] - 1;

        hd_core::Index source{}, destination{};
        source[direction] = lower ? first : end - 1;
        for (int j = 0; j < mesh.ntotal(transverse); ++j)
        {
            destination[transverse] = j;
            // Corners copy an active corner directly, independent of fill order.
            source[transverse] = std::clamp(j, transverse_first, transverse_last);
            for (int i = begin_ghost; i < end_ghost; ++i)
            {
                destination[direction] = i;
                mesh.at(destination) = mesh.at(source);
            }
        }
    }
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize_expansion_2D;
    hooks.userbc = expansion_outflow_boundary;
    return hooks;
}
