#include "user_config.h"  // Load this problem's settings before framework headers.
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <stdexcept>

static_assert(hd_core::NDIM == 1, "The Sod tube example requires HD_NDIM=1");

namespace
{
    /**
     * @brief Initialize the active cells of the standard Sod shock tube.
     * @param[in,out] mesh Allocated one-dimensional mesh; ghosts are filled by userbc.
     * @param[in] cfg Supplies the gas gamma; the discontinuity is at x=0.5.
     */
    void initialize_sod(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        for (int i = 0; i < mesh.nactive[0]; ++i)
        {
            const hd_core::Index index{i + mesh.nghost[0]};
            const bool left = mesh.cell_center(0, index) < 0.5;
            hd_core::Primitives W{};
            W[hd_core::RHO] = left ? 1.0 : 0.125;
            W[hd_core::VEL(0)] = 0.0;
            W[hd_core::PRE] = left ? 1.0 : 0.1;
            mesh.at(index) = hd_core::Cell(W, cfg.gamma);
        }
    }

    /**
     * @brief Fill every ghost layer by copying the closest active cell (zero gradient).
     * @param[in,out] mesh Copies both U and gamma; active cells remain unchanged.
     * @param[in] time Current time; this boundary is time independent.
     * @param[in] direction Must be 0 for this one-dimensional problem.
     * @param[in] side Lower or upper end of the tube.
     */
    void sod_boundary(hd_core::MeshBlock& mesh, double /*time*/, int direction,
                      boundary::BoundarySide side)
    {
        if (direction != 0)
            throw std::invalid_argument("Sod tube only has direction 0");
        const int first = mesh.nghost[0];
        const int end = first + mesh.nactive[0];
        if (side == boundary::BoundarySide::Lower)
        {
            for (int i = 0; i < first; ++i)
                mesh.at({i}) = mesh.at({first});
        }
        else if (side == boundary::BoundarySide::Upper)
        {
            for (int i = end; i < mesh.ntotal(0); ++i)
                mesh.at({i}) = mesh.at({end - 1});
        }
        else
        {
            throw std::invalid_argument("Invalid Sod boundary side");
        }
    }
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize_sod;
    hooks.userbc = sod_boundary;
    return hooks;
}
