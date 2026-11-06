#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <stdexcept>

static_assert(hd_core::NDIM == 1, "Spherical Sod requires HD_NDIM=1");

namespace
{
    void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        if (mesh.geometry != hd_core::Geometry::Spherical)
            throw std::invalid_argument("Spherical Sod requires spherical geometry");

        // A spherical discontinuity at r=0.5, with zero initial velocity.
        for (int i = 0; i < mesh.nactive[0]; ++i)
        {
            const hd_core::Index index{i + mesh.nghost[0]};
            const bool inner = mesh.cell_center(0, index) < 0.5;
            hd_core::Primitives W{};
            W[hd_core::RHO] = inner ? 1.0 : 0.125;
            W[hd_core::PRE] = inner ? 1.0 : 0.1;
            mesh.at(index) = hd_core::Cell(W, cfg.gamma);
        }
    }

    void outer_boundary(hd_core::MeshBlock& mesh, double, int direction,
                        boundary::BoundarySide side)
    {
        if (direction != 0 || side != boundary::BoundarySide::Upper)
            throw std::invalid_argument("Only the outer radial boundary uses userbc");
        const int end = mesh.nghost[0] + mesh.nactive[0];
        for (int i = end; i < mesh.ntotal(0); ++i)
            mesh.at({i}) = mesh.at({end - 1});
    }
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    hooks.userbc = outer_boundary;
    return hooks;
}
