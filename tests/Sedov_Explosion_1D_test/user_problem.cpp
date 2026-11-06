#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <stdexcept>
#include <iostream>

static_assert(hd_core::NDIM == 1, "Spherical Sod requires HD_NDIM=1");

namespace
{
    void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        double r_inj = 1/64.0;
        double V_inj = 0.0;
        if (mesh.geometry != hd_core::Geometry::Spherical)
            throw std::invalid_argument("Spherical Sod requires spherical geometry");

        // A spherical discontinuity at r=0.5, with zero initial velocity.
        for (int i = 0; i < mesh.nactive[0]; ++i)
        {
            const hd_core::Index index{i + mesh.nghost[0]};
            const auto flat = mesh.flat_index(index);
            const bool inner = mesh.cell_center(0, index) < r_inj;
            if(inner)
            {
                V_inj += mesh.volume_cell[flat];
            }
        }
        const double E_0 = 1;
        const double gamma = cfg.gamma;
        std::cout<<V_inj<<std::endl;
        if (V_inj <= 0.0)
            throw std::invalid_argument("Injection region contains no active cells");
        for (int i = 0; i < mesh.nactive[0]; ++i)
        {
            const hd_core::Index index{i + mesh.nghost[0]};
            const auto flat = mesh.flat_index(index);
            const bool inner = mesh.cell_center(0, index) < r_inj;
            hd_core::Conserved U{};
            U[hd_core::RHO] = 1.0;
            U[hd_core::ENERGY] =
                (inner ? E_0 / V_inj : 0.0) + 1e-5 / (gamma - 1.0);
            auto W = U.to_Primtivies(gamma);
            mesh.at(index).U = U;
            mesh.at(index).gamma = gamma;
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
