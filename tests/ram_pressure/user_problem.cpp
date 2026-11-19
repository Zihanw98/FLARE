#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <algorithm>
#include <stdexcept>

static_assert(hd_core::NDIM == 2, "ram_pressure requires a Cartesian 2D mesh");

namespace
{
    // Same bulk-parameter convention as the user's AMRVAC ram_pressure case.
    // chi = rho_inside/rho_outside; the circle is ten times less dense.
    // wind_mach uses the ambient (outside) sound speed. Code units throughout.
    constexpr double bubble_chi = 0.1;
    constexpr double bubble_radius = 1.0;
    constexpr double bubble_x = 0.0;
    constexpr double bubble_y = 0.0;
    constexpr double wind_density = 1.0;
    constexpr double wind_speed = 1.0;
    constexpr double wind_mach = 0.8;

    hd_core::Primitives wind_state(double gamma)
    {
        hd_core::Primitives w{};
        w[hd_core::RHO] = wind_density;
        const double cs = wind_speed / wind_mach;
        w[hd_core::PRE] = wind_density * cs*cs / gamma;
        w[hd_core::VEL(0)] = wind_speed;
        return w;
    }

    void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        if (cfg.geometry != "Cartesian" || cfg.constant_gravity_enabled ||
            cfg.spherical_gravity_enabled || cfg.spherical_self_gravity_enabled ||
            cfg.constant_heating_enabled)
            throw std::invalid_argument("ram_pressure requires Cartesian geometry without gravity or heating");
        if (bubble_x - bubble_radius <= cfg.lower[0] || bubble_x + bubble_radius >= cfg.upper[0] ||
            bubble_y - bubble_radius <= cfg.lower[1] || bubble_y + bubble_radius >= cfg.upper[1])
            throw std::invalid_argument("The entire cloud must lie inside the domain");

        for (int j = 0; j < mesh.nactive[1]; ++j)
            for (int i = 0; i < mesh.nactive[0]; ++i)
            {
                const hd_core::Index index{i + mesh.nghost[0], j + mesh.nghost[1]};
                const double dx = mesh.cell_center(0, index) - bubble_x;
                const double dy = mesh.cell_center(1, index) - bubble_y;
                auto w = wind_state(cfg.gamma);
                if (dx*dx + dy*dy < bubble_radius*bubble_radius)
                {
                    w[hd_core::RHO] = bubble_chi * wind_density;
                    w[hd_core::VEL(0)] = 0.0;
                }
                // Construct total energy separately from the common pressure:
                // E = p/(gamma-1) + rho*|v|^2/2, including wind kinetic energy.
                mesh.at(index) = hd_core::Cell(w, cfg.gamma);
            }
    }

    void fill_boundary(hd_core::MeshBlock& mesh, double, int direction,
                       boundary::BoundarySide side)
    {
        const bool lower = side == boundary::BoundarySide::Lower;
        const int first = mesh.nghost[direction];
        const int end = first + mesh.nactive[direction];
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            const auto index = mesh.coordinates(flat);
            if (!(lower ? index[direction] < first : index[direction] >= end))
                continue;
            auto source = index;
            for (int d = 0; d < hd_core::NDIM; ++d)
                source[d] = std::clamp(index[d], mesh.nghost[d],
                                       mesh.nghost[d] + mesh.nactive[d] - 1);
            // Maintain the wind across the full left boundary, including corners
            // subsequently visited by the y-boundary calls. Else copy the nearest
            // active cell for zero-gradient outflow.
            if (index[0] < mesh.nghost[0])
                mesh.cells[flat] = hd_core::Cell(wind_state(mesh.at(source).gamma), mesh.at(source).gamma);
            else
                mesh.cells[flat] = mesh.at(source);
        }
    }
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    hooks.userbc = fill_boundary;
    return hooks;
}
