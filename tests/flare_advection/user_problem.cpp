/**
 * @file user_problem.cpp
 * @brief Translate a cold, FLARE-shaped gas patch across periodic boundaries.
 * @details The background is exact vacuum: rho=p=0, with zero stored velocity.
 * Gas has rho=1, p=1e-8 and velocity (1,0.5). Positive-pressure gas expands into
 * vacuum, so rigid translation is a cold-gas approximation, not an exact Euler
 * solution. With gamma=1.4, the initial vacuum-front expansion speed relative
 * to the bulk is 2*c/(gamma-1) ~= 5.92e-4, well below one cell over t=1.
 *
 * In this directory: perl ../../build.pl; make; make test; make animate.
 */
#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <array>
#include <cmath>
#include <string_view>

static_assert(hd_core::NDIM == 2, "FLARE advection requires HD_NDIM=2");

namespace
{
    constexpr double gas_density = 1.0;
    constexpr double gas_pressure = 1e-8;
    constexpr double velocity_x = 1.0;
    constexpr double velocity_y = 0.5;

    // Five 5x7 letters; rows are written from top to bottom.
    constexpr std::array<std::array<std::string_view, 7>, 5> letters{{
        {{"11111", "10000", "10000", "11110", "10000", "10000", "10000"}},
        {{"10000", "10000", "10000", "10000", "10000", "10000", "11111"}},
        {{"01110", "10001", "10001", "11111", "10001", "10001", "10001"}},
        {{"11110", "10001", "10001", "11110", "10100", "10010", "10001"}},
        {{"11111", "10000", "10000", "11110", "10000", "10000", "11111"}}
    }};

    /** @brief Test a cell center in normalized domain coordinates against FLARE. */
    bool inside_letters(double x, double y)
    {
        // A 29x7 bitmap fills 87% of the width and 42% of the height.
        // These pixels are square on the default 2:1 physical domain.
        constexpr double left = 0.065, bottom = 0.29;
        constexpr double pixel_x = 0.03, pixel_y = 0.06;
        const int column = static_cast<int>(std::floor((x-left)/pixel_x));
        const int row_from_bottom = static_cast<int>(std::floor((y-bottom)/pixel_y));
        if (column < 0 || column >= 29 || row_from_bottom < 0 || row_from_bottom >= 7)
            return false;
        const int letter = column / 6;
        const int pixel = column % 6;
        return pixel < 5 && letters[letter][6-row_from_bottom][pixel] == '1';
    }

    /**
     * @brief Fill active cells; the framework subsequently fills periodic ghosts.
     * @param[in,out] mesh Allocated mesh whose active states are initialized.
     * @param[in] cfg Domain bounds and gas gamma.
     */
    void initialize_flare(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        for (int j=0; j<mesh.nactive[1]; ++j)
        {
            for (int i=0; i<mesh.nactive[0]; ++i)
            {
                const hd_core::Index index{i+mesh.nghost[0], j+mesh.nghost[1]};
                const double x = (mesh.cell_center(0,index)-cfg.lower[0]) / mesh.active_length(0);
                const double y = (mesh.cell_center(1,index)-cfg.lower[1]) / mesh.active_length(1);
                hd_core::Primitives w{};
                if (inside_letters(x,y))
                {
                    w[hd_core::RHO] = gas_density;
                    w[hd_core::PRE] = gas_pressure;
                    w[hd_core::VEL(0)] = velocity_x;
                    w[hd_core::VEL(1)] = velocity_y;
                }
                mesh.at(index) = hd_core::Cell(w,cfg.gamma);
            }
        }
    }
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize_flare;
    return hooks;
}
