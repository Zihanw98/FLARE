#pragma once
#include "core/mesh.h"
#include "config/config.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace test_case
{
inline void no_gravity(const config::SimulationConfig& cfg)
{
    if (cfg.constant_gravity_enabled || std::abs(cfg.gamma - 1.4) > 1e-14)
        throw std::invalid_argument("This benchmark requires gamma=1.4 and gravity disabled");
}

inline hd_core::Primitives shock_downstream(double mach, double gamma=1.4)
{
    hd_core::Primitives w{};
    const double ratio = (gamma+1)*mach*mach / ((gamma-1)*mach*mach+2);
    w[hd_core::RHO] = ratio;
    w[hd_core::PRE] = 1 + 2*gamma/(gamma+1)*(mach*mach-1);
    w[hd_core::VEL(0)] = mach*std::sqrt(gamma)*(1-1/ratio);
    return w;
}

// Fixed postshock inflow at x_lower; zero-gradient outflow elsewhere.
// Transverse coordinates are clamped so corner filling is order independent.
inline void shock_boundary(hd_core::MeshBlock& mesh, int direction,
                           boundary::BoundarySide side, double mach)
{
    const bool lower = side == boundary::BoundarySide::Lower;
    const int first = mesh.nghost[direction], end = first + mesh.nactive[direction];
    for (std::size_t flat=0; flat<mesh.cells.size(); ++flat)
    {
        const auto index = mesh.coordinates(flat);
        if (!(lower ? index[direction]<first : index[direction]>=end)) continue;
        if (direction==0 && lower)
            mesh.cells[flat] = hd_core::Cell(shock_downstream(mach), 1.4);
        else
        {
            auto source=index;
            for (int d=0; d<hd_core::NDIM; ++d)
                source[d]=std::clamp(index[d], mesh.nghost[d], mesh.nghost[d]+mesh.nactive[d]-1);
            mesh.cells[flat]=mesh.at(source);
        }
    }
}

struct Sod { double rho, u, p; };
inline Sod sod(double coordinate, double time)
{
    if (time==0) return coordinate<0 ? Sod{1,0,1} : Sod{.125,0,.1};
    constexpr double gamma=1.4;
    // Standard Sod star state, from its scalar pressure matching equation.
    constexpr double ps=.30313017805064685, us=.9274526200489499;
    const double cl=std::sqrt(gamma), cr=std::sqrt(gamma*.1/.125);
    const double cs=cl*std::pow(ps,(gamma-1)/(2*gamma));
    const double shock=cr*std::sqrt((gamma+1)/(2*gamma)*ps/.1+(gamma-1)/(2*gamma));
    const double speed=coordinate/time, b=(gamma-1)/(gamma+1);
    if (speed<=-cl) return {1,0,1};
    if (speed<us-cs)
    {
        const double u=2/(gamma+1)*(cl+speed);
        const double c=2/(gamma+1)*(cl-(gamma-1)*speed/2);
        return {std::pow(c/cl,2/(gamma-1)),u,std::pow(c/cl,2*gamma/(gamma-1))};
    }
    if (speed<=us) return {std::pow(ps,1/gamma),us,ps};
    if (speed<shock) return {.125*(ps/.1+b)/(b*ps/.1+1),us,ps};
    return {.125,0,.1};
}
}
