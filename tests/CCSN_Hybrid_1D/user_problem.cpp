#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include "microphysics/eos/nuclear.h"
#include "physcon/constants.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

static_assert(hd_core::NDIM == 1 && hd_core::NCOMP == 0);

namespace
{
void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    using namespace hd_core;
    using namespace hd_microphysics;
    if (mesh.geometry != Geometry::Spherical || cfg.lower[0] != 0.0
        || cfg.upper[0] > 1e10 || cfg.upper[0] < 3e8 || cfg.nactive[0] < 128
        || std::abs(cfg.gamma - 4.0/3.0) > 1e-12 || eos_type != Equation_of_State::Nuclear
        || !cfg.spherical_gravity_enabled || !cfg.spherical_self_gravity_enabled
        || cfg.spherical_central_mass != 0.0 || cfg.riemann_solver != "hll"
        || cfg.constant_gravity_enabled || cfg.constant_heating_enabled
        || mesh.bc_lower[0] != boundary::BoundaryCondition::Reflective
        || mesh.bc_upper[0] != boundary::BoundaryCondition::User)
        throw std::invalid_argument("CCSN toy requires Nuclear/HLL, gamma=4/3, spherical self gravity, no central mass/heating, reflective/outflow boundaries and 3e8<=rmax<=1e10 cm");

    // Fine at r=0, gradually stretched outward; no change to framework mesh rules.
    constexpr double stretch_radius = 2e6;
    const double stretch = std::asinh(cfg.upper[0] / stretch_radius);
    for (int i = 0; i <= mesh.ntotal(0); ++i)
        mesh.face[0][i] = stretch_radius * std::sinh(stretch
            * (i - mesh.nghost[0]) / mesh.nactive[0]);
    mesh.allocate();
    mesh.build_geometry();

    // n=3 Lane-Emden density profile for K_structure. Remove pressure support
    // by initializing with the smaller nuclear_K_cgs, without an energy bomb.
    constexpr double rho_c = 1e12; // Late-collapse toy core, still 200x below nuclear density.
    constexpr double K_structure = 5e14;
    constexpr double thermal_fraction = 0.01; // P_th / P_cold initially.
    constexpr double dxi = 1e-3;
    using Pair = std::array<double, 2>; // theta, -xi^2 theta'. Local initial data only.
    std::vector<Pair> solution{{1,0}, {1-dxi*dxi/6, dxi*dxi*dxi/3}};
    const auto derivative = [](double x, Pair y) -> Pair {
        return {-y[1]/(x*x), x*x*std::pow(std::max(0.0,y[0]),3)};
    };
    const auto add = [](Pair y, Pair k, double h) -> Pair {
        return {y[0]+h*k[0], y[1]+h*k[1]};
    };
    for (int i = 1; i < 10000 && solution.back()[0] > 0; ++i)
    {
        const double x = i*dxi;
        const auto y = solution.back();
        const auto k1 = derivative(x,y);
        const auto k2 = derivative(x+dxi/2,add(y,k1,dxi/2));
        const auto k3 = derivative(x+dxi/2,add(y,k2,dxi/2));
        const auto k4 = derivative(x+dxi,add(y,k3,dxi));
        solution.push_back({y[0]+dxi*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,
                            y[1]+dxi*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6});
    }
    if (solution.back()[0] > 0) throw std::runtime_error("Lane-Emden surface not found");
    const double a = std::sqrt(K_structure/(physcon::pi*physcon::G)) * std::pow(rho_c,-1.0/3);
    const auto density = [&](double r) {
        const double q = r / (a*dxi);
        double theta = 0;
        if (q < solution.size()-1)
        {
            const auto j = static_cast<std::size_t>(q);
            theta = std::max(0.0, solution[j][0]*(1-(q-j))+solution[j+1][0]*(q-j));
        }
        const double atmosphere = 1e3 / (1+std::pow(r/1e8,4));
        return std::max(rho_c*theta*theta*theta, atmosphere);
    };

    std::filesystem::create_directories(cfg.output_directory);
    std::ofstream grid(cfg.output_directory+"/grid.csv");
    grid << std::setprecision(17) << "r_left,r_right,volume\n";
    double mass = 0;
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        auto& cell = mesh.at_active({i});
        cell.gamma = cfg.gamma;
        cell.mean_molecular_weight = 1.0;
        const auto [l,r] = mesh.cell_face_pos(0, i+mesh.nghost[0]);
        double integral = 0;
        constexpr int nq = 16;
        for (int j = 0; j <= nq; ++j)
        {
            const double x = l+(r-l)*j/nq;
            integral += (j==0 || j==nq ? 1 : j%2 ? 4 : 2)*density(x)*x*x;
        }
        const double rho = integral / (nq*(l*l+l*r+r*r));
        double pc, uc;
        nuclear_cold_state(rho,cfg.gamma,pc,uc);
        Primitives w{};
        w[RHO] = rho;
        w[PRE] = pc*(1+thermal_fraction);
        cell.Set_Conserved_From_Primitives(w);
        const auto flat = mesh.flat_index({i},true);
        const double volume = mesh.volume_cell[flat];
        mass += volume*rho;
        grid << l << ',' << r << ',' << volume << '\n';
    }
    if (!grid) throw std::runtime_error("Cannot write toy grid");
    std::ofstream metadata(cfg.output_directory+"/model.json");
    metadata << std::setprecision(17)
        << "{\n\"rho_nuclear\":" << nuclear_density_cgs
        << ",\n\"gamma_low\":" << cfg.gamma
        << ",\n\"gamma_stiff\":" << nuclear_gamma_stiff
        << ",\n\"gamma_thermal\":" << nuclear_gamma_thermal
        << ",\n\"K_cold\":" << nuclear_K_cgs
        << ",\n\"K_structure\":" << K_structure
        << ",\n\"initial_thermal_fraction\":" << thermal_fraction
        << ",\n\"initial_central_density\":" << rho_c
        << ",\n\"initial_mass\":" << mass
        << ",\n\"initial_core_radius\":" << a*(solution.size()-1)*dxi
        << ",\n\"outer_radius\":" << cfg.upper[0] << "\n}\n";
    if (!metadata) throw std::runtime_error("Cannot write toy model metadata");
    std::cout << "Hybrid collapse: M=" << mass/1.98847e33 << " Msun, R_core="
              << a*(solution.size()-1)*dxi << " cm, rho_c=" << rho_c
              << ", initial pressure support=" << nuclear_K_cgs*(1+thermal_fraction)/K_structure
              << " of the n=3 reference.\n";
}

void outer_boundary(hd_core::MeshBlock& mesh, double, int direction,
                    boundary::BoundarySide side)
{
    if (direction != 0 || side != boundary::BoundarySide::Upper)
        throw std::invalid_argument("Toy user boundary is outer radial outflow only");
    const int end = mesh.nghost[0]+mesh.nactive[0];
    for (int i = end; i < mesh.ntotal(0); ++i)
        mesh.at({i}) = mesh.at({end-1});
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    hooks.userbc = outer_boundary;
    return hooks;
}
