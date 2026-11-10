#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include "physcon/constants.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

static_assert(hd_core::NDIM == 1, "Stellar equilibrium requires spherical 1D");

namespace
{
// IAU 2015 B3 nominal R_sun and GM_sun, converted to cgs.
constexpr double solar_radius = 6.957e10;
constexpr double solar_gm = 1.3271244e26;
// Stellar parameters in the unchanged solar mass and radius units.
constexpr double stellar_mass = 10.0;
constexpr double stellar_radius = 1000.0;
constexpr double injection_energy_erg = 1.6e51;
static_assert(stellar_mass > 0.0 && stellar_radius > 0.0);
// Scale the reference atmosphere with the star: rho ~ M/R^3, p/rho ~ M/R.
constexpr double atmosphere_density =
    1e-8 * stellar_mass / (stellar_radius * stellar_radius * stellar_radius);
constexpr double atmosphere_p_over_rho = 0.02 * stellar_mass / stellar_radius;

struct EmdenState { double theta, mass; }; // mass = -xi^2 theta'.
EmdenState add(EmdenState a, EmdenState b, double h)
{
    return {a.theta + h*b.theta, a.mass + h*b.mass};
}
EmdenState derivative(double x, EmdenState y)
{
    return {-y.mass/(x*x), x*x*std::pow(std::max(y.theta, 0.0), 3)};
}

struct Model
{
    static constexpr double dxi = 1e-4;
    std::vector<EmdenState> emden;
    std::vector<double> pressure;
    std::vector<double> atmosphere_mass;
    double xi1 = 0.0, mu1 = 0.0, rho_c = 0.0, dr = 0.0;

    void build(double outer_radius)
    {
        emden = {{1.0, 0.0},
                 {1.0-dxi*dxi/6.0, dxi*dxi*dxi/3.0}};
        for (int i=1; i<100000; ++i)
        {
            const double x=i*dxi;
            const auto y=emden.back();
            const auto k1=derivative(x,y);
            const auto k2=derivative(x+dxi/2,add(y,k1,dxi/2));
            const auto k3=derivative(x+dxi/2,add(y,k2,dxi/2));
            const auto k4=derivative(x+dxi,add(y,k3,dxi));
            const EmdenState next{
                y.theta+dxi*(k1.theta+2*k2.theta+2*k3.theta+k4.theta)/6,
                y.mass+dxi*(k1.mass+2*k2.mass+2*k3.mass+k4.mass)/6};
            emden.push_back(next);
            if (next.theta<=0)
            {
                const double f=y.theta/(y.theta-next.theta);
                xi1=x+f*dxi;
                mu1=y.mass+f*(next.mass-y.mass);
                break;
            }
        }
        if (xi1<=0 || mu1<=0)
            throw std::runtime_error("Lane-Emden surface not found");
        rho_c=stellar_mass/std::pow(stellar_radius,3)
            *std::pow(xi1,3)/(4*physcon::pi*mu1);

        // Integrate dp/dr = -rho G M/r^2 inward, including the tenuous atmosphere.
        // Boundary pressure is applied beyond the outermost ghost cell.
        constexpr int samples=65536;
        dr=outer_radius/samples;
        atmosphere_mass.assign(samples+1,0.0);
        const auto integrand=[&](double r) { return 4*physcon::pi*r*r*atmosphere(r); };
        for(int i=0; i<samples; ++i)
        {
            const double r=i*dr;
            atmosphere_mass[i+1]=atmosphere_mass[i]+dr/6.0*
                (integrand(r)+4*integrand(r+dr/2)+integrand(r+dr));
        }
        pressure.assign(samples+1, atmosphere(outer_radius)*atmosphere_p_over_rho);
        for(int i=samples-1; i>=0; --i)
        {
            const double r=i*dr;
            pressure[i]=pressure[i+1]+dr/6.0*
                (weight(r)+4*weight(r+dr/2)+weight(r+dr));
        }
    }

    EmdenState at(double r) const
    {
        if(r>=stellar_radius) return {0.0,mu1};
        const double x=std::max(r,0.0)/stellar_radius*xi1;
        if(x<dxi)
            return {1-x*x/6+std::pow(x,4)/40,
                    x*x*x/3-std::pow(x,5)/10};
        const double q=x/dxi;
        const auto i=static_cast<std::size_t>(q);
        const double f=q-i;
        return {std::max(0.0,emden[i].theta*(1-f)+emden[i+1].theta*f),
                emden[i].mass*(1-f)+emden[i+1].mass*f};
    }
    double density(double r) const
    {
        return rho_c*std::pow(at(r).theta,3)+atmosphere(r);
    }
    double atmosphere(double r) const
    {
        // Outside the stellar surface, isothermal atmosphere in the stellar monopole.
        // Its tiny self gravity is included by the subsequent pressure integration.
        return atmosphere_density*std::exp(physcon::G*stellar_mass*
            (1.0/std::max(r,stellar_radius)-1.0/stellar_radius)/atmosphere_p_over_rho);
    }
    double interpolate(const std::vector<double>& values, double r) const
    {
        const double q=std::clamp(r/dr,0.0,static_cast<double>(values.size()-1));
        const auto i=std::min(static_cast<std::size_t>(q),values.size()-2);
        const double f=q-i;
        return values[i]*(1-f)+values[i+1]*f;
    }
    double mass(double r) const
    {
        return stellar_mass*at(r).mass/mu1+interpolate(atmosphere_mass,r);
    }
    double weight(double r) const
    {
        return r==0.0 ? 0.0 : density(r)*physcon::G*mass(r)/(r*r);
    }
    double p(double r) const
    {
        return interpolate(pressure,r);
    }
    hd_core::Cell average(double a, double b, double gamma) const
    {
        // Simpson quadrature of conserved densities with the spherical measure.
        constexpr int intervals=32;
        double rho_int=0.0, p_int=0.0;
        const double h=(b-a)/intervals;
        for(int q=0; q<=intervals; ++q)
        {
            const double r=a+q*h;
            const double w=(q==0 || q==intervals) ? 1.0 : (q%2 ? 4.0 : 2.0);
            rho_int+=w*density(r)*r*r;
            p_int+=w*p(r)*r*r;
        }
        const double radial_volume=(b-a)*(a*a+a*b+b*b)/3.0;
        hd_core::Primitives W{};
        W[hd_core::RHO]=rho_int*h/(3*radial_volume);
        W[hd_core::PRE]=p_int*h/(3*radial_volume);
        return hd_core::Cell(W,gamma);
    }
} model;

void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
{
    // Keep solar units while stellar_mass and stellar_radius scale the profile.
    units.init_unitsystem_MLT(solar_gm/physcon::G_, solar_radius,
        std::sqrt(std::pow(solar_radius,3)/solar_gm),1.0);
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (mesh.geometry!=hd_core::Geometry::Spherical || cfg.lower[0]!=0.0 ||
        cfg.upper[0]<1.05*stellar_radius || !cfg.spherical_gravity_enabled ||
        !cfg.spherical_self_gravity_enabled || !std::isfinite(cfg.spherical_central_mass) ||
        cfg.spherical_central_mass<=0.0 ||
        cfg.spherical_central_mass>=stellar_mass ||
        cfg.constant_gravity_enabled || cfg.constant_heating_enabled ||
        mesh.bc_lower[0]!=boundary::BoundaryCondition::Reflective ||
        mesh.bc_upper[0]!=boundary::BoundaryCondition::User ||
        std::abs(cfg.gamma-5.0/3.0)>1e-12)
        throw std::invalid_argument("Mass-cut model requires x_lower=0 (automatic cut), x_upper>=1.05 stellar radii, 0<central_mass<stellar_mass, self gravity, gamma=5/3, reflective/user boundaries, and no additional heating or force");

    model.build(mesh.face[0].back());

    // Invert the original enclosed mass, including its tiny atmosphere.
    // The excised mass is represented once, by the configured central point mass.
    double lo=0.0, hi=stellar_radius;
    for(int iteration=0; iteration<100; ++iteration)
    {
        const double mid=0.5*(lo+hi);
        if(model.mass(mid)<cfg.spherical_central_mass) lo=mid;
        else hi=mid;
    }
    const double cut_radius=0.5*(lo+hi);
    auto shell_cfg=cfg;
    shell_cfg.lower[0]=cut_radius;
    mesh.build_mesh(shell_cfg);
    for (int i=mesh.nghost[0]; i<mesh.ntotal(0); ++i)
    {
        const auto [a,b]=mesh.cell_face_pos(0,i);
        mesh.at({i})=model.average(a,b,cfg.gamma);
    }

    std::filesystem::create_directories(cfg.output_directory);
    // Save the unheated envelope and exact shell volumes for energy/mass checks.
    std::ofstream profile(std::filesystem::path(cfg.output_directory)/"initial_envelope.csv");
    profile<<std::setprecision(17)<<"r_left,r_right,volume,rho,p,energy_density_before\n";
    double envelope_mass=0.0, envelope_energy=0.0;
    for(int i=0; i<mesh.nactive[0]; ++i)
    {
        const auto flat=mesh.flat_index({i},true);
        const auto [a,b]=mesh.cell_face_pos(0,i+mesh.nghost[0]);
        const auto& cell=mesh.cells[flat];
        const auto W=cell.U.to_Primtivies(cell.gamma);
        const double volume=mesh.volume_cell[flat];
        envelope_mass+=W[hd_core::RHO]*volume;
        envelope_energy+=cell.U[hd_core::ENERGY]*volume;
        profile<<a<<','<<b<<','<<volume<<','<<W[hd_core::RHO]<<','
               <<W[hd_core::PRE]<<','<<cell.U[hd_core::ENERGY]<<'\n';
    }
    profile.close();
    if(!profile) throw std::runtime_error("Cannot write initial envelope profile");

    // Instantaneous thermal bomb in exactly the first active shell.
    // U stores energy per volume; neither density nor momentum is changed.
    const auto injection_flat=mesh.flat_index({0},true);
    const double injected_energy=injection_energy_erg/physcon::units.unit_energy;
    mesh.cells[injection_flat].U[hd_core::ENERGY]+=
        injected_energy/mesh.volume_cell[injection_flat];
    mesh.cells[injection_flat].U[hd_core::INTERNAL_ENERGY] +=
        injected_energy/mesh.volume_cell[injection_flat];
    mesh.validate_active();

    std::ofstream info(std::filesystem::path(cfg.output_directory)/"model.json");
    info<<std::setprecision(17)
        <<"{\n  \"mass_unit_g\": "<<physcon::units.unit_mass
        <<",\n  \"length_unit_cm\": "<<physcon::units.unit_length
        <<",\n  \"time_unit_s\": "<<physcon::units.unit_time
        <<",\n  \"velocity_unit_cm_s\": "<<physcon::units.unit_velocity
        <<",\n  \"density_unit_g_cm3\": "<<physcon::units.unit_density
        <<",\n  \"pressure_unit_dyn_cm2\": "<<physcon::units.unit_pressure
        <<",\n  \"energy_unit_erg\": "<<physcon::units.unit_energy
        <<",\n  \"G_code\": "<<physcon::G
        <<",\n  \"stellar_mass_code\": "<<stellar_mass
        <<",\n  \"stellar_radius_code\": "<<stellar_radius
        <<",\n  \"core_mass_code\": "<<cfg.spherical_central_mass
        <<",\n  \"cut_radius_code\": "<<cut_radius
        <<",\n  \"enclosed_mass_at_cut_code\": "<<model.mass(cut_radius)
        <<",\n  \"outer_radius_code\": "<<cfg.upper[0]
        <<",\n  \"envelope_mass_code\": "<<envelope_mass
        <<",\n  \"envelope_energy_before_code\": "<<envelope_energy
        <<",\n  \"injection_energy_erg\": "<<injection_energy_erg
        <<",\n  \"injection_energy_code\": "<<injected_energy
        <<",\n  \"injection_cell_count\": 1"
        <<",\n  \"injection_outer_radius_code\": "
        <<mesh.cell_face_pos(0,mesh.nghost[0]).second
        <<",\n  \"injection_mass_code\": "
        <<mesh.cells[injection_flat].U[hd_core::RHO]*mesh.volume_cell[injection_flat]
        <<",\n  \"dynamical_time_code\": "
        <<std::sqrt(std::pow(stellar_radius,3)/(physcon::G*stellar_mass))
        <<",\n  \"xi1\": "<<model.xi1<<",\n  \"mu1\": "<<model.mu1
        <<",\n  \"central_density_code\": "<<model.rho_c
        <<",\n  \"central_pressure_code\": "<<model.p(0)
        <<",\n  \"atmosphere_density_code\": "<<atmosphere_density<<"\n}\n";
    info.close();
    if(!info) throw std::runtime_error("Cannot write stellar model metadata");
    std::cout<<"Stellar n=3 polytrope: M="<<stellar_mass
             <<", R="<<stellar_radius<<", G="<<physcon::G
             <<", time unit="<<physcon::units.unit_time<<" s\n";
    std::cout<<"Excised core mass="<<cfg.spherical_central_mass
             <<", cut radius="<<cut_radius<<", envelope mass="<<envelope_mass
             <<", injected energy="<<injection_energy_erg<<" erg\n";
}

void outer_boundary(hd_core::MeshBlock& mesh, double, int direction,
                    boundary::BoundarySide side)
{
    if(direction!=0 || side!=boundary::BoundarySide::Upper)
        throw std::invalid_argument("Stellar user boundary is the outer radial atmosphere");
    const int end=mesh.nghost[0]+mesh.nactive[0];
    for(int i=end; i<mesh.ntotal(0); ++i)
    {
        // Zero-gradient outflow, with inward velocity suppressed.
        auto W=mesh.at({end-1}).U.to_Primtivies(mesh.at({end-1}).gamma);
        W[hd_core::VEL(0)]=std::max(0.0,W[hd_core::VEL(0)]);
        mesh.at({i})=hd_core::Cell(W,mesh.at({end-1}).gamma);
    }
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize_units=initialize_units;
    hooks.initialize=initialize;
    hooks.userbc=outer_boundary;
    return hooks;
}
