// Standalone diagnostic main; uses the production solver without changing it.
// The runner copies and parameterizes user_problem.cpp in a separate directory.
#include "user_problem.cpp"
#include "hydro/riemann/factory.h"
#include "hydro/spatial.h"
#include "time/integrator.h"
#include "boundary/boundary.h"
#include "source/source.h"
#include "io/output.h"
#include <array>
#include <limits>

struct Diagnostics
{
    double mass=0, kinetic=0, internal=0, gravity=0;
    double mass_flux=0, energy_flux=0, vmax=0, vmax_r=0, vmax_mass=0;
    double fast_mass=0, fast_kinetic=0;
    double min_rho=std::numeric_limits<double>::infinity();
    double min_p=std::numeric_limits<double>::infinity();
    std::array<double,4> mass_velocity{};
    std::array<double,2> front{};
};

Diagnostics measure(const hd_core::MeshBlock& mesh,
                    const config::SimulationConfig& cfg,
                    const std::vector<double>& initial_pressure)
{
    Diagnostics d;
    double enclosed=cfg.spherical_central_mass;
    std::vector<std::pair<double,double>> speeds;
    const double c_code=2.99792458e10/physcon::units.unit_velocity;
    const double minimum_front_speed=1e7/physcon::units.unit_velocity;
    for(int i=0; i<mesh.nactive[0]; ++i)
    {
        const auto flat=mesh.flat_index({i},true);
        const auto& cell=mesh.cells[flat];
        const auto W=cell.U.to_Primtivies(cell.gamma);
        const auto [a,b]=mesh.cell_face_pos(0,i+mesh.nghost[0]);
        const double rho=W[hd_core::RHO], p=W[hd_core::PRE], v=W[hd_core::VEL(0)];
        const double dm=rho*mesh.volume_cell[flat], kinetic=0.5*dm*v*v;
        if(!(rho>0 && p>0 && std::isfinite(v)))
            throw std::runtime_error("Nonpositive or nonfinite diagnostic state");
        d.mass+=dm;
        d.kinetic+=kinetic;
        d.internal+=cell.U[hd_core::ENERGY]*mesh.volume_cell[flat]-kinetic;
        d.min_rho=std::min(d.min_rho,rho);
        d.min_p=std::min(d.min_p,p);
        if(v>d.vmax) { d.vmax=v; d.vmax_r=(a+b)/2; d.vmax_mass=dm; }
        if(v>0.1*c_code) { d.fast_mass+=dm; d.fast_kinetic+=kinetic; }
        speeds.emplace_back(v,dm);
        for(int j=0;j<2;++j)
            if(v>minimum_front_speed && p>initial_pressure[i]*(j==0?10:100))
                d.front[j]=b;
        // Exact 3-point Gauss integral of -G M(<r) dm/r for constant shell rho.
        // Includes point-core gravity and gas self-gravity once.
        constexpr std::array<double,3> nodes{-0.7745966692414834,0,0.7745966692414834};
        constexpr std::array<double,3> weights{5.0/9,8.0/9,5.0/9};
        for(int q=0;q<3;++q)
        {
            const double r=(a+b)/2+(b-a)/2*nodes[q];
            const double mr=enclosed+4*physcon::pi*rho/3*(r-a)*(r*r+r*a+a*a);
            d.gravity-=physcon::G*4*physcon::pi*rho*r*mr*(b-a)/2*weights[q];
        }
        enclosed+=dm;
    }
    std::sort(speeds.begin(),speeds.end(),std::greater<>());
    constexpr std::array<double,4> targets{0.001,0.01,0.1,1.0};
    double cumulative=0;
    std::size_t target=0;
    for(const auto [v,dm]:speeds)
    {
        cumulative+=dm;
        while(target<targets.size() && cumulative>=targets[target])
            d.mass_velocity[target++]=v;
    }
    // Undefined mass quantiles are NaN, e.g. after most of the gas exits.
    while(target<targets.size())
        d.mass_velocity[target++]=std::numeric_limits<double>::quiet_NaN();
    const auto flat=mesh.flat_index({mesh.nactive[0]-1},true);
    const auto& cell=mesh.cells[flat];
    const auto W=cell.U.to_Primtivies(cell.gamma);
    const double radius=cfg.upper[0], area=4*physcon::pi*radius*radius;
    const double v=std::max(0.0,W[hd_core::VEL(0)]);
    d.mass_flux=area*W[hd_core::RHO]*v;
    const double potential=-physcon::G*enclosed/radius;
    d.energy_flux=area*(cell.U[hd_core::ENERGY]+W[hd_core::PRE])*v
                  +potential*d.mass_flux;
    return d;
}

int main(int argc,char** argv)
{
    try
    {
        if(argc!=5) throw std::invalid_argument("usage: flare parameters atmosphere_multiplier heating_mass energy_erg");
        atmosphere_density*=std::stod(argv[2]);
        validation_heating_mass=std::stod(argv[3]);
        injection_energy_erg=std::stod(argv[4]);
        auto cfg=config::ReadParameters(argv[1]);
        cfg.output_step_every=1;
        hd_core::MeshBlock mesh;
        mesh.build_mesh(cfg);
        auto solver=hydro::make_riemann_solver(cfg.riemann_solver,cfg.riemann_tolerance,cfg.riemann_max_iterations);
        const auto spatial=hydro::parse_spatial_scheme("Godunov");
        const auto integrator=hd_time::parse_time_integrator(cfg.time_integrator);
        auto hooks=hd_user::make_user_hooks();
        hooks.initialize(mesh,cfg);
        boundary::FillBoundary(mesh,0,hooks);
        auto sources=hd_source::make_sources(cfg,mesh,hooks);
        std::vector<double> initial_pressure;
        for(int i=0;i<mesh.nactive[0];++i)
        {
            const auto [a,b]=mesh.cell_face_pos(0,i+mesh.nghost[0]);
            const auto cell=model.average(a,b,cfg.gamma);
            initial_pressure.push_back(cell.U.to_Primtivies(cell.gamma)[hd_core::PRE]);
        }
        std::ofstream log(std::filesystem::path(cfg.output_directory)/"diagnostics.csv");
        log<<std::setprecision(17)
           <<"step,time,mass,kinetic,internal,gravity,mass_flux,energy_flux,escaped_mass,escaped_energy,mass_residual,energy_residual,vmax,vmax_r,vmax_cell_mass,mass_gt_01c,kinetic_gt_01c,v_0001,v_001,v_01,v_1,front10,front100,min_rho,min_p\n";
        double time=0, escaped_mass=0, escaped_energy=0;
        std::size_t step=0;
        auto previous=measure(mesh,cfg,initial_pressure);
        const double mass0=previous.mass;
        const double energy0=previous.kinetic+previous.internal+previous.gravity;
        const auto write=[&](const Diagnostics& d)
        {
            log<<step<<','<<time<<','<<d.mass<<','<<d.kinetic<<','<<d.internal<<','<<d.gravity
               <<','<<d.mass_flux<<','<<d.energy_flux<<','<<escaped_mass<<','<<escaped_energy
               <<','<<d.mass+escaped_mass-mass0
               <<','<<d.kinetic+d.internal+d.gravity+escaped_energy-energy0
               <<','<<d.vmax<<','<<d.vmax_r<<','<<d.vmax_mass<<','<<d.fast_mass<<','<<d.fast_kinetic;
            for(auto v:d.mass_velocity) log<<','<<v;
            for(auto r:d.front) log<<','<<r;
            log<<','<<d.min_rho<<','<<d.min_p<<'\n';
        };
        write(previous);
        io::Output(mesh,cfg,hooks,step,time,0);
        double dt_next=std::numeric_limits<double>::infinity();
        // Matching physical times permit resolution comparisons of conserved fields.
        for(const double checkpoint:{10.,20.,30.,40.,50.,75.,100.})
        {
            const double stop=std::min(checkpoint,cfg.t_end);
            double dt=0;
            while(time<stop)
            {
                if(step>=cfg.maxsteps) throw std::runtime_error("Maximum steps exceeded");
                const double remaining=stop-time;
                dt=hd_time::advance_one_step(mesh,time,remaining,integrator,spatial,*solver,sources,cfg,hooks,cfg.cfl,dt_next);
                if(!(dt>0 && dt<=remaining && std::isfinite(dt) && time+dt>time))
                    throw std::runtime_error("Invalid timestep");
                time=dt==remaining ? stop : time+dt;
                ++step;
                const auto current=measure(mesh,cfg,initial_pressure);
                // Endpoint trapezoids approximate boundary transport; CFL/2 checks
                // distinguish this diagnostic error from conservation errors.
                escaped_mass+=dt/2*(previous.mass_flux+current.mass_flux);
                escaped_energy+=dt/2*(previous.energy_flux+current.energy_flux);
                write(current);
                previous=current;
            }
            io::Output(mesh,cfg,hooks,step,time,dt);
            log.flush();
            if(time>=cfg.t_end) break;
        }
        if(time!=cfg.t_end || !log) throw std::runtime_error("Incomplete diagnostics");
    }
    catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
