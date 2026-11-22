#include "test_support.h"
#include "hydro/riemann/exact.h"
#include "hydro/flux.h"
#include "microphysics/composition/species.h"
#include "microphysics/eos/eos.h"
#include "physcon/constants.h"
#include "config/config.h"
#include "time/integrator.h"
#include "user/user_api.h"
#include <filesystem>
#include <fstream>
#include <limits>
using namespace hd_core;
using hd_microphysics::Equation_of_State;
using hydro::ExactPressureMethod;
using hydro::ExactSolver;
constexpr auto rad_model = Equation_of_State::Ionized_Radiative;

Primitives make_state(double rho, double pressure, double u=0, double ion=0)
{
    Primitives w{};
    w[RHO]=rho; w[PRE]=pressure; w[VEL(0)]=u;
    if constexpr (NCOMP==2) { w[COMP(0)]=1-ion; w[COMP(1)]=ion; }
    return w;
}
void flux_checks()
{
    hd_microphysics::eos_type=rad_model;
    ExactSolver solver(1e-10,256), bisect(1e-10,256,ExactPressureMethod::Bisection);
    double speed,velocity;
    for(double ratio:{1e-6,1.,1e6}) {
        double T=1e4, R=physcon::R_g*1.35, pr=physcon::a_Rad*std::pow(T,4)/3;
        double rho=pr/(ratio*R*T),p=rho*R*T+pr;
        auto w=make_state(rho,p,2e5,0.35); auto before=w;
        Cell ctx; ctx.gamma=5./3;
        auto U=ctx.to_Conserved(w,rad_model);
        for(int d=0;d<NDIM;d++) {
            auto F=solver.flux(w,w,d,ctx.gamma,ctx.gamma,speed,velocity);
            auto expected=hydro::physical_flux(w,U,d);
            for(int k=0;k<NVAR;k++) CHECK(near(F[k],expected[k],1e-10));
        }
        CHECK(w.data==before.data);
        auto other=make_state(3*rho,p,2e5,0.75);
        auto F=solver.flux(w,other,0,ctx.gamma,1.4,speed,velocity);
        auto expected=hydro::physical_flux(w,U,0);
        for(int k=0;k<NVAR;k++) CHECK(near(F[k],expected[k],1e-10));
        w[VEL(0)]=other[VEL(0)]=0;
        F=solver.flux(w,other,0,ctx.gamma,1.4,speed,velocity);
        CHECK(F[RHO]==0 && F[ENERGY]==0 && F[MOMDEN(0)]==p);
    }
    // gamma=4/3: thermal+radiation energy is exactly 3P, giving an independent
    // analytic gamma-law oracle at ANY gas/radiation pressure ratio.
    for(auto [l,r]:{std::pair{make_state(1e-10,1e4),make_state(1e-11,1e3)},
                   std::pair{make_state(1e-10,1e4,-1e7),make_state(1e-10,1e4,1e7)},
                   std::pair{make_state(1e-10,1e4),Primitives{}},
                   std::pair{Primitives{},make_state(1e-10,1e4)},
                   std::pair{make_state(1e-10,1e4,-1e8),make_state(1e-10,1e4,1e8)}}) {
        hd_microphysics::eos_type=Equation_of_State::IdealGas_Default;
        auto expected=solver.flux(l,r,0,4./3,4./3,speed,velocity);
        double expected_speed=speed,expected_velocity=velocity;
        hd_microphysics::eos_type=rad_model;
        double hint=std::numeric_limits<double>::quiet_NaN();
        auto F=solver.flux_with_pressure_hint(l,r,0,4./3,4./3,hint,speed,velocity);
        CHECK(near(speed,expected_speed,1e-9)); CHECK(near(velocity,expected_velocity,1e-9));
        for(int k=0;k<NVAR;k++) CHECK(near(F[k],expected[k],2e-8));
        auto B=bisect.flux(l,r,0,4./3,4./3,speed,velocity);
        for(int k=0;k<NVAR;k++) CHECK(near(F[k],B[k],2e-8));
    }
    // Different pressure and composition: both pressure solvers must agree.
    for(double u:{-1e8,0.,1e8}) {
        auto l=make_state(1e-10,1e4,u,0.35),r=make_state(2e-11,1e2,-u,0.8);
        auto F=solver.flux(l,r,0,5./3,1.4,speed,velocity); double a=speed,v=velocity;
        auto B=bisect.flux(l,r,0,5./3,1.4,speed,velocity);
        CHECK(near(a,speed,2e-8));CHECK(near(v,velocity,2e-8));
        for(int k=0;k<NVAR;k++) CHECK(near(F[k],B[k],1e-7));
    }
    // Frozen ionization energy changes only the advected energy, not the waves.
    auto l=make_state(1e-10,1e4,0,.35),r=make_state(2e-11,1e3,0,.35);
    hd_microphysics::eos_type=Equation_of_State::IdealGas_Mixed;
    auto M=solver.flux(l,r,0,5./3,5./3,speed,velocity);double a=speed,v=velocity;
    hd_microphysics::eos_type=Equation_of_State::Ionized;
    auto I=solver.flux(l,r,0,5./3,5./3,speed,velocity);
    CHECK(speed==a && velocity==v && M[RHO]==I[RHO]);
    const double ei=.35/physcon::m_u*hd_microphysics::species_table.species_active[1].energy_eV*physcon::eV;
    CHECK(near(I[ENERGY]-M[ENERGY],ei*M[RHO],1e-9));
    CHECK(near(I[INTERNAL_ENERGY]-M[INTERNAL_ENERGY],ei*M[RHO],1e-9));
    std::cout<<"uniform/contact, gas-law oracle, vacuum, Newton/bisection and ion-energy flux passed\n";
}
void integration_check()
{
    hd_microphysics::eos_type=rad_model;
    config::SimulationConfig cfg; cfg.reconstructor="PCM";cfg.lts_max_levels=2;
    MeshBlock initial;
    for(int d=0;d<NDIM;d++) {
        initial.nactive[d]=4;initial.nghost[d]=2;
        initial.bc_lower[d]=initial.bc_upper[d]=boundary::BoundaryCondition::Periodic;
        initial.set_linear_uniform_faces(d,0.,1.);
    }
    initial.allocate();initial.build_geometry();
    for(std::size_t flat=0;flat<initial.cells.size();flat++) {
        auto index=initial.coordinates(flat);
        bool left=index[0]<initial.nghost[0]+2;
        auto w=make_state(left?1e-10:2e-11,left?1e4:1e3,0.,left?.35:.75);
        auto& cell=initial.cells[flat];cell.gamma=5./3;
        cell.Set_Conserved_From_Primitives(w);
        hd_microphysics::refresh_thermodynamics(cell);
    }
    ExactSolver solver(1e-9,256);
    for(auto advance:{hd_time::advance_euler,hd_time::advance_lts_euler,hd_time::advance_lts_rk2}) {
        auto mesh=initial;
        auto total=[&](const MeshBlock& m) {
            Conserved U{};
            for(auto flat:m.active_to_storage) U+=m.cells[flat].U*m.volume_cell[flat];
            return U;
        };
        auto before=total(mesh);
        double dt=advance(mesh,0.,1e-10,hydro::SpatialScheme::Godunov,solver,{},cfg,{},.4,nullptr);
        CHECK(dt>0 && dt<=1e-10);
        auto after=total(mesh);
        for(int k:{RHO,ENERGY,COMP(0),COMP(1)})
            CHECK(std::abs((after[k]-before[k])/before[k])<2e-10);
        for(auto flat:mesh.active_to_storage) {
            mesh.cells[flat].U.validate();
            hd_microphysics::validate_cell_composition(mesh.cells[flat]);
            CHECK(mesh.cells[flat].temperature>0);
        }
    }
    std::cout<<"Euler/LTS Euler/LTS RK2 periodic radiation step conserves mass, total energy and species\n";
}
int main(int argc, char** argv)
{
    return run_test("Frozen-composition exact EoS regression", [&] {
        if constexpr (NCOMP != 2)
        {
            std::cout << "Radiation tests require HD_NCOMP=2 (H I/H II); skipped\n";
            return;
        }
        CHECK(argc == 2);
        const auto directory = std::filesystem::absolute(argv[1]);
        std::filesystem::create_directories(directory);
        config::SimulationConfig cfg;
        cfg.composition_activate_file = (directory / "active.txt").string();
        std::ofstream(cfg.composition_activate_file) << "h1 0 1\n";
        hd_microphysics::initialize_species_table(cfg);
        flux_checks();
        integration_check();
    });
}
