#include "test_support.h"
#include "core/cell.h"
#include "microphysics/eos/nuclear.h"
#include "physcon/unit.h"
#include "physcon/constants.h"
#include "hydro/riemann/factory.h"
#include "hydro/flux.h"
#include <limits>

namespace
{
using namespace hd_core;
using namespace hd_microphysics;
constexpr auto model = Equation_of_State::Nuclear;

void check_eos()
{
    const double density_unit = physcon::units.unit_density;
    double pl, ul, pr, ur;
    const double rn = nuclear_density_cgs/density_unit;
    nuclear_cold_state(rn*(1-1e-9),4.0/3,pl,ul);
    nuclear_cold_state(rn*(1+1e-9),4.0/3,pr,ur);
    CHECK(near(pl,pr,5e-9));
    CHECK(near(ul,ur,5e-9));
    for (double rho_cgs : {1e10,1e14,2e14,3e14,1e15})
    for (double temperature_kelvin : {0.0,1e9,1e11})
    {
        Cell c;
        c.gamma = 4.0/3;
        c.U[RHO] = rho_cgs/density_unit;
        c.U[MOMDEN(0)] = c.U[RHO]*1e8/physcon::units.unit_velocity;
        c.temperature = temperature_kelvin/physcon::units.unit_temperature;
        const double temperature = c.temperature;
        double p = 0;
        evaluate_eos(c,EoSInput::RhoT,model,p);
        const double energy = c.U.internal_energy();
        const double pressure_cgs = p*physcon::units.unit_pressure;
        const double pc_cgs = nuclear_K_cgs*std::pow(nuclear_density_cgs,4.0/3)
            *std::pow(rho_cgs/nuclear_density_cgs,rho_cgs>nuclear_density_cgs ? nuclear_gamma_stiff : 4.0/3);
        CHECK(near(pressure_cgs,pc_cgs+rho_cgs*(physcon::k_B_/physcon::m_u_)*temperature_kelvin,1e-10));
        const auto conserved = c.U.data;
        evaluate_eos(c,EoSInput::RhoE,model,p);
        CHECK(c.U.data == conserved);
        CHECK(std::abs(c.temperature-temperature) <= 1e-11*std::max(temperature,1e11/physcon::units.unit_temperature));
        CHECK(near(c.Get_Pressure_From_Updated_State(model),p));
        const double p_input = p;
        evaluate_eos(c,EoSInput::RhoP,model,p);
        CHECK(p == p_input);
        CHECK(near(c.U.internal_energy(),energy));

        // Independent isentropic finite difference: P_th scales as rho^gamma_th.
        // Skip the kink, where the two one-sided derivatives intentionally differ.
        if (rho_cgs == nuclear_density_cgs) continue;
        const double rho = c.U[RHO], h = 1e-5;
        double pc, uc;
        nuclear_cold_state(rho,c.gamma,pc,uc);
        const double pt = p-pc;
        nuclear_cold_state(rho*(1+h),c.gamma,pr,ur);
        nuclear_cold_state(rho*(1-h),c.gamma,pl,ul);
        const double derivative = (pr-pl+pt*(std::pow(1+h,nuclear_gamma_thermal)
            -std::pow(1-h,nuclear_gamma_thermal)))/(2*h*rho);
        CHECK(near(std::pow(c.get_sound_speed(model),2),derivative,1e-8));
        CHECK(near((ur/(rho*(1+h))-ul/(rho*(1-h)))/(2*h*rho),pc/(rho*rho),1e-8));

        Cell invalid = c;
        invalid.U[INTERNAL_ENERGY] = 0.99*uc;
        invalid.U[ENERGY] = 0.99*uc+invalid.U.kinetic_energy();
        expect_throw<std::invalid_argument>([&] { evaluate_eos(invalid,EoSInput::RhoE,model,p); });
        double bad_p = 0.99*pc;
        expect_throw<std::invalid_argument>([&] { evaluate_eos(c,EoSInput::RhoP,model,bad_p); });
    }
}

void check_flux(const std::string& name)
{
    eos_type = model;
    const auto solver = hydro::make_riemann_solver(name,1e-10,100);
    Cell c;
    c.gamma = 4.0/3;
    c.U[RHO] = 3e14;
    c.temperature = 1e11;
    c.U[MOMDEN(0)] = c.U[RHO]*1e8;
    double p;
    evaluate_eos(c,EoSInput::RhoT,model,p);
    auto w = c.Get_Primitives(model);
    double speed, velocity;
    const auto flux = solver->flux(w,w,0,c.gamma,c.gamma,speed,velocity);
    const auto expected = hydro::physical_flux(w,c.U,0);
    for (std::size_t k=0;k<flux.data.size();++k) CHECK(near(flux[k],expected[k]));
    CHECK(near(velocity,w[VEL(0)]));
    CHECK(near(speed,std::abs(w[VEL(0)])+c.get_sound_speed(model)));
    auto reflected = w;
    reflected[VEL(0)] = -w[VEL(0)];
    const auto wall = solver->flux(reflected,w,0,c.gamma,c.gamma,speed,velocity);
    CHECK(std::abs(wall[RHO]) < 1e-12*std::abs(expected[RHO]));
    CHECK(velocity == 0);
    CHECK(near(wall[ENERGY],0));
    const auto collision = solver->flux(w,reflected,0,c.gamma,c.gamma,speed,velocity);
    CHECK(collision[MOMDEN(0)] > p);
    CHECK(velocity == 0.0 && collision[RHO] == 0.0);
    eos_type = Equation_of_State::IdealGas_Default;
}

void check_exact_vacuum(const hydro::RiemannSolver& solver)
{
    const double g = 4.0/3.0, stiff = nuclear_gamma_stiff;
    const double rn = nuclear_density_cgs / physcon::units.unit_density;
    Primitives gas{}, vacuum{}, face{};
    gas[RHO] = 2.0*rn;
    face[RHO] = rn;
    double uc;
    nuclear_cold_state(gas[RHO], g, gas[PRE], uc);
    nuclear_cold_state(rn, g, face[PRE], uc);

    // Choose the incident velocity so the face lies between the two sonic
    // characteristics at the kink. Integrate the two cold polytropes analytically.
    const double c0 = std::sqrt(stiff*gas[PRE]/gas[RHO]);
    const double c_high = std::sqrt(stiff*face[PRE]/rn);
    const double c_low = std::sqrt(g*face[PRE]/rn);
    const double gain_to_kink = 2.0*(c0-c_high)/(stiff-1.0);
    const double gain_to_vacuum = gain_to_kink + 2.0*c_low/(g-1.0);
    const int d = NDIM-1;
    double speed, velocity, hint = gas[PRE];
    for (double sign : {1.0, -1.0})
    {
        face[VEL(d)] = sign*0.5*(c_low+c_high);
        gas[VEL(d)] = face[VEL(d)] - sign*gain_to_kink;
        const auto expected = hydro::physical_flux(face, d, g);
        const auto one_sided = solver.flux_with_pressure_hint(
            sign > 0 ? gas : vacuum, sign > 0 ? vacuum : gas, d, g, g, hint, speed, velocity);
        for (int k=0; k<NVAR; ++k) CHECK(near(one_sided[k], expected[k], 2e-9));
        CHECK(near(velocity, face[VEL(d)], 2e-9));
        CHECK(near(speed, std::max(std::abs(gas[VEL(d)]-sign*c0),
                                   std::abs(gas[VEL(d)]+sign*gain_to_vacuum)), 2e-9));
        CHECK(hint == 0.0);

        // A second gas receding beyond the vacuum edge gives the same face state.
        auto other = gas;
        other[VEL(d)] = gas[VEL(d)] + sign*3.0*gain_to_vacuum;
        const auto generated = solver.flux(sign > 0 ? gas : other, sign > 0 ? other : gas,
                                           d, g, g, speed, velocity);
        for (int k=0; k<NVAR; ++k) CHECK(near(generated[k], expected[k], 2e-9));
        CHECK(near(velocity, face[VEL(d)], 2e-9));
    }
    auto left = gas, right = gas;
    left[VEL(d)] = -2.0*gain_to_vacuum;
    right[VEL(d)] = 2.0*gain_to_vacuum;
    const auto central = solver.flux(left, right, d, g, g, speed, velocity);
    CHECK(central.data == Conserved{}.data && velocity == 0.0);
    CHECK(near(speed, right[VEL(d)]+c0));
    const auto empty = solver.flux(vacuum, vacuum, d, g, g, speed, velocity);
    CHECK(empty.data == Conserved{}.data && speed == 0.0 && velocity == 0.0);
}
}

int main()
{
    return run_test("nuclear",[] {
        physcon::units.init_unitsystem_MLT(1,1,1,1);
        check_eos();
        for (const char* name : {"hll", "exact", "exact_bisection"}) check_flux(name);
        // The same physical model must survive non-CGS density/energy/T units.
        physcon::units.init_unitsystem_RHOLV(1e12,1e6,1e9,1e10);
        check_eos();
        eos_type = model;
        for (const char* name : {"exact", "exact_bisection"})
            check_exact_vacuum(*hydro::make_riemann_solver(name,1e-10,100));
        eos_type = Equation_of_State::IdealGas_Default;
        physcon::units.init_unitsystem_MLT(1,1,1,1);
    });
}
