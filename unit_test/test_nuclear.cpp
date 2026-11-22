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

void check_flux()
{
    eos_type = model;
    expect_throw<std::invalid_argument>([] { hydro::make_riemann_solver("exact",1e-10,100); });
    const auto solver = hydro::make_riemann_solver("hll",1e-10,100);
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
    eos_type = Equation_of_State::IdealGas_Default;
}
}

int main()
{
    return run_test("nuclear",[] {
        physcon::units.init_unitsystem_MLT(1,1,1,1);
        check_eos();
        check_flux();
        // The same physical model must survive non-CGS density/energy/T units.
        physcon::units.init_unitsystem_RHOLV(1e12,1e6,1e9,1e10);
        check_eos();
        physcon::units.init_unitsystem_MLT(1,1,1,1);
    });
}
