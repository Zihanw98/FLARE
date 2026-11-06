#include "test_support.h"
#include "physcon/unit.h"
#include "physcon/constants.h"

int main()
{
    return run_test("physical constants across unit changes", [] {
        physcon::units.init_unitsystem_MLT(2.0,3.0,5.0,7.0);
        const auto& u=physcon::units;
        CHECK(near(physcon::G/physcon::G_, 2.0*25.0/27.0));
        CHECK(near(physcon::c/physcon::c_, 5.0/3.0));
        CHECK(near(physcon::m_u/physcon::m_u_, 0.5));
        // E = k_B T: conversion back to cgs must recover the original constant.
        CHECK(near(physcon::k_B*u.unit_energy/u.unit_temperature/physcon::k_B_,1.0));
        // F = sigma T^4, F has dimensions energy/(time area).
        CHECK(near(physcon::sigma_SB*u.unit_energy/
            (u.unit_time*u.unit_length*u.unit_length*std::pow(u.unit_temperature,4))
            /physcon::sigma_SB_,1.0));
        CHECK(near(physcon::a_Rad*u.unit_energy/
            (std::pow(u.unit_length,3)*std::pow(u.unit_temperature,4))
            /(4*physcon::sigma_SB_/physcon::c_),1.0));
        physcon::units.init_unitsystem_RHOLV(2.0,3.0,4.0,7.0);
        CHECK(near(u.unit_density,2.0));
        CHECK(near(u.unit_velocity,4.0));
        physcon::units.init_unitsystem_MLT(1.0,1.0,1.0,1.0);
        CHECK(physcon::G==physcon::G_);
        CHECK(physcon::k_B==physcon::k_B_);
        CHECK(physcon::sigma_SB==physcon::sigma_SB_);
    });
}
