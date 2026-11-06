#include "unit.h"
#include "constants.h"
#include <cmath>
namespace physcon
{
    Unit units{};

    void Unit::init_unitsystem_MLT(double M0, double L0, double T0, double Temp0)
    {
        unit_mass = M0;
        unit_length = L0;
        unit_time = T0;
        unit_temperature = Temp0;

        unit_density = unit_mass / std::pow(unit_length,3);

        unit_velocity = unit_length / unit_time;

        unit_acceleration = unit_velocity / unit_time;

        unit_force = unit_mass * unit_acceleration;

        unit_energy = unit_force * unit_length;

        unit_pressure = unit_energy / std::pow(unit_length,3);
        reset_constant();

    }
    void Unit::init_unitsystem_RHOLV(double Rho0, double L0, double V0, double Temp0)
    {
        init_unitsystem_MLT(Rho0*std::pow(L0,3),L0,L0/V0,Temp0);
    }
    void Unit::reset_constant()
    {
        c = c_ / (unit_length/unit_time);
        G = G_ /  (std::pow(unit_length,3) * (1/unit_mass) * std::pow(unit_time,-2));
        m_u = m_u_ / unit_mass;
        k_B = k_B_ * unit_temperature / unit_energy;
        sigma_SB = sigma_SB_ * unit_time * std::pow(unit_length,2)
                 * std::pow(unit_temperature,4) / unit_energy;
        a_Rad = 4.0 * sigma_SB / c;
    }
}
