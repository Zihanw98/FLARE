#pragma once
namespace physcon
{
    struct Unit
    {
        double unit_mass = 1.0;         //g
        double unit_length = 1.0;       //cm
        double unit_time = 1.0;         //s

        double unit_density = 1.0;      //g cm^-3
        double unit_velocity = 1.0;     // cm s^-1
        double unit_acceleration = 1.0; //cm s^-2
        double unit_force = 1.0;        //dyn
        double unit_energy = 1.0;       //  g cm^2 s^-2
        double unit_pressure = 1.0;     //erg cm^-3


        //independent
        double unit_temperature = 1.0;  //K
        
        void init_unitsystem_MLT(double M0, double L0, double T0, double Temp0);

        void init_unitsystem_RHOLV(double Rho0, double L0, double V0, double Temp0);

        // C_phys is constant, C_phys = C_code * [units], so C_code = C_phys * [units_0] ^(-1), units_0 -> 100 units_0 => C_code -> C_code * 0.01 
        //NEED A DETAILED DOXYGEN
        void reset_constant();
    };
    extern Unit units;
}