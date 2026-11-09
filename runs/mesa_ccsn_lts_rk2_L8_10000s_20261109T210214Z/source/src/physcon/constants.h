#pragma once
#include <cmath>
/**
 * @file constants.h
 * @brief Shared physical constants.
 */
namespace physcon
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double N_A = 6.02214076e23;    

    constexpr double c_ = 29979245800;
    
    inline double c = c_;
    
    constexpr double G_ = 6.67430e-8;
    
    inline double G = G_;

    constexpr double m_u_ = 1.66053906892e-24; 

    inline double m_u = m_u_;

    constexpr double k_B_ = 1.380649e-16;      //Boltzmann Constant, erg/K

    inline double k_B = k_B_;

    constexpr double sigma_SB_ = 5.670374419e-5;//Stefan-Bolzmann Constant, erg cm^-2 K^-4

    inline double sigma_SB = sigma_SB_;

    inline double a_Rad = 4.0 * sigma_SB_ / c_;

    
}