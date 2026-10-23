#pragma once
#include "user_config.h"
#include <array>
#include <cstddef>
#include <stdexcept>
namespace hd
{

inline constexpr int NDIM  = HD_NDIM;
inline constexpr int NCOMP = HD_NCOMP;
inline constexpr int NUSER = HD_NUSER;

//Total Number of Variables: for conservatives and primitives...
inline constexpr int NVAR = 2 + NDIM + NCOMP + NUSER;

static_assert(NDIM >= 1 && NDIM <= 3);
static_assert(NCOMP >= 0 && NUSER >= 0);

//Primitives component index
inline constexpr int RHO = 0//both primitives and conservatives
inline constexpr int PRE = NDIM + 1; 
constexpr VEL(int d)
{
    if(d<0 || d>=NDIM)
    {
        throw std::out_of_range("Invalid Dimension");
    }
    return 1+d;
}
//Conservatives component index
inline constexpr int ENERGY = PRE;
constexpr MOMDEN(int d)
{
    if(d<0 || d>=NDIM)
    {
        throw std::out_of_range("Invalid Dimension");
    }
    return 1+d;
}
//Hydrodynamics component number
inline constexpr int NHYDRO = NDIM + 2;


}