#pragma once
#include "user_config.h"
#include "hdmath.h"
#include <array>
#include <cstddef>
#include <stdexcept>
namespace hd_core
{

    inline constexpr int NDIM  = HD_NDIM;
    inline constexpr int NCOMP = HD_NCOMP;
    inline constexpr int NUSER = HD_NUSER;
    using VECDIM = std::array<double, NDIM>;

    //Total Number of Variables: for conservatives and primitives...
    inline constexpr int NVAR = 2 + NDIM + NCOMP + NUSER;

    static_assert(NDIM >= 1 && NDIM <= 3);
    static_assert(NCOMP >= 0 && NUSER >= 0);

    //Primitives component index
    inline constexpr int RHO = 0;//both primitives and conservatives
    inline constexpr int PRE = NDIM + 1; 
    constexpr int VEL(int d)
    {
        if(d<0 || d>=NDIM)
        {
            throw std::out_of_range("Invalid Dimension");
        }
        return 1+d;
    }
    //Conservatives component index
    inline constexpr int ENERGY = PRE;
    constexpr int MOMDEN(int d)
    {
        if(d<0 || d>=NDIM)
        {
            throw std::out_of_range("Invalid Dimension");
        }
        return 1+d;
    }
    //Hydrodynamics component number
    inline constexpr int NHYDRO = NDIM + 2;
    inline constexpr int COMP_BEGIN = NHYDRO;
    inline constexpr int USER_BEGIN = COMP_BEGIN + NCOMP;

    constexpr int COMP(int i)
    {
        if(i<0 || i>=NCOMP)
        {
            throw std::out_of_range("Invalid Dimension");        
        }
        return COMP_BEGIN + i;
    }
    constexpr int CONUSR(int i)
    {
        if(i<0 || i>=NUSER)
        {
            throw std::out_of_range("Invalid Dimension");        
        }
        return USER_BEGIN + i;
    }
    struct Conserved;
    struct Primitives
    {
        std::array<double,NVAR> data{};
        double& operator[](int i)
        {
            return data[i];
        }
        const double& operator[](int i) const
        {
            return data[i];
        }
        Conserved to_Conserved(double gamma) const;
        hdmath::Vector Velocity() const;
    };
    struct Conserved
    {
        std::array<double,NVAR> data{};
        double& operator[](int i)
        {
            return data[i];
        }
        const double& operator[](int i) const
        {
            return data[i];
        }
        Primitives to_Primtivies(double gamma) const;
        hdmath::Vector Velocity() const;
    };
    inline hdmath::Vector Primitives::Velocity() const
    {
        std::vector<double> part(data.begin() + 1, data.begin() + NDIM);
        return hdmath::Vector(part);
    }
    inline hdmath::Vector Conserved::Velocity() const
    {
        std::vector<double> part(data.begin() + 1, data.begin() + NDIM);
        auto momden = hdmath::Vector(part);
        return hdmath::Vector(momden/data[RHO]);
    }
    inline Conserved Primitives::to_Conserved(double gamma) const
    {
        Conserved W;
        W[RHO] = data[RHO];
        double E_k = 0;
        for(int i=0;i<NDIM;i++)
        {
            W[MOMDEN(i)] = W[RHO] * data[VEL(i)];
        }
        W[ENERGY] = data[PRE] / (gamma - 1);
    }
    inline Primitives Conserved::to_Primtivies(double gamma) const
    {
        Primitives W;
    }
}