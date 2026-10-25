#pragma once
#include "user_config.h"
#include "hdmath.h"
#include <array>
#include <cstddef>
#include <cmath>
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
        // Constant-gamma ideal gas; vacuum states are not supported.
        void validate(double gamma) const;
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
        void validate(double gamma) const;
    };
    inline void Primitives::validate(double gamma) const
    {
        if (!std::isfinite(gamma) || gamma <= 1.0)
            throw std::invalid_argument("Primitives: gamma must be finite and > 1");

        // Extra components must be finite, but their physical bounds depend
        // on the user's definition (not every extra is a mass fraction).
        for (double value : data)
            if (!std::isfinite(value))
                throw std::invalid_argument("Primitives: state contains NaN or infinity");

        if (data[RHO] <= 0.0)
            throw std::invalid_argument("Primitives: density must be positive");
        if (data[PRE] <= 0.0)
            throw std::invalid_argument("Primitives: pressure must be positive");
    }

    inline void Conserved::validate(double gamma) const
    {
        if (!std::isfinite(gamma) || gamma <= 1.0)
            throw std::invalid_argument("Conserved: gamma must be finite and > 1");

        for (double value : data)
            if (!std::isfinite(value))
                throw std::invalid_argument("Conserved: state contains NaN or infinity");

        const double rho = data[RHO];
        if (rho <= 0.0)
            throw std::invalid_argument("Conserved: density must be positive");

        double kinetic = 0.0;
        for (int d = 0; d < NDIM; ++d)
        {
            const double momentum = data[MOMDEN(d)];
            const double velocity = momentum / rho;
            if (!std::isfinite(velocity))
                throw std::invalid_argument("Conserved: derived velocity is not finite");
            // Avoid forming momentum*momentum before dividing by rho.
            kinetic += (0.5 * momentum) * velocity;
        }
        if (!std::isfinite(kinetic))
            throw std::invalid_argument("Conserved: kinetic energy is not finite");

        const double internal = data[ENERGY] - kinetic;
        if (!std::isfinite(internal) || internal <= 0.0)
            throw std::invalid_argument("Conserved: internal energy must be finite and positive");

        const double pressure = (gamma - 1.0) * internal;
        if (!std::isfinite(pressure) || pressure <= 0.0)
            throw std::invalid_argument("Conserved: derived pressure must be finite and positive");
    }

    inline hdmath::Vector Primitives::Velocity() const
    {
        std::vector<double> part(data.begin() + 1, data.begin() + NDIM + 1);
        return hdmath::Vector(part);
    }
    inline hdmath::Vector Conserved::Velocity() const
    {
        std::vector<double> part(data.begin() + 1, data.begin() + NDIM +1 );
        auto momden = hdmath::Vector(part);
        return hdmath::Vector(momden/data[RHO]);
    }
    inline Conserved Primitives::to_Conserved(double gamma) const
    {
        Conserved W;
        W[RHO] = data[RHO];
        auto vel = Velocity();
        for(int i=0;i<NDIM;i++)
        {
            W[MOMDEN(i)] = W[RHO] * data[VEL(i)];
        }
        W[ENERGY] = data[PRE] / (gamma - 1) + 0.5 * W[RHO] * vel.squared_norm();
        return W;
    }
    inline Primitives Conserved::to_Primtivies(double gamma) const
    {
        Primitives W;
        W[RHO] = data[RHO];
        auto vel = Velocity();
        for(int i=0;i<NDIM;i++)
        {
            W[VEL(i)] = data[MOMDEN(i)] / W[RHO];
        }
        W[PRE] = (gamma - 1) * (data[ENERGY] - 0.5 * W[RHO] * vel.squared_norm());
        return W;
    }
}
