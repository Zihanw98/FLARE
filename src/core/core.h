/**
 * @file core.h
 * @brief Physical states and mesh storage.
 */
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
    using Index = std::array<int, NDIM>;
    inline Index Index_add_offset(const Index& index,const Index& offset)
    {
        Index ret = index;
        for(int i=0;i<NDIM;i++)
        {
            ret[i] += offset[i];
        }
        return ret;
    }
    //Total Number of Variables: for conservatives and primitives...
    inline constexpr int NVAR = 2 + NDIM + NCOMP + NUSER;

    static_assert(NDIM >= 1 && NDIM <= 3);
    static_assert(NCOMP >= 0 && NUSER >= 0);

    //Primitives component index
    inline constexpr int RHO = 0;//both primitives and conservatives
    inline constexpr int PRE = NDIM + 1; 
    /**
     * @brief Return the storage index of a velocity component.
     * @param[in] d Zero-based component index in [0, NDIM).
     * @return The array offset 1 + d.
     * @throws std::out_of_range If the component index is outside its valid range.
     */
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
    /**
     * @brief Return the storage index of a momentum-density component.
     * @param[in] d Zero-based component index in [0, NDIM).
     * @return The array offset 1 + d.
     * @throws std::out_of_range If the component index is outside its valid range.
     */
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

    /**
     * @brief Return the storage index of a composition component.
     * @param[in] i Zero-based component index in [0, NCOMP).
     * @return The array offset COMP_BEGIN + i.
     * @throws std::out_of_range If the component index is outside its valid range.
     */
    constexpr int COMP(int i)
    {
        if(i<0 || i>=NCOMP)
        {
            throw std::out_of_range("Invalid Dimension");        
        }
        return COMP_BEGIN + i;
    }
    /**
     * @brief Return the storage index of a user-defined conserved component.
     * @param[in] i Zero-based component index in [0, NUSER).
     * @return The array offset USER_BEGIN + i.
     * @throws std::out_of_range If the component index is outside its valid range.
     */
    constexpr int CONUSR(int i)
    {
        if(i<0 || i>=NUSER)
        {
            throw std::out_of_range("Invalid Dimension");        
        }
        return USER_BEGIN + i;
    }
    struct Primitives;
    struct Conserved;
    

    /** @brief Primitive ideal-gas variables with compile-time component storage. */
    struct Primitives
    {
        std::array<double,NVAR> data{};
        /**
         * @brief Access a state component without bounds checking.
         * @param[in] i Component index in [0, NVAR).
         * @return A mutable reference to the stored component.
         * @pre The index is valid; an out-of-range access is undefined behavior.
         */
        double& operator[](int i)
        {
            return data[i];
        }
        /**
         * @brief Access a state component without bounds checking.
         * @param[in] i Component index in [0, NVAR).
         * @return A read-only reference to the stored component.
         * @pre The index is valid; an out-of-range access is undefined behavior.
         */
        const double& operator[](int i) const
        {
            return data[i];
        }
        /**
         * @brief Convert primitive hydrodynamic components to conserved variables.
         * @param[in] gamma Finite ideal-gas adiabatic index greater than one.
         * @return Density, momentum densities, and total energy density; extra components remain zero.
         * @pre The input is a valid primitive state for gamma.
         * @note Zero density returns an all-zero state before validation. Only hydrodynamic components are converted.
         * @throws std::invalid_argument If the resulting non-vacuum state fails validation.
         */
        Conserved to_Conserved(double gamma) const;
        /**
         * @brief Copy the stored primitive velocity components into a vector.
         * @return An NDIM-component velocity vector.
         * @note Does not validate density, pressure, or velocity values.
         */
        hdmath::Vector Velocity() const;
        // Ideal gas with a finite gamma > 1; zero-density vacuum is supported.
        /**
         * @brief Validate primitive variables for the ideal-gas equation of state.
         * @param[in] gamma Finite adiabatic index greater than one.
         * @throws std::invalid_argument If any component is non-finite, density is negative, gas pressure is not positive, or vacuum has nonzero pressure or velocity.
         * @note Zero-density vacuum is allowed when pressure and all velocities are zero; finite extra components are not otherwise constrained.
         */
        void validate(double gamma) const;
    };


    /** @brief Conserved variables, fluxes, or RHS values with componentwise arithmetic. */
    struct Conserved
    {
        std::array<double,NVAR> data{};
        /**
         * @brief Access a state component without bounds checking.
         * @param[in] i Component index in [0, NVAR).
         * @return A mutable reference to the stored component.
         * @pre The index is valid; an out-of-range access is undefined behavior.
         */
        double& operator[](int i)
        {
            return data[i];
        }
        /**
         * @brief Access a state component without bounds checking.
         * @param[in] i Component index in [0, NVAR).
         * @return A read-only reference to the stored component.
         * @pre The index is valid; an out-of-range access is undefined behavior.
         */
        const double& operator[](int i) const
        {
            return data[i];
        }

        // Component-wise arithmetic also applies to fluxes and RHS values,
        // which need not satisfy physical-state positivity constraints.
        /**
         * @brief Add another object componentwise in place.
         * @param[in] other Right-hand Conserved operand.
         * @return A reference to this object.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved& operator+=(const Conserved& other)
        {
            for (int k = 0; k < NVAR; ++k)
                data[k] += other.data[k];
            return *this;
        }

        /**
         * @brief Subtract another object componentwise in place.
         * @param[in] other Right-hand Conserved operand.
         * @return A reference to this object.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved& operator-=(const Conserved& other)
        {
            for (int k = 0; k < NVAR; ++k)
                data[k] -= other.data[k];
            return *this;
        }

        /**
         * @brief Add two objects componentwise.
         * @param[in] other Right-hand Conserved operand.
         * @return A new Conserved containing the result.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved operator+(const Conserved& other) const
        {
            Conserved result = *this;
            result += other;
            return result;
        }

        /**
         * @brief Subtract another object componentwise.
         * @param[in] other Right-hand Conserved operand.
         * @return A new Conserved containing the result.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved operator-(const Conserved& other) const
        {
            Conserved result = *this;
            result -= other;
            return result;
        }

        /**
         * @brief Negate every component.
         * @return A new Conserved containing the result.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved operator-() const
        {
            Conserved result = *this;
            for (double& value : result.data)
                value = -value;
            return result;
        }

        /**
         * @brief Scale every component in place.
         * @param[in] scalar Scalar multiplier.
         * @return A reference to this object.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved& operator*=(double scalar)
        {
            for (double& value : data)
                value *= scalar;
            return *this;
        }

        /**
         * @brief Return a scaled copy.
         * @param[in] scalar Scalar multiplier.
         * @return A new Conserved containing the result.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved operator*(double scalar) const
        {
            Conserved result = *this;
            result *= scalar;
            return result;
        }

        /**
         * @brief Divide every component in place.
         * @param[in] scalar Scalar divisor; must not be zero.
         * @return A reference to this object.
         * @throws std::domain_error If scalar is zero.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved& operator/=(double scalar)
        {
            if (scalar == 0.0)
                throw std::domain_error("Conserved: division by zero");

            // Divide directly to avoid overflowing the reciprocal of a tiny scalar.
            for (double& value : data)
                value /= scalar;
            return *this;
        }

        /**
         * @brief Return a copy divided by a scalar.
         * @param[in] scalar Scalar divisor; must not be zero.
         * @return A new Conserved containing the result.
         * @throws std::domain_error If scalar is zero.
         * @note Operates on all NVAR components without physical-state validation; also usable for fluxes and RHS values.
         */
        Conserved operator/(double scalar) const
        {
            Conserved result = *this;
            result /= scalar;
            return result;
        }

        /**
         * @brief Convert conserved hydrodynamic components to primitive variables.
         * @param[in] gamma Finite ideal-gas adiabatic index greater than one.
         * @return Density, velocities, and pressure; extra components remain zero.
         * @pre The input is a valid conserved state for gamma.
         * @note Zero density returns an all-zero primitive state before validating momentum or energy. The existing API spelling to_Primtivies is retained.
         * @throws std::invalid_argument If velocity extraction or validation of the converted gas state fails.
         */
        Primitives to_Primtivies(double gamma) const;
        /**
         * @brief Compute velocity from momentum density and mass density.
         * @return An NDIM-component velocity vector, or a zero vector for a valid vacuum.
         * @throws std::invalid_argument If density is non-finite or negative, or a zero-density state has nonzero momentum or energy.
         * @note This accessor is not a substitute for full Conserved::validate() on non-vacuum input.
         */
        hdmath::Vector Velocity() const;
        /**
         * @brief Validate conserved variables and the derived ideal-gas pressure.
         * @param[in] gamma Finite adiabatic index greater than one.
         * @throws std::invalid_argument If components are non-finite, density is negative, derived gas quantities are invalid, or vacuum has nonzero momentum or energy.
         * @note Positive-density gas requires finite positive internal energy and pressure. Zero-density vacuum is permitted.
         */
        void validate(double gamma) const;
    };

    /**
     * @brief Multiply a conserved state by a scalar on the left.
     * @param[in] scalar Scalar multiplier.
     * @param[in] state Operand to scale without modifying it.
     * @return A new object containing the componentwise product.
     */
    inline Conserved operator*(double scalar, const Conserved& state)
    {
        return state * scalar;
    }

    inline void Primitives::validate(double gamma) const
    {
        if (!std::isfinite(gamma) || gamma <= 1.0)
            throw std::invalid_argument("Invalid gamma");

        for (double value : data)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument("Non-finite primitive state");
        }

        if (data[RHO] == 0.0)
        {
            if (data[PRE] != 0.0)
                throw std::invalid_argument(
                    "Vacuum must have zero pressure"
                );

            for (int d = 0; d < NDIM; ++d)
            {
                if (data[VEL(d)] != 0.0)
                    throw std::invalid_argument(
                        "Vacuum velocity must be stored as zero"
                    );
            }

            return;
        }

        if (data[RHO] < 0.0)
            throw std::invalid_argument("Negative density");

        if (data[PRE] <= 0.0)
            throw std::invalid_argument("Gas pressure must be positive");
    }

    inline void Conserved::validate(double gamma) const
    {
        if (!std::isfinite(gamma) || gamma <= 1.0)
            throw std::invalid_argument("Conserved: gamma must be finite and > 1");

        for (double value : data)
            if (!std::isfinite(value))
                throw std::invalid_argument("Conserved: state contains NaN or infinity");

        const double rho = data[RHO];
        if (rho == 0.0)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                if (data[MOMDEN(d)] != 0.0)
                    throw std::invalid_argument(
                        "Vacuum must have zero momentum"
                    );
            }

            if (data[ENERGY] != 0.0)
                throw std::invalid_argument(
                    "Vacuum must have zero energy"
                );

            return;
        }

        if (rho < 0.0)
            throw std::invalid_argument("Negative density");

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
        const double rho = data[RHO];

        if (!std::isfinite(rho) || rho < 0.0)
            throw std::invalid_argument("Invalid density");

        if (rho == 0.0)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                if (data[MOMDEN(d)] != 0.0)
                    throw std::invalid_argument(
                        "Vacuum must have zero momentum"
                    );
            }

            if (data[ENERGY] != 0.0)
                throw std::invalid_argument(
                    "Vacuum must have zero energy"
                );

            return hdmath::Vector(NDIM, 0.0);
        }

        std::vector<double> part(
            data.begin() + 1,
            data.begin() + NDIM + 1
        );
        
        return hdmath::Vector(part) / rho;
    }

    inline Conserved Primitives::to_Conserved(double gamma) const
    {
        if (data[RHO] == 0.0)
        {
            return Conserved{};
        }
        Conserved W;
        W[RHO] = data[RHO];
        double speed_squared = 0.0;
        for(int i=0;i<NDIM;i++)
        {
            speed_squared += data[VEL(i)] * data[VEL(i)];
            W[MOMDEN(i)] = W[RHO] * data[VEL(i)];
        }
        W[ENERGY] = data[PRE] / (gamma - 1) + 0.5 * W[RHO] * speed_squared;
        W.validate(gamma);
        return W;
    }
    inline Primitives Conserved::to_Primtivies(double gamma) const
    {
        if (data[RHO] == 0.0)
        {
            return Primitives{};
        }
        Primitives W;
        W[RHO] = data[RHO];
        if (!std::isfinite(data[RHO]) || data[RHO] < 0.0)
            throw std::invalid_argument("Invalid density");
        double speed_squared = 0.0;
        for(int i=0;i<NDIM;i++)
        {
            W[VEL(i)] = data[MOMDEN(i)] / W[RHO];
            speed_squared += W[VEL(i)] * W[VEL(i)];
        }
        W[PRE] = (gamma - 1) * (data[ENERGY] - 0.5 * W[RHO] * speed_squared);
        W.validate(gamma);
        return W;
    }
}
