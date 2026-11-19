/**
 * @file hdmath.h
 * @brief Mathematical and string utilities.
 */
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace hdmath
{

    /**
     * @brief Solve f(x) = 0 by Newton-Raphson using a caller-supplied derivative.
     * @param[in] function Callable double(double); may capture equation parameters.
     * @param[in] derivative Callable returning df/dx for the same equation.
     * @param[in] initial_guess Finite starting value in the function's domain.
     * @param[in] residual_tolerance Positive absolute tolerance in units of f.
     * Normalize f and df/dx by the same physical scale for a relative residual.
     * @param[in] max_iterations Positive maximum number of Newton updates.
     * @return x for which abs(f(x)) <= residual_tolerance.
     * @throws std::invalid_argument For invalid iteration settings.
     * @throws std::runtime_error For nonfinite arithmetic, a zero derivative,
     * stagnation above the residual tolerance, or exhaustion of the update limit.
     * @note This is unbracketed Newton iteration; the caller chooses an appropriate
     * starting value. Callback exceptions propagate. No cell or EoS state is owned
     * or updated by this function; any callback side effects belong to the caller.
     */
    template<class Function, class Derivative>
    double newton_raphson(Function&& function, Derivative&& derivative,
                         double initial_guess, double residual_tolerance = 1e-10,
                         int max_iterations = 100)
    {
        if (!std::isfinite(initial_guess) || !std::isfinite(residual_tolerance) ||
            residual_tolerance <= 0.0 || max_iterations <= 0)
            throw std::invalid_argument("newton_raphson: invalid initial guess, tolerance or iteration limit");

        double x = initial_guess;
        for (int iteration = 0; ; ++iteration)
        {
            const double residual = function(x);
            if (!std::isfinite(residual))
                throw std::runtime_error("newton_raphson: nonfinite residual");
            if (std::abs(residual) <= residual_tolerance)
                return x;
            if (iteration == max_iterations)
                break;

            const double slope = derivative(x);
            if (!std::isfinite(slope) || slope == 0.0)
                throw std::runtime_error("newton_raphson: derivative must be finite and nonzero");
            const double correction = residual / slope;
            const double next = x - correction;
            if (!std::isfinite(correction) || !std::isfinite(next))
                throw std::runtime_error("newton_raphson: nonfinite Newton update");
            // A small step alone does not establish that the equation is solved.
            if (next == x)
                throw std::runtime_error("newton_raphson: stagnated before residual convergence");
            x = next;
        }
        throw std::runtime_error("newton_raphson: maximum iterations exceeded");
    }

    // A dynamically sized mathematical vector of doubles.
    /** @brief Dynamically sized mathematical vector with checked component access. */
    class Vector
    {
    public:
        /**
         * @brief Construct an empty mathematical vector.
         */
        Vector() = default;
        /**
         * @brief Construct a vector with uniformly initialized components.
         * @param[in] size Number of components.
         * @param[in] value Initial value of every component.
         */
        explicit Vector(std::size_t size, double value = 0.0);
        /**
         * @brief Construct a vector from an initializer list.
         * @param[in] values Components in their desired order.
         */
        Vector(std::initializer_list<double> values);
        /**
         * @brief Construct a vector by taking ownership of a by-value storage argument.
         * @param[in] values Component storage moved into the vector; an lvalue argument is copied first.
         */
        explicit Vector(std::vector<double> values);

        // Templates must be defined in the header for arbitrary array lengths.
        /**
         * @brief Copy a fixed-size array into a dynamic mathematical vector.
         * @param[in] values Array containing the components in order.
         * @tparam N Number of array entries.
         */
        template<std::size_t N>
        explicit Vector(const std::array<double, N>& values)
            : data_(values.begin(), values.end()) {}

        /**
         * @brief Return the number of vector components.
         * @return The current component count.
         */
        std::size_t size() const noexcept;
        /**
         * @brief Check whether the vector has no components.
         * @return True if size() is zero; otherwise false.
         */
        bool empty() const noexcept;

        // Unlike std::vector::operator[], both accessors check bounds.
        /**
         * @brief Access a vector component with bounds checking.
         * @param[in] index Zero-based component offset.
         * @return A mutable reference to the selected value.
         * @throws std::out_of_range If index is not less than size().
         */
        double& operator[](std::size_t index);
        /**
         * @brief Access a vector component with bounds checking.
         * @param[in] index Zero-based component offset.
         * @return A read-only reference to the selected value.
         * @throws std::out_of_range If index is not less than size().
         */
        const double& operator[](std::size_t index) const;

        /**
         * @brief Add another object componentwise in place.
         * @param[in] other Right-hand Vector operand.
         * @return A reference to this object.
         * @throws std::invalid_argument If vector lengths differ.
         */
        Vector& operator+=(const Vector& other);
        /**
         * @brief Subtract another object componentwise in place.
         * @param[in] other Right-hand Vector operand.
         * @return A reference to this object.
         * @throws std::invalid_argument If vector lengths differ.
         */
        Vector& operator-=(const Vector& other);
        /**
         * @brief Add two objects componentwise.
         * @param[in] other Right-hand Vector operand.
         * @return A new Vector containing the result.
         * @throws std::invalid_argument If vector lengths differ.
         */
        Vector operator+(const Vector& other) const;
        /**
         * @brief Subtract another object componentwise.
         * @param[in] other Right-hand Vector operand.
         * @return A new Vector containing the result.
         * @throws std::invalid_argument If vector lengths differ.
         */
        Vector operator-(const Vector& other) const;
        /**
         * @brief Negate every component.
         * @return A new Vector containing the result.
         */
        Vector operator-() const;
        /**
         * @brief Scale every component in place.
         * @param[in] scalar Scalar multiplier.
         * @return A reference to this object.
         */
        Vector& operator*=(double scalar);
        /**
         * @brief Return a scaled copy.
         * @param[in] scalar Scalar multiplier.
         * @return A new Vector containing the result.
         */
        Vector operator*(double scalar) const;
        // Division by zero throws std::domain_error; no scalar / vector overload.
        /**
         * @brief Divide every component in place.
         * @param[in] scalar Scalar divisor; must not be zero.
         * @return A reference to this object.
         * @throws std::domain_error If scalar is zero.
         */
        Vector& operator/=(double scalar);
        /**
         * @brief Return a copy divided by a scalar.
         * @param[in] scalar Scalar divisor; must not be zero.
         * @return A new Vector containing the result.
         * @throws std::domain_error If scalar is zero.
         */
        Vector operator/(double scalar) const;

        /**
         * @brief Compute the Euclidean inner product with another vector.
         * @param[in] other Vector of the same length.
         * @return The sum of pairwise component products.
         * @throws std::invalid_argument If vector lengths differ.
         */
        double dot(const Vector& other) const;
        /**
         * @brief Compute the squared Euclidean norm.
         * @return The sum of squared components, or zero for an empty vector.
         * @note Direct multiplication may overflow or underflow; norm() uses iterative hypot instead.
         */
        double squared_norm() const; // Sum of squares, not elementwise squaring.
        /**
         * @brief Compute the Euclidean norm using iterative hypot accumulation.
         * @return The nonnegative vector length, or zero for an empty vector.
         */
        double norm() const;
        /**
         * @brief Compute the three-dimensional cross product.
         * @param[in] other Right-hand three-component vector.
         * @return A new vector equal to this vector crossed with other.
         * @throws std::invalid_argument If either vector does not have exactly three components.
         */
        Vector cross(const Vector& other) const; // Both operands must have size 3.
        /**
         * @brief Copy a half-open range of vector components.
         * @param[in] start Inclusive first offset.
         * @param[in] stop Exclusive final offset.
         * @return A new vector containing entries in [start, stop); empty when start == stop.
         * @throws std::out_of_range Unless start <= stop <= size().
         */
        Vector slice(std::size_t start, std::size_t stop) const;

    private:
        /**
         * @brief Require matching vector lengths before a binary operation.
         * @param[in] other Vector whose length is checked against this vector.
         * @throws std::invalid_argument If vector lengths differ.
         */
        void require_same_size(const Vector& other) const;
        std::vector<double> data_;
    };

    // Scalar on the left: 2.0 * v.
    /**
     * @brief Multiply a vector by a scalar on the left.
     * @param[in] scalar Scalar multiplier.
     * @param[in] vector Operand to scale without modifying it.
     * @return A new object containing the componentwise product.
     */
    Vector operator*(double scalar, const Vector& vector);

    double interp_linear(const std::vector<double>& value, const std::vector<double>& grid,
                         const std::vector<double>& interp_pos, int ndim);
} // namespace hdmath
