/**
 * @file hdmath.cpp
 * @brief Mathematical and string utilities.
 */
#include "hdmath.h"

#include <cmath>
#include <stdexcept>
#include <vector>
#include <array>
#include <utility>

namespace hdmath
{

Vector::Vector(std::size_t size, double value) : data_(size, value) {}

Vector::Vector(std::initializer_list<double> values) : data_(values) {}

Vector::Vector(std::vector<double> values) : data_(std::move(values)) {}

std::size_t Vector::size() const noexcept { return data_.size(); }
bool Vector::empty() const noexcept { return data_.empty(); }

double& Vector::operator[](std::size_t index) { return data_.at(index); }
const double& Vector::operator[](std::size_t index) const { return data_.at(index); }

void Vector::require_same_size(const Vector& other) const
{
    if (size() != other.size())
        throw std::invalid_argument("hdmath::Vector: vector lengths must match");
}

Vector& Vector::operator+=(const Vector& other)
{
    require_same_size(other);
    for (std::size_t i = 0; i < size(); ++i)
        data_[i] += other.data_[i];
    return *this;
}

Vector& Vector::operator-=(const Vector& other)
{
    require_same_size(other);
    for (std::size_t i = 0; i < size(); ++i)
        data_[i] -= other.data_[i];
    return *this;
}

Vector Vector::operator+(const Vector& other) const
{
    Vector result = *this;
    result += other;
    return result;
}

Vector Vector::operator-(const Vector& other) const
{
    Vector result = *this;
    result -= other;
    return result;
}

Vector Vector::operator-() const
{
    Vector result = *this;
    for (double& value : result.data_)
        value = -value;
    return result;
}

Vector& Vector::operator*=(double scalar)
{
    for (double& value : data_)
        value *= scalar;
    return *this;
}

Vector Vector::operator*(double scalar) const
{
    Vector result = *this;
    result *= scalar;
    return result;
}

Vector& Vector::operator/=(double scalar)
{
    if (scalar == 0.0)
        throw std::domain_error("hdmath::Vector: division by zero");

    // Divide directly: computing 1/scalar first may overflow for tiny scalars.
    for (double& value : data_)
        value /= scalar;
    return *this;
}

Vector Vector::operator/(double scalar) const
{
    Vector result = *this;
    result /= scalar;
    return result;
}

Vector operator*(double scalar, const Vector& vector)
{
    return vector * scalar;
}

double Vector::dot(const Vector& other) const
{
    require_same_size(other);
    double result = 0.0;
    for (std::size_t i = 0; i < size(); ++i)
        result += data_[i] * other.data_[i];
    return result;
}

double Vector::squared_norm() const { return dot(*this); }

double Vector::norm() const
{
    // Avoid unnecessary overflow/underflow from explicitly squaring components.
    double result = 0.0;
    for (double value : data_)
        result = std::hypot(result, value);
    return result;
}

Vector Vector::cross(const Vector& other) const
{
    if (size() != 3 || other.size() != 3)
        throw std::invalid_argument("hdmath::Vector: cross product requires two 3D vectors");

    return {
        data_[1] * other.data_[2] - data_[2] * other.data_[1],
        data_[2] * other.data_[0] - data_[0] * other.data_[2],
        data_[0] * other.data_[1] - data_[1] * other.data_[0]
    };
}
Vector Vector::slice(std::size_t start, std::size_t stop) const
{
    if (start > stop || stop > size())
        throw std::out_of_range("Invalid slice range");

    Vector result(stop - start);

    for (std::size_t i = start; i < stop; ++i)
        result[i - start] = data_[i];

    return result;
}
    double interp_linear(const std::vector<double>& value,
                        const std::vector<double>& grid,
                        const std::vector<double>& interp_pos,
                        int ndim)
    {
        if (ndim < 1 || ndim > 3)
            throw std::invalid_argument("interp_linear: ndim must be 1, 2, or 3");

        const auto dim = static_cast<std::size_t>(ndim);
        const std::size_t npoints = std::size_t{1} << dim;

        if (value.size() != npoints || grid.size() != 2 * dim ||
            interp_pos.size() != dim)
            throw std::invalid_argument("interp_linear: invalid input size");

        // t[d] 是查询点沿第 d 个方向从低端到高端的相对位置。
        std::array<double, 3> t{};
        for (std::size_t d = 0; d < dim; ++d)
        {
            const double lo = grid[2 * d];
            const double hi = grid[2 * d + 1];
            const double width = hi - lo;
            const double pos = interp_pos[d];

            if (!std::isfinite(width) || width <= 0.0)
                throw std::invalid_argument("interp_linear: invalid grid interval");
            if (!std::isfinite(pos) || pos < lo || pos > hi)
                throw std::out_of_range("interp_linear: point outside grid cell");

            t[d] = (pos - lo) / width;
        }

        double result = 0.0;
        for (std::size_t corner = 0; corner < npoints; ++corner)
        {
            double W = 1.0;
            for (std::size_t d = 0; d < dim; ++d)
            {
                const bool high = ((corner >> d) & 1u) != 0;
                W *= high ? t[d] : (1.0 - t[d]);
            }
            result += W * value[corner];
        }
        return result;
    }
} // namespace hdmath
