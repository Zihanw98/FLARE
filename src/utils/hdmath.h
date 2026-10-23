#pragma once

#include <array>
#include <cstddef>
#include <initializer_list>
#include <vector>

namespace hdmath
{

// A dynamically sized mathematical vector of doubles.
class Vector
{
public:
    Vector() = default;
    explicit Vector(std::size_t size, double value = 0.0);
    Vector(std::initializer_list<double> values);
    explicit Vector(std::vector<double> values);

    // Templates must be defined in the header for arbitrary array lengths.
    template<std::size_t N>
    explicit Vector(const std::array<double, N>& values)
        : data_(values.begin(), values.end()) {}

    std::size_t size() const noexcept;
    bool empty() const noexcept;

    // Unlike std::vector::operator[], both accessors check bounds.
    double& operator[](std::size_t index);
    const double& operator[](std::size_t index) const;

    Vector& operator+=(const Vector& other);
    Vector& operator-=(const Vector& other);
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator-() const;
    Vector& operator*=(double scalar);
    Vector operator*(double scalar) const;
    // Division by zero throws std::domain_error; no scalar / vector overload.
    Vector& operator/=(double scalar);
    Vector operator/(double scalar) const;

    double dot(const Vector& other) const;
    double squared_norm() const; // Sum of squares, not elementwise squaring.
    double norm() const;
    Vector cross(const Vector& other) const; // Both operands must have size 3.
    Vector slice(std::size_t start, std::size_t stop) const;

private:
    void require_same_size(const Vector& other) const;
    std::vector<double> data_;
};

// Scalar on the left: 2.0 * v.
Vector operator*(double scalar, const Vector& vector);

} // namespace hdmath
