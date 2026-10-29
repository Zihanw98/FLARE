#include "mesh.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hd_core
{
namespace
{
void validate_axis(const std::vector<double>& faces, int count)
{
    if (faces.size() != static_cast<std::size_t>(count) + 1)
        throw std::invalid_argument("Face count must equal total cells + 1");

    for (std::size_t j = 0; j < faces.size(); ++j)
    {
        if (!std::isfinite(faces[j]))
            throw std::invalid_argument("Face coordinate must be finite");
        if (j > 0)
        {
            const double width = faces[j] - faces[j - 1];
            if (!std::isfinite(width) || width <= 0.0)
                throw std::invalid_argument("Faces must have finite positive separation");
        }
    }
}
}

int MeshBlock::ntotal(int d) const
{
    if (d < 0 || d >= NDIM)
        throw std::out_of_range("Invalid mesh direction");

    if (nactive[d] <= 0 || nghost[d] < 0)
        throw std::invalid_argument("Invalid mesh size");
    if (nghost[d] > (std::numeric_limits<int>::max() - nactive[d]) / 2)
        throw std::overflow_error("Mesh size overflow");

    return nactive[d] + 2 * nghost[d];
}

void MeshBlock::set_uniform_faces(int d, double lower, double upper)
{
    const int count = ntotal(d); // Checks direction, sizes, and integer overflow.
    if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower)
        throw std::invalid_argument("Uniform grid bounds must be finite and increasing");

    const double length = upper - lower;
    const double width = length / nactive[d];
    if (!std::isfinite(length) || !std::isfinite(width) || width <= 0.0)
        throw std::invalid_argument("Invalid uniform grid spacing");

    std::vector<double> coordinates(static_cast<std::size_t>(count) + 1);
    const auto first = static_cast<std::size_t>(nghost[d]);
    const auto last = first + static_cast<std::size_t>(nactive[d]);
    for (std::size_t j = 0; j < coordinates.size(); ++j)
    {
        if (j == first)
            coordinates[j] = lower;
        else if (j == last)
            coordinates[j] = upper;
        else
            coordinates[j] = std::fma(static_cast<double>(j) - nghost[d], width, lower);
    }
    // Detect overflow in ghost coordinates or collapsed faces due to rounding.
    validate_axis(coordinates, count);
    face[d].swap(coordinates); // Leave previous geometry unchanged on failure.
}

double MeshBlock::cell_width(int d, int i) const
{
    const int count = ntotal(d);
    if (i < 0 || i >= count)
        throw std::out_of_range("Invalid mesh cell index");
    if (face[d].size() != static_cast<std::size_t>(count) + 1)
        throw std::invalid_argument("Face count must equal total cells + 1");

    const auto j = static_cast<std::size_t>(i);
    const double lower = face[d][j];
    const double upper = face[d][j + 1];
    const double width = upper - lower;
    if (!std::isfinite(lower) || !std::isfinite(upper) ||
        !std::isfinite(width) || width <= 0.0)
        throw std::invalid_argument("Cell faces must have finite positive separation");
    return width;
}

double MeshBlock::cell_center(int d, int i) const
{
    const double width = cell_width(d, i);
    return face[d][static_cast<std::size_t>(i)] + 0.5 * width;
}

double MeshBlock::active_length(int d) const
{
    const int count = ntotal(d);
    validate_axis(face[d], count);
    const auto first = static_cast<std::size_t>(nghost[d]);
    const auto last = first + static_cast<std::size_t>(nactive[d]);
    const double length = face[d][last] - face[d][first];
    if (!std::isfinite(length) || length <= 0.0)
        throw std::invalid_argument("Active domain length must be finite and positive");
    return length;
}

void MeshBlock::validate_geometry() const
{
    for (int d = 0; d < NDIM; ++d)
        (void)active_length(d); // Validates all faces, including ghost faces.
}

void MeshBlock::allocate()
{
    validate_geometry();
    std::size_t count = 1;

    for (int d = 0; d < NDIM; ++d)
    {
        const auto n = static_cast<std::size_t>(ntotal(d));

        if (count > cells.max_size() / n)
            throw std::length_error("Too many mesh cells");

        count *= n;
    }

    cells.resize(count);
}

}
