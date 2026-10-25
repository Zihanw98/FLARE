#include "mesh.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hd_core
{

int MeshBlock::total(int d) const
{
    if (d < 0 || d >= NDIM)
        throw std::out_of_range("Invalid mesh direction");

    if (nactive[d] <= 0 || nghost[d] < 0)
        throw std::invalid_argument("Invalid mesh size");
    if (nghost[d] > (std::numeric_limits<int>::max() - nactive[d]) / 2)
        throw std::overflow_error("Mesh size overflow");

    return nactive[d] + 2 * nghost[d];
}

double MeshBlock::spacing(int d) const
{
    if (d < 0 || d >= NDIM)
        throw std::out_of_range("Invalid mesh direction");

    if (nactive[d] <= 0)
        throw std::invalid_argument("Active cell count must be positive");
    if (!std::isfinite(length[d]) || length[d] <= 0.0)
        throw std::invalid_argument("Domain length must be finite and positive");
    const double result = length[d] / nactive[d];
    if (!std::isfinite(result) || result <= 0.0)
        throw std::invalid_argument("Invalid derived cell spacing");
    return result;
}

void MeshBlock::allocate()
{
    std::size_t count = 1;

    for (int d = 0; d < NDIM; ++d)
    {
        (void)spacing(d); // Also reject underflow to a zero cell width.
        if (nactive[d] <= 0 || nghost[d] < 0)
            throw std::invalid_argument("Invalid mesh size");

        if (!std::isfinite(length[d]) || length[d] <= 0.0)
            throw std::invalid_argument("Invalid domain length");

        if (nghost[d] >
            (std::numeric_limits<int>::max() - nactive[d]) / 2)
            throw std::overflow_error("Mesh size overflow");

        const auto n = static_cast<std::size_t>(total(d));

        if (count > cells.max_size() / n)
            throw std::length_error("Too many mesh cells");

        count *= n;
    }

    cells.resize(count);
}

}
