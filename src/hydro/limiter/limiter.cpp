/**
 * @file limiter.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "limiter.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hydro
{
    Limiter parse_limiter(const std::string& name)
    {
        if (utils::iequals(name, "none")) return Limiter::none;
        if (utils::iequals(name, "minmod")) return Limiter::minmod;
        if (utils::iequals(name, "mc")) return Limiter::MC;
        if (utils::iequals(name, "vanleer") || utils::iequals(name, "van_leer"))
            return Limiter::van_Leer;
        if (utils::iequals(name, "mcbeta") || utils::iequals(name, "mc_beta"))
            return Limiter::MCbeta;
        if (utils::iequals(name, "superbee")) return Limiter::superbee;
        throw std::invalid_argument("Unknown limiter: " + name
            + "; expected none, minmod, mc, vanleer, mcbeta, or superbee");
    }

    double limited_slope(double left, double right, Limiter limiter, double beta)
    {
        if (!std::isfinite(left) || !std::isfinite(right))
            throw std::invalid_argument("Limiter slopes must be finite");
        switch (limiter)
        {
            case Limiter::none:
                // Unlimited centered slope, including opposite signs or a zero side.
                // Halve first to avoid overflowing left + right.
                return 0.5 * left + 0.5 * right;
            case Limiter::minmod:
            case Limiter::MC:
            case Limiter::van_Leer:
            case Limiter::superbee:
                break;
            case Limiter::MCbeta:
                if (!std::isfinite(beta) || beta < 1.0 || beta > 2.0)
                    throw std::invalid_argument("limiter_beta must be finite and in [1, 2]");
                break;
            default:
                throw std::invalid_argument("Invalid limiter enum");
        }
        // Compare signs directly: multiplying slopes can overflow or underflow.
        if (!((left > 0.0 && right > 0.0) || (left < 0.0 && right < 0.0)))
            return 0.0;

        const double small = std::min(std::abs(left), std::abs(right));
        const double large = std::max(std::abs(left), std::abs(right));
        const auto capped_product = [small](double factor, double cap) {
            return small >= cap / factor ? cap : factor * small;
        };
        double magnitude = small;
        switch (limiter)
        {
            case Limiter::none: // Already returned above.
                break;
            case Limiter::minmod:
                break;
            case Limiter::MC:
                magnitude = capped_product(2.0, small + 0.5 * (large - small));
                break;
            case Limiter::van_Leer:
                magnitude = small / (0.5 + 0.5 * (small / large));
                break;
            case Limiter::MCbeta:
                magnitude = capped_product(beta, small + 0.5 * (large - small));
                break;
            case Limiter::superbee:
                magnitude = capped_product(2.0, large);
                break;
        }
        return std::copysign(magnitude, left);
    }
}
