#pragma once
#include <string>
namespace hydro
{
    enum struct Limiter
    {
        none,
        minmod,
        MC,
        van_Leer,
        MCbeta,
        superbee
    };

    // Parse once during setup; names are case-insensitive.
    Limiter parse_limiter(const std::string& name);

    // Left/right derivatives dq/dx on a uniform grid. Inputs must be finite.
    // none returns the unlimited centered slope, not zero (PCM).
    // beta is used only by MCbeta and must be in [1, 2].
    double limited_slope(double left, double right, Limiter limiter,
                         double beta = 1.5);
}
