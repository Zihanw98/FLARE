#pragma once
#include "core.h"
namespace hd_core
{
    struct Cell
    {
        Conserved U;
        double gamma = 5.0/3;  
        Cell() = default;
        explicit Cell(const Conserved& state, double g = 5.0/3)
            : U(state), gamma(g)
        {
            U.validate(gamma);
        }
        explicit Cell(const Primitives& state, double g = 5.0/3)
            : U(state.to_Conserved(g)), gamma(g)
        {
        }
    };
}
