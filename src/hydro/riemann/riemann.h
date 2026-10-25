#pragma once
#include "core.h"
namespace hydro
{
    class RiemannSolver
    {
        public:
            virtual ~RiemannSolver() = default;
            virtual hd_core::Conserved flux(
                const hd_core::Primitives& left,
                const hd_core::Primitives& right,
                int direction,
                double gamma_L,
                double gamma_R
                //later, we will add an EoS parameter here to support different EoS
            ) const = 0;
    };
}