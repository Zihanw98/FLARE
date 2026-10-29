#pragma once
#include "riemann.h"
namespace hydro
{
    struct InterfaceState
    {   
        //Those need to be recorded when sampling interface...
        hd_core::Primitives W;
        double gamma;
    };
    class ExactSolver final: public RiemannSolver
    {
        public:
            explicit ExactSolver(
                double tolerance = 1e-8,
                int max_iterations = 100
            );
            hd_core::Conserved flux(
                const hd_core::Primitives& left,
                const hd_core::Primitives& right,
                int direction,
                double gamma_L,
                double gamma_R
            ) const override;
    
        private:
        double tolerance_;
        int max_iterations_;
    };
}