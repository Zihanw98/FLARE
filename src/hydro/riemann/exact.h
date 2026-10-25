#pragma once
#include "riemann.h"
namespace hydro
{
    class ExactSolver final: public RiemannSolver
    {
        public:
            explicit ExactSolver(
                double tolerance = 1e-8,
                int max_iterations = 100
            );
    
        private:
        double tolerance_;
        int max_iterations_;
    };
}