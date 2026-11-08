/**
 * @file exact.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "riemann.h"
namespace hydro
{
    /** @brief Sampled primitive interface state and its originating adiabatic index. */
    struct InterfaceState
    {   
        //Those need to be recorded when sampling interface...
        hd_core::Primitives W;
        double gamma;
    };
    // Bisection is retained as a numerical/performance reference for tests.
    enum class ExactPressureMethod { SafeguardedNewton, Bisection };
    /** @brief Exact ideal-gas Riemann solver with gas/vacuum interface handling. */
    class ExactSolver final: public RiemannSolver
    {
        public:
            /**
             * @brief Construct an exact ideal-gas Riemann solver.
             * @param[in] tolerance Finite positive relative pressure-bracket tolerance.
             * @param[in] max_iterations Positive maximum number of pressure iterations.
             * @param[in] method Safeguarded Newton (default) or reference bisection.
             * @throws std::invalid_argument If tolerance, iteration count or method is invalid.
             */
            explicit ExactSolver(
                double tolerance = 1e-8,
                int max_iterations = 100,
                ExactPressureMethod method = ExactPressureMethod::SafeguardedNewton
            );
            /**
             * @brief Solve the ideal-gas Riemann problem and return the stationary-face flux.
             * @param[in] left Primitive state immediately to the left of the face.
             * @param[in] right Primitive state immediately to the right of the face.
             * @param[in] direction Face-normal axis in [0, hd_core::NDIM).
             * @param[in] gamma_L Left ideal-gas adiabatic index greater than one.
             * @param[in] gamma_R Right ideal-gas adiabatic index greater than one.
             * @return Conserved flux in the positive axis direction; the exact solver returns zero when the sampled interface is vacuum.
             * @pre Each input is a valid gas or vacuum state for its adiabatic index.
             * @post Every successful call updates wave_timestep.face_speed, including zero-flux cases.
             * @note The shared wave-speed record is not thread-safe. The exact implementation leaves extra-component fluxes zero.
             * @throws std::out_of_range If direction is invalid.
             * @throws std::invalid_argument If an input state or adiabatic index is invalid.
             * @throws std::runtime_error If star pressure fails to converge or wave-speed bounds are invalid.
             */
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
        ExactPressureMethod method_;
    };
}
