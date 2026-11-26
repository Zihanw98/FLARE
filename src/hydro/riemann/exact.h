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
    // Both methods use the same wave functions; bisection is a root-finder reference.
    enum class ExactPressureMethod { SafeguardedNewton, Bisection };
    /**
     * @brief Riemann solver for the configured EoS, with gas/vacuum interfaces.
     * @note Default/Mixed/Ionized retain the original analytic wave/sampling branches.
     * Ionized_Radiative uses Hugoniot for shocks and RK4 integration for
     * rarefactions, assuming a single gas/radiation temperature. X, mu and each side's gas gamma
     * are frozen across its waves; Saha and radiation transport are not solved.
     * Constant specific ionization energy cancels from the wave relations and
     * is restored by the existing EoS conversion when building the energy flux.
     * Nuclear uses the low-density cold exponent supplied on each side and
     * samples isentropic expansion across the nuclear-density transition and into vacuum.
     */
    class ExactSolver final: public RiemannSolver
    {
        public:
            /**
             * @brief Construct a Riemann solver; each flux call reads the global EoS type.
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
             * @brief Solve the frozen-composition Riemann problem and return the stationary-face flux.
             * @param[in] left Primitive state immediately to the left of the face.
             * @param[in] right Primitive state immediately to the right of the face.
             * @param[in] direction Face-normal axis in [0, hd_core::NDIM).
             * @param[in] gamma_L Left gas adiabatic index (low-density cold exponent for Nuclear).
             * @param[in] gamma_R Right gas adiabatic index (low-density cold exponent for Nuclear).
             * @param[out] face_speed Maximum absolute outer wave speed, for CFL control.
             * @param[out] face_velocity Signed sampled normal fluid velocity, for pressure work.
             * @return Conserved flux in the positive axis direction; the exact solver returns zero when the sampled interface is vacuum.
             * @pre Each input is a valid gas or vacuum state. Composition models
             * require an initialized species table and the face mass fractions.
             * For Ionized_Radiative, PRE includes gas and radiation pressure.
             * @post Every successful call writes face_speed and face_velocity, including zero-flux cases.
             * @note Output references must refer to distinct caller-owned values.
             * @throws std::out_of_range If direction is invalid.
             * @throws std::invalid_argument If an input state or adiabatic index is invalid.
             * @throws std::runtime_error If a pressure/density/temperature solve,
             * rarefaction integral or vacuum limit fails to converge, or a wave
             * state exceeds floating-point range. Numerical radiation-wave
             * accuracy has a double-precision floor; very near-isothermal gamma
             * can require a vacuum integration beyond the supported log range.
             */
            hd_core::Conserved flux(
                const hd_core::Primitives& left,
                const hd_core::Primitives& right,
                int direction,
                double gamma_L,
                double gamma_R,
                double& face_speed,
                double& face_velocity
            ) const override;

            /** Same exact solve, optionally starting Newton from the previous
             * star pressure on this face. Invalid/out-of-bracket hints fall back
             * to the ordinary initial guess. Bisection ignores the input hint.
             * The hint is updated only on success; vacuum writes zero. */
            hd_core::Conserved flux_with_pressure_hint(
                const hd_core::Primitives& left, const hd_core::Primitives& right,
                int direction, double gamma_L, double gamma_R, double& pressure_hint,
                double& face_speed, double& face_velocity) const override;
    
        private:
        hd_core::Conserved flux_impl(
            const hd_core::Primitives& left, const hd_core::Primitives& right,
            int direction, double gamma_L, double gamma_R, double* pressure_hint,
            double& face_speed, double& face_velocity) const;
        double tolerance_;
        int max_iterations_;
        ExactPressureMethod method_;
    };
}
