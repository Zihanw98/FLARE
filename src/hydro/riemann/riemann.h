/**
 * @file riemann.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"
#include <array>
#include <vector>
namespace hydro
{
    /** @brief Polymorphic interface for numerical fluxes and interface wave-speed reporting. */
    class RiemannSolver
    {
        public:
            /**
             * @brief Destroy a solver safely through the polymorphic base interface.
             */
            virtual ~RiemannSolver() = default;
            /**
             * @brief Compute a numerical flux and return caller-owned interface measurements.
             * @param[in] left Primitive state immediately to the left of the face.
             * @param[in] right Primitive state immediately to the right of the face.
             * @param[in] direction Face-normal axis in [0, hd_core::NDIM).
             * @param[in] gamma_L Left ideal-gas adiabatic index greater than one.
             * @param[in] gamma_R Right ideal-gas adiabatic index greater than one.
             * @param[out] face_speed Maximum absolute outer wave speed, for CFL control.
             * @param[out] face_velocity Signed sampled normal fluid velocity, for pressure work.
             * @return Conserved flux in the positive axis direction; the exact solver returns zero when the sampled interface is vacuum.
             * @pre Each input is a valid gas or vacuum state for its adiabatic index.
             * @post Every successful call writes face_speed and face_velocity, including zero-flux cases.
             * @note Output references must refer to distinct caller-owned values.
             */
            virtual hd_core::Conserved flux(
                const hd_core::Primitives& left,
                const hd_core::Primitives& right,
                int direction,
                double gamma_L,
                double gamma_R,
                double& face_speed,
                double& face_velocity
                //later, we will add an EoS parameter here to support different EoS
            ) const = 0;

            /** Optional per-face pressure guess, owned by the caller. On success
             * replace it with the solved star pressure (zero if unavailable).
             * Other solvers retain their ordinary flux implementation. */
            virtual hd_core::Conserved flux_with_pressure_hint(
                const hd_core::Primitives& left, const hd_core::Primitives& right,
                int direction, double gamma_L, double gamma_R, double& pressure_hint,
                double& face_speed, double& face_velocity) const
            {
                const auto result = flux(left, right, direction, gamma_L, gamma_R, face_speed, face_velocity);
                pressure_hint = 0.0;
                return result;
            }
    };
}
