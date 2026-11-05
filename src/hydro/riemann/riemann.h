/**
 * @file riemann.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"
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
            // 实现必须通过 wave_timestep.h 更新 wave_timestep.face_speed，
            // 包括返回零通量的分支；网格计算紧接着读取此值。
            /**
             * @brief Compute a numerical flux and record the interface propagation speed.
             * @param[in] left Primitive state immediately to the left of the face.
             * @param[in] right Primitive state immediately to the right of the face.
             * @param[in] direction Face-normal axis in [0, hd_core::NDIM).
             * @param[in] gamma_L Left ideal-gas adiabatic index greater than one.
             * @param[in] gamma_R Right ideal-gas adiabatic index greater than one.
             * @return Conserved flux in the positive axis direction; the exact solver returns zero when the sampled interface is vacuum.
             * @pre Each input is a valid gas or vacuum state for its adiabatic index.
             * @post Every successful call updates wave_timestep.face_speed, including zero-flux cases.
             * @note The shared wave-speed record is not thread-safe. The exact implementation leaves extra-component fluxes zero.
             */
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
