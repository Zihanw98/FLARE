#pragma once
#include "core.h"
namespace hydro
{
    class RiemannSolver
    {
        public:
            virtual ~RiemannSolver() = default;
            // 实现必须通过 wave_timestep.h 更新 wave_timestep.face_speed，
            // 包括返回零通量的分支；网格计算紧接着读取此值。
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
