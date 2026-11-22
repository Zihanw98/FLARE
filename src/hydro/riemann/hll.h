#pragma once
#include "riemann.h"

namespace hydro
{
    // Two-wave flux for EoS models without an analytic exact Riemann solution.
    // Wave speeds use the selected EoS. No star-pressure iteration is performed.
    class HLLSolver final : public RiemannSolver
    {
    public:
        hd_core::Conserved flux(const hd_core::Primitives& left,
            const hd_core::Primitives& right, int direction,
            double gamma_L, double gamma_R,
            double& face_speed, double& face_velocity) const override;
    };
}
