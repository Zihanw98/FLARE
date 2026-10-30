#include "flux.h"


namespace hydro
{
    hd_core::Conserved physical_flux(
        const hd_core::Primitives& W,
        int direction,
        double gamma
    )
    {
        if (direction < 0 || direction >= hd_core::NDIM)
            throw std::out_of_range("Invalid flux direction");
        W.validate(gamma);
        hd_core::Conserved F{};
        if(W[hd_core::RHO]==0.0)
        {
            return F;
        }
        const hd_core::Conserved U = W.to_Conserved(gamma);
        const double v_n = W[hd_core::VEL(direction)];
        const double p = W[hd_core::PRE];
        F[hd_core::RHO] = U[hd_core::RHO] * v_n;
        for(int d=0;d < hd_core::NDIM; d++)
        {
            F[hd_core::MOMDEN(d)] = U[hd_core::MOMDEN(d)] * v_n + (d==direction? p : 0 ); 
        }
        F[hd_core::ENERGY] = (U[hd_core::ENERGY] + p) * v_n;
        return F;
    }
}