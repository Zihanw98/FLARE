/**
 * @file flux.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "flux.h"
#include "core/cell.h"

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
        if(W[hd_core::RHO]==0.0)
        {
            return hd_core::Conserved{};
        }
        hd_core::Cell context;
        context.gamma = gamma;

        const hd_core::Conserved U = context.to_Conserved(W);
        return physical_flux(W, U, direction);
    }

    hd_core::Conserved physical_flux(
        const hd_core::Primitives& W,
        const hd_core::Conserved& U,
        int direction
    )
    {
        if (direction < 0 || direction >= hd_core::NDIM)
            throw std::out_of_range("Invalid flux direction");
        U.validate();
        hd_core::Conserved F{};
        if (U[hd_core::RHO] == 0.0) return F;

        const double v_n = W[hd_core::VEL(direction)];
        const double p = W[hd_core::PRE];
        F[hd_core::RHO] = U[hd_core::RHO] * v_n;
        for(int d=0;d < hd_core::NDIM; d++)
        {
            F[hd_core::MOMDEN(d)] = U[hd_core::MOMDEN(d)] * v_n + (d==direction? p : 0 ); 
        }
        F[hd_core::INTERNAL_ENERGY] = U[hd_core::INTERNAL_ENERGY] * v_n;
        F[hd_core::ENERGY] = (U[hd_core::ENERGY] + p) * v_n;
        for (int k = 0; k < hd_core::NCOMP; ++k)
        {
            const int idx = hd_core::COMP(k);
            F[idx] = F[hd_core::RHO] * W[idx];
        }
        return F;
    }
}
