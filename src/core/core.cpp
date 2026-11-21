#include "cell.h"

namespace hd_core
{
void Primitives::update_internal_energy(double gamma)
{
    Cell context;
    context.gamma = gamma;
    const auto U = context.to_Conserved(*this);
    data[INTERNAL_ENERGY] = data[RHO] == 0.0 ? 0.0 : U[INTERNAL_ENERGY] / data[RHO];
}
}
