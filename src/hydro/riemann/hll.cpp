#include "hll.h"
#include "core/cell.h"
#include "hydro/flux.h"
#include "microphysics/eos/eos.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hydro
{
hd_core::Conserved HLLSolver::flux(const hd_core::Primitives& left,
    const hd_core::Primitives& right, int direction, double gamma_L, double gamma_R,
    double& face_speed, double& face_velocity) const
{
    using namespace hd_core;
    if (direction < 0 || direction >= NDIM)
        throw std::out_of_range("Invalid HLL flux direction");
    if (left[RHO] <= 0.0 || right[RHO] <= 0.0)
        throw std::invalid_argument("HLL requires positive-density states on both sides");
    Cell L, R;
    L.gamma = gamma_L;
    R.gamma = gamma_R;
    L.Set_Conserved_From_Primitives(left);
    R.Set_Conserved_From_Primitives(right);
    hd_microphysics::refresh_thermodynamics(L);
    hd_microphysics::refresh_thermodynamics(R);
    const double vl = left[VEL(direction)], vr = right[VEL(direction)];
    const double cl = L.get_sound_speed(), cr = R.get_sound_speed();
    const double sl = std::min(vl - cl, vr - cr);
    const double sr = std::max(vl + cl, vr + cr);
    face_speed = std::max(std::abs(sl), std::abs(sr));
    const auto fl = physical_flux(left, L.U, direction);
    const auto fr = physical_flux(right, R.U, direction);
    if (sl >= 0.0) { face_velocity = vl; return fl; }
    if (sr <= 0.0) { face_velocity = vr; return fr; }
    const auto middle = (sr * R.U - sl * L.U + fl - fr) / (sr - sl);
    if (!std::isfinite(middle[RHO]) || middle[RHO] <= 0.0)
        throw std::invalid_argument("HLL intermediate density is not positive");
    face_velocity = middle[MOMDEN(direction)] / middle[RHO];
    return (sr * fl - sl * fr + (sl * sr) * (R.U - L.U)) / (sr - sl);
}
}
