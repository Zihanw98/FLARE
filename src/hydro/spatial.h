/**
 * @file spatial.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"
#include "star_pressure_cache.h"
#include "reconstruction.h"
#include "limiter.h"
#include "riemann/riemann.h"
#include <string>
#include <vector>

namespace hd_core
{
    struct MeshBlock;
}

namespace hydro
{
class RiemannSolver;

/** @brief Available hydrodynamic spatial operators. */
enum class SpatialScheme
{
    Godunov ///< Finite-volume Godunov flux divergence.
};
/**
 * @brief Parse the name of an available spatial discretization.
 * @param[in] name Case-insensitive spatial scheme name; currently godunov.
 * @return The selected SpatialScheme enum.
 * @throws std::invalid_argument If the name is unknown or unavailable.
 */
SpatialScheme parse_spatial_scheme(const std::string& name);
// 填充 rhs，同时通过 max_inverse_dt 更新本次波速限制。
/**
 * @brief Initialize the RHS and apply the selected hydrodynamic spatial operator.
 * @param[in] mesh Allocated mesh with valid geometry and ghost states for this stage.
 * @param[in] scheme Spatial discretization to dispatch; currently Godunov.
 * @param[in] solver Riemann solver that returns face speed and velocity through output references on every face.
 * @param[out] rhs Resized and zeroed array; active entries receive the flux divergence.
 * @param[in] reconstructor PCM or PLM face reconstruction.
 * @param[in] limiter Slope limiter used by PLM; ignored by PCM.
 * @param[in] limiter_beta MCbeta parameter in [1, 2] when that limiter is used.
 * @pre Ghost states are filled at the current stage time; PCM needs at least one ghost layer and PLM at least two.
 * @throws std::invalid_argument If mesh/storage or a scheme argument is invalid.
 * @throws std::runtime_error If the solver reports invalid wave speeds or reconstruction fails.
 * @note Resets and updates the shared max_inverse_dt. This shared state is intended for serial mesh evaluation; solver and geometry exceptions propagate.
 */
void compute_hydro_rhs(
    const hd_core::MeshBlock& mesh,
    SpatialScheme scheme,
    const RiemannSolver& solver,
    std::vector<hd_core::Conserved>& rhs,
    const Reconstructor& reconstructor = Reconstructor::PCM,
    const Limiter& limiter = Limiter::none,
    const double limiter_beta = 1.5,
    hydro::StarPressureCache* star_pressure_cache = nullptr
);

}
