/**
 * @file godunov.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"
#include "star_pressure_cache.h"
#include "reconstruction.h"
#include "limiter.h"
#include "riemann/riemann.h"
#include <vector>

namespace hd_core
{
    struct MeshBlock;
}
namespace hydro
{
    class RiemannSolver;
    // ghost 状态必须预先填好，rhs 必须已按 mesh.cells.size() 分配。
    // 重置并更新 wave_timestep 中的本次 CFL 限制。
    /**
     * @brief Compute the Godunov flux-divergence RHS and record the CFL rate.
     * @param[in] mesh Allocated mesh with valid geometry and ghost states for this stage.
     * @param[in] solver Riemann solver that writes wave_timestep.face_speed on every face.
     * @param[out] rhs Preallocated array with mesh.cells.size() entries; only active entries are overwritten.
     * @param[in] reconstructor PCM or PLM face reconstruction.
     * @param[in] limiter Slope limiter used by PLM; ignored by PCM.
     * @param[in] limiter_beta MCbeta parameter in [1, 2] when that limiter is used.
     * @pre Ghost states are filled at the current stage time; PCM needs at least one ghost layer and PLM at least two.
     * @throws std::invalid_argument If mesh/storage or a scheme argument is invalid.
     * @throws std::runtime_error If the solver reports invalid wave speeds or reconstruction fails.
     * @note Resets and updates the shared wave_timestep record. This shared state is intended for serial mesh evaluation; solver and geometry exceptions propagate.
     */
    void compute_godunov_rhs(
        const hd_core::MeshBlock& mesh,
        const RiemannSolver& solver,
        std::vector<hd_core::Conserved>& rhs,
        const Reconstructor& reconstructor = Reconstructor::PCM,
        const Limiter& limiter = Limiter::none,
        const double limiter_beta = 1.5,
        hydro::StarPressureCache* star_pressure_cache = nullptr
    );
}
