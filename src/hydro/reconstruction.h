/**
 * @file reconstruction.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"
#include "mesh.h"
#include "limiter.h"

namespace hydro
{
    /** @brief Piecewise-constant or piecewise-linear face reconstruction. */
    enum struct Reconstructor
    {
        PCM, ///< Use cell-centered primitive states.
        PLM, ///< Linearly reconstruct primitive states with optional slope limiting.
    };
    /**
     * @brief Select PLM by name, otherwise fall back to PCM.
     * @param[in] name Reconstruction name compared to PLM without case sensitivity.
     * @return Reconstructor::PLM for PLM; Reconstructor::PCM for every other string.
     * @note Unknown names are not rejected by the current implementation.
     */
    Reconstructor parse_reconstructor(const std::string& name);
    /**
     * @brief Reconstruct left and right primitive states at a cell interface.
     * @param[in] reconstruct Reconstruction method to dispatch.
     * @param[in] mesh Mesh containing initialized cell states and face geometry.
     * @param[in] direction Axis normal to the interface, in [0, hd_core::NDIM).
     * @param[in] right_index Storage coordinates of the cell immediately to the right of the interface.
     * @param[in] limiter Slope limiter; unused by PCM.
     * @param[in] limiter_beta MCbeta parameter in [1, 2] when used; unused by PCM.
     * @return Pair of left and right interface states, in that order.
     * @pre direction is valid and all stencil cells, including ghosts, are initialized; the direction is not checked before indexing.
     * @note PLM uses two ghost layers and falls back independently to the cell-centered state if an extrapolated side fails physical validation.
     * @throws std::runtime_error If the selected reconstruction is invalid or PLM has fewer than two ghost layers.
     */
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct(const Reconstructor& reconstruct, const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta=1.5);
    /**
     * @brief Use adjacent cell-centered primitive states at an interface.
     * @param[in] mesh Mesh containing initialized cell states and face geometry.
     * @param[in] direction Axis normal to the interface, in [0, hd_core::NDIM).
     * @param[in] right_index Storage coordinates of the cell immediately to the right of the interface.
     * @param[in] limiter Slope limiter; unused by PCM.
     * @param[in] limiter_beta MCbeta parameter in [1, 2] when used; unused by PCM.
     * @return Pair of left and right interface states, in that order.
     * @pre direction is valid and all stencil cells, including ghosts, are initialized; the direction is not checked before indexing.
     * @note Needs the adjacent left and right cells; does not use limiter settings.
     * @throws std::out_of_range If an adjacent cell is outside mesh storage.
     */
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PCM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta=1.5);
    /**
     * @brief Reconstruct primitive states with limited piecewise-linear slopes.
     * @param[in] mesh Mesh containing initialized cell states and face geometry.
     * @param[in] direction Axis normal to the interface, in [0, hd_core::NDIM).
     * @param[in] right_index Storage coordinates of the cell immediately to the right of the interface.
     * @param[in] limiter Slope limiter; unused by PCM.
     * @param[in] limiter_beta MCbeta parameter in [1, 2] when used; unused by PCM.
     * @return Pair of left and right interface states, in that order.
     * @pre direction is valid and all stencil cells, including ghosts, are initialized; the direction is not checked before indexing.
     * @note PLM uses two ghost layers and falls back independently to the cell-centered state if an extrapolated side fails physical validation.
     * @throws std::runtime_error If the selected reconstruction is invalid or PLM has fewer than two ghost layers.
     */
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PLM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta=1.5);
}