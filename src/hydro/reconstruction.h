/**
 * @file reconstruction.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"
#include "mesh.h"
#include <cstdint>
#include "limiter.h"
#include <array>
#include <vector>

namespace hydro
{
    /** @brief Piecewise-constant or piecewise-linear face reconstruction. */
    enum struct Reconstructor
    {
        PCM, ///< Use cell-centered primitive states.
        PLM, ///< Linearly reconstruct primitive states with optional slope limiting.
        PLM_corrected, ///< PLM with a shared positive density/pressure slope factor.
    };
    /** Per-stage lazy primitive/slope cache. The mesh must remain unchanged
     * between reset() calls; reset after every boundary fill/stage/retry.
     * Only stencil cells are converted (unused multidimensional ghost corners
     * need not contain valid states). Direct reconstruction APIs remain available.
     */
    class ReconstructionCache
    {
    public:
        ReconstructionCache(const hd_core::MeshBlock& mesh, Reconstructor method,
                            Limiter limiter, double beta);
        ReconstructionCache(const ReconstructionCache&) = delete;
        ReconstructionCache& operator=(const ReconstructionCache&) = delete;
        /** Rebind a new immutable stage while retaining allocation capacity. */
        void reset(const hd_core::MeshBlock& mesh, Reconstructor method,
                   Limiter limiter, double beta);
        std::pair<hd_core::Primitives, hd_core::Primitives> face(
            int direction, const hd_core::Index& right_index);
    private:
        // cell_flat indexes mesh.cells; axis_index is the storage coordinate along direction.
        const hd_core::Primitives& primitive(std::size_t cell_flat);
        const hd_core::Primitives& slope(int direction, std::size_t cell_flat, int axis_index);
        const hd_core::MeshBlock* mesh_;
        Reconstructor method_;
        Limiter limiter_;
        double beta_;
        std::vector<hd_core::Primitives> primitives_;
        std::vector<std::uint64_t> primitive_epoch_;
        std::array<std::vector<hd_core::Primitives>, hd_core::NDIM> slopes_;
        std::array<std::vector<std::uint64_t>, hd_core::NDIM> slope_epoch_;
        std::uint64_t epoch_ = 0;
        std::array<std::vector<double>, hd_core::NDIM> centers_, widths_;
        // Change in the flattened storage offset when moving one cell along each axis.
        std::array<std::size_t, hd_core::NDIM> strides_{};
    };
    /**
     * @brief Select PLM or PLM_corrected by name, otherwise fall back to PCM.
     * @param[in] name Reconstruction name compared without case sensitivity.
     * @return Matching PLM variant, or PCM for every other string.
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
    /** @brief PLM with one slope factor per cell and direction shared by all
     * hydrodynamic primitives and both faces. For positive cell centers, density
     * and pressure retain at least 1e-10 of their center value (up to roundoff).
     * This constrains face states, not the positivity of a subsequent time step.
     * Validation fallback remains a safeguard for vacuum/numerical exceptions.
     */
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PLM_corrected(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta=1.5);
}
