/**
 * @file mesh.h
 * @brief Structured mesh geometry and cell indexing, including ghost regions.
 * @note Storage coordinates include ghost offsets; active coordinates start at zero.
 *       Flat offsets address the complete cells array, with axis 0 varying fastest.
 *       Independent face coordinates on each axis define a tensor-product grid.
 * @see hd_core::MeshBlock
 */
#pragma once
#include "core.h"
#include "cell.h"
#include "facefield.h"
#include "boundary.h"
#include <array>
#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace config { struct SimulationConfig; }
/// @brief Physical state and structured-mesh data types.
namespace hd_core
{
    enum struct Geometry
    {
        Cartesian,
        Cylindrical,
        Spherical
    };

    /**
     * @brief Store mesh geometry, boundary conditions, and all cell states.
     * @pre Set nactive, nghost, and faces before calling allocate().
     * @note Filling ghost states does not update face coordinates. Changing an existing
     *       grid requires reinitialization or conservative remapping of physical data.
     */
    struct MeshBlock
    {

        Geometry geometry = Geometry::Cartesian;
        Index nactive{}; ///< Positive active-cell counts on each axis.
        Index nghost{};  ///< Nonnegative ghost-cell counts per side on each axis.

        /**
         * @brief Return the cell count on one axis, including both ghost regions.
         * @param[in] d Axis index in [0, NDIM).
         * @return nactive[d] + 2 * nghost[d].
         * @throws std::out_of_range If the axis is invalid.
         * @throws std::invalid_argument If nactive is nonpositive or nghost is negative.
         * @throws std::overflow_error If the total exceeds the int range.
         */
        int ntotal(int d) const;

        /**
         * @brief Return the cell counts on all axes, including ghost regions.
         * @return An Index whose components are ntotal(d).
         * @note Propagates the per-axis validation errors from ntotal(int).
         */
        Index ntotal() const;
        
        /**
         * @brief Face coordinates on each axis, including both ghost regions.
         * @note face[d] has ntotal(d) + 1 finite, strictly increasing entries.
         *       Active faces run from nghost[d] through nghost[d] + nactive[d], inclusive.
         */
        std::array<std::vector<double>,NDIM> face;

        // Full storage including ghosts, axis 0 fastest. Direction d has shape
        // ntotal(), with one extra entry on axis d; not the hydro active-face layout.
        // at(I) selects the lower d-face of storage cell I; at(I + e_d)
        // selects its upper face. Unlike hydro fluxes, I includes ghost offsets.
        std::array<FaceField<double>,NDIM> face_area;

        // Same layout as cells. Cartesian: unit omitted dimensions; spherical:
        // full 4*pi shells, with signed-radius ghosts extended using r^2 dr.
        std::vector<double> volume_cell;

        std::array<boundary::BoundaryCondition,NDIM> bc_upper; ///< Upper boundary condition on each axis.
        std::array<boundary::BoundaryCondition,NDIM> bc_lower; ///< Lower boundary condition on each axis.

        /**
         * @brief Create uniformly spaced faces on one axis, including ghost faces.
         * @param[in] d Axis index in [0, NDIM).
         * @param[in] lower Finite lower boundary of the active domain.
         * @param[in] upper Finite upper boundary, greater than lower.
         * @pre nactive and nghost have been set for this axis.
         * @post On success, face[d] contains ntotal(d) + 1 finite, increasing coordinates.
         * @throws std::invalid_argument If bounds or generated spacings are invalid.
         * @note Preserves previous faces on failure; does not allocate cells or initialize physical states. Axis/count errors from ntotal() propagate.
         */
        void set_linear_uniform_faces(int d, double lower, double upper);
        void set_log_uniform_faces(int d,double lower, double upper);
        /**
         * @brief Return the cell width along a selected axis.
         * @param[in] direction Axis index in [0, NDIM).
         * @param[in] axis_index Cell storage index along this axis, including the ghost offset.
         * @return The positive distance between adjacent faces.
         * @throws std::out_of_range If the axis or selected cell index is invalid.
         * @throws std::invalid_argument If mesh counts, face count, or selected face spacing are invalid.
         * @note Uses storage coordinates, including ghost cells; the index is not a flattened cell offset.
         */
        double cell_width(int direction, int axis_index) const;

        /**
         * @brief Return the cell center along a selected axis.
         * @param[in] direction Axis index in [0, NDIM).
         * @param[in] axis_index Cell storage index along this axis, including the ghost offset.
         * @return The coordinate midpoint of the two adjacent faces.
         * @throws std::out_of_range If the axis or selected cell index is invalid.
         * @throws std::invalid_argument If mesh counts, face count, or selected face spacing are invalid.
         * @note Uses storage coordinates, including ghost cells; the index is not a flattened cell offset.
         */
        double cell_center(int direction, int axis_index) const;

        double cell_volume_center(int direction, int axis_index) const;
        
        std::pair<double,double> cell_face_pos(int direction, int axis_index) const;        
        /**
         * @brief Return the cell width along a selected axis.
         * @param[in] direction Axis index in [0, NDIM).
         * @param[in] storage_index Multidimensional storage coordinates; only the selected axis component is used.
         * @return The positive distance between adjacent faces.
         * @throws std::out_of_range If the axis or selected cell index is invalid.
         * @throws std::invalid_argument If mesh counts, face count, or selected face spacing are invalid.
         * @note Uses storage coordinates, including ghost cells; the index is not a flattened cell offset.
         */
        double cell_width(int direction, const Index& storage_index) const;

        /**
         * @brief Return the cell center along a selected axis.
         * @param[in] direction Axis index in [0, NDIM).
         * @param[in] storage_index Multidimensional storage coordinates; only the selected axis component is used.
         * @return The coordinate midpoint of the two adjacent faces.
         * @throws std::out_of_range If the axis or selected cell index is invalid.
         * @throws std::invalid_argument If mesh counts, face count, or selected face spacing are invalid.
         * @note Uses storage coordinates, including ghost cells; the index is not a flattened cell offset.
         */
        double cell_center(int direction, const Index& storage_index) const;

        double cell_volume_center(int direction, const Index& storage_index) const;

        std::pair<double,double> cell_face_pos(int direction, const Index& storage_index) const;        
        /**
         * @brief Measure the active domain length on one axis.
         * @param[in] d Axis index in [0, NDIM).
         * @return The upper active face coordinate minus the lower active face coordinate.
         * @throws std::invalid_argument If faces or the active length are invalid.
         * @note Validates every face on this axis, including ghost faces; ntotal() errors propagate.
         */
        double active_length(int d) const;

        /**
         * @brief Validate dimensions, coordinates and supported coordinate-system ranges.
         * @throws std::invalid_argument If coordinates are invalid or geometry is unsupported.
         * @note Checks all ghost coordinates for monotonicity, but radial nonnegativity
         *       only in the active domain. Does not require caches or allocated cells;
         *       build_geometry() validates computed areas and volumes before storing them.
         */
        void validate_geometry() const;

        /**
         * @brief Validate every active conserved state using its cell adiabatic index.
         * @pre Cell storage matches the current mesh dimensions.
         * @throws std::invalid_argument If dimensions or an active physical state are invalid.
         * @throws std::logic_error If storage does not match the mesh dimensions.
         * @note Skips ghost states and leaves all data unchanged; geometry is checked separately by validate_geometry().
         */
        void validate_active() const;

        /**
         * @brief Validate geometry and resize cell storage to include all ghost cells.
         * @pre nactive, nghost, and faces are initialized on every axis.
         * @throws std::length_error If the cell count exceeds the container capacity.
         * @note Preserves existing elements when resizing and default-constructs new cells. Does not fill initial states; geometry validation errors propagate.
         */
        void allocate();

        std::vector<Cell> cells; ///< All cells including ghosts, with axis 0 varying fastest.

        /**
         * @brief Access a cell using storage coordinates.
         * @param[in] coordinates Exactly NDIM storage coordinates, including ghost offsets.
         * @return A mutable reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note An initializer list must have exactly NDIM entries; otherwise std::invalid_argument is thrown.
         */
        Cell& at(std::initializer_list<int> coordinates);

        /**
         * @brief Read a cell using storage coordinates.
         * @param[in] coordinates Exactly NDIM storage coordinates, including ghost offsets.
         * @return A read-only reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note An initializer list must have exactly NDIM entries; otherwise std::invalid_argument is thrown.
         */
        const Cell& at(std::initializer_list<int> coordinates) const;

        /**
         * @brief Access an active cell using local active coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, nactive[d]); ghost offsets are added internally.
         * @return A mutable reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note An initializer list must have exactly NDIM entries; otherwise std::invalid_argument is thrown.
         */
        Cell& at_active(std::initializer_list<int> coordinates);

        /**
         * @brief Read an active cell using local active coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, nactive[d]); ghost offsets are added internally.
         * @return A read-only reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note An initializer list must have exactly NDIM entries; otherwise std::invalid_argument is thrown.
         */
        const Cell& at_active(std::initializer_list<int> coordinates) const;

        /**
         * @brief Convert a flattened cell offset to storage coordinates.
         * @param[in] flat Offset in cells, including ghost cells.
         * @return Storage coordinates with axis 0 varying fastest.
         * @throws std::out_of_range If flat is outside cells.
         * @throws std::logic_error If storage and dimensions disagree.
         * @throws std::length_error If dimensions exceed the container capacity.
         */
        Index coordinates(std::size_t flat) const;

        /**
         * @brief Convert a flattened cell offset to coordinates relative to the active origin.
         * @param[in] flat Offset in the full cells array, including ghost cells.
         * @return coordinates(flat) with nghost subtracted on every axis.
         * @note Ghost offsets are accepted and return coordinates outside the active range; errors from coordinates() propagate.
         */
        Index coordinates_active(std::size_t flat) const;

        /**
         * @brief Access a cell using storage coordinates.
         * @param[in] coordinates Exactly NDIM storage coordinates, including ghost offsets.
         * @return A mutable reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note The coordinate tuple is checked by flat_index().
         */
        Cell& at(const Index& coordinates);

        /**
         * @brief Read a cell using storage coordinates.
         * @param[in] coordinates Exactly NDIM storage coordinates, including ghost offsets.
         * @return A read-only reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note The coordinate tuple is checked by flat_index().
         */
        const Cell& at(const Index& coordinates) const;

        /**
         * @brief Access an active cell using local active coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, nactive[d]); ghost offsets are added internally.
         * @return A mutable reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note The coordinate tuple is checked by flat_index().
         */
        Cell& at_active(const Index& coordinates);

        /**
         * @brief Read an active cell using local active coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, nactive[d]); ghost offsets are added internally.
         * @return A read-only reference to the selected cell.
         * @pre Cell storage matches the mesh dimensions.
         * @throws std::out_of_range If any coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @note The coordinate tuple is checked by flat_index().
         */
        const Cell& at_active(const Index& coordinates) const;

        /**
         * @brief Flatten cell coordinates with axis 0 varying fastest.
         * @param[in] coordinates Storage coordinates, or local active coordinates when active is true.
         * @param[in] active Whether to check active bounds and add ghost offsets; false uses storage bounds.
         * @return Offset in the full cells array, including ghost cells.
         * @pre Cell storage matches the current dimensions.
         * @throws std::out_of_range If a coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @throws std::length_error If dimensions exceed container capacity.
         * @throws std::invalid_argument If the coordinate list does not contain exactly NDIM entries.
         * @note Count validation errors from ntotal() also propagate.
         */
        std::size_t flat_index(
            std::initializer_list<int> coordinates,
            bool active = false
        ) const;

        /**
         * @brief Flatten cell coordinates with axis 0 varying fastest.
         * @param[in] coordinates Storage coordinates, or local active coordinates when active is true.
         * @param[in] active Whether to check active bounds and add ghost offsets; false uses storage bounds.
         * @return Offset in the full cells array, including ghost cells.
         * @pre Cell storage matches the current dimensions.
         * @throws std::out_of_range If a coordinate is outside the selected region.
         * @throws std::logic_error If storage and dimensions disagree.
         * @throws std::length_error If dimensions exceed container capacity.
         * @note Count validation errors from ntotal() also propagate.
         */
        std::size_t flat_index(
            const Index& coordinates,
            bool active = false
        ) const;

        void build_mesh(config::SimulationConfig& cfg);

        /**
         * @brief Allocate and recompute volume_cell and face_area from current faces.
         * @note Supports Cartesian 1D/2D/3D and spherical 1D. Does not modify cells
         *       or require them to be allocated. Call again after changing geometry,
         *       dimensions or faces. Failure preserves the previous caches.
         */
        void build_geometry();

        bool is_active(size_t flat) const;


    };
    Geometry parse_geometry(const std::string& name);
}
