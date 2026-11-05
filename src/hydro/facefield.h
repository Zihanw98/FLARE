/**
 * @file facefield.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "core.h"
#include <cstddef>
#include <initializer_list>
#include <vector>

namespace hydro
{
    /** @brief Dense face-centered storage with axis 0 varying fastest. */
    struct FaceField
    {

        hd_core::Index shape{};

        std::vector<hd_core::Conserved> faces;

        /**
         * @brief Allocate a zero-initialized face field with the requested dimensions.
         * @param[in] dimensions Positive face counts on every axis.
         * @post shape equals dimensions and faces contains their product in zero-valued entries.
         * @throws std::invalid_argument If any dimension is nonpositive.
         * @throws std::length_error If the total exceeds container capacity.
         * @note Replaces all old face data on success and preserves it on allocation failure.
         */
        void allocate(
            const hd_core::Index& dimensions
        );

        // Face-field local coordinates; axis 0 varies fastest.
        /**
         * @brief Convert a flattened face offset to local face-field coordinates.
         * @param[in] flat Offset in faces.
         * @return Local coordinates with axis 0 varying fastest.
         * @throws std::out_of_range If flat is outside faces.
         * @throws std::logic_error If shape is invalid or inconsistent with storage.
         * @throws std::length_error If dimensions exceed container capacity.
         */
        hd_core::Index coordinates(std::size_t flat) const;

        /**
         * @brief Access a face using local face-field coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]); no mesh ghost offset is added.
         * @return A mutable reference to the stored face value.
         * @pre The field has been allocated with a valid shape.
         * @note Bounds and storage errors from flat_index() propagate.
         */
        hd_core::Conserved& at(
            std::initializer_list<int> coordinates
        );

        /**
         * @brief Access a face using local face-field coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]); no mesh ghost offset is added.
         * @return A read-only reference to the stored face value.
         * @pre The field has been allocated with a valid shape.
         * @note Bounds and storage errors from flat_index() propagate.
         */
        const hd_core::Conserved& at(
            std::initializer_list<int> coordinates
        ) const;

        /**
         * @brief Access a face using local face-field coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]); no mesh ghost offset is added.
         * @return A mutable reference to the stored face value.
         * @pre The field has been allocated with a valid shape.
         * @note Bounds and storage errors from flat_index() propagate.
         */
        hd_core::Conserved& at(const hd_core::Index& coordinates);
        /**
         * @brief Access a face using local face-field coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]); no mesh ghost offset is added.
         * @return A read-only reference to the stored face value.
         * @pre The field has been allocated with a valid shape.
         * @note Bounds and storage errors from flat_index() propagate.
         */
        const hd_core::Conserved& at(const hd_core::Index& coordinates) const;

        // Convert face-field local coordinates to a faces storage index.
        /**
         * @brief Flatten local face-field coordinates with axis 0 varying fastest.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]).
         * @return Offset into faces.
         * @throws std::out_of_range If any coordinate is outside the field.
         * @throws std::logic_error If shape is invalid or inconsistent with storage.
         * @throws std::length_error If dimensions exceed container capacity.
         * @throws std::invalid_argument If the initializer list has the wrong number of entries.
         */
        std::size_t flat_index(
            std::initializer_list<int> coordinates
        ) const;

        /**
         * @brief Flatten local face-field coordinates with axis 0 varying fastest.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]).
         * @return Offset into faces.
         * @throws std::out_of_range If any coordinate is outside the field.
         * @throws std::logic_error If shape is invalid or inconsistent with storage.
         * @throws std::length_error If dimensions exceed container capacity.
         */
        std::size_t flat_index(const hd_core::Index& coordinates) const;
    };
}
