/**
 * @file facefield.h
 * @brief Generic dense face-centered storage and indexing.
 */
#pragma once
#include "core.h"
#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace hd_core
{
    /**
     * @brief Dense face-centered storage with axis 0 varying fastest.
     * @tparam T Value-initializable face value (e.g. double or Conserved).
     * @note The caller supplies the complete face shape, including the extra
     *       normal face. No direction or ghost offset is inferred here.
     *       Very Important: For a field normal to d, index I names the lower face of cell I;
     *       increment I[d] to select that cell's upper face.
     */
    template<typename T>
    struct FaceField
    {

        hd_core::Index shape{};

        std::vector<T> faces;

        /**
         * @brief Allocate a value-initialized face field with the requested dimensions.
         * @param[in] dimensions Positive face counts on every axis.
         * @post shape equals dimensions and faces contains their product in value-initialized entries.
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
        T& at(
            std::initializer_list<int> coordinates
        );

        /**
         * @brief Access a face using local face-field coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]); no mesh ghost offset is added.
         * @return A read-only reference to the stored face value.
         * @pre The field has been allocated with a valid shape.
         * @note Bounds and storage errors from flat_index() propagate.
         */
        const T& at(
            std::initializer_list<int> coordinates
        ) const;

        /**
         * @brief Access a face using local face-field coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]); no mesh ghost offset is added.
         * @return A mutable reference to the stored face value.
         * @pre The field has been allocated with a valid shape.
         * @note Bounds and storage errors from flat_index() propagate.
         */
        T& at(const hd_core::Index& coordinates);
        /**
         * @brief Access a face using local face-field coordinates.
         * @param[in] coordinates Exactly NDIM coordinates in [0, shape[d]); no mesh ghost offset is added.
         * @return A read-only reference to the stored face value.
         * @pre The field has been allocated with a valid shape.
         * @note Bounds and storage errors from flat_index() propagate.
         */
        const T& at(const hd_core::Index& coordinates) const;

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

// Template definitions must be visible to every caller.
namespace hd_core
{

    template<typename T>
    void FaceField<T>::allocate(
        const hd_core::Index& dimensions
    )
    {
        std::size_t total = 1;

        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            if (dimensions[d] <= 0)
            {
                throw std::invalid_argument(
                    "FaceField dimensions must be positive"
                );
            }

            const auto n =
                static_cast<std::size_t>(dimensions[d]);

            if (total > faces.max_size() / n)
            {
                throw std::length_error("FaceField is too large");
            }

            total *= n;
        }

        // Value-initialize every face; arithmetic values and Conserved start at zero.
        // 如果分配失败，原来的 shape 和 faces 保持不变。
        std::vector<T> new_faces(total);

        faces.swap(new_faces);
        shape = dimensions;
    }

    template<typename T>
    hd_core::Index FaceField<T>::coordinates(std::size_t flat) const
    {
        if (flat >= faces.size())
            throw std::out_of_range("Face index out of range");

        hd_core::Index index{};
        std::size_t total = 1;

        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            if (shape[d] <= 0)
                throw std::logic_error("FaceField has invalid or uninitialized shape");

            const auto n = static_cast<std::size_t>(shape[d]);
            if (total > faces.max_size() / n)
                throw std::length_error("FaceField is too large");

            total *= n;
            index[d] = static_cast<int>(flat % n);
            flat /= n;
        }

        if (faces.size() != total)
            throw std::logic_error("FaceField storage does not match its shape");

        return index;
    }

    template<typename T>
    std::size_t FaceField<T>::flat_index(
        std::initializer_list<int> coordinates
    ) const
    {
        if (coordinates.size()
            != static_cast<std::size_t>(hd_core::NDIM))
        {
            throw std::invalid_argument(
                "Number of coordinates must equal NDIM"
            );
        }

        hd_core::Index index{};
        std::copy(coordinates.begin(), coordinates.end(), index.begin());
        return flat_index(index);
    }

    template<typename T>
    std::size_t FaceField<T>::flat_index(const hd_core::Index& coordinates) const
    {
        std::size_t offset = 0;
        std::size_t stride = 1;
        int d = 0;

        for (int coordinate : coordinates)
        {
            if (shape[d] <= 0)
            {
                throw std::logic_error(
                    "FaceField has invalid or uninitialized shape"
                );
            }

            if (coordinate < 0 || coordinate >= shape[d])
            {
                throw std::out_of_range(
                    "FaceField coordinate out of range"
                );
            }

            const auto n = static_cast<std::size_t>(shape[d]);

            if (stride > faces.max_size() / n)
            {
                throw std::length_error("FaceField is too large");
            }

            offset += static_cast<std::size_t>(coordinate) * stride;
            stride *= n;

            ++d;
        }

        if (faces.size() != stride)
        {
            throw std::logic_error(
                "FaceField storage does not match its shape"
            );
        }

        return offset;
    }

    template<typename T>
    T& FaceField<T>::at(
        std::initializer_list<int> coordinates
    )
    {
        return faces.at(flat_index(coordinates));
    }

    template<typename T>
    const T& FaceField<T>::at(
        std::initializer_list<int> coordinates
    ) const
    {
        return faces.at(flat_index(coordinates));
    }

    template<typename T>
    T& FaceField<T>::at(const hd_core::Index& coordinates)
    {
        return faces.at(flat_index(coordinates));
    }

    template<typename T>
    const T& FaceField<T>::at(const hd_core::Index& coordinates) const
    {
        return faces.at(flat_index(coordinates));
    }

} // namespace hd_core
