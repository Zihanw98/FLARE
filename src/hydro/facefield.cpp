#include "facefield.h"

#include <algorithm>
#include <stdexcept>

namespace hydro
{

    void FaceField::allocate(
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

        // 创建新的数组，所有通量分量初始化为零。
        // 如果分配失败，原来的 shape 和 faces 保持不变。
        std::vector<hd_core::Conserved> new_faces(total);

        faces.swap(new_faces);
        shape = dimensions;
    }

    hd_core::Index FaceField::coordinates(std::size_t flat) const
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

    std::size_t FaceField::flat_index(
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

    std::size_t FaceField::flat_index(const hd_core::Index& coordinates) const
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

    hd_core::Conserved& FaceField::at(
        std::initializer_list<int> coordinates
    )
    {
        return faces.at(flat_index(coordinates));
    }

    const hd_core::Conserved& FaceField::at(
        std::initializer_list<int> coordinates
    ) const
    {
        return faces.at(flat_index(coordinates));
    }

    hd_core::Conserved& FaceField::at(const hd_core::Index& coordinates)
    {
        return faces.at(flat_index(coordinates));
    }

    const hd_core::Conserved& FaceField::at(const hd_core::Index& coordinates) const
    {
        return faces.at(flat_index(coordinates));
    }

} // namespace hydro
