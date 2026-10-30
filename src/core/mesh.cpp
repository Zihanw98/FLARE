#include "mesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hd_core
{
namespace
{
    void validate_axis(const std::vector<double>& faces, int count)
    {
        if (faces.size() != static_cast<std::size_t>(count) + 1)
            throw std::invalid_argument("Face count must equal total cells + 1");

        for (std::size_t j = 0; j < faces.size(); ++j)
        {
            if (!std::isfinite(faces[j]))
                throw std::invalid_argument("Face coordinate must be finite");
            if (j > 0)
            {
                const double width = faces[j] - faces[j - 1];
                if (!std::isfinite(width) || width <= 0.0)
                    throw std::invalid_argument("Faces must have finite positive separation");
            }
        }
    }
    }

    int MeshBlock::ntotal(int d) const
    {
        if (d < 0 || d >= NDIM)
            throw std::out_of_range("Invalid mesh direction");

        if (nactive[d] <= 0 || nghost[d] < 0)
            throw std::invalid_argument("Invalid mesh size");
        if (nghost[d] > (std::numeric_limits<int>::max() - nactive[d]) / 2)
            throw std::overflow_error("Mesh size overflow");

        return nactive[d] + 2 * nghost[d];
    }
    Index MeshBlock::ntotal() const
    {
        Index result{};

        for (int d = 0; d < NDIM; ++d)
        {
            result[d] = ntotal(d);
        }

        return result;  
    }
    
    void MeshBlock::set_uniform_faces(int d, double lower, double upper)
    {
        const int count = ntotal(d); // Checks direction, sizes, and integer overflow.
        if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower)
            throw std::invalid_argument("Uniform grid bounds must be finite and increasing");

        const double length = upper - lower;
        const double width = length / nactive[d];
        if (!std::isfinite(length) || !std::isfinite(width) || width <= 0.0)
            throw std::invalid_argument("Invalid uniform grid spacing");

        std::vector<double> coordinates(static_cast<std::size_t>(count) + 1);
        const auto first = static_cast<std::size_t>(nghost[d]);
        const auto last = first + static_cast<std::size_t>(nactive[d]);
        for (std::size_t j = 0; j < coordinates.size(); ++j)
        {
            if (j == first)
                coordinates[j] = lower;
            else if (j == last)
                coordinates[j] = upper;
            else
                coordinates[j] = std::fma(static_cast<double>(j) - nghost[d], width, lower);
        }
        // Detect overflow in ghost coordinates or collapsed faces due to rounding.
        validate_axis(coordinates, count);
        face[d].swap(coordinates); // Leave previous geometry unchanged on failure.
    }

    double MeshBlock::cell_width(int direction, int axis_index) const
    {
        // 只检查所选方向的单元数（含 ghost），不是整个网格的单元总数。
        const int axis_cell_count = ntotal(direction);
        if (axis_index < 0 || axis_index >= axis_cell_count)
            throw std::out_of_range("Cell index along the selected axis is out of range");
        if (face[direction].size() != static_cast<std::size_t>(axis_cell_count) + 1)
            throw std::invalid_argument("Face count must equal total cells + 1");

        // face[direction] 是该方向的一维面坐标数组，允许非均匀间距。
        // 沿该方向的第 axis_index 个单元由相邻两个面夹住。
        const auto lower_face_index = static_cast<std::size_t>(axis_index);
        const double lower_face_position = face[direction][lower_face_index];
        const double upper_face_position = face[direction][lower_face_index + 1];
        const double width = upper_face_position - lower_face_position;
        if (!std::isfinite(lower_face_position) || !std::isfinite(upper_face_position) ||
            !std::isfinite(width) || width <= 0.0)
            throw std::invalid_argument("Cell faces must have finite positive separation");
        return width;
    }

    double MeshBlock::cell_center(int direction, int axis_index) const
    {
        // cell_width() 同时检查方向、方向内下标及两侧面坐标是否合法。
        const double width = cell_width(direction, axis_index);
        const auto lower_face_index = static_cast<std::size_t>(axis_index);
        const double lower_face_position = face[direction][lower_face_index];
        // 返回所选方向上的坐标中点。
        return lower_face_position + 0.5 * width;
    }

    double MeshBlock::cell_width(int direction, const Index& storage_index) const
    {
        if (direction < 0 || direction >= NDIM)
            throw std::out_of_range("Invalid mesh direction");
        // 张量积网格只需要该方向的坐标下标，无需转成 flat。
        return cell_width(direction, storage_index[direction]);
    }

    double MeshBlock::cell_center(int direction, const Index& storage_index) const
    {
        if (direction < 0 || direction >= NDIM)
            throw std::out_of_range("Invalid mesh direction");
        return cell_center(direction, storage_index[direction]);
    }

    double MeshBlock::active_length(int d) const
    {
        const int count = ntotal(d);
        validate_axis(face[d], count);
        const auto first = static_cast<std::size_t>(nghost[d]);
        const auto last = first + static_cast<std::size_t>(nactive[d]);
        const double length = face[d][last] - face[d][first];
        if (!std::isfinite(length) || length <= 0.0)
            throw std::invalid_argument("Active domain length must be finite and positive");
        return length;
    }

    void MeshBlock::validate_geometry() const
    {
        for (int d = 0; d < NDIM; ++d)
            (void)active_length(d); // Validates all faces, including ghost faces.
    }

    void MeshBlock::allocate()
    {
        validate_geometry();
        std::size_t count = 1;

        for (int d = 0; d < NDIM; ++d)
        {
            const auto n = static_cast<std::size_t>(ntotal(d));

            if (count > cells.max_size() / n)
                throw std::length_error("Too many mesh cells");

            count *= n;
        }

        cells.resize(count);
    }
    Index MeshBlock::coordinates(std::size_t flat) const
    {
        if (flat >= cells.size())
            throw std::out_of_range("Cell index out of range");

        Index index{};
        std::size_t total = 1;

        for (int d = 0; d < NDIM; ++d)
        {
            const auto n = static_cast<std::size_t>(ntotal(d));
            if (total > cells.max_size() / n)
                throw std::length_error("Mesh is too large");

            total *= n;
            index[d] = static_cast<int>(flat % n);
            flat /= n;
        }

        if (cells.size() != total)
            throw std::logic_error("Cell storage does not match mesh dimensions");

        return index;
    }
    Index MeshBlock::coordinates_active(std::size_t flat) const
    {
        auto coordinates = MeshBlock::coordinates(flat);
        for(int i=0;i<NDIM;i++)
        {
            coordinates[i] -= nghost[i];
        }
        return coordinates;
    }

    std::size_t MeshBlock::flat_index(
        std::initializer_list<int> coordinates,
        bool active
    ) const
    {
        // 必须恰好提供 NDIM 个坐标。
        if (coordinates.size() != static_cast<std::size_t>(NDIM))
        {
            throw std::invalid_argument(
                "Number of coordinates must equal NDIM"
            );
        }

        Index index{};
        std::copy(coordinates.begin(), coordinates.end(), index.begin());
        return flat_index(index, active);
    }

    std::size_t MeshBlock::flat_index(
        const Index& coordinates,
        bool active
    ) const
    {
        std::size_t offset = 0;
        std::size_t stride = 1;
        int d = 0;

        for (int coordinate : coordinates)
        {
            const int n = ntotal(d);  // 包括 ghost cells
            const int limit = active ? nactive[d] : n;

            if (coordinate < 0 || coordinate >= limit)
            {
                throw std::out_of_range(
                    active ? "Active cell coordinate out of range"
                           : "Cell coordinate out of range"
                );
            }

            const auto count = static_cast<std::size_t>(n);

            // 防止总格数乘法溢出或超过容器容量限制。
            if (stride > cells.max_size() / count)
            {
                throw std::length_error("Mesh is too large");
            }

            // ntotal() has checked that the sizes and their sum fit in int.
            const int storage_coordinate =
                active ? coordinate + nghost[d] : coordinate;
            offset += static_cast<std::size_t>(storage_coordinate) * stride;
            stride *= count;
            ++d;
        }

        // 检查是否已按当前网格尺寸分配 cells。
        if (cells.size() != stride)
        {
            throw std::logic_error(
                "Cell storage does not match mesh dimensions"
            );
        }

        return offset;
    }

    Cell& MeshBlock::at(std::initializer_list<int> coordinates)
    {
        return cells.at(flat_index(coordinates));
    }

    const Cell& MeshBlock::at(
        std::initializer_list<int> coordinates
    ) const
    {
        return cells.at(flat_index(coordinates));
    }

    Cell& MeshBlock::at(const Index& coordinates)
    {
        return cells.at(flat_index(coordinates));
    }

    const Cell& MeshBlock::at(const Index& coordinates) const
    {
        return cells.at(flat_index(coordinates));
    }

    Cell& MeshBlock::at_active(std::initializer_list<int> coordinates)
    {
        return cells.at(flat_index(coordinates, true));
    }

    const Cell& MeshBlock::at_active(
        std::initializer_list<int> coordinates
    ) const
    {
        return cells.at(flat_index(coordinates, true));
    }

    Cell& MeshBlock::at_active(const Index& coordinates)
    {
        return cells.at(flat_index(coordinates, true));
    }

    const Cell& MeshBlock::at_active(const Index& coordinates) const
    {
        return cells.at(flat_index(coordinates, true));
    }
}
