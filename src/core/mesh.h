#pragma once
#include "core.h"
#include "cell.h"
#include <array>
#include <cstddef>
#include <initializer_list>
#include <vector>
namespace hd_core
{
    struct MeshBlock
    {
        Index nactive{};
        Index nghost{};
        
        int ntotal(int d) const;
        Index ntotal() const;
        
        // Tensor-product coordinates, including ghost faces.
        // face[d].size() == total(d) + 1.
        // Active faces run from nghost[d] through nghost[d] + nactive[d].
        std::array<std::vector<double>,NDIM> face;

        // Set nactive/nghost first. Bounds describe the ACTIVE region;
        // ghost faces are extended with the same uniform spacing.
        void set_uniform_faces(int d, double lower, double upper);

        // axis_index 是沿 direction 的单元下标，包含 ghost 偏移，
        // 不是 cells 数组的 flat 下标；单元两侧面下标为 axis_index 和 axis_index + 1。
        double cell_width(int direction, int axis_index) const;
        // 返回该方向两侧面坐标的中点，不是一般曲线网格的体积质心。
        double cell_center(int direction, int axis_index) const;
        // storage_index 是包含 ghost 的多维存储坐标，只使用所选方向的分量。
        double cell_width(int direction, const Index& storage_index) const;
        double cell_center(int direction, const Index& storage_index) const;
        double active_length(int d) const;

        void validate_geometry() const;
        // Validates geometry, then allocates cells. Does not set physical states.
        // After allocation, do not change sizes/coordinates without reinitializing
        // or conservatively remapping the physical data.
        void allocate();

        std::vector<Cell> cells;

        Cell& at(std::initializer_list<int> coordinates);
        const Cell& at(std::initializer_list<int> coordinates) const;

        // Active-local indices: 0 <= coordinates[d] < nactive[d].
        // Ghost offsets are added internally; ghost cells cannot be selected.
        Cell& at_active(std::initializer_list<int> coordinates);
        const Cell& at_active(std::initializer_list<int> coordinates) const;


        // Storage coordinates including ghost cells; axis 0 varies fastest.
        Index coordinates(std::size_t flat) const;
        // Input is a cells storage index, not an active-only flat index.
        // Ghost cells map outside [0, nactive[d]).
        Index coordinates_active(std::size_t flat) const;
        Cell& at(const Index& coordinates);
        const Cell& at(const Index& coordinates) const;

        Cell& at_active(const Index& coordinates);
        const Cell& at_active(const Index& coordinates) const;
        // Always returns a cells storage index. With active=true, input
        // coordinates are active-local and ghost offsets are added internally.
        std::size_t flat_index(
            std::initializer_list<int> coordinates,
            bool active = false
        ) const;

        std::size_t flat_index(
            const Index& coordinates,
            bool active = false
        ) const;
 
    };
}
