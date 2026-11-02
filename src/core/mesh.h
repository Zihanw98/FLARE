/**
 * @file mesh.h
 * @brief 定义带 ghost 单元的结构化网格块及其坐标、索引访问接口。
 * @details
 * 网格维度由编译期常量 hd_core::NDIM 决定。每个方向分别保存一维面坐标，
 * 可以表示均匀或非均匀的张量积网格。
 *
 * 本文件使用三种索引：
 * - 存储坐标：包含 ghost 偏移，每个分量位于 [0, ntotal(d))。
 * - 活跃区局部坐标：活跃区从 0 开始，不包含 ghost 偏移。
 * - flat 下标：cells 中的一维存储位置，方向 0 变化最快，包含 ghost 单元。
 *
 * @see hd_core::MeshBlock
 */
#pragma once
#include "core.h"
#include "cell.h"
#include "boundary.h"
#include <array>
#include <cstddef>
#include <initializer_list>
#include <vector>
/// @brief 核心物理状态量与网格数据类型。
namespace hd_core
{
    /**
     * @brief 保存网格几何、边界类型和全部单元状态。
     * @details 先设置 nactive、nghost 和各方向 face，再调用 allocate()。
     *          面坐标独立于单元状态；填充 ghost 状态不会自动改变面坐标。
     * @note 修改已分配网格的尺寸或几何后，需要重新初始化或守恒重映射物理数据。
     */
    struct MeshBlock
    {
        Index nactive{}; ///< 各方向的活跃单元数，必须大于 0。
        Index nghost{};  ///< 各方向每一侧的 ghost 层数，必须非负。
        
        /**
         * @brief 返回某方向的总单元数，包括两侧 ghost。
         * @param[in] d 方向编号，满足 0 <= d < NDIM。
         * @return nactive[d] + 2 * nghost[d]。
         * @throws std::out_of_range 方向编号越界。
         * @throws std::invalid_argument 活跃单元数或 ghost 层数非法。
         * @throws std::overflow_error 总单元数超出 int 范围。
         */
        int ntotal(int d) const;

        /**
         * @brief 返回各方向的总单元数，包括两侧 ghost。
         * @return 各分量分别为 ntotal(d) 的 Index。
         * @see ntotal(int) const
         */
        Index ntotal() const;
        
        /**
         * @brief 各方向的面坐标，包含 ghost 区域的面。
         * @details face[d].size() 必须等于 ntotal(d) + 1。
         *          活跃区域的面下标为 nghost[d] 到 nghost[d] + nactive[d]，含两端。
         *          所有面坐标必须有限并严格递增。
         */
        std::array<std::vector<double>,NDIM> face;
        std::array<boundary::BoundaryCondition,NDIM> bc_upper; ///< 各方向上边界的边界条件类型。
        std::array<boundary::BoundaryCondition,NDIM> bc_lower; ///< 各方向下边界的边界条件类型。

        /**
         * @brief 为一个方向生成均匀面坐标，并向两侧延伸 ghost 面。
         * @param[in] d 方向编号，满足 0 <= d < NDIM。
         * @param[in] lower 活跃区域的下边界坐标。
         * @param[in] upper 活跃区域的上边界坐标，必须大于 lower。
         * @pre 已设置该方向的 nactive 和 nghost。
         * @note 单元宽度为 (upper - lower) / nactive[d]，ghost 使用相同宽度。
         *       本函数不分配 cells，不填充物理状态；失败时保留原有面坐标。
         * @throws std::invalid_argument 边界非有限或非递增，或生成的面间距不是有限正值。
         * @see ntotal(int) const
         */
        void set_uniform_faces(int d, double lower, double upper);

        /**
         * @brief 返回单元沿指定方向的宽度。
         * @param[in] direction 方向编号，满足 0 <= direction < NDIM。
         * @param[in] axis_index 该方向的存储下标，包含 ghost 偏移。
         * @return face[direction][axis_index + 1] - face[direction][axis_index]。
         * @note axis_index 是单个方向的下标，不是 cells 的 flat 下标。
         * @throws std::out_of_range 方向或该方向的单元下标越界。
         * @throws std::invalid_argument 网格尺寸、面数量或面间距非法。
         */
        double cell_width(int direction, int axis_index) const;

        /**
         * @brief 返回单元沿指定方向的中心坐标。
         * @param[in] direction 方向编号，满足 0 <= direction < NDIM。
         * @param[in] axis_index 该方向的存储下标，包含 ghost 偏移。
         * @return 两侧面坐标的中点。
         * @note ghost 使用相同公式，其中心可以位于活跃区域以外。
         *       这里计算坐标中点，不是一般曲线网格的体积质心。
         * @see cell_width(int, int) const
         */
        double cell_center(int direction, int axis_index) const;

        /**
         * @brief 使用多维存储坐标取得单元宽度。
         * @param[in] direction 要查询的方向。
         * @param[in] storage_index 包含 ghost 偏移的多维坐标，仅使用所选方向的分量。
         * @return cell_width(direction, storage_index[direction])。
         * @see cell_width(int, int) const
         */
        double cell_width(int direction, const Index& storage_index) const;

        /**
         * @brief 使用多维存储坐标取得单元中心。
         * @param[in] direction 要查询的方向。
         * @param[in] storage_index 包含 ghost 偏移的多维坐标，仅使用所选方向的分量。
         * @return cell_center(direction, storage_index[direction])。
         * @see cell_center(int, int) const
         */
        double cell_center(int direction, const Index& storage_index) const;

        /**
         * @brief 返回指定方向活跃区域的长度。
         * @param[in] d 方向编号，满足 0 <= d < NDIM。
         * @return 活跃区域两端面的坐标差，不包含 ghost 宽度。
         * @throws std::invalid_argument 该方向的面数量、面坐标或区域长度非法。
         * @see ntotal(int) const
         */
        double active_length(int d) const;

        /**
         * @brief 检查所有方向的尺寸及面坐标，包括 ghost 区域。
         * @throws std::invalid_argument 网格尺寸、面数量或面坐标不合法。
         * @see active_length(int) const
         */
        void validate_geometry() const;

        /**
         * @brief 验证几何并将 cells 调整到包含 ghost 的总单元数。
         * @pre 已设置各方向的 nactive、nghost 和 face。
         * @note 本函数不填充用户物理初态；resize 保留已有元素，新元素使用 Cell 默认值。
         * @throws std::length_error 总单元数超过 cells 的容量限制。
         * @see validate_geometry()
         */
        void allocate();

        std::vector<Cell> cells; ///< 全部存储单元，包含 ghost；方向 0 的下标变化最快。

        /**
         * @brief 按存储坐标取得可修改的单元引用。
         * @param[in] coordinates 恰好 NDIM 个坐标分量，包含 ghost 偏移。
         * @return 对应单元的引用。
         * @see flat_index(std::initializer_list<int>, bool) const
         */
        Cell& at(std::initializer_list<int> coordinates);

        /**
         * @brief 按存储坐标取得只读的单元引用。
         * @param[in] coordinates 恰好 NDIM 个坐标分量，包含 ghost 偏移。
         * @return 对应单元的只读引用。
         * @see flat_index(std::initializer_list<int>, bool) const
         */
        const Cell& at(std::initializer_list<int> coordinates) const;

        /**
         * @brief 按活跃区局部坐标取得可修改的单元引用。
         * @param[in] coordinates 恰好 NDIM 个分量，每个分量满足 0 <= coordinates[d] < nactive[d]。
         * @return 对应活跃单元的引用，内部自动添加 ghost 偏移。
         * @note 此接口不能选中 ghost 单元。
         * @see flat_index(std::initializer_list<int>, bool) const
         */
        Cell& at_active(std::initializer_list<int> coordinates);

        /**
         * @brief 按活跃区局部坐标取得只读的单元引用。
         * @param[in] coordinates 恰好 NDIM 个合法的活跃区坐标分量。
         * @return 对应活跃单元的只读引用，内部自动添加 ghost 偏移。
         * @see flat_index(std::initializer_list<int>, bool) const
         */
        const Cell& at_active(std::initializer_list<int> coordinates) const;

        /**
         * @brief 将 cells 的 flat 下标转换为包含 ghost 偏移的多维存储坐标。
         * @param[in] flat cells 的存储下标，满足 flat < cells.size()。
         * @return 对应的 Index，方向 0 变化最快。
         * @throws std::out_of_range flat 超出存储范围。
         * @throws std::logic_error cells 大小与当前网格尺寸不一致。
         */
        Index coordinates(std::size_t flat) const;

        /**
         * @brief 将 cells 的 flat 下标转换为以活跃区为原点的多维坐标。
         * @param[in] flat 包含 ghost 单元的存储下标，不是活跃区单独展平后的下标。
         * @return coordinates(flat) 的每个分量减去对应的 nghost。
         * @note 输入可以选中 ghost，此时返回值至少一个分量不在 [0, nactive[d]) 内。
         * @see coordinates(std::size_t) const
         */
        Index coordinates_active(std::size_t flat) const;

        /**
         * @brief 按 Index 存储坐标取得可修改的单元引用。
         * @param[in] coordinates 包含 ghost 偏移的多维存储坐标。
         * @return 对应单元的引用。
         * @see flat_index(const Index&, bool) const
         */
        Cell& at(const Index& coordinates);

        /**
         * @brief 按 Index 存储坐标取得只读的单元引用。
         * @param[in] coordinates 包含 ghost 偏移的多维存储坐标。
         * @return 对应单元的只读引用。
         * @see flat_index(const Index&, bool) const
         */
        const Cell& at(const Index& coordinates) const;

        /**
         * @brief 按 Index 活跃区坐标取得可修改的单元引用。
         * @param[in] coordinates 活跃区局部坐标，内部自动添加 ghost 偏移。
         * @return 对应活跃单元的引用。
         * @see flat_index(const Index&, bool) const
         */
        Cell& at_active(const Index& coordinates);

        /**
         * @brief 按 Index 活跃区坐标取得只读的单元引用。
         * @param[in] coordinates 活跃区局部坐标，内部自动添加 ghost 偏移。
         * @return 对应活跃单元的只读引用。
         * @see flat_index(const Index&, bool) const
         */
        const Cell& at_active(const Index& coordinates) const;

        /**
         * @brief 将坐标列表转换为 cells 的存储下标。
         * @param[in] coordinates 恰好 NDIM 个坐标分量。
         * @param[in] active 为 true 时输入是活跃区局部坐标；否则输入是存储坐标。
         * @return 包含 ghost 的 cells 存储下标，active=true 时内部自动添加 ghost 偏移。
         * @throws std::invalid_argument 坐标分量个数不等于 NDIM。
         * @see flat_index(const Index&, bool) const
         */
        std::size_t flat_index(
            std::initializer_list<int> coordinates,
            bool active = false
        ) const;

        /**
         * @brief 将 Index 转换为 cells 的存储下标，方向 0 变化最快。
         * @param[in] coordinates 存储坐标或活跃区局部坐标，由 active 决定。
         * @param[in] active 为 true 时检查活跃区范围并添加 ghost 偏移，默认 false。
         * @return 始终返回完整 cells 数组中的下标。
         * @pre cells 已按当前网格尺寸分配。
         * @throws std::out_of_range 输入坐标超出所选区域的范围。
         * @throws std::logic_error cells 大小与当前网格尺寸不一致。
         * @throws std::length_error 网格大小超过 cells 的容量限制。
         * @see ntotal(int) const
         */
        std::size_t flat_index(
            const Index& coordinates,
            bool active = false
        ) const;
 
    };
}
