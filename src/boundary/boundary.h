#pragma once
#include <string>

namespace hd_core { struct MeshBlock; }
namespace hd_user { struct UserHooks; }
namespace boundary
{
    enum class BoundarySide
    {
        Upper,
        Lower
    };
    enum class BoundaryCondition
    {
        Periodic,
        Outflow,
        Reflective,
        Cont,
        User
    };
    BoundaryCondition parse_bc(const std::string& bc);
    /**
     * @brief 填充指定方向、指定端的边界，优先调用用户 hook。
     * @param[in,out] mesh 复制完整 Cell（U 和 gamma）到指定端的所有 ghost 层。
     * @param[in] direction 坐标方向，范围为 [0, NDIM)。
     * @param[in] side 当前方向的 Lower 或 Upper 端。
     * @param[in] time 当前物理时间；周期边界不依赖时间。
     * @param[in] hooks 若 userbc 非空，则调用它并立即返回，覆盖该端的内置处理。
     * @note 横向覆盖所有存储单元。逐方向填完两端后，全周期网格角区也被更新。
     *       没有 userbc 时使用 mesh 中的边界类型，目前内置实现只有周期边界。
     * @throws std::invalid_argument 方向、端点、时间或周期配对不合法。
     * @throws std::runtime_error User 边界未注册 hook，或选择了尚未实现的内置边界。
     */
    void ApplyBoundaryCondition(
        hd_core::MeshBlock& mesh,
        int direction,
        BoundarySide side,
        double time,
        const hd_user::UserHooks& hooks
    );
}
