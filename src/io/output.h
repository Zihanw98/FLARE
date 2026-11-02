#pragma once
#include <cstddef>

namespace hd_core { struct MeshBlock; }
namespace config { struct SimulationConfig; }
namespace hd_user { struct UserHooks; }

namespace io
{
    /**
     * @brief 打印当前步信息，并按步数间隔保存网格 CSV。
     * @param[in,out] mesh 输出所有存储单元（包括 ghost）；内置输出只读取网格。
     *                    现有 useroutput 接口接收可修改引用，用户应只做输出。
     * @param[in] cfg 输出目录、文件后缀、step_every、t_end 和 maxsteps。
     * @param[in] hooks CSV 成功保存后调用可选的 useroutput(mesh, cfg)。
     * @param[in] step 已完成的步数；0 表示初始状态。
     * @param[in] time 当前物理时间。
     * @param[in] dt 刚完成一步使用的步长；初始状态传 0。
     * @return 无返回值。
     * @note 每次调用打印活跃单元的密度/压力范围、总质量和总能量。
     *       step=0、step=1、达到 t_end/maxsteps 或 output_step_every 的倍数时写 CSV。
     *       x/y/z 为单元中心；缺失方向的坐标和速度填 0，is_ghost 标记 ghost。
     *       文件前 3 行为 # step=...、# time=...、# dt=...，记录该快照的
     *       步数、物理时间和本步步长；初始快照三者均为 0，第 4 行为 CSV 列名。
     *       按 mesh.cells 的存储顺序输出，不重填 ghost，不改变其存储状态。
     * @throws std::invalid_argument 输出间隔、时间参数或网格状态不合法。
     * @throws std::runtime_error 文件无法打开或写入失败。
     * @throws std::filesystem::filesystem_error 输出目录无法创建。
     */
    void Output(
        hd_core::MeshBlock& mesh,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        std::size_t step,
        double time,
        double dt
    );
}
