# 本轮修复验证记录

基准为本轮修复前的工作区：保留此前的 exact 和 implicit 改动，与本轮修复后的版本比较。
旧算例使用原网格、物理初值和终止时间；减少快照频率，输出放在临时目录。均匀加热另用相同初值增加隐式运行（dt=0.04，容差 1e-10）。

| 旧算例 | 方法 | 终止时间 | 步数（前后相同） | 最终状态最大差异 |
| --- | --- | --- | --- | --- |
| sod_tube | Euler | 0.2 | 436 | 0 |
| User_Heating_1D | RK2 | 2 | 208 | 0 |
| User_Heating_1D | Implicit | 2 | 50 | 2.72e-11 |
| Periodic_Acoustic_2D | Euler | 0.707106781187 | 227 | 0 |

差异按每个字段的 `max(1, max|基准值|)` 缩放；速度接近零时相当于报告绝对误差。比较 rho、各速度分量、P、内能密度和总能量密度。三个显式运行上述字段完全相同，隐式差异小于求解容差。

| 新算例 | 方法 | 步数 | 最大积分误差 |
| --- | --- | --- | --- |
| mixed_contact | rk2 | 493 | 2.22e-16 |
| mixed_contact | lts_rk2 | 124 | 2.22e-16 |
| ionized_contact | rk2 | 493 | 2.22e-16 |
| ionized_contact | lts_rk2 | 124 | 4.44e-16 |
| ionized_heating | euler | 3 | 1.11e-16 |
| ionized_heating | rk2 | 4 | 1.11e-16 |
| ionized_heating | lts_euler | 3 | 1.11e-16 |
| ionized_heating | lts_rk2 | 3 | 1.11e-16 |
| ionized_heating | implicit | 5 | 1.74e-12 |

积分误差包含质量、两种组分质量和总能量；加热时先扣除解析的加热贡献。另检查了 T、mu、电子数密度、电离能、组分和、接触压力与速度，以及显式加热的首步内能时间尺度。两个接触算例的组分 L1 误差约 0.0124，满足 PCM 的 0.025 阈值。

单元回归：1D 无组分完整默认测试集通过（exact_eos 的辐射部分按设计跳过）；2D 无组分的 state/cell_conversion/eos/implicit/integrator/rk2/source/lts_rk2/lts_euler 通过；1D 和 2D 的 HD_NCOMP=2 下 eos/cell_conversion/implicit/exact_eos 通过。

辐射 exact 回归覆盖均匀流、移动/静止接触、gamma=4/3 的解析参照、激波/稀疏波/真空、Newton 与二分压力解、冻结电离能输运，以及周期网格的一步守恒。
