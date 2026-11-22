# 验证记录

2026-11-22，使用本目录的初态和 `parameters.txt`，全程未启用组分。
默认 768 单元，另把 `nxa` 改为 1536、输出间隔改为 20 步运行加密检查；其他物理与数值设置相同。

| 诊断 | 768 单元（默认动画） | 1536 单元 |
| --- | ---: | ---: |
| 初始质量 / \(M_\odot\) | 1.4852049535 | 1.4852049532 |
| 中心首次越过核密度 / ms | 9.2313 | 8.5848 |
| 中心密度最大值对应时刻 / ms | 9.5930 | 9.0476 |
| 最大中心密度 / \(\mathrm{g\,cm^{-3}}\) | \(2.3954\times10^{14}\) | \(3.1752\times10^{14}\) |
| 诊断首次确认外传激波 / ms | 9.7811 | 9.0783 |
| 10.4 ms 的激波位置 / km | 65.01 | 107.18 |
| 最后快照的压强跃变 | 77.37 | 85.91 |
| 最后快照的密度跃变 | 5.480 | 4.950 |
| 最后快照的上游激波参考系 Mach 数 | 7.322 | 8.386 |
| 最大质量相对误差 | \(8.88\times10^{-16}\) | \(1.33\times10^{-15}\) |

这两组都通过 `analyze.py` 的核密度穿越、中心反弹、热能非负、EoS 一致性、质量守恒、向外传播与压缩跃变检查。
激波诊断在压缩前沿两侧取数个单元，因此跃变数值同时含有背景径向梯度和数值展宽。

PCM/HLL 的数值耗散影响坍缩和热压增长：这里峰值密度相差约 33%，反弹时刻相差约 5.7%。因此 **目前是定性反弹演示，不是已收敛的 CCSN 物理结果**。默认动画保留完整模拟快照中的状态，反弹附近降低播放速度。

## 单元与回归测试

- 新 `nuclear` 测试：RhoT/RhoP/RhoE 往返；RhoE 不改 U；核密度处冷压力和冷能量连续；(de_c/d\rho=P_c/\rho^2)；声速与独立等熵有限差分相符；CGS/非 CGS 单位一致；非法负热能拒绝；HLL 均匀态、反射边界态；Nuclear/exact 错误配对拒绝。
- `NDIM=1, HD_NCOMP=0`：`nuclear eos cell_conversion godunov integrator rk2 source wave_timestep implicit`，9 组通过。
- `NDIM=1, HD_NCOMP=2`：`nuclear eos cell_conversion exact_eos`，4 组通过，保留之前 Mixed/Ionized/辐射分支回归覆盖。
- 后续只增加了 factory 中 Nuclear/exact 的报错提示，再次运行 `nuclear eos cell_conversion` 通过。
- Python 脚本语法检查及 `git diff --check` 通过。

既有 `exact`、引力、双能量同步、重构和积分器实现没有因本次 Nuclear 模型继续改动。
配套 HLL 是独立可选求解器，只有配置 `riemann_solver=hll` 才会使用。

## 产物

默认运行数据、`validation.json`、`history.csv` 和动画在本目录的 `output_bounce/` 中。加密检查的数据在本次机器的 `/tmp/flare-nuclear-fine-output/`，此临时路径不是长期交付文件；关键结果已记录在上表。
