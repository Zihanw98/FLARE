# ram pressure（2D）

二维低密度圆泡处在高密度长条背景中，圆内密度为外部的十分之一，背景风沿 `+x` 方向流动；
圆内初始静止，全域初始压强相同。圆团随后自由受流体作用而运动、变形。
外部密度归一化为 1，圆内为 0.1，即内外密度比为 1:10。

沿用本地 AMRVAC `tests/project/ram_pressure/mod_usr.t` 的参数约定：

```text
rho_circle = bubble_chi * wind_density
p0 = wind_density * (wind_speed / wind_mach)^2 / gamma
```

AMRVAC 源文件位于 `/home/u0190250/Research/amrvac/tests/project/ram_pressure/mod_usr.t`。
这里设置低密度圆泡 `bubble_chi=rho_inside/rho_outside=0.1`，与 AMRVAC 源文件中的默认值一致。
沿用 `wind_density=1`、`wind_speed=1`、`wind_mach=0.8`、`bubble_radius=1`、
圆心 `(0,0)` 的约定。采用 FLARE 自身的求解器与均匀网格，是同类初值的独立算例。

默认计算域为 `[-6,58] × [-4,4]`，网格 `768 × 96`，所有数值均为无量纲代码单位。
沿风向长度由 32 加倍到 64，延长下游区域；网格间距仍为 `dx=dy=1/12`。
理想气体 `gamma=5/3`，共同初始压强 `p0=0.9375`。

| 区域 | 密度 rho | vx | vy | 初始压强 p |
| --- | ---: | ---: | ---: | ---: |
| 圆内 `x²+y²<1` | 0.1 | 0 | 0 | 0.9375 |
| 圆外背景 | 1 | 1 | 0 | 0.9375 |

按单元中心判断圆内外。密度和速度为锐利间断，不加界面平滑或额外扰动。
内能密度均为 `p0/(gamma-1)=1.40625`；总能量密度圆内为 `1.40625`、背景为 `1.90625`。
相同组分下，圆内初始温度为背景的十倍。

左侧持续输入背景风；右侧和上下两侧采用零梯度出流。
四侧均用 `user` 边界回调实现，左侧角单元也保持风的状态。
使用固定实验室参考系，无重力、加热或冷却；不启用 frame tracking。
数值格式为 exact Riemann solver、PLM/MC 重构和 RK2 时间推进。
等压仅针对初始时刻，压力会随流场演化。

外部声速 `c_s,out=sqrt(gamma*p0/rho_out)=1.25`，因此背景风速
`v_w=0.8*c_s,out=1`；Mach 数明确以外部声速计算。圆内初始速度仍为零。
动压尺度 `rho_w*v_w²=1`，默认 `t_end=40`，即 `40 R/v_w`。

在仓库根目录运行：

```bash
make -C tests/ram_pressure -j4 test
make -C tests/ram_pressure animate
```

`test` 编译、运行并检查初始圆形密度分布、圆内静止、背景风速、外部 Mach 数 0.8、全域等压，
以及初始内能和总能量、初末边界、所有保存时刻的正密度/压力和完整结束时间。
检查器使用 `output/run.log` 中本次运行的文件清单，避免混入旧输出。
`make check` 只检查已有输出。
此前高密度圆团的运行结果保存在 `output_dense_chi10/`，当前低密度版本写入 `output/`。
此前长度 32、运行至 `t=12` 的低密度结果及其参数保存在 `output_light_L32_t12/`。
此前长度 64、运行至 `t=20` 的低密度结果及其参数保存在 `output_light_L64_t20/`。

- 初值图：`output/diagnostics/initial.png`
- 初末密度对比：`output/diagnostics/comparison.png`
- 检查摘要：`output/diagnostics/summary.json`
- 密度动画：`output/animations/ram_pressure.gif`

动画每帧时长相同，标题标出实际模拟时间。
Python 后处理依赖 numpy、matplotlib、Pillow。

网格、计算区域、结束时间、输出频率和数值格式在 `parameters.txt` 中调整；
密度比、风速、Mach 数、圆心和半径在 `user_problem.cpp` 顶部的命名常数中调整，
修改后重新编译。共同压强由 Mach 数自动计算。
`check.py` 验证上述默认物理状态；改变这些常数时也需更新其期望值。
保持 `nxa:nya=8:1` 可得到正方形网格；默认半径约 12 格，适合演示。
更高分辨率可取 `1536 × 192`，并使用新的输出目录保留不同运行。

这是二维圆形截面的流体算例。输出中没有云物质示踪量，
密度图用于观察形态，不能据此精确测定云物质质量或混合率。
