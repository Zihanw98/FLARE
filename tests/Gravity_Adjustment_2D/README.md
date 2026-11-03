# 均匀气体的重力调整

二维理想气体初始为 `rho=p=1`、速度为零，重力向下 `g_y=-0.1`。
x 周期、y 上下反射壁面；默认网格 32×64，运行到 t=4。
气体先向下加速，然后在壁面作用下产生压缩、膨胀和反射声波。
不假设最终静止，也不把 `p=rho*g*h` 当作可压缩气体的精确终态。

在项目根目录运行：

```bash
make -C tests/Gravity_Adjustment_2D -j4 test
make -C tests/Gravity_Adjustment_2D check
```

`test` 编译、运行、检查；`check` 只分析已有的 `run.log` 及其中列出的快照。
Python 检查脚本依赖 numpy、matplotlib。
输出保存在本算例的 `output/`；日志为 `run.log`，诊断为
`output/diagnostics/metrics.csv` 和 `diagnostics.png`。
普通 `make run` 不生成检查所需的 run.log，请用 `make test`。

检查包括：初态、首步内部气体的 `v_y=g_y*dt`、压力响应、正密度/正压力、
水平对称性、质量守恒，以及气体能量加重力势能的相对漂移小于 0.1%。
这个能量阈值是粗网格回归界限，不是精确守恒的宣称。
诊断图给出初末竖直剖面、速度幅值、上下压力差和总能量漂移。
检查只覆盖保存的快照；时间或网格修改后应重新运行，不能混用旧日志。
