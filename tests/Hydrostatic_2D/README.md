# 静水平衡保持

二维稳定等温分层，无初始速度或扰动。以底部为 y=0、`p/rho=1`、
底部 `rho=p=1`，解析静态解为：

```text
rho(y) = exp(g_y*y)
p(y)   = rho(y)
dp/dy  = rho*g_y
```

初始化使用该指数分布的精确单元平均值。等温仅描述初始平衡态，
演化仍使用绝热理想气体 Euler 方程，没有冷却或强制恒温。
选择稳定分层是为了避免把对流不稳定性混入静水平衡误差。
x 周期、y 上下使用与 RTI 相同的内置反射壁面；32×64，t=2。
本测试检查通量、重力源项和边界的组合，不单独隔离任何一个模块。

在项目根目录运行：

```bash
make -C tests/Hydrostatic_2D -j4 test
make -C tests/Hydrostatic_2D check
make -C tests/Hydrostatic_2D -j4 diagnose
```

`test` 编译、运行并严格检查所有保存快照中的最大马赫数、密度相对变化、
压力相对变化均不超过 `1e-10`。**当前普通 Godunov 通量、中心重力源项和
镜像压力边界不保持静水平衡，预期此检查失败。** 这不代表初始化公式错误，
也不是一般一阶格式精度测试；它专门暴露缺失的 well-balanced 性质。
`check` 对已有日志重复严格检查。`diagnose` 运行相同算例，仍检查初态、
正性、水平对称性、质量守恒和总能量漂移，但对平衡漂移仅报告数值，
不将其当作通过静水平衡测试。

Python 依赖 numpy、matplotlib。输出在本算例的 `output/`，日志为 `run.log`，
诊断在 `output/diagnostics/metrics.csv` 和 `diagnostics.png`。
普通 `make run` 不生成检查所需的 run.log，请用 `make test` 或 `make diagnose`。
指标包含竖直压力梯度残差；该诊断使用剖面的有限差分，不等同于求解器 RHS。
总能量统计包含重力势能，0.1% 漂移阈值只作粗网格回归界限。
后续修改边界或平衡格式时，可用本测试比较误差随网格和时间的变化。
