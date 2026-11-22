# H/He 扩散与守恒演示

二维 Cartesian 周期盒，48×48 个有效 cell，边长 `1e8 cm`，演化至 `20 s`。
初始中央半径 `0.2 L` 的圆形区域富 H，外部富 He，过渡宽度 `0.025 L`；
两侧各保留 `1e-4` 的微量另一组分。初始 `rho=1e-10 g/cm³`、
`P=1 dyn/cm²`，速度为零。使用 `IdealGas_Mixed`、`gamma=5/3`，
两个组分是中性 H I 和 He I。

初始压力和密度均匀，因此温度会随平均分子量变化。
没有剪切扰动；关闭扩散时，这个静止接触结构应保持不动。
默认使用 `LTS_RK2`，`lts_max_levels=3` 是允许的层数上限，
并不保证这个均匀网格、常数扩散系数算例实际出现不同时间层。

## 扩散怎样接入

本次仅新增这个算例目录，复用已经接好的 `SourceKind::Transport`。
`user_problem.cpp` 中的 `diffusion_flux()` 在同一阶段的相邻 cell 间计算：

\[
J_H=-\rho_f D\frac{X_{H,R}-X_{H,L}}{x_R-x_L},\qquad J_{He}=-J_H,
\]

其中 `D=1e13 cm²/s`，`rho_f` 取两侧算术平均。
`J_H` 和 `J_He` 是相对于流体平均速度的质量通量，单位 `g/(cm² s)`。
这两个通量之和为零，不给总质量、动量额外添加通量。

扩散同时输运组分的焓：

\[
F_E=F_u=(h_H-h_{He})J_H,\qquad
h_k=\frac{\gamma}{\gamma-1}\frac{R_gT_{\mathrm{donor},k}}{A_k}.
\]

每种组分携带自己出发侧的温度，`A_H=1`、`A_He=4`。
两种组分方向相反，所以温度也分别取相反侧；这是焓通量的迎风离散，
避免用热邻居的平均温度计算冷 cell 的流出能量。
总能量和辅助内能的通量都显式填写；积分器负责面积、体积、时间积分和两侧相反的更新。
所以每种组分可以在 cell 之间移动，但周期盒内各组分总质量、总质量和总能量应守恒。
焓输运会改变局部压力并激发流动，因此这不是固定密度热方程的解析解测试。
辅助内能还会与动能交换，不能要求它的全域积分单独守恒。

这是常数 Fick 系数的数值演示，没有实现真实 H/He 碰撞系数表、热扩散、
压强扩散、Saha 或辐射。步长估计在组分扩散率基础上乘 `8`，
为本算例的焓耦合留出余量；不是任意多组分模型的通用稳定性公式。

## 运行和绘图

以下命令在 FLARE 仓库根目录执行：

```bash
make -j4 PROBLEM_DIR=tests/Diffusion_HHe_2D
(cd tests/Diffusion_HHe_2D && ./flare -par parameters.txt)
python3 tests/Diffusion_HHe_2D/check.py
python3 animate.py tests/Diffusion_HHe_2D/output diffusion_hhe \
  --composition --thermodynamics --cgs --fps 10 --dpi 80 \
  --output tests/Diffusion_HHe_2D/output/animations/diffusion_hhe.gif
```

动画显示二维分布和中央切线：H 向外扩散，He 向内混合，并同步展示 `rho`、
速度、压力、`mu` 和 `T`。全部绘图由仓库现有的 `animate.py` 完成。
绘图需要 NumPy、Matplotlib 和 Pillow。
若 Python 环境安装了 Plotly，也可把输出后缀改为 `.html`，生成可拖动时间的交互版。

`check.py` 排除 ghost，以 `rho*X*cell_volume` 积分每个组分的质量，
检查每个已保存时刻的质量、组分质量、总能量、动量和 `sum(X)`，
同时核对输出的 `mu`、`T`，检查双能量修正没有掩盖能量误差。
它还验证 H 峰值下降、原圆形区域以外的 H 增加和质量加权混合指标
`sum(rho*V*4*X_H*X_He)/sum(rho*V)` 增加。
结果保存为 `output/history.csv` 和 `output/validation.json`。
二维积分的质量与能量分别是单位厚度下的 `g/cm` 和 `erg/cm`。

## 关闭扩散的对照

不修改任何源文件，用独立构建目录和可执行文件把扩散系数设为零：

```bash
make -j4 PROBLEM_DIR=tests/Diffusion_HHe_2D \
  BUILD_DIR=tests/Diffusion_HHe_2D/.build/control \
  TARGET=tests/Diffusion_HHe_2D/output_control/flare \
  CPPFLAGS=-DHHE_DIFFUSIVITY_CGS=0
(cd tests/Diffusion_HHe_2D && ./output_control/flare -par parameters_control.txt)
python3 tests/Diffusion_HHe_2D/check.py tests/Diffusion_HHe_2D/output_control \
  --parameters tests/Diffusion_HHe_2D/parameters_control.txt --control
```

两个参数文件只有输出目录不同；关闭扩散后不注册这个 Transport 源项。
对照检查额外要求组分分布和速度保持初值，避免把流体数值混合误认为物理扩散。
更改参数或重复运行前，建议使用新的输出目录，避免混入以前的 CSV。

## 可以改哪里

- `parameters.txt`：网格、终止时间、积分器、输出频率。可把积分器改成 `rk2` 或 `lts_euler`。
- `user_problem.cpp`：初始圆的半径、过渡宽度、密度、压力和扩散系数。
- `diffusion_flux()`：质量通量、焓通量和相应时间步限制。
- `check.py`：若改变圆半径、组分或演化时间，需要同步调整对应的演示判据。

本算例验证传输行为与守恒。扩散算子的解析衰减、跨层共享通量和重试守恒，
已有 `unit_test/test_transport.cpp` 中的独立测试；这里不把外观或守恒误差
当成空间收敛阶、真实扩散系数或跨层误差的证明。
