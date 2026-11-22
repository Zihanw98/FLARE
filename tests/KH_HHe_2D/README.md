# H/He Kelvin–Helmholtz 混合演示

这是独立的新算例，不依赖旧 KH 算例的输出。使用二维笛卡尔网格、两种中性组分 H I/He I，以及 `IdealGas_Mixed`。`HD_NDIM=2, HD_NCOMP=2, HD_NUSER=0`。

## 初始状态

- 区域：`[0,1e8 cm]^2`，128×128 活动网格，四边周期边界。
- 中间 `0.25L<y<0.75L` 为 H I；上下两侧为 He I。
- 初始 rho=1e-10 g/cm³、P=100 dyn/cm²，gamma=5/3。
- 中间流体向右，上下两侧向左，速度约 ±5e5 cm/s；界面处速度在厚度约 0.02L 内平滑过渡。
- 加入振幅 1e4 cm/s 的竖直速度扰动，沿 x 有两个波长，扰动集中在两条组分边界附近。
- 演化到 300 s。两侧速度差为 1e6 cm/s，因此 L/Δv=100 s。
- exact、PLM/MC、RK2、CFL=0.35。

初始压力保持平衡，温度随组分变化：纯 H 的 mu=1、T≈12027 K，纯 He 的 mu=4、T≈48109 K。由 `Set_Conserved_From_Primitives` 的 RhoP 模式建立能量；共用初始化流程填充 ghost 并更新温度、数密度等。

初始条件写在 `user_problem.cpp` 顶部常量和 initialize 中；网格、终止时间、数值方法和输出设置在 `parameters.txt`。初始 H 带的边缘采用单元覆盖比例，默认 128² 下边界正好落在网格面上。这个算例直接初始化，不需要生成 profile。

## 运行和检查

从项目根目录运行：

```bash
make -j4 PROBLEM_DIR=tests/KH_HHe_2D
make PROBLEM_DIR=tests/KH_HHe_2D run
python3 tests/KH_HHe_2D/check.py
```

重复运行或更改配置时请更换 `output_directory`，避免不同运行的快照混在一起。`check.py OTHER_OUTPUT --parameters OTHER_PARAMETERS` 可检查另一份配置的输出。

检查器验证所有快照的各组分质量、总质量、总能量、总动量、组分范围/总和，以及 mu、T 的独立公式。还检查竖直流动增长、组分边界沿 x 的变形，以及 H 进入初始 He 区域的质量比例，用来确认演示确实发生了卷入。后续压力和密度会随 KH 流动变化，不要求它们始终均匀。

结果保存为 `output/validation.json` 和 `output/history.csv`。二维质量、能量按第三方向单位长度积分。`mixed_mass_fraction` 是处于 `0.1<X_H<0.9` 的单元质量占比；`mixing_measure` 是质量加权的 `4*X_H*X_He`。

## 动画

使用项目通用绘图脚本，不需要新的绘图实现：

```bash
python3 animate.py tests/KH_HHe_2D/output kh_hhe \
  --composition --thermodynamics --cgs --slice-y 2.5e7 \
  --fps 12 --dpi 80 \
  --output tests/KH_HHe_2D/output/animations/kh_hhe.gif
```

左侧显示 rho、vx、vy、P、H/He 质量分数、mu 和 T；右侧取初始下方剪切层附近的切面。要在 rho/P 图上叠加速度流线，可再加 `--streamlines`。

## 如何理解混合

动画展示剪切不稳定性把两种组分卷入涡旋、拉成条带。没有启用 Saha、辐射、引力、加热、显式黏性或组分物理扩散。组分作为固定电离态随流体运动；网格内的中间质量分数包含有限体积平均和数值扩散。混合指标不是分子扩散系数或已经收敛的物理混合率。

这是可视化和守恒检查算例，未以网格收敛或解析增长率验证 KH 的定量动力学。更改初始参数后，检查器中用于确认明显卷入的阈值未必适用于新的结束时刻。

已完成的 64²、128² 运行结果见 `VALIDATION.md`。
