# 50×50 二维 H / He 压力不平衡 profile 算例

均匀 Cartesian 网格，区域 `[0,1]×[0,1]`，每边两层 ghost，四面周期边界。
框架按 `parameters.txt` 生成网格，再从 `profile.dat` 读取 2500 个活跃单元的初值。
`user_problem.cpp` 中的初始化 hook 保留读入状态，可在这里接入后续 Saha 初始化。

沿 x 方向在区域中点 `x=0.5` 分界；每侧 25×50 个单元，沿 y 不变：

| 量 | 左侧：纯 H-1，`x<0.5` | 右侧：纯 He-4，`x>=0.5` |
| --- | --- | --- |
| 密度 `rho` | 1 | 4 |
| 压力 `pre` | 1 | 0.1 |
| 速度 `vx, vy, vz` | 0, 0, 0 | 0, 0, 0 |
| `gamma` | 5/3 | 5/3 |
| 中性氢质量分数 `HI` | 1 | 0 |
| 电离氢质量分数 `HII` | 0，由读取器补零 | 0，由读取器补零 |
| 中性氦质量分数 `He_I` | 0 | 1 |

按核质量 `A*m_u` 近似，初始核数密度为 `n=rho/(A*m_u)`，因此
`rho_He=4*rho_H` 保证两侧数密度相同。初始全为中性原子，粒子数密度也相同。
压力特意不平衡：`P_H/P_He=10`，按理想气体关系 `P=n*k_B*T`，
对应 `T_H/T_He=10`。温度不作为 profile 输入；当前由密度和压力指定初态。

这里的 `HI` 就是 H I（中性氢）；当前物种表没有 `H_I` 这个列名。
`composition.txt` 激活 `h1 0 1` 和 `he4 0 0`，因此 `HD_NCOMP=3`。
profile 只提供中性态列，不需要添加全零的 HII 列。当前模板使用固定 gamma
流体求解，没有执行 EoS/Saha 初始化。

## 运行

从项目根目录运行：

```bash
make -C tests/Hydrogen_Profile_2D -j4 run
```

`make -C tests/Hydrogen_Profile_2D -j4 test` 也是运行此算例。
默认运行到 `t=10`，CSV 保存在本目录 `output/` 下，前缀为 `hhe_profile`；文件中 `is_ghost=0`
表示活跃单元，`comp_0=HI`、`comp_1=HII`、`comp_2=He_I`，物种名称也写在 CSV 注释头中。
压力差会驱动流动。x 方向仍为周期边界，因此除了 `x=0.5`，域的首尾相接处
也有一个 H / He 界面；较晚时刻两处界面产生的波会相互作用。要观察早期演化，
可把 `t_end` 改为 `0.1`。

## 动画

运行后，生成流体量以及 HI、HII、He_I 质量分数的动画：

```bash
make -C tests/Hydrogen_Profile_2D animate
```

算例的 `animate.py` 从参数文件读取输出目录和前缀，调用根目录 `animate.py`
的 `--composition` 选项。组分名称从 CSV 的 `comp_k_name` 读取，质量分数色标和
剖面纵轴固定为 `[0,1]`，全零的 HII 也会显示。GIF 和末帧 PNG 默认写入
`output/animations/`。`--cgs` 不改变质量分数，`--log` 只影响密度和压力。

也可直接使用通用脚本或输出交互式 HTML：

```bash
python3 animate.py tests/Hydrogen_Profile_2D/output hhe_profile --composition
python3 tests/Hydrogen_Profile_2D/animate.py --plotly
```

只生成动画不会重跑模拟。若重新运行时改变了组分或输出步数，先把旧 CSV 移到其他
目录，避免不同运行的快照混在一起；脚本会拒绝组分、网格或时间序列不一致的输入。
此次改用 `hhe_profile` 前缀，与之前 `hydrogen_profile` 的输出区分。

## 修改初值

可以直接编辑 `profile.dat`，普通编译和运行不会重新生成它。格式为：

```text
rho vx vy vz pre gamma HI He_I
1 0 0 0 1 1.6666666666666667 1 0
4 0 0 0 0.10000000000000001 1.6666666666666667 0 1
```

上面两行分别展示左右两侧的状态。实际文件首行为列名，后面恰好 2500 行，不含 ghost。
按先 x 后 y 的顺序排列：第一行数据为 `(i,j)=(0,0)`，第 50 行为 `(49,0)`，
第 51 行为 `(0,1)`。零基索引 `(i,j)` 对应文件物理行号 `2 + j*50 + i`。
2D 输入仍需保留 `vz` 列，默认填零。

批量修改时，编辑 `generate_profile.py` 顶部的 `RHO_H`、`PRESSURE_H`、
`PRESSURE_HE` 或 `initial_state(x, y, x_mid)`，再执行：

```bash
make -C tests/Hydrogen_Profile_2D inputs
```

此命令会覆盖 `profile.dat`。脚本从 `parameters.txt` 读取网格大小、区域和 gamma，
分界位置自动取 x 区域的中点。
修改网格大小后也需要同步更新 profile 行数。直接编辑 profile 时，每行的 gamma
才是该单元使用的值；仅修改参数文件不会改写已有 profile。

输入中的密度、压力、速度均采用代码单位。当前单位缩放为
`M0=1 g, L0=1 cm, t0=1 s, Temp0=1 K`，可在 `user_problem.cpp` 修改。
若更改单位，profile 中的数值不会自动转换。

调整氢氦比例时保持 `HI+He_I=1`。若要手动指定部分电离，可在 profile 增加
`HII` 列，并保证 `HI+HII+He_I=1`。
若通过生成脚本添加组分列，需要同时修改返回字典和 `columns`。
增加其他核素时，同步修改 `composition.txt`、`user_config.h` 中的 `HD_NCOMP`
以及 profile 的组分列。
