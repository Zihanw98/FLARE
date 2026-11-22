# 二维组分输运和数值扩散

这个算例通过 profile 初始化，测试组分守恒、接触压力平衡、数值扩散，以及输出的平均分子量和温度。没有物理扩散、Saha、引力或加热。

编译配置为 `HD_NDIM=2, HD_NCOMP=2, HD_NUSER=0`。默认使用 exact、PLM/MC 和 RK2。

| `--case` | 组分和分布 | 速度 | 目的 |
| --- | --- | --- | --- |
| `stationary` | 圆内 H I、圆外 He I | 0 | 静止接触是否保持、是否产生虚假压力波 |
| `smooth` | H I/He I 平滑周期团块 | vx=vy=1e6 cm/s | 定量误差、网格收敛、峰值衰减 |
| `sharp` | 圆内 H I、圆外 He I | vx=vy=1e6 cm/s | 界面展宽、方向误差、组分越界 |
| `ionized` | 圆内 H I、圆外 H II | vx=vy=1e6 cm/s | 固定电离态随流输运，检查电离能和压力平衡 |

区域为 `[0,1e8 cm]^2`，四边周期边界；rho=1e-10 g/cm³，P=100 dyn/cm²，gamma=5/3。t_end=100 s，运动版本恰好平移一个周期。圆心在 `(0.35L,0.5L)`，半径 `0.15L`。

平滑版本采用 `X0=0.1+0.8*g(x-xc)*g(y-yc)`，其中 `g(s)=cos(pi*s/L)^8`，`X1=1-X0`。这是严格周期的光滑函数。profile 保存其解析单元平均值；圆形版本保存圆与单元的重叠面积比例，所以切过圆边界的单元初始就是混合组分。

输入只有 rho、v、P、gamma、X，读取器通过 RhoP 建立能量。温度随 mu 改变，不强制温度均匀。Mixed 下 `mu=1/(X_H+X_He/4)`；Ionized 下 `mu=1/(X_HI+2*X_HII)`，内能还包括 H II 的电离能。没有反应导致的电离态变化。

## 编译、生成和运行

以下命令在项目根目录执行：

```bash
make -j4 PROBLEM_DIR=tests/Composition_Advection_2D
python3 tests/Composition_Advection_2D/generate.py --case sharp --n 128
tests/Composition_Advection_2D/flare \
  -par tests/Composition_Advection_2D/output_sharp_128/parameters.txt
python3 tests/Composition_Advection_2D/check.py \
  tests/Composition_Advection_2D/output_sharp_128
```

生成目录内有 `profile.dat`、`composition.txt`、`parameters.txt` 和供验证使用的 `model.json`。快照在该目录的 `output/` 下。生成器不会覆盖已有输入；需要新运行时用 `--directory` 指定新目录。运行模拟前也应使用新的输出目录，避免与旧快照混合。

`parameters.txt`（算例根目录）是生成器模板，不附带大型 profile 文件。修改初始分布可编辑 `generate.py` 或生成后的 `profile.dat`；修改解析模型后，应同步更新 `model.json` 和验证公式，不能再把原解析解用于自定义问题。

其他版本只需改变 `--case`；不同网格只改 `--n`，不用重新编译。比较 PCM：

```bash
python3 tests/Composition_Advection_2D/generate.py \
  --case smooth --n 64 --reconstructor PCM \
  --directory tests/Composition_Advection_2D/output_smooth_pcm_64
tests/Composition_Advection_2D/flare \
  -par tests/Composition_Advection_2D/output_smooth_pcm_64/parameters.txt
```

也可用 `--solver hll` 比较接触扩散。HLL 的静止组分界面可能展宽；检查器只对 exact 静止版本要求接近舍入误差。移动算例的 L1 阈值是这些中等分辨率测试的回归范围，不是任意分辨率/求解器的精度保证。

## CSV 输出和检查

框架 CSV 保留所有原列的位置，在末尾追加：

- `mean_molecular_weight`：无量纲 mu；
- `temperature`：代码温度，乘文件头 `unit_temperature` 得到 K。本算例基准单位全为 1，因此该列直接是 K。

已有的 `comp_0/comp_1` 是质量分数 X，不是 rho*X；物种名称、A 和 charge 仍放在文件头。ghost 行同样包含这些量，可用 `is_ghost=0` 筛选活动单元。输出在临时 Cell 上从当前 U 恢复热力学量，不修改网格 U 或缓存。真空的 T、mu 无定义，写 0 作为标记。

`check.py` 检查所有已保存时刻：各物种质量、总质量和总能量、X 的取值和总和、rho/P/v 平衡、mu/T 的独立公式、电离能与内能关系，以及解析平移后的单元平均组分。中间时刻也参与检查，避免“根本没有输运，但一个周期后与初态一致”的误通过。

结果写入 `validation.json` 和 `history.csv`。二维积分的质量/能量按第三方向单位长度计算。`composition_L1` 是全域平均绝对质量分数误差，不是相对误差；不连续界面的 Linf 不用于要求二阶收敛。`mixed_area_fraction` 统计 `0.1<X0<0.9`，主要用于尖锐界面版本；初始切割单元本身已有混合，应对比初始值和末值。`mixing_measure=mean(X0*(1-X0))` 反映展宽，但不代表物理混合速率。

## 使用通用脚本绘图

```bash
python3 animate.py tests/Composition_Advection_2D/output_sharp_128/output blob \
  --composition --thermodynamics --cgs --fps 12 \
  --output tests/Composition_Advection_2D/output_sharp_128/animations/composition_thermodynamics.gif
```

显示 rho、vx、vy、P、两个组分、mu、T 的二维图和固定切面。切面默认 y=L/2；团块穿过周期边界时会出现在两侧。`--thermodynamics` 对旧 CSV 是可选的，缺少新列的旧输出仍能按原命令绘制。`--log` 仍只影响 rho 和 P；mu/T 保持线性。

复现结果见 `VALIDATION.md`。
