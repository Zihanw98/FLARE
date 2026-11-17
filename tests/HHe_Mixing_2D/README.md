# 二维 H/He 混合：周期 Kelvin–Helmholtz 双剪切层

本算例使用 FLARE 的组分表、守恒通量、PLM/MC 重构、RK2 时间推进和标准 CSV 输出。
编译配置为 `HD_NDIM=2, HD_NCOMP=2, HD_NUSER=0`，激活中性 `h1`、`he4`
（物种名 `HI`、`He_I`）。初始化通过物种表查找索引，不假定输入文件的排列顺序。
原始变量保存质量分数 `X_k`，守恒变量保存 `rho X_k`。

## 运行

在 FLARE 根目录执行：

```bash
make -C tests/HHe_Mixing_2D -j4 test
make -C tests/HHe_Mixing_2D animate
```

`test` 会构建、运行、记录 `output/run.log`，并检查所有输出快照。
`animate` 生成两组 GIF 和首末帧 PNG。需要 Python、NumPy、Matplotlib、Pillow。
如需同时生成 MP4（另需 ffmpeg/libx264）：

```bash
python3 tests/HHe_Mixing_2D/animate.py --format both
```

较快的 32×32 短测试和只检查已有结果：

```bash
make -C tests/HHe_Mixing_2D smoke
make -C tests/HHe_Mixing_2D check
python3 tests/HHe_Mixing_2D/animate.py \
  --parameters tests/HHe_Mixing_2D/output/smoke/parameters.txt
```

`smoke` 在 `output/smoke/` 运行到 `t=0.15`，不覆盖完整算例的输出。
使用 `test` 创建绘图所需的运行清单；框架原有的裸 `make run` 不创建 `run.log`。
输出清单只引用当前运行的快照，所以先前运行遗留的其他 CSV 不会混入。
重新运行相同配置会覆盖该输出目录下的同名快照、日志和诊断。

## 初值与物理解释

默认单位正方形、128×128 活跃单元、两层 ghost、四面周期边界、`gamma=5/3`。
采用精确 Riemann 求解器、PLM、MC、RK2、CFL=0.35，运行到 `t=3`。
设 `x,y` 是相对于区域大小归一化后的坐标：

```text
f(y) = [1 - tanh(cos(2 pi y)/0.15)]/2
rho  = 1 + f
p    = 2.5
vx   = f - 0.5
vy   = 0.02 sin(4 pi x) exp[-(cos(2 pi y)/0.3)^2]
X_He = f
X_H  = 1 - f
```

中心带富 He、上下两带富 H；界面位于 `y≈0.25,0.75`，沿 x 有两个扰动波长。
余弦形式使初始层及其导数周期连续。初态已含一条有限宽度的混合过渡带。
该无黏 Euler 算例不包含显式分子扩散、反应或电离，质量分数随流体平流。
`4 X_H X_He` 为单元平均混合指标：纯组分接近 0，等质量混合为 1；
它包含有限体积数值扩散的影响，不能解释成物理扩散系数或湍流混合率。

## 温度和单位

当前 `src/microphysics/eos/eos.cpp` 中的混合 EOS 仍是占位实现，
`Cell.temperature` 没有随流体更新。本算例在 Python 中按中性单原子理想气体计算温度：

```text
mu = 1 / (X_H/1 + X_He/4)
T[K] = (p_code/rho_code) (unit_length/unit_time)^2 mu m_u/k_B
```

采用质量数近似 `m_i=A_i m_u`，常数与 `src/physcon/constants.h` 一致。
这是一项后处理诊断，不调用未实现的 EOS，也不计算电离平衡。
组分不改变固定的 gamma，但会影响由粒子数推导的温度。

单位由 `initialize_units` 设置，并写入各 CSV：

| 量 | 一个代码单位对应 |
|---|---|
| 长度 | `1e18 cm` |
| 密度 | `1e-24 g/cm^3` |
| 速度 | `1e5 cm/s` |
| 时间 | `1e13 s` |
| 压力、能量密度 | `1e-14 erg/cm^3` |
| 温度 | `1 K` |

图中除温度外均标为代码单位，温度明确标为 K。CSV 保留原有的 `comp_0,comp_1`
列，并增加 `comp_k_name`、`comp_k_A`、`comp_k_charge` 注释头。
Python 根据这些元数据识别 H/He；即使物种索引改变，也不会颠倒颜色或温度权重。

## 输出文件

所有文件均在本算例目录下：

- `output/hhe_mixing_XXXXXXXX.csv`：含 ghost 标记、流体状态和组分质量分数。
- `output/run.log`：完整步进日志及当前运行的 CSV 清单。
- `output/diagnostics/metrics.csv`：各快照的质量、H/He 质量、能量、质量分数范围和混合指标。
- `output/diagnostics/summary.json`：检查结果和守恒误差摘要。
- `output/animations/hhe_mixing.gif`：H、He 质量分数和混合指标三个面板。
- `output/animations/hydro_fields.gif`：密度、压力、温度、速度大小、vx、vy、内能密度、总能量密度八个面板。
- 同目录的 `*_first.png`、`*_last.png`；`--format both` 还生成同名 `.mp4`。

色标在整段动画中固定，速度分量色标以零为中心，组分及混合指标范围固定为 `[0,1]`。
绘图复用根目录 `animate.py` 的 CSV 读取、网格检查和坐标边界函数，
以及 `animate3d.py` 的运行清单读取。ghost 不参与绘图或积分。
默认最多选取 60 帧并保留首末帧，`--max-frames 120 --stride 1 --fps 10 --dpi 100`
可调整输出。每帧对应一个选中快照，播放间隔并非等物理时间。

## 检查与实测

检查初态解析表达式、网格、全时段有限性与正密度/压力、`0≤X_k≤1`、
`X_H+X_He=1`、总质量/H 质量/He 质量/总能量守恒、组分移动、扰动增长及完整终止时间。
完整运行允许的相对守恒误差为 `2e-9`，质量分数误差容限为 `5e-10`。

2026-11-17 默认运行：5307 步，91 个快照，终止时间 3.0。
H、He、总质量及总能量的最大相对漂移均不超过 `2.23e-16`；
`max|X_H+X_He-1|=3.89e-14`，质量分数范围约 `[1.626e-6,0.999998374]`。
质量加权混合指标从 `0.0964393` 增至 `0.3181859`，横向动能峰值为初始值的约 995 倍。
这些是当前网格的回归指标，不是解析 KH 解的收敛证明。

Python 回归测试（包括合成数据、打乱行序、物种索引交换、温度单位、坏数据拒绝、GIF 帧数）：

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
make -C unit_test -j4 test
make -C unit_test -j4 check NDIM=2 CPPFLAGS=-DHD_NCOMP=3 \
  BUILD_DIR="$PWD/unit_test/build/ndim2_species" \
  TESTS='state godunov reconstruction optimization output_boundary species eos lts_rk2'
```

组分构建会执行原有的三组分重构/缓存检查。回归修正同时更新了两项旧预期：
CSV 测试按编译时组分数检查列和物种元数据；LTS-RK2 多维测试由原来的“必须抛错”
改为与全局 RK2 等价、守恒及发生输运的检查。
