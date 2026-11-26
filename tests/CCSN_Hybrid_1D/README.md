# Nuclear hybrid EoS：1D 球对称坍缩反弹演示

这个算例从一个晚期坍缩的 toy core 开始，展示中心密度越过核密度、EoS 变硬、核心反弹和激波向外传播。使用已有的牛顿自引力和球坐标源项，`HD_NDIM=1`、`HD_NCOMP=0`。没有能量注入，也没有调用 Saha。

2026-11-23 的能量复核发现：现有 1536 单元、15 ms 运行的总能量（含自引力势能）
增加了约 `1.2444e52 erg`，可由 PCM/HLL 质量通量与单元中心引力做功的不一致重建。
因此，已有质量守恒和反弹形态检查不构成能量守恒或物理可信度验证，外层接近光速的
结果不能用于物理预测。数值和复现方法见 [VALIDATION.md](VALIDATION.md) 的能量复核。

## EoS

沿用已有的 `Equation_of_State::Nuclear` 枚举。参数暂时集中在
`src/microphysics/eos/nuclear.h`，没有扩展运行配置：

| 量 | 当前设置 |
| --- | --- |
| 核密度 `nuclear_density_cgs` | (2\times10^{14}\ \mathrm{g\,cm^{-3}}) |
| 低密度冷分支指数 \(\Gamma_1\) | `cell.gamma`，本算例为 (4/3) |
| 高密度冷分支指数 \(\Gamma_2\) | 2.5 |
| 热分支指数 \(\Gamma_{\rm th}\) | 1.5，与 `cell.gamma` 分开 |
| \(K_1\) | (4\times10^{14})，CGS；\(\Gamma_1=4/3\) 时单位为 \(\mathrm{dyn\,cm^2\,g^{-4/3}}\) |

令 \(\rho_n\) 为核密度，\(K_2=K_1\rho_n^{\Gamma_1-\Gamma_2}\)：

\[
P_c=\begin{cases}K_1\rho^{\Gamma_1},&\rho\le\rho_n,\\K_2\rho^{\Gamma_2},&\rho>\rho_n.\end{cases}
\]

冷比内能为

\[
e_c=\begin{cases}
\dfrac{K_1\rho^{\Gamma_1-1}}{\Gamma_1-1},&\rho\le\rho_n,\\
\dfrac{K_2\rho^{\Gamma_2-1}}{\Gamma_2-1}
+\dfrac{P_n}{\rho_n}\left(\dfrac1{\Gamma_1-1}-\dfrac1{\Gamma_2-1}\right),&\rho>\rho_n.
\end{cases}
\]

所以阈值处 (P_c,e_c) 连续，且满足 (de_c/d\rho=P_c/\rho^2)。声速在阈值处有意保留导数跳变。

\[
u=\rho e_c+u_{\rm th},\qquad
P=P_c+(\Gamma_{\rm th}-1)u_{\rm th},\qquad
c_s^2=\frac{\Gamma_cP_c+\Gamma_{\rm th}P_{\rm th}}{\rho}.
\]

三种输入统一通过 `evaluate_eos`：

- `RhoE`：从双能量规则选出的内能密度扣除冷内能，得到热能、温度和压强；不改 `U`。
- `RhoP`：先取 (P-P_c)，得到热能，然后建立两套能量。
- `RhoT`：采用 (P_{\rm th}=\rho R_g T/\mu) 和 \(c_V=R_g/[\mu(\Gamma_{\rm th}-1)]\)，再加冷能量。

这里的温度是为接口定义的 **toy 热温度**，本例固定 \(\mu=1\)，不代表真实简并核物质的温度。允许完全冷的 (T=0)。明显的负热能、(P<P_c) 会报错；只容忍冷能量量级 (10^{-12}) 内的减法舍入，不改写守恒量补能量。

这类冷多方分支加热分支的模型可参考 [Einstein Toolkit EOS_Omni](https://www.einsteintoolkit.org/thornguide/EinsteinEOS/EOS_Omni/documentation.html)；核密度变硬的坍缩模型见 [MPA core-collapse simulations](https://wwwmpa.mpa-garching.mpg.de/rel_hydro/axi_core_collapse/)。本算例采用牛顿引力，不含中微子、光致解离、核反应、GR 或真实核物质表；只演示反弹激波形成，不用于预测真实 CCSN 爆炸或中子星结构。

## 初态和计算域

- 以 \(K_{\rm structure}=5\times10^{14}\)、\(n=3\) 的 Lane–Emden 解给出密度，中心密度 (10^{12}\ \mathrm{g\,cm^{-3}})，约为核密度的 1/200。
- 核质量约 (1.485\,M_\odot)，初始核半径约 337 km。初始速度为零，初始热压为冷压的 1%。演化用较小的 \(K_1\)，总压力约为参考平衡压力的 80.8%，因此自引力驱动坍缩。
- 活动区域 (0\le r\le10^{10}\ \mathrm{cm}\)。768 个非均匀球壳，中心单元宽约 0.240 km，外部逐渐拉宽。中心反射，外边界由算例自己的回调复制末端状态。
- 框架已有的几何接口重新建立体积和面面积。没有修改网格、引力、边界或积分器的实现。
- 默认参数使用 `hll` 通量、PCM 和 RK2。将参数文件中的 `riemann_solver` 设为 `exact` 或 `exact_bisection` 可使用 Nuclear 精确黎曼求解器。两者分别用带保护的牛顿法和二分法求星区压力。
- 运行到 10.4 ms，刚好覆盖反弹后的早期外传激波。动画放大内部 180–500 km；计算域仍延伸到 (10^{10}\ \mathrm{cm}\)。

## 运行和看动画

```bash
cd tests/CCSN_Hybrid_1D
make -j4
./flare -par parameters.txt > bounce.log
python3 analyze.py
python3 animate.py
```

重新运行时使用空的输出目录或新的 `output_directory`，避免混入上次的快照。

输出位于 `output_bounce/`：

- `bounce_*.csv`：标准框架快照。
- `grid.csv`、`model.json`：球壳体积和物理参数。
- `history.csv`、`validation.json`：中心密度历史与激波诊断。
- `animations/nuclear_bounce.gif`：密度、速度、热压占比和中心密度历史四联动画。
- `animations/nuclear_bounce_last.png`：最后一帧。

动画只使用保存的数值状态，不对状态插值，并在反弹附近放慢播放。
`analyze.py` 检查核密度穿越、中心密度峰值后的下降、质量守恒、热能非负，以及外传压缩激波的密度/压力跃变和上游激波参考系 Mach 数。激波位置按速度压缩前沿估计，跃变采样覆盖数个单元，不作为精确 Rankine–Hugoniot 拟合。

## 高分辨率、延长运行及阈值对比

`parameters_hires.txt` 使用 1536 个径向单元，中心单元宽约 0.120 km，运行到 15 ms。
每 10 步保存一次；加密后物理时间采样也更密。结果写入独立的 `output_bounce_hires/`。

```bash
# 在本算例目录中运行；输出目录必须为空或使用新的目录名。
make -j4
./flare -par parameters_hires.txt > bounce_hires.log
python3 analyze.py --output output_bounce_hires --t-end 0.015

# 回到仓库根目录，用通用脚本生成核心近景。
cd ../..
python3 animate.py tests/CCSN_Hybrid_1D/output_bounce_hires bounce \
  --plot spherical --cgs --log --rmax 1e7 \
  --radial-grid tests/CCSN_Hybrid_1D/output_bounce_hires/grid.csv \
  --density-threshold 2e14 --map-pixels 801 --dpi 200 \
  --stride 4 --fps 12 \
  --output tests/CCSN_Hybrid_1D/output_bounce_hires/animations/bounce_threshold_100km.mp4
```

将 `--rmax` 改为 `1.5e8`，输出文件名改为 `bounce_threshold_1500km.mp4`，即可查看完整外传激波。
两部视频均为 2400×1800 像素、328 帧、约 27.3 秒；球面对称映射采用 801×801 个采样点。
每帧对应真实保存状态，播放等帧间隔，不等于物理时间匀速播放。

`--density-threshold 2e14` 将 60% 的色标分配给核密度的 0.5–1.5 倍区间；
低密度为蓝色，阈值为浅黄色，高于阈值转为橙红色。分段对数色标保留完整密度范围，
所有帧共用同一尺度，色条标注真实密度，径向剖面增加阈值虚线。此选项要求 `--log`，
阈值单位与显示单位相同（这里为 CGS）；它只改变显示，不改变物理模型或快照数据。
`--map-pixels` 控制球面图采样，`--dpi` 控制导出像素，两者与模拟网格数分别设置。

同一绘图命令可将输出扩展名改为 `.html` 生成可暂停和逐帧查看的离线动画（需要 Plotly）。
本次另生成 `bounce_threshold_100km.html`（401×401 映射以减轻浏览器负担）、
`bounce_threshold_100km.gif`（1600×1200、328 帧）、`bounce_threshold_peak_last.png`
和同一峰值状态的新旧色标对比图 `density_threshold_comparison.png`。

### 原 colormap、密度下限和对数半径

以下命令在仓库根目录执行，使用原来的 viridis 对数密度色标：

```bash
python3 animate.py tests/CCSN_Hybrid_1D/output_bounce_hires bounce \
  --plot spherical --cgs --log --log-x \
  --radial-grid tests/CCSN_Hybrid_1D/output_bounce_hires/grid.csv \
  --rmax 1e7 --profile-range 1e4 1e7 \
  --density-min 1e8 --map-pixels 801 --dpi 200 \
  --stride 4 --fps 12 \
  --output tests/CCSN_Hybrid_1D/output_bounce_hires/animations/bounce_viridis_cut1e8_logr_100km.mp4
```

| 选项 | 含义（使用 `--cgs` 时） |
| --- | --- |
| `--log-x` / `--log-radius` | 右侧所有径向剖面的半径使用对数轴；左侧由 `--log-map-radius` 独立控制 |
| `--log-map-radius` | 左侧圆图也按对数半径展开，保持角度和圆形球壳 |
| `--map-linthresh 1e4` | 左侧中心 0.1 km 内保留线性半径，处理 r=0；需要 `--log-map-radius` |
| `--profile-range 1e4 1e7` / `--xlim 1e4 1e7` | 右侧半径显示 0.1–100 km；对数轴要求两个端点严格为正 |
| `--rmax 1e7` | 左侧球面图显示至半径 100 km；未指定剖面范围时也控制剖面上限 |
| `--density-min 1e8` | 密度色条和剖面纵轴从 10⁸ g cm⁻³ 开始，较低密度使用色条底端颜色 |
| `--density-max 4e14` | 可选：指定密度色条和剖面纵轴上限，较高密度使用色条顶端颜色 |

未指定 `--profile-range` 时，对数半径的默认下限取第一个正半径单元中心，避免对零取对数。
不使用 `--cgs` 时，所有范围参数都按代码单位解释。范围与对数轴选项同时支持 MP4、GIF、PNG 和 HTML；
原始快照、径向坐标和 HTML 悬停读数均保留真实物理值。密度上下限不能与 `--density-threshold` 同用。

1500 km 全景使用 `--rmax 1.5e8 --profile-range 1e4 1.5e8`，并换用相应输出文件名。
剖面范围可以独立于左侧图设置，例如左侧放大 100 km、右侧查看 0.1–1500 km。

左右同时采用对数半径的演化动画：

```bash
python3 animate.py tests/CCSN_Hybrid_1D/output_bounce_hires bounce \
  --plot spherical --cgs --log --log-x \
  --log-map-radius --map-linthresh 1e4 \
  --radial-grid tests/CCSN_Hybrid_1D/output_bounce_hires/grid.csv \
  --rmax 1e7 --profile-range 1e4 1e7 \
  --density-min 1e8 --map-pixels 801 --dpi 200 \
  --stride 4 --fps 12 \
  --output tests/CCSN_Hybrid_1D/output_bounce_hires/animations/bounce_viridis_cut1e8_logmap_100km.mp4
```

`--log-map-radius` 要求 `--plot spherical`。未指定 `--map-linthresh` 时，
中心线性区半径取第一个正半径单元中心。其单位与绘图一致（`--cgs` 时为 cm），且必须小于左图外半径。
显示半径采用 `s=r/r0`（r≤r0）和 `s=1+ln(r/r0)`（r>r0），其中 r0 为中心线性区半径。
映射保持圆形和角度；对数区的相邻十倍半径占据相同径向距离，图上面积不代表物理面积。
虚线圆环和沿轴刻度对应真实半径；HTML 悬停显示采样球壳中心的真实半径与物理量。
MP4、GIF 和 HTML 共用同一球壳映射，原始数据不变。1500 km 版本将两个范围上限改为 `1.5e8`。

已生成 1500 km GIF：`output_bounce_hires/animations/bounce_viridis_cut1e8_logmap_1500km.gif`，
分辨率 2400×1800，328 帧，约 27.33 秒，循环播放；其物理时间与视图范围和同名 MP4 一致。

把速度色标固定为 −0.1c 到 +0.1c（超范围使用两端颜色，速度剖面保留完整数值）：

```bash
python3 animate.py tests/CCSN_Hybrid_1D/output_bounce_hires bounce \
  --plot spherical --cgs --log --log-x \
  --log-map-radius --map-linthresh 1e4 \
  --radial-grid tests/CCSN_Hybrid_1D/output_bounce_hires/grid.csv \
  --rmax 1.5e8 --profile-range 1e4 1.5e8 \
  --density-min 1e8 --velocity-min=-2.99792458e9 --velocity-max=2.99792458e9 \
  --map-pixels 801 --dpi 200 --stride 4 --fps 12 \
  --output tests/CCSN_Hybrid_1D/output_bounce_hires/animations/bounce_viridis_cut1e8_logmap_1500km_vlim0p1c.gif
```

EoS 单元测试：

```bash
make -f unit_test/Makefile check NDIM=1 TESTS=nuclear
```

该命令在仓库根目录运行，覆盖三个输入模式、冷阈值连续性、冷能量热力学导数、等熵声速、非 CGS 单位、非法热能和 HLL 均匀/反射状态。
