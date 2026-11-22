# Nuclear hybrid EoS：1D 球对称坍缩反弹演示

这个算例从一个晚期坍缩的 toy core 开始，展示中心密度越过核密度、EoS 变硬、核心反弹和激波向外传播。使用已有的牛顿自引力和球坐标源项，`HD_NDIM=1`、`HD_NCOMP=0`。没有能量注入，也没有调用 Saha。

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
- 使用新加的 `hll` 通量、PCM 和 RK2。`exact` 目前不支持 Nuclear，错误配对会给出明确提示。HLL 当前要求左右密度均为正，本例有稀薄大气。
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

EoS 单元测试：

```bash
make -f unit_test/Makefile check NDIM=1 TESTS=nuclear
```

该命令在仓库根目录运行，覆盖三个输入模式、冷阈值连续性、冷能量热力学导数、等熵声速、非 CGS 单位、非法热能和 HLL 均匀/反射状态。
