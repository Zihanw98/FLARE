# MESA CCSN：独立隐式积分器

沿用 `CCSN_MESA_day13` 的单位、初始化和边界代码，三个接口完全不改。
切除约 4.00126 M_sun，外侧补到 1800 个代码长度，补层密度 1e-15 g/cm³；
第一活动壳注入 1.6e51 erg。保留导入网格，不生成多方球。

`parameters.txt` / `parameters_adaptive.txt` 使用此前合并后的 3892 行 profile，
切除后 2573 壳 + 100 外壳；`parameters_native.txt` 使用原始半径修复 profile，
切除后 2599 壳 + 100 外壳。空间离散均为 PLM_corrected / MC / exact，
exact 容差设为 1e-12。时间单位 1592.8584656366752 秒。

## 运行

`parameters.txt` 默认使用合并网格和自适应隐式；等价配置为 `parameters_adaptive.txt`。
`parameters_native.txt` 使用原始网格和自适应隐式；`parameters_fixed.txt` 保留固定步长验证。

```ini
# 下列设置已写入自适应算例的参数文件。
time_integrator = implicit
dt_implicit = 6.2780216922807006e-05  # 0.1 物理秒，初始尝试及上限
implicit_adaptive = true
implicit_tolerance = 1e-6
implicit_newton_max_iterations = 60
implicit_report_file = output_config_adaptive/solver.csv
```

```sh
make -C tests/CCSN_MESA_implicit_standalone -j4
cd tests/CCSN_MESA_implicit_standalone

# 合并网格，目标 10 物理秒；无需额外隐式命令行参数。
./flare -par parameters.txt > config_adaptive.log 2>&1

# 原始网格，同样的步长上限和目标时间。
./flare -par parameters_native.txt > config_native.log 2>&1

# 固定 0.1 秒：参数文件中 implicit_adaptive=false，失败直接退出。
./flare -par parameters_fixed.txt > config_fixed.log 2>&1
```

新运行分别写入 `output_config_adaptive`、`output_config_native`、`output_config_fixed`，
此前命令行版的 `output_adaptive`、`output_native`、`output_fixed` 保留供对比。

## 2026-11-09 独立命令行版验证结果

| 配置 | 接受 / 拒绝尝试 | 到达物理时间 | 接受步长范围（秒） | 整次运行墙钟时间 | 峰值 RSS |
| --- | --- | --- | --- | --- | --- |
| 合并网格，固定 0.1 秒 | 17 / 1 | 1.7 秒，下一步失败 | 0.1 | 3.31 秒 | 11.71 MiB |
| 合并网格，自适应 | 125 / 16 | 10 秒，完成 | 7.9191e-4 – 0.1 | 20.94 秒 | 11.74 MiB |
| 原始网格，自适应 | 190 / 56 | 10 秒，完成 | 3.0136e-5 – 0.1 | 56.58 秒 | 11.51 MiB |

这些是当次实测，期间有并行编译和另一个算例运行，不作为严格性能基准。
**固定 0.1 秒尚未跑通；成功的是明确启用缩步的自适应模式。**
失败尝试均不推进时间。没有放宽默认 `1e-6` 残差容差，也没有换空间格式。

独立读取每一步快照重算 `U_new-U_old-dt*RHS_new`：合并网格全部 125 步、
原始网格全部 190 步通过；最大缩放残差分别为 `9.83105e-7`、`9.47844e-7`。
所有保存状态密度、压力为正，网格与导入 profile 对应一致，注入能量为 1.6e51 erg。
最大相对质量变化分别为 `1.32e-9` 和 `9.57e-11`。

```sh
python3 check.py --parameters parameters_adaptive.txt
make check-residual CONFIG=parameters_adaptive.txt
python3 check.py --parameters parameters_native.txt
make check-residual CONFIG=parameters_native.txt
```

`check.py` 检查 profile 映射、初态、注入与保存状态，质量相对变化阈值 1e-6；
`verification/check_residual.cpp` 独立检查每一步守恒残差，不调用隐式求解器的
残差/变量实现。检查器输出的 `max_scaled_residual` 已除以容差 1e-6，应小于约 1。
这是爆炸初期的收敛验证，尚未验证时间离散精度或完整超新星演化。

原显式兼容性：Sod 400 单元到 t=0.2 的 46 个 CSV，以及 KH2D 200×200 到
t=0.01 的 23 个 CSV，与加入模块之前逐字节一致。17 组单元测试 × 3 维通过。
完整模块说明见 `docs/implicit.md`。

## 接入配置文件与统一积分器后的复验

改为 `time_integrator=implicit`、参数文件控制步长，并经 `advance_one_step`
调用后，合并网格与原始网格再次完成 10 秒，分别接受 125 / 190 步，拒绝
16 / 56 次尝试。包括初态的全部 126 / 191 个快照与此前独立命令行版逐字节一致；
每步独立守恒残差、正性和 profile 映射检查通过。

本次整次墙钟时间为 20.12 / 54.68 秒，峰值 RSS 为 12064 / 12060 KiB。
固定 0.1 秒仍在接受 17 步后失败，没有静默转为自适应。
旧 Sod 和 KH2D 的 46 / 23 个输出文件也逐字节一致。
17 组单元测试分别在 NDIM=1、2、3 下通过；新增配置解析、自动缩步、
成功增步、跨步历史和末步截断检查。二维、三维目前仅有小网格验证，尚无大规模性能基准。

本次接入不修改 `UserHooks`、`SourceTerm`、重力实现与隐式模块的私有线性代数、
变量变换内核。工作区同时存在的 LTS 和压力缓存命名改动属于其他工作，保留原样。
详细复验记录位于仓库根目录 `output_implicit_config_20261109/validation.json`。
