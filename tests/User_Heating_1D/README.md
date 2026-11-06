# 用户自定义恒定加热源项

一维 Cartesian 周期区域 `[0, 1]` 中，初始气体均匀：`rho=1`、`p=1`、速度为零。
用户源项以恒定的单位质量加热率 `q` 注入能量：

```text
dE/dt = rho*q
rho(t) = 1
v(t) = 0
p(t) = 1 + (gamma-1)*q*t
```

由于没有压力梯度，气体保持静止，密度不变，压力线性增长。
默认 `gamma=1.4`、`q=0.5`、`t_end=2`，最终压力应为 `1.4`，
单位体积总能量从 `2.5` 增长到 `3.5`。

## 运行和验证

在项目根目录运行：

```bash
make -C tests/User_Heating_1D -j4 test
```

默认结果保存在本目录的 `output/` 和 `run.log`；
`output/diagnostics/metrics.csv` 保存解析压力、最大压力误差、质量和总能量。
检查脚本只依赖 Python 标准库。

`test` 还会在临时目录中用 Euler、RK2 分别测试 `q=0`、`0.5`、`1.5`，
核对全部保存快照中的 active/ghost 状态、质量和能量收支，并检查源项步长限制。
`q=1.5` 的验证特意缩小能量变化比例，使首步由加热限制决定。
临时算例不会改动默认输出。另检查负加热率、零步长比例和未知参数名被拒绝。
同时比较内置 `constant_heating` 与用户实现：Euler/RK2、`q=0.5/500`、
相同的 `0.1` 时间尺度系数，检查全部 CSV 和运行日志逐字节一致。

只运行模拟或检查已有结果：

```bash
make -C tests/User_Heating_1D run
make -C tests/User_Heating_1D check
```

`run` 把日志输出到终端；`check` 需要 `test` 保存的 `run.log`。
更改参数后请重新运行 `test`，避免使用旧输出。

生成包含解析解对照的 GIF 动画（需要 Matplotlib 和 Pillow）：

```bash
make -C tests/User_Heating_1D animate
```

动画保存在 `output/animations/heating_analytic.gif`，左侧展示压力空间剖面，
右侧展示压力随时间增长；使用真实保存帧，不插值，按物理时间间隔播放。
末帧另存为 `heating_analytic_last.png`。

## 自定义参数

编辑 `heating.txt`，无需重新编译：

```text
specific_heating_rate = 0.5
max_fractional_energy_change = 0.1
```

- `specific_heating_rate`：每单位质量、每单位时间增加的能量，必须非负；零表示不加热。
- `max_fractional_energy_change`：允许单步增加的内能比例，范围 `(0,1]`。
  源项给出 `dt <= fraction * internal_energy / (rho*q)`，最终步长仍与流体 CFL 限制取最小值。

这是算例自己的参数文件，由 `user_problem.cpp` 读取、校验，未向 `SimulationConfig`
增加字段。文件路径相对于运行目录；上述 Makefile 命令会在本算例目录运行。
直接调用可执行文件时，也应在包含 `heating.txt` 的目录运行。

## 用户源项的实现位置

全部实现都在 `user_problem.cpp` 中：

1. `HeatingParams` 定义源项参数，`read_heating_params()` 读取并校验。
2. `heating()` 向 active 单元的能量 RHS 累加变化率，不乘 `dt`；必要时收紧步长。
3. `make_heating()` 用 lambda 按值捕获参数，返回 `SourceTerm`。
4. `register_sources()` 把源项加入列表。
5. `make_user_hooks()` 设置 `hooks.register_sources`。

参数在创建源项时读取一次，回调没有会随调用次数变化的物理状态，因此适用于 RK2 重试。
以后新增其他源项，可以在同一个注册函数中继续 `push_back`，不必修改框架的 `make_sources()`。

## 对照内置加热

框架也提供 `make_constant_heating(double heating_rate)`，只捕获一个 double。
在 `parameters.txt` 中设置以下两项即可启用：

```text
constant_heating_enabled = true
heating_rate = 0.5
```

内置源项用固定的 `0.1` 系数限制加热时间尺度。
在本算例中切换到内置加热时，将 `heating.txt` 的 `specific_heating_rate` 改为 `0`，
避免两个源项同时加热。`make test`/`check` 支持用户加热、内置加热以及两者叠加，
解析压力按总加热率计算。`test` 还会在隔离目录中自动运行两种实现的对照测试。
