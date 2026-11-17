# 独立功能测试

每组测试有自己的 `main()`，生成独立可执行文件，不启动 FLARE 模拟主循环，
不需要用户问题文件，也不依赖 Sod 运行结果。使用 C++17 编译器和 GNU Make，
无需额外 C++ 测试框架。

在项目根目录运行全部测试（默认分别编译和检查 1D、2D、3D）：

```bash
make -C unit_test -j4
```

只运行某组、某个维度：

```bash
make -C unit_test -j4 TESTS=wave_timestep DIMENSIONS=1
make -C unit_test -j4 TESTS="mesh facefield godunov" DIMENSIONS="2 3"
```

| TESTS 名称 | 检查内容 |
| --- | --- |
| `optimization` | 转换无堆分配及原公式逐位一致；非均匀网格所有 limiter 的缓存/直接重构一致；随机 exact 通量和波速对比、压力初值复用与无效初值回退、静止气体解析通量、近真空、失败统计 |
| `reconstruction` | PLM_corrected 在非均匀网格上的密度/压力正性、共享斜率系数、阈值连续性、低密度、各方向缓存/直接接口一致及未受限状态逐位兼容 |
| `state` | 原始量/守恒量转换、动能与压力、算术运算、真空、非法状态、各方向物理通量 |
| `eos` | EoS 名称解析、函数指针分派、非法选择及未实现物理函数的明确报错 |
| `species` | 物种表、单态索引、按 `(Z,A)` 查询所有已激活电离态索引（按电荷排序，含负离子） |
| `initialization` | 普通 profile 与 MESA 网格、流体状态读取及错误回滚 |
| `profile_composition` | 专用 `HD_NCOMP=10` 测试：基态输入、可选电离态、同位素别名、MESA 核素映射、质量分数校验与回滚 |
| `lts_euler` | PCM/PLM/PLM_corrected 多层子循环、16:1 步长、同一子步时刻的 PLM 邻域预测、跨层与周期边界守恒、Euler 等价性、界面压力历史与异常回滚、减少界面求解、球对称几何源项、源项限制、激波跨界、重试回滚及维度限制 |
| `lts_rk2` | 独立的多层 Heun 子循环、单层/全细层 RK2 等价性、固定网格上的二阶时间收敛（含非局部时变源项）、跨层/周期守恒、球对称自引力、激波、界面求解次数、阶段源项限制及状态/压力缓存回滚 |
| `mesh` | flat/Index、active/ghost 偏移、const 访问、非均匀间距、越界与损坏的几何 |
| `facefield` | 各方向面数量、flat/Index、const 访问、清零、错误分配和越界 |
| `godunov` | 每个面的左右单元和 gamma、ghost 偏移、散度负号、非均匀网格宽度、面波速归约 |
| `limiter` | 五种斜率限制器、none 中心斜率、线性保持与极值归零、极端数值、MC-beta 端点、配置选择及非法参数 |
| `wave_timestep` | 已有波速回归测试：强激波、真空、CFL、Euler、失败时状态保持不变 |
| `source` | 捕获参数的生命周期与实例独立性、重力/加热 RHS 和 Euler/RK2 更新、加热内能时间尺度、零加热和真空、自定义源项步长限制、球坐标几何源项 |
| `units` | 非平凡质量、长度、时间、温度单位下的常量换算、跨编译单元共享常量、恢复 cgs |
| `output_boundary` | 已有边界/输出测试：多层周期 ghost、角点、用户 hook、参数读取、CSV 和输出频率 |

后两组复用 `tests/test_wave_timestep.cpp` 和 `tests/test_output_boundary.cpp`，
避免复制测试代码。之前的 `tests/run_wave_timestep.sh` 也转到这个统一入口。
每组是独立程序，其中波速/积分器、边界/输出组同时包含相关模块的集成回归检查。

编译产物在 `unit_test/build/ndim1/`、`ndim2/`、`ndim3/`。
框架编译为测试专用静态库，排除 `src/main.cpp`；每个测试只链接它实际引用的对象。
所有维度使用框架默认配置（`NCOMP=NUSER=0`），不读取 Sod 的 `user_config.h`。
旧测试的 `assert` 强制开启；断言或检查失败时返回非零状态并停止本次运行。
参数和 CSV 测试在一次性临时目录中执行，不改动用户输出。

`profile_composition` 不属于默认的无组分测试；它激活 H⁻/H I/H II、He-3 和 He-4
的三个电离态以及中性 C-12，共十个组分。使用独立构建目录运行：

```bash
make -C unit_test -j4 check NDIM=1 CPPFLAGS=-DHD_NCOMP=10 \
  BUILD_DIR="$PWD/unit_test/build/ndim1_profile_composition" \
  TESTS='species profile_composition'
```

将 `NDIM` 和构建目录中的 `ndim1` 改为 2 或 3 可检查多维普通 profile；
MESA 网格构建只在 1D 下执行。

## 组分 profile 输入

先设置 `enable_composition=true`，在 `composition_activate_file` 中激活所需核素和
电荷范围，`HD_NCOMP` 仍须等于所有激活状态数。例如：

```text
h1 -1 1
he4 0 2
```

对应 `HD_NCOMP=6`，但普通 profile 只需给出两个中性态列：

```text
rho vx vy vz pre gamma HI He_I
1   0  0  0  1   1.6666666666666667 0.7 0.3
```

也可将 `HI He_I` 换成核素名 `h1 he4`（核素名不区分大小写）。同一中性态的两种
名称不能同时出现。省略的电离态质量分数设为零；若显式提供 `HII`、`He_II` 等列，
则保留其值并一起参与总和校验。每个已激活核素的中性态列仍然必需，即使其值为零。

MESA 保持原始六行表头格式，通过 `name_base` 匹配 `h1`、`he3`、`he4`、`c12` 等列。
同位素总质量分数全部写入其 `charge=0` 状态。若 MESA 同时含有 `h1` 和 `prot`，
两者相加存入 H I；仅有 `prot` 时也可作为氢输入。元素总量列
`x_mass_fraction_H`、`y_mass_fraction_He` 等不参与累加，以免重复计数。
内侧填充区复制来源壳层的组分；外侧无来源行的填充区仍由用户初始化。

两个读取入口都会拒绝重复/歧义组分列、负值、超过 1 或非有限值，以及
`|sum(X)-1| >= 1e-6` 的行。总和只包括映射到已激活物种的输入；若漏掉有显著丰度的
核素会报错，不会自动归一化。读取阶段不调用 Saha，也不重新计算电离能或 EoS。
