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
