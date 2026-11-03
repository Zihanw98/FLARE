# 周期声波（2D）

单位正方形，四面周期边界，`64×64` 活跃单元。理想气体 gamma=1.4，无重力。
背景 rho=1、p=1/gamma，因此声速 c=1。波矢 `k=(2π,2π)`，沿对角线传播；
密度扰动幅度 epsilon=1e-4。初始化：

```
δ = epsilon sin(k·x) sinc(kx Δx/2) sinc(ky Δy/2)
rho = 1 + δ
p   = 1/gamma + δ
v   = δ k/|k|
```

这里 `sinc(a)=sin(a)/a`，对正弦扰动做单元平均。原始变量的扰动平均是精确的；
转换成守恒变量的非线性乘积在 O(epsilon²) 上不是严格的解析单元平均。
结束时间 `1/sqrt(2)` 为一个线性声波周期。

从仓库根目录：

```bash
make -C tests/Periodic_Acoustic_2D -j4 test
make -C tests/Periodic_Acoustic_2D animate
make -C tests/Periodic_Acoustic_2D converge
```

检查所有保存时刻的 Fourier 相位（防止“一周期末不动也通过”）、振幅、密度 L1 误差，
以及周期域质量、总能量守恒。解析参考是线性化 Euler 解，不是有限振幅 Euler 的严格解；
epsilon 很小，当前网格上的误差主要是数值耗散。`converge` 保留默认输出，另跑 `128×128`，
要求误差比低于 0.8，结果在 `output/convergence/convergence.json`。

默认实测：64² 一周期振幅保留 0.7866，128² 保留 0.8869；密度 L1 从 1.357e-5
降为 7.196e-6，观测收敛阶约 0.915，符合当前一阶格式的表现。

输出清单在 `output/run.log`，指标在 `output/diagnostics/`，动画在 `output/animations/`。
依赖 numpy、matplotlib、Pillow；MP4 另需 ffmpeg/libx264。`make check` 只分析已有输出。
使用 `test` 运行以更新清单；框架裸 `make run` 不写检查器所需的新运行清单。

参考：[Euler 方程及声学特征速度](https://www.clawpack.org/riemann_book/html/Euler.html)。
