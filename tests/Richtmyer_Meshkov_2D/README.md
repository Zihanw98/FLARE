# Richtmyer–Meshkov 不稳定性（2D）

区域 `[0,2]×[0,1]`，160×80 网格，无重力，gamma=1.4。
初始压力 p=1、速度为零，轻重气体通过平滑界面连接：

```
x_interface(y) = 0.7 + 0.02 cos(2π y)
rho(x,y) = 2 + tanh((x-x_interface)/0.015)   # 左侧1，右侧3
```

从 x=0.35 向右传播一个 Mach 2 平面激波。其后状态由 Rankine–Hugoniot 关系给出：
rho=8/3、p=4.5、vx=1.25 sqrt(1.4)≈1.479、vy=0。
x 下边界持续输入该激波后状态，x 上边界零梯度出流，y 两侧周期。
运行到 t=0.6。激波到达起伏接触层后产生涡量，界面先受压，再出现扰动增长。

```bash
make -C tests/Richtmyer_Meshkov_2D -j4 test
make -C tests/Richtmyer_Meshkov_2D animate
```

检查正密度/压力、初值、后期涡量以及冲击后的界面扰动增长。
**不能简单要求振幅从初始值单调增长**：激波首先压缩起伏。
检查器在预计冲击到达后的窗口找压缩振幅，再要求后期增长且末段拟合斜率为正。
默认结果从初始约 0.020 压到约 0.012，再增长到约 0.019；这个时间段主要是初期增长，
不是已发展出巨大蘑菇头的混合层。

当前框架没有完整的被动组分输运，因此界面使用 `p/rho^gamma=0.6` 的等值线作为近似位置。
轻气体初始该值为1，重气体约0.215，足以在此算例初期区分两边；但它不是物质示踪量，
也会受激波产熵与数值扩散影响。输出的振幅、增长率是诊断代理，不能用作严格理论增长率
或混合质量的验证。单模平滑接触层、低分辨率与一阶耗散都会影响卷曲程度。

质量/能量经 x 边界流入流出，不做闭域守恒断言。`make check` 分析已有输出；
分析与动画只读取 `output/run.log` 清单，请使用 `make test` 更新清单。
输出在 `output/diagnostics/` 与 `output/animations/`。
依赖 numpy、matplotlib、Pillow；MP4 另需 ffmpeg/libx264。

参考：[MPI-AMRVAC 的 RM 测试介绍](https://dev.amrvac.org/md_tests_2hd_2Richtmyer__Meshkov_2Richtmeyer__Meshkov.html)。
此处参数是为 FLARE 构造的轻量案例，不声称复现该项目的全部配置。
