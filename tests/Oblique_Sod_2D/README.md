# 斜向 Sod（2D）

单位正方形，128² 网格，无重力，gamma=1.4，t_end=0.15。
定义法向坐标 `s=(2(x-0.5)+(y-0.5))/sqrt(5)`，传播法向与 x 轴夹角为 26.565°。
`s<0` 为 `(rho,p,v)=(1,1,0)`，`s>=0` 为 `(0.125,0.1,0)`。
注意这个角度描述波的传播方向，不是间断线本身相对 x 轴的角度。

所有边界使用随时间变化的平面 Sod 解析解；零梯度边界会破坏斜向一维解，
因此本算例不使用普通 outflow。初值和边界以 8×8 子单元求积平均守恒变量。
这是求积近似，不声称在跨间断单元获得严格解析平均。

```bash
make -C tests/Oblique_Sod_2D -j4 test
make -C tests/Oblique_Sod_2D animate
```

独立 Python 标量压力求根给出 Sod 参考，沿 `s` 采样并同样做守恒量求积。
检查密度、法向速度、压力 L1 和理论上为零的切向速度 RMS。
默认实测依次约为 0.01393、0.02435、0.01259、0.006715。
图中散点投影到波法向，与解析剖面对照；这能显示数值扩散和方向误差。

边界有质量与能量流动，不要求全域质量/能量恒定。该算例检查斜向传播的误差，
没有把 26.565° 和其他角度之间的差异直接判为失败。
若改变分辨率，可复用初始化；改变 gamma 会明确拒绝运行，因参考为 gamma=1.4。

`make check` 分析已有输出。所有分析只使用 `output/run.log` 中的文件；
`test` 更新清单，裸 `make run` 不更新。指标、末帧、动画分别在
`output/diagnostics/`、`output/animations/`。
Python 依赖 numpy、matplotlib、Pillow；MP4 另需 ffmpeg/libx264。

参考：[Euler Riemann 解](https://www.clawpack.org/riemann_book/html/Euler.html)。
