# 3D 球形高压区膨胀

现有 Godunov + Euler + 理想气体求解器的简单三维集成测试，不修改框架。
单位立方体内初始密度均为 1、速度为 0，中心半径 0.12 的球内压力为 5，
外部压力为 1，gamma=1.4，无重力。初值按单元中心采样，因此低分辨率球面有阶梯形状。
默认 32×32×32 活跃单元，六面反射边界，运行至 t=0.08。
高压球膨胀，向外推动压缩波，中心气体变稀。默认时间用于观察早期膨胀；
继续延长时间后会出现反射壁影响。

这是有限半径、有限背景压力的球形膨胀，不是可直接套用 Sedov 点爆炸解析解的测试。
检查离散反射/坐标置换对称性也不等于证明任意方向的球对称精度。

从 FLARE 根目录运行：

```bash
make -C tests/Spherical_Blast_3D -j4 test
make -C tests/Spherical_Blast_3D animate
```

`test` 编译、运行并检查初值、正密度/压力、质量/总气体能量守恒、
三方向反射和坐标置换对称性，以及实际向外膨胀。检查器要求等格数立方区域。
`check` 只分析已有输出。`run` 是框架原有的裸运行目标，不生成检查器使用的运行清单；
推荐使用 `test`，它会生成 `output/run.log`，让检查器和动画排除旧运行残留 CSV。
记录写在 `output/diagnostics/metrics.csv`。

## 独立 3D 动画脚本

根目录的 `animate3d.py` 不依赖现有 `animate.py`，适用于 FLARE 的三维 CSV。
依赖 numpy、matplotlib、Pillow；MP4 另需 ffmpeg。默认导出 GIF 和末帧 PNG。
左侧是真实坐标下的三张彩色正交切面，右侧是对应的 XY、XZ、YZ 二维切片；
这是剖切图，不是等值面或体渲染。

```bash
python3 animate3d.py tests/Spherical_Blast_3D/output mesh_blast3d --field p
python3 animate3d.py tests/Spherical_Blast_3D/output mesh_blast3d --field rho
python3 animate3d.py tests/Spherical_Blast_3D/output mesh_blast3d --field speed --stride 2
python3 animate3d.py tests/Spherical_Blast_3D/output mesh_blast3d --field p --output /tmp/blast3d.mp4
```

支持 `rho/p/vx/vy/vz/speed`，`--slice-x/y/z` 指定切片位置，默认选域中心附近的单元。
32 格时中心在两格之间，图上会标注实际采用的坐标。
`--vmin/--vmax` 指定色标范围；否则用所选帧全局范围，动画不会逐帧自动缩放。
`--stride` 对已保存文件抽帧，始末帧始终保留。播放是等帧时长，物理时刻显示在标题中。
默认读取输出目录的 `run.log` 清单；也可用 `--log` 指定日志。
无日志时按文件名发现快照并检查时间递增，无法排除所有同时间段旧数据，须使用独立目录。

更细的演示可把 `nxa/nya/nza` 都改为 48 或 64，并同时换 `output_directory`。
三维加密一倍会增加 8 倍单元数，时间步通常也变小，运行及输出开销明显上升。
`make animate` 使用默认目录和前缀；改输出配置后请直接调用上面的脚本并传入新目录/前缀。

其他值得后续实现的算例：激波撞密度球（团块压缩与剥离）、
Richtmyer–Meshkov（激波穿过起伏密度界面）、斜向 Sod（方向误差）、
周期小振幅声波（相位和耗散误差）。它们可复用当前流体物理，但需要各自的初始化。

## 立体表面动画与 Blender 模型

`render3d.py` 提取真实三角网格等值面，并进行光照着色及缓慢旋转拍摄。
默认固定 `p=1.12`，颜色表示表面局部密度；移除前方四分之一区域，露出内外表面。
切口保持开放，没有虚构的端盖。等压面表示阈值区域的边界，并非流体物质表面，
也不能直接当作精确的激波位置。粗网格的棱角仍然保留。

```bash
make -C tests/Spherical_Blast_3D render
# 或从 FLARE 根目录调用，--full 可显示完整外观：
python3 render3d.py tests/Spherical_Blast_3D/output mesh_blast3d --full --output-dir /tmp/blast-full
```

默认输出至 `output/animations/surface3d/`：

- `blast_surface.gif`、`blast_surface.mp4`（有 ffmpeg/libx264 时）、末帧 PNG。
- `surface_final.obj` 与 `surface_final_cutaway.obj`：可在 Blender 中导入的最终完整/剖开网格。
- `--obj-sequence` 还输出 `meshes/frame_*.obj`，对应每个原始输出时刻；
  `manifest.json` 记录文件、物理时间及步数，方便在 Blender 中组织网格切换动画。

默认把原始快照在物理时间上均匀插值为 72 帧，12 fps 播放。先线性插值压力和密度，
再提取等值面；插值仅用于显示，不是额外的流体求解。色标、阈值和取景范围固定。
`--frames`、`--fps` 可调整时长，`--preview` 只导出末帧和模型，便于快速检查。
无需 Blender 即可渲染；OBJ 是可导入的几何资产，不是包含材质/灯光的 `.blend` 场景。
表面使用 marching tetrahedra：对相邻单元中心构成的立方体作一致的六四面体剖分，
在边上插值得到交点。方法参考 [Paul Bourke](https://paulbourke.net/geometry/polygonise/)。
