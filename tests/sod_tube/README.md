# Sod 激波管

也可以直接从本问题目录运行 Perl 编译入口，不需要操作框架目录：

```bash
perl ../../build.pl -j4
./flare -par parameters.txt > run.log
python3 animate.py
```

`build.pl` 默认以运行命令时的当前目录作为问题目录，找到该目录的用户源码和
`user_config.h`，调用框架 Makefile，直接在该目录生成 `flare`。
中间文件以隐藏的 `.flare-*` 文件存放在同一目录，不创建 `build/` 子目录。
它只负责编译，重复执行会增量构建。用户自己的问题目录**不需要 Makefile**：

```bash
cd /path/to/my_problem
perl /path/to/FLARE/build.pl -j4
./flare -par parameters.txt > run.log
```

问题目录需提供实现 `hd_user::make_user_hooks()` 的 `.cpp` 文件；顶层多个
`.cpp` 会一起编译。可选的本地 `user_config.h` 优先于框架默认文件。
运行时参数放在 `parameters.txt` 中。
还可以从其他位置使用 `--problem /path/to/my_problem`，或加 `--show-config`
只检查选中的配置和源码。脚本使用系统 Perl、GNU Make 和 C++17 编译器。

在本目录运行：

```bash
make -j4 test
python3 -m pip install -r requirements.txt  # 仅首次缺少绘图依赖时需要
make animate
```

`make test` 编译并运行一维 Sod 算例，检查数值结果，保存 CSV 和本次运行的
`run.log`。`make animate` **只读取已有输出**，不会重新运行模拟。

生成的文件在 `output/animations/`：

- `sod_colormap.gif`：密度、x 方向速度、压力的二维色图演化。
- `sod_plot.gif`：这三个量的一维曲线演化。
- `sod_colormap_last.png`、`sod_plot_last.png`：末帧静态图。

色图是把一维数据沿纵向重复显示，纵向没有物理演化，不代表二维模拟。
所有图都排除 ghost 单元，色标及曲线坐标范围在整段动画中固定。
每个已保存快照对应一帧；由于时间步自适应，匀速播放快照不代表匀速物理时间。
CSV 开头有三行注释，分别记录已完成步数、当前物理时间和刚完成一步使用的步长。
初始快照三者均为 0；第 4 行开始是原来的 CSV 表格，不在每个单元重复记录时间。

```text
# step=0
# time=0
# dt=0
x,y,z,gamma,rho,vx,vy,vz,p,is_ghost
```

Pandas 可用 `pd.read_csv(path, comment="#")` 读取表格，也可以指定 `skiprows=3`。
动画通过 `snapshot_io.py` 单独读取注释头，使用完整精度时间，无需日志。
旧版 CSV 没有时间信息时自动尝试 `run.log`，两者都没有则只标注步数。
显式使用 `--log` 可只读取该日志列出的快照，避免混入旧运行残留的 CSV。

可以调整速度、抽帧或只生成一种动画：

```bash
python3 animate.py --fps 6 --kind plot
python3 animate.py --stride 2 --output output/my_animation
python3 animate.py --input /path/to/output --log /path/to/run.log
make animate ANIMATION_ARGS="--fps 3 --kind colormap"
```

当前 `parameters.txt` 每 10 步保存一次。如需更细的演化过程，可进一步减小
`step_every`，再运行 `make test` 和 `make animate`。改变播放速度不会增加数值快照。
脚本当前用于固定、均匀的一维 Sod 网格，默认匹配 `mesh_sod_*.csv`；
可以用 `--pattern` 修改匹配规则，或用 `--no-log` 明确读取全部匹配快照。
