# Sod 激波管

在问题目录中运行 Perl 脚本生成 Makefile，然后单独编译：

```bash
perl ../../build.pl
make -j4
./flare -par parameters.txt > run.log
python3 animate.py
```

`build.pl` 不接受参数，只在当前问题目录生成一个引用 FLARE 框架编译规则的
Makefile，然后退出。它不调用 Make 或编译器，也不生成 `.build/` 或 `flare`。
运行脚本前必须在问题目录准备 `user_config.h`。生成的 Makefile 固定使用这份
配置，每个源文件编译前都会通过完整路径加载它；缺失时脚本或 Make 会报错，
不会回退到框架或其他 testcase 的配置。若希望使用默认值，可以先把框架的
`src/user/user_config.h` 模板复制到问题目录。
首次及以后编译都使用 `make` 或 `make -j4`。已有的 Makefile 会保留，因此本
Sod 目录里的 `test` 和 `animate` 目标也会保留，直接运行 `make` 即可。
中间文件（`.o`、`.d` 和编译设置）集中保存在问题目录的隐藏子目录 `.build/` 中，
可执行文件仍为外面的 `./flare`。重复执行 `make` 会自动复用 `.build/`：没有变化时不重新
编译；源文件或头文件变化时，只重新编译受影响的对象，再链接可执行文件。
修改编译选项或 `user_config.h` 时会重新编译，以保持各模块的配置一致。
用户自己的问题目录**不需要手写 Makefile**：

```bash
cd /path/to/my_problem
perl /path/to/FLARE/build.pl  # 首次：只生成 Makefile
make -j4                    # 编译；后续修改源码后也用这条命令
./flare -par parameters.txt > run.log
```

生成的 Makefile 通过绝对路径引用 FLARE 框架 Makefile；如果移动了框架目录，
需要更新其中的 `include` 路径。脚本不会覆盖已有的 Makefile，也不会启动模拟。

问题目录需提供实现 `hd_user::make_user_hooks()` 的 `.cpp` 文件；顶层多个
`.cpp` 会一起编译。编译配置必须来自该问题目录的 `user_config.h`。
为方便编辑器解析，在用户 `.cpp` 的最前面、所有框架头文件之前写
`#include "user_config.h"`。带引号的包含会先查找当前源文件所在目录，使编辑器
解析各 testcase 的用户源码时也能先读取各自的配置。共用的框架 `.cpp` 仍需要
通过编辑器的编译配置选择问题目录；Makefile 编译时会为所有源码统一设置该目录。
运行时参数放在 `parameters.txt` 中。
生成 Makefile 后，可以运行 `make show-config` 检查选中的配置和源码。
生成步骤只需要 Perl；编译步骤需要 GNU Make 和 C++17 编译器。

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
