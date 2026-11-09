# 后向 Euler 积分器

求解内核位于 `src/time/implicit.cpp`，与 Euler、RK2 并列，
统一接入原 `TimeIntegrator` / `advance_one_step` 分派。
公共框架只增加配置字段、方法选择和工作区传递；`UserHooks`、`SourceTerm`、
`evaluate_RHS` 的返回类型与参数均不改变，显式积分、重构和 Riemann 算法不改变。
模块通过现有接口计算几何、重力、加热、用户源项和边界，不需要 dependency hook。

## 使用

全部设置写入参数文件，正常运行 `./flare -par parameters.txt`：

```ini
time_integrator = implicit
dt_implicit = 0.01
implicit_adaptive = true
```

`dt_implicit` 与 `t_end` 一样采用代码时间单位；只有时间单位为 1 秒时才等于物理秒。
MESA 算例的单位约为 1592.858 秒，所以 0.1 物理秒对应 `6.2780216922807006e-05`。
此前的 `--implicit-*` 命令行参数已迁移至参数文件。

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `dt_implicit` | 隐式模式必填 | 初始尝试步长及自适应上限；显式算例不需要设置 |
| `implicit_adaptive` | `true` | 非线性迭代失败后减半重试；成功且 Newton ≤ 8 次则下步增长 1.5 倍 |
| `implicit_tolerance` | `1e-6` | 最大局部缩放守恒残差容差 |
| `implicit_newton_max_iterations` | `60` | 每次尝试的 Newton 上限 |
| `implicit_report_file` | 空 | 每次成功/失败尝试、迭代数、残差和分段耗时的 CSV 路径 |
| `implicit_verbose` | `false` | 输出每次 Newton 的残差 |

`implicit_adaptive=false` 时为固定步长，只允许最后一步截到 `t_end`。
失败保持调用前的全部单元状态。自适应模式每步最多尝试 16 次；接受的实际 dt
返回主循环，`dt_next` 同步记录下一步建议。自适应只管理非线性收敛，没有时间
截断误差控制；用户源项回调异常及违反源项硬步长限制仍会报错退出。

主程序为隐式模式创建一个 `hd_implicit::Integrator` 工作区，传到
`advance_one_step` 新增的最后一个可选参数。隐式模式要求传入此工作区，以复用
染色、稀疏存储、报告文件和步长历史；原显式调用省略该参数即可，接口返回值不变。
工作区的设置在构造时从配置复制，一次运行中保持不变。

## 求解与优化

仍求解守恒方程 `R = U - U_old - dt*RHS(U, t_new) = 0`。
Newton 坐标为 `d(log rho), dv/c_ref, d(log p)`；行尺度分别为局部密度、
`rho*c_ref`、总能量，取新旧状态尺度的较大值。附加守恒分量保持线性更新。
低密度区没有统一的绝对 `1e-8` 扰动。增量采用 `expm1`，减少近零扰动的消减误差。
要求正密度、正压力，暂不支持真空。

- **矩阵自由的完整 Jacobian 作用**：`dU/dq` 解析计算，RHS 使用中心方向差分。
  每次重新填真实边界、计算全部源项，自引力与任意非局部用户依赖均包含其中。
  无须在公共接口新增 dependency hook；也不根据源项名称猜测耦合。
- **coloring 只组装近似预条件矩阵**：用现有流体算子和球坐标几何源项，固定
  ghost 状态，并省略物理源项导数。这些近似只影响线性迭代速度。
  最终残差与矩阵自由作用始终使用完整方程，不冻结实际引力/边界。
- 预条件矩阵为 CSR，按实际 PCM/PLM 局部 stencil 建图和贪心染色，适用于
  Cartesian 1/2/3D 与框架支持的球坐标 1D。拓扑和颜色跨时间步复用。
- ILU(0) + GMRES(30)，最多 120 次线性迭代，不存储稠密全局 Jacobian。
  GMRES 检查实际算子残差；ILU 失败时使用对角预条件。
- 线性相对容差为 0.01，避免早期 Newton 过度求解。线搜索采用固定行尺度的
  L2 残差下降，但最终接受仍要求每个分量的最大缩放残差满足用户容差。
- 候选网格和扰动网格保留容量；一步成功后才提交。已有 performance 仍可统计
  RHS/重构/Riemann 等耗时，模块自己的 CSV 记录预条件矩阵及线性求解时间。

局部预条件矩阵不含周期绕回、用户边界和非局部源项耦合。这些项如果主导刚性，
可能增加 Krylov 次数。PLM/MC 切换依然不光滑，不能保证任意大步长收敛。
源回调必须无外部副作用；现有源项报告的有限步长上限仍检查，不会静默忽略。
此模块不改变空间重构，也不自动换用 PCM、HLL 或修改网格。

## 文件职责

| 文件 | 职责 |
| --- | --- |
| `implicit.h` | 读取配置的工作区入口与隐藏实现的声明 |
| `implicit.cpp` | 残差、Newton、方向差分、染色预条件、自适应步长和报告 |
| `implicit_variables.h` | 对数增量、局部尺度及 `dU/dq`，仅模块内部使用 |
| `implicit_linear.h/.cpp` | 私有 CSR、ILU(0) 和 GMRES，不向公共框架添加类型 |

以上代码文件均直接位于 `src/time/`，没有单独的 `implicit/` 子目录。
`main.cpp` 只通过通用 `-par` 入口读取参数文件、创建所选积分器的工作区；
没有 `--implicit-*` 参数解析。文件整理不改变求解算法、容差或自适应规则。

回归测试：`make -C unit_test -j4 TESTS=implicit_standalone`。
覆盖解析非局部冷却、`1e-18` 密度、失败回滚、非局部用户边界、Cartesian
1/2/3D 和球坐标自引力，并独立重算已接受状态的守恒残差。
