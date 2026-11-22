# 多层守恒 RK2 局部时间步（`lts_rk2`）

[English version](lts_rk2_explanation_en.md)

本文配合 [一维参考代码](lts_rk2.cpp) 阅读，说明 LTS-RK2 如何分配时间层级、构造阶段状态、递归推进，并通过共享面通量积分保持跨层输运守恒。[当前多维实现](../src/time/lts_rk2.cpp) 沿用相同思路；涉及多维索引、邻域遍历和 CFL 检查的差异在相应位置说明。

阅读主线是：**选择宏步 → 分配单元与面的层级 → 准备起点阶段 → 递归推进细层 → 结算本层 → 同步或重试**。

## 1. 时间层级与基本记号

设宏步从物理时刻 $t_0$ 开始，长度为 $H$。层级 $\ell$ 的局部步长为

$$
h_\ell = \frac{H}{2^\ell}.
$$

层级越高，步长越短；下文将其称为“细层”，将步长较长的层级称为“粗层”或“慢层”。一次宏步内层级固定，每次重新尝试宏步时重新分配。

| 代码变量 | 含义 |
| --- | --- |
| `macro_step_` | 宏步长 $H$ |
| `depth` | 当前递归层级 $\ell$ |
| `fraction` | 当前局部步起点在宏步内的位置 $f$，对应时刻 $t_0+fH$ |
| `span` | 当前局部步占宏步的比例 $2^{-\ell}$ |
| `local_step` | 当前局部步长 $h_\ell=H\,\mathrm{span}$ |
| `end_fraction` | 局部步终点位置 $f+2^{-\ell}$ |
| `max_depth_` | 允许的最大层级编号，等于 `lts_max_levels - 1` |
| `deepest_` | 本次分层实际使用的最高层级 |

层级从 0 开始编号。因此，`lts_max_levels = 8` 表示允许层级 0–7，最大步长比为 $2^7=128:1$。只有 `max_depth_ = 8`（共 9 个层级）时，最大步长比才是 $256:1$。

本文用 $U_i$ 表示单元的守恒状态密度，用 $F_f$ 表示面上的数值通量。两者之间还需要面积、体积和时间积分，不能直接把“面通量”视为“单元的时间导数”。

## 2. 工作数组与索引

### 2.1 起点、试算状态与阶段视图

| 数组 | 保存的内容 |
| --- | --- |
| `start_` | 宏步起点的完整网格快照，重试时保持不变 |
| `state_` | 各单元最近一次完成局部 RK2 步后的状态；不同单元可能对应不同时刻 |
| `view_` | 为一次阶段计算构造的、对应同一物理时刻的近似网格状态；在 `build_stage_view()` 中复制并预测更新 |
| `base_[i]` | 单元 $i$ 当前局部步的起点状态 $U_i^0$ |
| `slope_[i]` | 该局部步的起点时间导数 $k_{1,i}$，用于线性预测 |
| `local_step_start_[i]` | 该局部步的起点分数位置 |
| `source_start_[i]` | 该局部步起点的源项 |
| `pressure_work_start_[i]` | 该局部步起点的内能压力功 |
| `pending_flux_integral_[i]` | 已累计给单元 $i$、尚未结算的净面通量积分 |

`state_` 保存的是状态，不是变化率；`slope_` 保存的是起点斜率，不是 RK2 两阶段平均后的斜率。线性预测只写入 `view_`，不会直接更新 `state_` 或正式网格。

每次宏步尝试开始时，`reset_attempt()` 执行 `state_.cells = start_.cells`。本层 Phase I 从 `state_` 取出局部步起点保存到 `base_`，并计算 `slope_`，此时不推进本层的 `state_`。递归进入细层后，细层会在自己的 Phase II 更新自己的 `state_`，所以不同层的状态可以对应不同时刻。

本层 Phase II 结算单元时，`auto& cell = state_.cells[flat]` 引用的就是试算数组中的单元；对 `cell.U` 的赋值将本层状态更新到局部步终点。整个宏步成功后才通过 `mesh_.cells.swap(state_.cells)` 提交到正式网格；若尝试失败，下次 `reset_attempt()` 会从 `start_` 重新开始。

`view_` 是给重构、Riemann 求解和源项计算使用的临时网格快照。`state_` 中不同层的状态可以处在不同时刻，`view_` 则把这些状态临时放到本次阶段的共同目标时刻。粗层尚未收齐细层的通量积分，不能提前把预测值当作已经完成的 RK2 更新写入 `state_`。

每次真正调用 `build_stage_view()`，都会先执行 `view_.cells = state_.cells`，复制整份 cell 数组，包括各单元当前的 `U` 和其他成员。随后才逐个判断 `needs_prediction`：不需要预测的单元保留复制来的 `U`；需要预测的单元则用 `base_ + slope_ * elapsed` 覆盖 `view_` 中的 `U`。起点只预测更粗层，终点预测当前层和更粗层，已经到达目标时刻的更细层不再预测。预测只发生在 `view_`，不会推进 `state_`。

面通量和源项由这个同步到目标时刻的 `view_` 计算。细层复用起点阶段时，跳过起点的 `build_stage_view()`，也就不会再次执行 `state_` 到 `view_` 的复制；它直接读取已经算好的 `flux_`、`speed_`、`face_velocity_` 和 `source_`，其中已经包含粗层在目标起点的预测状态。粗层的 `state_` 可以继续停留在自己的局部步起点。

### 2.2 面量与星区压力缓存

`flux_`、`speed_` 和 `face_velocity_` 分别保存面通量、CFL 波速和面速度。`first_flux_` 保存一个面在其所属层级当前局部步起点的通量，避免后续计算覆盖它。

`flux_` 是 `LTSRK2Step` 的成员数组，所有递归调用共用同一份，不会在进入下一层时重新创建，也不会自动随 `state_` 或 `view_` 的变化而更新。实际求新通量发生在 `solve_faces(stage, faces)` 内：对传入集合中的每个面重构左右状态，调用 Riemann 求解器，再将返回值赋给 `flux_[face]`。不在本次 `faces` 集合中的条目保持原值。

| 执行位置 | 对 `flux_` 的写入 | 结果对应的状态 |
| --- | --- | --- |
| `run()` 的 `solve_faces(start_, all_faces_)` | 更新全部独立面，再保存到 `initial_flux_` | 宏步共同起点的 `start_` |
| 每次尝试的 `reset_attempt()` | `flux_ = initial_flux_`，整份恢复，不重新求解 | 同一宏步起点；失败尝试中的覆盖不会保留 |
| 非复用步起点的 `solve_faces(view_, faces_needed_by_level_[depth])` | 更新所有 `owner >= depth` 的面 | 以当前更新态和粗层预测态构造的起点 `view_` |
| 本层终点的 `solve_faces(view_, faces_owned_by_level_[depth])` | 只更新 `owner == depth` 的面 | 本层及粗层预测、细层已完成更新的终点 `view_` |

`build_stage_view()` 只准备状态，不写 `flux_`；`evaluate_sources()` 写源项数组，也不写 `flux_`。`first_flux_[face] = flux_[face]` 则是从工作通量复制一份起点值到另一个数组，不会修改右侧的 `flux_`。因此要判断能否复用，应追踪这个面最近一次由哪个 `solve_faces()` 写入，而不是仅看目前 `state_` 或 `view_` 的时刻。

`hints_` 保存上一次黎曼求解得到的星区压力 $p_*$，作为下一次求解的迭代初猜；`0` 表示没有可用的压力初猜。多维访问方式为

```cpp
hints_[d][face_flat]
```

其中 `d` 是法向方向，`face_flat` 是该方向 active 面场的展平索引，包含切向位置。计算时先减去各方向的 ghost 偏移，再按面场尺寸展平；法向尺寸为 `nactive[d] + 1`，其余方向为对应的 `nactive`。

例如，二维 active 网格为 $2\times2$ 时，x 向面场尺寸为 $3\times2$。若 `(i, j)` 为减去 ghost 偏移后的面坐标，则

```cpp
face_flat = i + 3 * j;
```

| 面片 | 面坐标 | `face_flat` | 缓存位置 |
| --- | --- | ---: | --- |
| 下方一行的面片 | `(1, 0)` | 1 | `hints_[0][1]` |
| 上方一行的面片 | `(1, 1)` | 4 | `hints_[0][4]` |

两个面片虽然具有相同的 x 坐标，但连接不同的单元，因此使用不同的缓存位置。

几套面索引需要区分：

| 对象 | 索引布局 |
| --- | --- |
| `hints_[d][face_flat]` | 按方向分组的 active 面场 |
| 常规 Godunov 的 `flux[d].faces[face_flat]` | 与压力缓存对应的 active 面场布局 |
| 多维 LTS 的 `flux_[face_id]`、`speed_[face_id]`、`face_velocity_[face_id]` | 所有方向的独立面统一编号，对应 `mesh.face_info[face_id]` |
| `mesh.face_area[d].faces[flat]` | 含 ghost 区域的完整面场 |

一维有 $n$ 个 active 单元时，非周期边界有 $n+1$ 个独立面，周期边界有 $n$ 个独立面。旧版一维面数组仍保留 $n+1$ 个槽位，但周期端点代表同一个面，只独立求解和累计一次。多维版本直接用 `mesh.face_info.size()` 得到独立面数。

压力缓存仍保留周期边界的两个端点槽位，而 `face_info` 将它们合并为一个独立面。因此，LTS 的 `face_id` 不能直接当作 `hints_` 的 `face_flat`；求解时需要转换索引，并同步周期端点的压力缓存。

## 3. `run()`：选择宏步并准备一次尝试

### 3.1 初始扫描

首先在同步起点求解所有独立面，并计算源项：

```cpp
solve_faces(start_, all_faces_);
initial_flux_ = flux_;
initial_speed_ = speed_;
initial_velocity_ = face_velocity_;
initial_hints_ = hints_;
const double source_limit = evaluate_sources(start_, time_, 0.0);
initial_source_ = source_;
```

这些结果既用于选步，也供最初的 RK 阶段复用；`initial_*` 则用于失败后的恢复。

在一维中，单元 $i$ 的 CFL 变化率和允许步长为

$$
r_i=\frac{\max(a_{i-1/2},a_{i+1/2})}{\Delta x_i},
\qquad
\Delta t_i^{\mathrm{CFL}}=
\begin{cases}
\mathrm{CFL}/r_i, & r_i>0,\\
\infty, & r_i=0.
\end{cases}
$$

这里 $a_f$ 对应 `speed_[f]`。当前多维版本将各坐标方向的贡献相加：

$$
r_i=\sum_d\frac{\max(a_{i,d,-},a_{i,d,+})}{\Delta x_{i,d}}.
$$

遍历所有 active 单元后，得到 `local_dt_[i]` 及其最小值 $\Delta t_{\min}$。

### 3.2 选择宏步

```cpp
macro_step_ = std::min({
    dt_limit_, std::ldexp(minimum_local_dt, max_depth_), source_limit});
```

即

$$
H=\min\left(\mathrm{dt\_limit},\;2^{\ell_{\max}}\Delta t_{\min},\;H_{\mathrm{source}}\right).
$$

三项分别对应外部给定的步长上限、允许层数能够容纳的流体步长，以及源项步长限制。`dt_limit_` 可用于保证不越过指定的结束时刻。

### 3.3 重置、推进与重试

每次尝试先检查 $H$ 是否能够推进时间，再调用 `reset_attempt()`：重新分配单元和面的层级，将试算状态和阶段缓存恢复到同步起点，并清空待结算的面通量积分。随后调用

```cpp
advance_level(0, 0.0, true);
```

若某个阶段返回 `false`，失败沿递归调用逐层返回，整个宏步作废；外层将 $H$ 减半，从同步起点重新尝试。所有层成功后，填充终点边界，再将 `state_` 和压力缓存提交给正式网格及外部缓存。超出重试次数或遇到直接抛出的异常时，程序报错。

## 4. 分配单元层级与面所有权

### 4.1 `assign_time_levels()`：先满足 CFL，再建立缓冲区

对每个单元，从层级 0 开始，逐级缩短步长，直到满足

$$
H/2^\ell\leq\Delta t_i^{\mathrm{CFL}}.
$$

随后进行两次层级提升：

1. **扩展同层缓冲区。** 根据该单元的 CFL 层级 `depth`，调用 `promote_neighbors(cell, depth)`，将中心单元和缓冲范围内的邻居提升到至少该层级。
2. **平滑层级过渡。** 从最高层向下处理，对当前处于 `depth` 层的单元调用 `promote_neighbors(cell, depth - 1)`，使其周围至少达到下一较粗层级。

两次操作都通过 `std::max` 更新层级，只提升，不降低已有层级。最终相邻单元的层级差不超过 1；缓冲宽度较大时，还会约束更远的邻居。

缓冲宽度取决于重构模板：

| 重构方法 | `level_buffer_` | 原因 |
| --- | ---: | --- |
| PCM | 1 | 一个面的黎曼问题需要左右两个单元的状态 |
| PLM / PLM_corrected | 2 | 除面两侧单元外，计算其空间斜率还需要外侧邻居 |

例如，计算单元 `i` 的右侧面时，要重构 `i` 和 `i+1`；后者的斜率还需要 `i+2`。

### 4.2 多维邻域的遍历方式

多维 `promote_neighbors()` 沿每个坐标轴的低侧和高侧分别向外行走，最多走 `level_buffer_` 格：

```cpp
for (int d = 0; d < hd_core::NDIM; ++d)
    for (int side : {-1, 1})
    {
        int neighbor = cell;
        for (int distance = 0; distance < level_buffer_; ++distance)
        {
            const auto& pair = mesh_.cell_faces[neighbor][d];
            neighbor = side < 0
                ? left_cell(static_cast<int>(pair.lower))
                : right_cell(static_cast<int>(pair.upper));
            if (neighbor < 0) break;
            level_[neighbor] = std::max(level_[neighbor], depth);
        }
    }
```

每处理一个新方向，`neighbor` 都从中心单元重新开始。周期边界通过面拓扑自动连接首尾；非周期边界外没有 active 邻居，返回负值后停止。

以三维中心单元 `(i,j,k)`、缓冲宽度 2 为例：

| 方向 | 第一个邻居 | 第二个邻居 |
| --- | --- | --- |
| x 低侧 / 高侧 | `(i±1,j,k)` | `(i±2,j,k)` |
| y 低侧 / 高侧 | `(i,j±1,k)` | `(i,j±2,k)` |
| z 低侧 / 高侧 | `(i,j,k±1)` | `(i,j,k±2)` |

单次调用直接访问的是三条坐标轴上的邻域，不包含斜对角单元。远离边界且邻居不重复时，宽度 1 覆盖中心加 6 个邻居，宽度 2 覆盖中心加 12 个邻居；多次提升和平滑后，整体被提升区域可以更大。

### 4.3 `assign_faces_to_levels()`：由细层负责积分

每个面由两侧 active 单元中更高的层级负责：

$$
\ell_f=\max(\ell_{\mathrm{left}},\ell_{\mathrm{right}}).
$$

在非周期边界，缺失的一侧用负单元编号表示，面的层级取有效一侧的层级。

| 层级表 | `depth` 对应的内容 |
| --- | --- |
| `cells_by_level_[depth]` | 层级恰好为 `depth` 的单元 |
| `faces_owned_by_level_[depth]` | 负责层级恰好为 `depth` 的面 |
| `faces_needed_by_level_[depth]` | 负责层级大于或等于 `depth` 的面 |

三个表均有 `deepest_ + 1` 个层级槽位。`owned` 决定在哪一层进行面通量积分；`needed` 用于构造本层起点斜率，并为共享起点的细层准备通量。

`needed` 可以包含不直接接触本层单元的面，因为后续细层也会复用这些起点结果。反之，负责层级低于 `depth` 的面不可能接触 `depth` 层单元，因此本层不需要它们。

## 5. `build_stage_view()`：把不同层级放到同一时刻

这个函数构造临时阶段状态，供后面的重构、面通量、源项和终点压力功计算使用。它不执行本层的守恒结算，也不把 `view_` 整体写回 `state_`。

细层求解面通量时，邻近慢单元可能尚未完成自己的局部步，不能直接使用它在 `state_` 中的旧时刻状态。程序为它构造线性预测：

$$
\widetilde U_i(t)=U_i^0+(t-t_i)k_{1,i}.
$$

其中 $t_i$ 是该单元当前局部步的起点。用分数时间表示，就是

```cpp
const double elapsed =
    (fraction - local_step_start_[cell_index]) * macro_step_;
cell.U = base_[cell_index] + slope_[cell_index] * elapsed;
```

例如，设 A、B、C 分别属于层级 0、1、2，C 即将从 $3H/4$ 推进到 $H$。此时 A 的局部步仍从 0 开始，B 的局部步已从 $H/2$ 开始，则

$$
\widetilde U_A(3H/4)=U_A(0)+\frac{3H}{4}k_{1,A}(0),
$$

$$
\widetilde U_B(3H/4)=U_B(H/2)+\frac{H}{4}k_{1,B}(H/2).
$$

这里省略了共同的物理起始时刻 $t_0$；层级分布仅作示意，不考虑缓冲区可能引起的进一步提升。C 已完成前一个细步，直接使用其 $3H/4$ 时刻的 RK2 更新状态。

函数先无条件执行 `view_.cells = state_.cells`。此时复制的还是各单元最近完成局部步后的状态，时间不一定相同；只复制还没有完成阶段视图。接下来按以下规则决定哪些单元保留复制值，哪些需要覆盖成预测值：

| 阶段 | 更粗层 `level < depth` | 当前层 `level == depth` | 更细层 `level > depth` |
| --- | --- | --- | --- |
| 起点，`endpoint=false` | 线性预测到起点 | 使用已更新状态 | 使用已更新状态 |
| 终点，`endpoint=true` | 线性预测到终点 | 用起点斜率预测终点 | 使用已到达终点的更新状态 |

对应判断为

```cpp
const bool needs_prediction =
    level_[cell_index] < depth || (endpoint && level_[cell_index] == depth);
```

`needs_prediction` 判断的是复制之后要不要继续预测，不是要不要复制：

- 为 `false` 时执行 `continue`，`view_` 中该单元的 `U` 就保持从 `state_` 复制来的值。更细层总是属于这一种；构建起点视图时，当前层也属于这一种。
- 为 `true` 时，将 `view_` 中该单元的 `U` 替换为 `base_[i] + slope_[i] * elapsed`。`elapsed` 从该单元自己的 `local_step_start_[i]` 算到本次目标 `fraction`，不是一律加一个细层步长。

等待结算的粗层、以及构建终点视图时的当前层，其 `state_.U` 仍是自己的局部步起点，与 `base_` 对应，所以可以理解为给复制来的起点状态加上斜率预测增量。不过代码明确从 `base_` 计算并赋值，不是在上一次 `view_` 的预测值上继续累加；每次调用都会重新复制、重新预测。

为什么细层可以直接使用 `state_`？递归安排保证：进入本层局部步时，本层及所有更细层共享起点；本层构造终点阶段之前，更细层已经完成覆盖该时间区间的全部子步。

预测后检查状态的物理有效性，并按目标物理时刻填充边界。当前多维实现还会更新 `view_` 的热力学缓存，因此“不再预测”指保留复制来的 `U`，不表示跳过后面的边界填充和热力学更新。由此得到的 `view_` 是同一时刻的近似状态，既可用于面重构，也可用于依赖整张网格的源项，如自引力。PLM 必须在这个阶段视图上重新计算空间斜率，不能混用不同时刻的邻居状态。

两个布尔参数要分清：`advance_level(..., true)` 中的 `true` 表示复用起点阶段；`build_stage_view(..., true)` 中的 `true` 表示构建终点视图，需要把当前层也纳入预测。`build_stage_view()` 本身没有“是否复用”的参数。

## 6. `advance_level()`：一次局部 RK2 步

调用

```cpp
advance_level(depth, fraction, reuse_start_stage);
```

将当前层单元从 $t_0+fH$ 推进到 $t_0+(f+2^{-\ell})H$，并通过递归把所有更细层推进到相同终点。流程分为起点准备、细层递归和本层结算三部分。

### 6.1 第一阶段：准备起点通量、源项和预测斜率

如果起点阶段不能复用，就重新构造阶段视图、求解所需面的通量并计算源项：

```cpp
if (!reuse_start_stage)
{
    if (!build_stage_view(depth, fraction, false)) return false;
    solve_faces(view_, faces_needed_by_level_[depth]);
    if (macro_step_ > evaluate_sources(
            view_, time_ + fraction * macro_step_, macro_step_)) return false;
}
```

随后检查波速，并保存**本层负责的面**的起点通量：

```cpp
if (!check_wave_speeds(faces_needed_by_level_[depth])) return false;
for (int face : faces_owned_by_level_[depth])
    first_flux_[face] = flux_[face];
```

接着对本层单元保存 `base_`、`local_step_start_` 和 `source_start_`，计算起点斜率。在一维中，

$$
k_{1,i}=\frac{A_{i-1/2}F_{i-1/2}^{(0)}-A_{i+1/2}F_{i+1/2}^{(0)}}{V_i}
+S_i^{(0)}+Q_i^{(0)}.
$$

$Q_i$ 是仅在独立内能分量上非零的压力功项；代码使用

$$
Q_{i,e}=-(\gamma_i-1)e_i\,
\frac{A_{i+1/2}v_{i+1/2}-A_{i-1/2}v_{i-1/2}}{V_i}.
$$

这里 $e_i$ 是独立内能密度，$v_f$ 是黎曼求解得到的面速度。多维版本对各方向的通量差和面积加权面速度差分别求和。

所得 $k_{1,i}$ 写入 `slope_[i]`，用于本层的 Euler 终点预测和细层所需的慢区预测。程序还会检查 $U_i^0+h_\ell k_{1,i}$ 是否物理有效。此时尚未将本层单元更新到终点。

这些保存、斜率计算和波速检查在复用步中也会执行。循环只遍历 `cells_by_level_[depth]`，因此 `base_` 和起点压力读取的是已经到达 `fraction` 的本层 `state_`；相邻粗层的影响通过由 `view_` 算好的面量和源项进入。本层仍需保存自己的 `first_flux_`、`source_start_`，并计算自己的 `pressure_work_start_` 和 `slope_`。

### 6.2 起点复用：复用的是阶段计算

`reuse_start_stage=true` 表示所需面通量、波速、面速度和源项已经由同一个起点阶段计算得到，可以直接使用。此时跳过整个起点 `if (!reuse_start_stage)` 块：不调用 `build_stage_view()`，不执行 `view_.cells = state_.cells`，也不重新求面量或源项。复用的是这些计算结果，不是“再复制一次 `state_` 就得到同步状态”。

| 复用数组 | 本层的用途 |
| --- | --- |
| `flux_` | 保存本层负责面的起点通量，并计算本层单元的通量散度 |
| `speed_` | 检查本层所需面的 CFL 波速限制 |
| `face_velocity_` | 计算本层起点压力功中的速度散度 |
| `source_` | 保存本层 `source_start_`，参与起点斜率和最终两阶段积分 |

- 宏步第一次调用复用 `run()` 的初始扫描。
- 每层的第一个细半步与父层共享起点，复用父层准备好的结果。
- 第二个细半步从新的时刻开始，必须重新构造阶段视图并求值。

粗层的 `state_` 尚未推进到这个起点时，缓存仍然可用，因为准备缓存时已经先在 `view_` 中预测了粗层。复用跳过的是重新构造起点视图、重新求面量和源项；并不要求全部 `state_` 同步，也不会跳过本层的起点斜率计算。

例如宏步从 0 开始，有 0、1、2 三层。细层完成前半个宏步后，进入 `advance_level(1, 0.5, false)`，其起点视图按以下顺序构建：

| 层级 | `state_` 保存的时刻 | 刚复制到 `view_` 时 | `needs_prediction` | 起点视图构建完成后 |
| --- | --- | --- | --- | --- |
| 0 | $0$ | $U_0(0)$ | `true`，更粗层 | $U_0(0)+(H/2)k_{1,0}(0)$，预测到 $H/2$ |
| 1 | $H/2$ | $U_1(H/2)$ | `false`，起点同层 | 保留复制值，已经在 $H/2$ |
| 2 | $H/2$ | $U_2(H/2)$ | `false`，更细层 | 保留复制值，已经在 $H/2$ |

层级 1 用这份 `view_` 计算 `owner >= 1` 的面量和源项，结果保存到 `flux_`、`speed_`、`face_velocity_` 和 `source_`；0 层的 `state_` 仍在 0。随后调用 `advance_level(2, 0.5, true)` 时，不再构建起点视图，也不会把 0 层的旧值重新复制过来参与起点求值，而是直接使用刚算好的 $H/2$ 阶段结果。2 层自己的 `base_` 和起点压力从它已到达 $H/2$ 的 `state_` 读取，然后计算自己的起点斜率。

起点复用不会取消终点视图：2 层随后调用 `build_stage_view(2, 0.75, true)` 时，仍会重新复制 `state_`，再把 0 层从 0、1 层从 $H/2$、当前 2 层从 $H/2$ 分别预测到 $3H/4$。预测只改 `view_`，2 层完成 Phase II 结算后才把自己的 `state_` 更新到 $3H/4$。

从父层准备起点阶段到调用第一个细半步之间，父层只保存局部步起点并计算斜率，没有推进单元或覆盖这些阶段缓存。父层的 `needed` 面集合又包含子层所需的面，因此缓存的时刻、状态和面覆盖范围都满足复用要求。

具体到 `flux_`，这里没有“把父层数组传给子层”的额外动作：父子调用访问的是同一个成员数组。父层先把自己的 `owned` 面复制到 `first_flux_`，再读取 `flux_` 计算本层 `slope_`；这两步都不改写 `flux_`，父层的终点 `solve_faces()` 还要等两个细半步全部返回后才会执行。所以子层进入 `reuse=true` 分支时，它需要的起点通量仍是刚才准备好的值。父层本身若也是复用步，沿递归链向前追溯，最近的写入来自共同起点的某次非复用求解，或者宏步最初的扫描。

用第 7 节的四单元、三层例子跟踪具体条目：`C0、C1、C2、C3` 的层级为 `0、1、2、1`，`f2` 连接 `C1` 和 `C2`，归层级 2 负责。层级 1 的 `needed` 为 `f1`–`f4`，其中包含层级 2 需要的 `f2、f3`。

| 执行顺序 | 本次写入哪些通量 | `flux_[f2]` 的含义或用途 |
| --- | --- | --- |
| 宏步初始扫描，随后进入 `(0,0,true)`、`(1,0,true)`、`(2,0,true)` | 初始扫描写全部面；这三个调用的起点都不写 | 三层直接读取同一份 0 时刻的起点通量 |
| 细层走完前两个细步，层级 1 完成 $H/2$ 终点求解 | 层级 1 此时只写自己的 `f1、f4` | `f2` 仍是层级 2 上一细步在 $H/2$ 的预测终点通量 |
| `advance_level(1,0.5,false)` 的起点求解 | 重算并覆盖 `f1、f2、f3、f4` | `f2` 变成以 `C1、C2` 的 $H/2$ 更新态构造的起点通量 |
| 层级 1 保存 `first_flux_`、计算 `slope_` | 不写 `flux_` | `f2` 的值保持不变 |
| `advance_level(2,0.5,true)` 的起点 | 不写 `flux_`，将 `f2、f3` 的现值保存到 `first_flux_` | 直接读取父层刚写入的 $H/2$ 起点通量，并计算层级 2 的斜率 |
| 层级 2 的 $3H/4$ 终点求解 | 覆盖自己的 `f2、f3` | `f2` 变成预测终点通量；用它与保存的 `first_flux_[f2]` 做梯形积分 |
| `advance_level(2,0.75,false)` 的起点求解 | 再次覆盖 `f2、f3` | 使用 `C2` 的 RK2 更新态和慢层预测态，得到新的 $3H/4$ 起点通量 |

这里“复用步不更新”只指它自己的起点求解被跳过，不表示整个递归期间 `flux_` 都不变。后面的细层求解、本层终点求解仍会按面覆盖它。`first_flux_` 保存起点值，正是为了在工作数组被终点值覆盖后，仍能完成两阶段积分；各层只保存自己 `owned` 的面，子层不会覆盖父层负责面的起点副本。

此外，上一个子步的终点通量通常不能直接作为下一个子步的起点通量：前者使用本层的 Euler 预测状态，后者需要本层经过 RK2 修正后的状态。即使时刻相同，状态也不同。

<a id="recursion-feature-1"></a>

### 6.3 递归：先让细层走完两个半步

```cpp
if (depth < deepest_)
{
    if (!advance_level(depth + 1, fraction, true)) return false;
    if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
}
```

第一次调用把所有更细层推进到 $f+\mathrm{span}/2$；第二次调用再把它们推进到 $f+\mathrm{span}$。每个子调用内部也执行相同的两个半步，因此越细的层越早到达当前调用的终点。

第一次调用的 `true` 依赖父层已经准备好的起点阶段；较粗层即使尚未推进，其预测状态也已用于该阶段的面量和源项计算。第二次调用的起点已经改变，而且细层已经完成 RK2 更新，所以使用 `false` 重建起点视图和缓存。

这保证了本层进入第二阶段时，所有更细层都已到达 `end_fraction`，而本层仍保留局部步起点状态。

即使当前层没有单元，也必须先完成细层递归，之后才能通过

```cpp
if (cells_by_level_[depth].empty()) return true;
```

跳过本层结算。

### 6.4 第二阶段：构造终点视图，只求解本层负责的面

```cpp
if (!build_stage_view(depth, end_fraction, true)) return false;
solve_faces(view_, faces_owned_by_level_[depth]);
if (!check_wave_speeds(faces_owned_by_level_[depth]) ||
    macro_step_ > evaluate_sources(
        view_, time_ + end_fraction * macro_step_, macro_step_)) return false;
```

无论起点是否复用，到了本层的这个终点阶段，都会重新调用 `build_stage_view()`。它先从 `state_` 复制所有 cell，更细层已经到达终点，保留复制来的 `U`；当前层和更粗层则用各自的 `base_ + slope_ * elapsed` 覆盖 `view_` 中的 `U`。源项也在这个终点视图上重新计算。

这里的 `build_stage_view()` 只更新 `view_`，当前层的 `state_` 仍保存本局部步的起点。后面的面循环也只累积通量积分；真正写入本层 `state_` 的操作在 6.6 节的单元结算循环中。

此时只求解 `owned` 面，因为更细层负责的面已经在细层子步中完成了通量积分，并将贡献写入相邻单元的累积数组。本层无需再次积分这些面，否则会重复计入同一份输运。

`flux_` 是按面更新的工作缓存，并不始终代表整张网格在某一共同状态下的通量。起点 `needed` 集合的求解和复用规则保证当前斜率所需的通量有效；终点则只更新本层负责的面。

### 6.5 共享面通量积分

对本层负责的面 $f$，计算梯形积分

$$
I_f=A_f\frac{h_\ell}{2}\left(F_f^{(0)}+F_f^{(1)}\right).
$$

随后把同一份积分以相反符号写入两侧 active 单元：

```cpp
if (left >= 0) pending_flux_integral_[left] -= integral;
if (right >= 0) pending_flux_integral_[right] += integral;
```

正向通量从左流向右，所以左单元减少、右单元增加。内部面的两份贡献在全域体积加权求和时相消；边界面只向存在的 active 单元记账。

若一个面连接粗、细两个单元，细单元每完成一个子步就结算自己的累积量，粗单元则保留同一面贡献，直到自己的较长局部步结束。例如，粗步为 $H$、细步为 $H/2$ 时，粗单元收到

$$
I_f^{[0,H]}=I_f^{[0,H/2]}+I_f^{[H/2,H]}.
$$

这些积分直接保存在共享的 `pending_flux_integral_` 中，并不需要在函数返回时逐层复制。递归只决定何时产生贡献、何时结算各单元。

### 6.6 结算本层单元

本层单元的所有相邻面贡献到齐后，更新

$$
U_i^{\mathrm{new}}=U_i^0
+\frac{J_i}{V_i}
+\frac{h_\ell}{2}\left(S_i^{(0)}+S_i^{(1)}\right)
+\frac{h_\ell}{2}\left(Q_i^{(0)}+Q_i^{(1)}\right),
$$

其中 $J_i=\mathrm{pending\_flux\_integral\_}[i]$ 是带符号的净面通量积分。$J_i/V_i$ 已经是状态增量，不能再乘一次时间步；压力功只加到独立内能分量。

源项和压力功均使用起点、终点的两阶段平均。单层情况下，这对应普通 Heun/RK2 更新；跨层面则使用细层子步积分之和。

上式实际写入的是 `state_` 中的本层单元。当前多维代码中的对应赋值为：

```cpp
const auto flat = storage_index(cell_index);
const double volume = mesh_.volume_cell[flat];
auto& cell = state_.cells[flat];
cell.U = base_[cell_index]
    + pending_flux_integral_[cell_index] / volume
    + (source_start_[cell_index] + source_[flat]) * (0.5 * local_step);
cell.U[hd_core::INTERNAL_ENERGY] +=
    0.5 * local_step * (pressure_work_start_[cell_index] + pressure_work_end);
```

其中 `pressure_work_end` 在累加前由终点 `view_` 的压力和终点面速度计算。由于 `cell` 是引用，赋值直接修改 `state_`；这不是把预测的 `view_` 整体复制回来。只遍历 `cells_by_level_[depth]`，因此当前层推进到 `end_fraction`，细层保留已经完成的更新，更粗层继续等待自己的 Phase II。

更新后检查物理有效性、执行双能量同步，并仅清空**本层单元**的 `pending_flux_integral_`。更粗层单元的累积量继续保留，等待它们自己的局部步结束。

共享面积分保证的是输运部分的离散守恒。源项可以改变全域守恒量，双能量同步也可能修正总能量；两者应与内部面的守恒抵消区分。

## 流程总结

0. 外部run函数各种各种初始化,包括初始化数组,时间步的分配,cell和face的层级分配,计算出最大步长macro_step(由所有cell中最小的可行dt * $2^{\mathrm{max_depth}}$),然后开始尝试循环。
1. run先执行solve_faces(start_,all_faces_),把所有独立面的起点通量写入flux_,再保存到initial_flux_。每次尝试在reset_attempt中用state_.cells = start_.cells恢复状态,用flux_ = initial_flux_恢复起点通量,然后呼叫advance_level(0,0,true)。因此最外层起点直接读取的是本次宏步已经求好并恢复的通量,不会使用失败尝试留下的终点值。
2. 每一层的advance_level有层数参数depth,起点参数fraction(分数步长)和是否复用reuse. 计算当前的分数步长span = $2^{-\mathrm{depth}}$和终点end_fraction = fraction + span.先检测是否为复用步,如果不是,调用 `build_stage_view(depth, fraction, false)`，在函数内复制RK2推进的cell状态结果state_到临时数组view_中,由于步长不一,每个cell在state_里的状态对应的时间是不一致的。复制完再逐个判断needs_prediction,起点时只有更低层为true,就用各cell自己的base_ + slope_ * ((fraction-local_step_start_)*macro_step_)覆盖view_中对应的U,预测到fraction处；不是在旧view_的预测值上继续加。由于递归特性,当前层和更高层级的cell必然已经推进到当前起点fraction处,因此不用修改从state_处复制的值,于是整个view_同步。然后用view_里的cell状态去求一次黎曼问题,求出所有owner为当前层级和更高层级的面上的通量flux_,因为只有当前/更高层级owner的面才可能是当前层级的cell的所属面,通量才会对这些cell有影响。

   如果是复用步,意味着这个起点所需的flux_、speed_、face_velocity_和source_已经准备好,所以跳过build_stage_view、solve_faces和evaluate_sources,连这次起点的state_到view_复制也不执行。更低层级的state_虽然可能还在更早的时刻,但它们在当前起点的预测状态已经通过之前的view_参与了这些计算。父层的needed求解已经覆盖子层所需的面,保存first_flux_和计算slope_又只读取flux_,所以进入第一个细半步时,需要的通量没有被改写。因此不用为了本层再预测和求解一遍,后面仍然要保存本层自己的起点并计算slope_。
3. 保存起点处由本层负责的各个面的通量flux_到first_flux_下。起点处各相邻面的通量乘上面积、按流入减流出求和后除以cell体积,得到不含源项和压力功的状态变化率；再加上源项和内能压力功,得到RK2第一部分的完整斜率$k_1$。最终积分时,源项和压力功都采用起点和终点值的平均。

 然后遍历本层的所有cell,从state_保存当前局部步的起点到base_,从source_保存起点源项到source_start_,再用已有的flux_、face_velocity_和本层起点压力计算pressure_work_start_和slope_。这些在复用步中也要执行。因为只遍历cells_by_level_[depth],所以读取的state_和压力都属于已经到达fraction的本层cell,不会把更低层级尚未推进的状态当作本层起点。Phase I只保存起点和斜率,本层state_要等自己的Phase II才更新。
4. 两次推进更高一层的递归,先执行第一次advance_level(depth+1,fraction,true),这是让更高一层从fraction走到fraction+span/2,因为父层已经准备好这个共同起点的面量和源项,且父层的needed集合包含细层所需的面,所以复用为true；更低层级的预测状态也已经包含在这些结果里。注意这是递归,这里会由下一层advance_level(depth+1)跳转到步骤2继续执行这个流程,直到最高一层计算完毕,再逐层返回到这一层；没有cell的层也要先完成更高层的递归,之后才跳过本层结算。因此,当执行完这一次advance_level返回后,所有高于depth的层级都推进到了fraction + span/2。

接着这一层继续执行第二次advance_level(depth+1,fraction + span/2,false),这是让更高一层从fraction + span/2 走到end_fraction,由于这次推进起点是本层的起点向后走半个本层的步长,因此本层无法与这一次推进共起点,所以设置reuse=false。即使前一个细步计算过这个时刻的终点通量,那也是用预测终点状态得到的；现在细层已经完成RK2更新,需要用更新后的细层状态和重新预测的慢层状态构造新的起点view_。同理是递归呼叫,更高一层跳回步骤2继续执行流程，因此只有再所有更高层全部推进到end_fraction后,才会执行本层接下来的所有逻辑。

注意以上过程没有讨论不合法的情形,如果出现合法性检查不通过的层,直接逐级return false,回到run函数里的尝试(attempt)循环中,将最大步长macro_step缩小一半再次尝试。

5.  接下来进入RK2的Phase II。要计算终点处的cell状态对时间的变化率(终点斜率),我们首先通过build_stage_view预测终点end_fraction处的cell状态view_。将state_再次复制到view_中,由于上述提到的递归特性,更高层级的cell已经推进到了终点处,所以更高层级的cell不用修改。而当前层级还没有推进到终点,所以当前层级以及更低层级的needs_prediction为true,用各自的base_ + slope_ * elapsed覆盖view_中的U,外推到同一终点；更高层级为false,只保留复制来的U。这一步无论起点reuse是否为true都要执行,build_stage_view的第三个true是endpoint参数。 此处只更新 `view_`，本层 `state_` 仍保持局部步起点，直到第8步结算时才更新。

6. 用view_中的cell状态计算solve_faces,这一次只计算owner为当前层级的face的终点通量,并覆盖共享flux_里对应的条目；后面的梯形积分用first_flux_保存的起点值和这里新写入的终点值。终点值基于预测状态,下一次非复用起点仍要用更新后的状态重新求解,不能直接当作下一步的起点值。owner更低的面不可能是本层cell的相邻面；owner更高的面已经在细层子步中完成积分,其贡献直接累计在共享的pending_flux_integral_里,等待对应cell结算。这里使用细层已累计的通量积分,和Phase I的reuse_start_stage复用起点面量、源项是两个不同的步骤。

7. 计算owner为这一层的面的通量积分,将起点first_flux_和终点flux_取平均,乘上面积和本层局部时间步,再把其负值加到pending_flux_integral_[left],正值加到pending_flux_integral_[right]。这些积分保存在共享数组里,直到对应cell的层级结算时,除以cell体积后作为状态增量使用,并不是随递归返回逐层复制。

8. 积分推进步，遍历 `cells_by_level_[depth]`，通过 `auto& cell = state_.cells[flat]` 引用本层cell。执行 `cell.U = ...`，从 `base_[cell_index]` 出发，加上 `pending_flux_integral_[cell_index]` 除以体积以及两阶段源项积分，再给独立内能加上两阶段压力功。这些赋值直接把本层 `state_` 推进到 `end_fraction`；更高层级保留递归中已经完成的状态，更低层级继续等待自己的 Phase II。最后清理掉pending_flux_integral_中属于这一层cell的通量。这是由于后续更低层的积分推进中,已经用不到了(一个cell不可能属于两个层级)，这样所有层级完成后这个数组也全部清零。在中间步的计算中, 如果depth层级的某一次advance_level完成了的话,那么depth层和更高层的cell对应的pending_flux_integral_都会清零，计算更低层(如 depth-1 层)时也不会影响到更低层的cell。
举例:   f1  C1  f2 C2  f3 C3 f4
C1属于第1层,C2属于第0层,C3属于第二层，f1的owner_level维1,f2的owner_level也为1,f3和f4的owner_level为2。在层级2的时候会计算owner level为2的面f3和f4的通量,并把累积通量积分加给相邻cell C2(f3)和C3(f3和f4)的累积通量积分,然后层级2的推进步会令这一层的cell C3用得到的f3和f4给的累积通量积分推进状态，然后清零C3得到的累积通量积分。而C2属于层级0，因此f3给C2的通量不会清零。然后退到层级1，会把面f1和f2的累积通量积分加给C1(f1)和C2(f1,f2)的累积通量积分,然后层级1再令C1用从f1得到的累积通量积分推进状态，然后清零C1得到的累积通量积分，但C2得到的通量仍然不会清零,这时C2得到的累积通量积分已经是面f2和(负的)面f3之和。最后退回层级0,没有属于层级0的面,然后推进结算层级0的cell C2，这时就能结算在高层级逐层传递下来的累积通量积分pending_flux_integral_[C2],把相邻面f2和f3贡献的通量积分之和加载C2的cell状态中，这样尽管C2的相邻面属于更高层级，我们依然正确结算了C2的所有相邻面的贡献,并且每个面的贡献只需计算一次,不重复计算浪费性能。最后全部清零。

9. 返回真值，完成这一步的推进。递归回退到上一级后会逐层完成结算，除非合法性检查不通过，那所有的层级会快速返回false，退回run函数中的attempt循环，尝试0.5倍的宏步。 当最外层 `advance_level(0, 0, true)` 成功返回、所有层都走完宏步后，`run()` 填充 `state_` 的终点边界，再执行 `mesh_.cells.swap(state_.cells)` 提交到正式网格；否则下一次 `reset_attempt()` 重新从 `start_` 开始。


## 7. 四个单元、三个层级的完整示例

设四个单元采用非周期边界，层级已经确定为 `[0, 1, 2, 1]`：

```text
边界 | f0 | C0 | f1 | C1 | f2 | C2 | f3 | C3 | f4 | 边界
层级：      0         1         2         1
步长：      H        H/2       H/4       H/2
```

这个层级分布用于展示递归顺序；实际 `assign_time_levels()` 可能因缓冲区要求继续提升邻居。

<a id="face-ownership"></a>

### 7.1 面的负责层级

| 面 | 相邻 active 单元 | 负责层级 |
| --- | --- | ---: |
| `f0` | `C0` | 0 |
| `f1` | `C0`、`C1` | 1 |
| `f2` | `C1`、`C2` | 2 |
| `f3` | `C2`、`C3` | 2 |
| `f4` | `C3` | 1 |

因此，各层的面集合为

| 层级 | `faces_owned_by_level_` | `faces_needed_by_level_` |
| --- | --- | --- |
| 0 | `{f0}` | `{f0, f1, f2, f3, f4}` |
| 1 | `{f1, f4}` | `{f1, f2, f3, f4}` |
| 2 | `{f2, f3}` | `{f2, f3}` |

<a id="initial-stage"></a>
<a id="level-0-start"></a>
<a id="level-1-first-step"></a>

### 7.2 在共同起点准备各层

为简洁起见，将宏步起点记为 0。初始扫描计算所有面的通量、波速、面速度，以及所有单元的源项。随后依次进入

```text
advance_level(0, 0.0, true)
└─ advance_level(1, 0.0, true)
   └─ advance_level(2, 0.0, true)
```

三层复用同一份起点阶段结果，但分别保存自己负责的面通量和单元预测数据：

| 层级 | 保存的起点面通量 | 保存起点状态与斜率的单元 |
| --- | --- | --- |
| 0 | `f0` | `C0` |
| 1 | `f1`、`f4` | `C1`、`C3` |
| 2 | `f2`、`f3` | `C2` |

例如，`C1` 的起点斜率由左面 `f1`、右面 `f2` 的通量差除以体积，再加上源项和内能压力功得到。虽然 `f2` 归层级 2 负责积分，层级 1 仍需要它的起点通量来构造 `C1` 的预测斜率。

<a id="level-2-first-step"></a>
<a id="level-2-second-step"></a>

### 7.3 `C2` 先走两个细步：$0\rightarrow H/4\rightarrow H/2$

第一个细步没有更细层可递归，直接构造 $H/4$ 的终点视图。此时四个单元都用各自的起点斜率预测到 $H/4$，求解 `f2`、`f3`，并计算

$$
I_f^{[0,H/4]}=A_f\frac{H}{8}\left(F_f^{(0)}+F_f^{(1)}\right),
\qquad f\in\{f2,f3\}.
$$

其中阶段上标表示当前细步的起点和预测终点。累积符号为

| 面 | 左单元贡献 | 右单元贡献 |
| --- | --- | --- |
| `f2` | `C1`：$-I_{f2}$ | `C2`：$+I_{f2}$ |
| `f3` | `C2`：$-I_{f3}$ | `C3`：$+I_{f3}$ |

`C2` 用自己的净积分、两阶段源项和压力功完成更新，并清空自己的累积量。`C1`、`C3` 尚未结束局部步，保留收到的贡献。

第二个细步调用

```cpp
advance_level(2, 0.25, false);
```

在 $H/4$ 的起点视图中，`C2` 使用刚完成的 RK2 更新状态，其余单元仍使用慢区预测。因此必须重新计算 `f2`、`f3` 的起点通量，不能复用上一步基于预测终点得到的通量。

随后重复预测、求解、积分和更新，将 `C2` 推进到 $H/2$。至此，`C1`、`C3` 已收齐来自细层面的两个 $H/4$ 积分，但尚未更新到 $H/2$。

<a id="level-1-first-commit"></a>

### 7.4 结算层级 1：`C1`、`C3` 到达 $H/2$

构造层级 1 的终点视图时，`C0`、`C1`、`C3` 预测到 $H/2$，`C2` 直接使用已更新的 $H/2$ 状态。接着求解本层负责的 `f1`、`f4`，计算它们在 $[0,H/2]$ 上的梯形积分。

记 $I_f^{[a,b]}$ 为该时间区间内已完成的面积分，则待结算的净积分为

$$
J_{C1}=I_{f1}^{[0,H/2]}
-I_{f2}^{[0,H/4]}-I_{f2}^{[H/4,H/2]},
$$

$$
J_{C3}=I_{f3}^{[0,H/4]}+I_{f3}^{[H/4,H/2]}
-I_{f4}^{[0,H/2]}.
$$

加上各自的源项和压力功后，`C1`、`C3` 更新到 $H/2$，并清空自己的累积量。`f1` 给 `C0` 的贡献 $-I_{f1}^{[0,H/2]}$ 继续保留。

<a id="level-1-second-step"></a>

### 7.5 后半个宏步：$H/2\rightarrow H$

调用 `advance_level(1, 0.5, false)` 时，`C1`、`C2`、`C3` 已同步到 $H/2$，只有 `C0` 需要从 0 预测到 $H/2$。层级 1 重新求解所需面 `f1`–`f4`，把新起点通量逐面写入共享的 `flux_`，覆盖此前留在数组里的终点阶段值。随后保存 `first_flux_[f1]、first_flux_[f4]` 并计算 `C1`、`C3` 的新起点斜率，这些操作只读取 `flux_`，不会改写刚算好的 `flux_[f2]、flux_[f3]`，所以它们可以直接供第一个细半步使用。

随后执行

```cpp
advance_level(2, 0.5, true);   // H/2 -> 3H/4
advance_level(2, 0.75, false); // 3H/4 -> H
```

第一个细步复用层级 1 刚准备的 `flux_`、`speed_`、`face_velocity_` 和 `source_`。此时 `C0` 的 `state_` 仍在 0，但上述缓存已经用包含 `C0` 在 $H/2$ 的预测状态的 `view_` 计算过。这里不会再次调用起点 `build_stage_view()`，连 `state_` 到 `view_` 的复制也跳过。`C2` 从自己的 $H/2$ 更新态保存 `base_` 并计算起点斜率，无需再次预测 `C0` 或重算起点面量和源项；到 $3H/4$ 的终点阶段仍须重新构建视图。

第二个细步从 $3H/4$ 开始，使用已更新的 `C2`，并在新 `view_` 中预测 `C0`、`C1` 和 `C3`，重新计算起点阶段。`C2` 到达 $H$ 后，层级 1 再构造终点视图，完成 `f1`、`f4` 的后半步积分，并将 `C1`、`C3` 更新到 $H$。

<a id="level-0-commit"></a>

### 7.6 最后结算 `C0`

此时 `C1`、`C2`、`C3` 已到达 $H$，`C0` 的试算状态仍在 0，但已收到 `f1` 在两个半步中的贡献。

层级 0 将 `C0` 预测到 $H$，其余单元直接使用已更新状态；求解 `f0` 的终点通量并完成其宏步积分。`C0` 的净积分为

$$
J_{C0}=I_{f0}^{[0,H]}
-I_{f1}^{[0,H/2]}-I_{f1}^{[H/2,H]}.
$$

加入 `C0` 的两阶段源项和压力功后，四个单元全部同步到 $H$。最终填充边界，并提交整个宏步。

<a id="execution-timeline"></a>

### 7.7 递归执行顺序

```text
准备 C0 在 0 时刻的斜率
│
├─ 准备 C1、C3 在 0 时刻的斜率
│  ├─ C2：0 -> H/4
│  └─ C2：H/4 -> H/2
├─ 结算 C1、C3 到 H/2
│
├─ 准备 C1、C3 在 H/2 时刻的新斜率
│  ├─ C2：H/2 -> 3H/4
│  └─ C2：3H/4 -> H
├─ 结算 C1、C3 到 H
│
└─ 结算 C0 到 H
```

## 8. 源项求值与步长检查

### 8.1 `evaluate_sources()`：分别检查，再检查合计

该函数读取阶段网格 `stage`，将各个源项的贡献累加到 `source_`，不会直接把源项加到 `stage` 的状态上。

每次调用都会清空源项工作区并使尺度缓存失效，然后依次：求值一个源项、补齐所需的内能源项、检查回调给出的步长上限、检查该源项的状态变化限制，再将其累加。全部求值后，还要对源项总和检查一次。

分别检查可以避免相抵源项掩盖过大的单项变化；合计检查可以发现多个同向源项叠加后超限的情况。

### 8.2 允许变化量与阶段内缓存

`scan_source_changes()` 使用硬编码的 10% 系数，为每个单元计算

$$
\mathcal A_\rho=0.1\rho,\qquad
\mathcal A_{m_d}=0.1\rho(|v_d|+c),\qquad
\mathcal A_e=0.1\frac{p}{\gamma-1},
\qquad c=\sqrt{\frac{\gamma p}{\rho}}.
$$

这些量分别限制密度、各方向动量密度和独立内能密度的源项变化。动量尺度包含声速，因此不是简单地限制为当前动量的 10%。

对于某个源项分量的变化率 $R$，检查条件为

$$
H|R|\leq\mathcal A.
$$

这里约束的是整个宏步 $H$，即使当前正在执行细层子步。它是基于阶段变化率的步长限制，并非对完整非线性演化变化量的严格保证。

同一次 `evaluate_sources()` 中，第一个源项填充密度和允许变化量缓存，后续源项及合计检查复用这些尺度。下一次阶段求值会重新填充，不能因为仍使用同一个 `view_` 对象就沿用旧尺度。

所有源项分量都要检查有限性；真空单元允许零源项，当前实现不支持通过非零源项从真空创造气体。

### 8.3 初始选步与后续试步采用不同返回方式

| 调用情形 | `source_change_limit()` 的行为 |
| --- | --- |
| 初始扫描，`trial_H = 0` | 对所有非零变化率取 $\min(\mathcal A/\lvert R\rvert)$，返回实际步长上限 |
| 后续阶段，`trial_H > 0` | 检查当前 $H$ 是否超限：通过返回 $\infty$，超限返回 0 |

初始上限参与宏步选择。后续阶段已经知道正在尝试的 $H$，只需要判断是否可接受，无需再计算精确的最小上限。源项回调自己给出的上限仍独立参与 `evaluate_sources()` 的最终限制。

后续阶段使用

```cpp
scan_source_changes(stage, rates, [&](double allowance, double rate) {
    exceeded |= source_change_exceeds(trial_H, twice_H, allowance, rate);
});
```

只要一个分量超限，`exceeded` 就保持为真。不过扫描仍继续，以免遗漏其他单元或源项的非法值。返回 0 后，上层拒绝本次宏步，并尝试减半后的 $H$。

### 8.4 避免不必要除法的快速判断

```cpp
bool source_change_exceeds(double H, double twice_H, double allowance, double rate)
{
    const double magnitude = std::abs(rate);
    if (twice_H * magnitude < allowance) return false;
    return rate != 0.0 && H > allowance / magnitude;
}
```

若 $2H|R|<\mathcal A$，当前步长显然没有超限，可以直接接受。否则回到原本的除法比较 $H>\mathcal A/|R|$；零变化率不限制步长。倍数 2 为快速判断留出余量，不改变 10% 的实际限制。

## 9. 从一维代码读到多维实现

原文的一维访问方式

```cpp
auto& cell = view_.cells[ng_ + cell_index];
```

在当前多维代码中对应

```cpp
const auto flat = storage_index(cell_index);
auto& cell = view_.cells[flat];
```

其中

```cpp
std::size_t LTSRK2Step::storage_index(int cell_id) const
{
    return mesh_.active_to_storage[cell_id];
}
```

`cell_index` 是 active 单元编号，`flat` 是包含 ghost 区域的完整存储编号。`base_`、`slope_` 和通量累积数组按 active 编号访问；`MeshBlock::cells`、`source_` 和单元体积按存储编号访问。多维中不能用一个固定 ghost 偏移代替这层映射。

其他主要对应关系为：

| 一维参考实现 | 当前多维实现 |
| --- | --- |
| 面 `i` 和 `i+1` 是单元 `i` 的左右面 | `mesh.cell_faces[i][d]` 给出方向 `d` 的低侧面和高侧面 |
| 使用一对面的通量差 | 对各方向的面积加权通量差求和 |
| 使用一对面速度计算压力功 | 对各方向的面积加权面速度差求和 |
| 用左右面波速和单元宽度计算 CFL 变化率 | 对各方向贡献求和；动态检查使用本次尝试中记录的面波速峰值 |
| 通过一维邻接关系扩展缓冲区 | 沿面拓扑在各坐标方向遍历邻居 |

波速检查使用相邻单元各自完整的局部步长，不能仅用当前面的细步长检查粗单元。若动态 CFL、源项限制或状态有效性检查失败，仍按同一套宏步回退逻辑处理。

<a id="core-principles"></a>

阅读代码时，把三个操作分开即可贯穿整个流程：**预测**为阶段计算提供同一时刻的状态；**面积分**为相邻单元提供共享的输运贡献；**局部结算**在贡献到齐后更新该层单元。递归安排它们的先后顺序，使所有层最终重新同步到宏步终点。
