# 多维 CFL 公式

本文记录 FLARE 的 Godunov 方法使用的多维 CFL 计算方式，并说明它与 LTS RK2 中的面波速、单元宽度和局部时间步的对应关系。多维部分针对当前的 Cartesian 对齐网格，允许各方向的网格间距不均匀。

## 1. 单个 cell 的 CFL 公式

对 active cell $i$，先定义第 $d$ 个方向上的代表波速：

$$
a_{i,d} = \max\left(a_{i,d}^{-},\;a_{i,d}^{+}\right).
$$

其中 $a_{i,d}^{-}$ 和 $a_{i,d}^{+}$ 分别是该 cell 在第 $d$ 个方向的低侧面、高侧面上，由黎曼求解器返回的非负波速尺度。

将各方向的“波速除以单元宽度”相加，得到该 cell 的 CFL 变化率：

$$
\boxed{
r_i = \sum_{d=0}^{D-1}
\frac{\max\left(a_{i,d}^{-},\;a_{i,d}^{+}\right)}{\Delta x_{i,d}}
}
$$

相应的局部 CFL 时间步上限为：

$$
\boxed{\Delta t_{i,\mathrm{CFL}} = \frac{C}{r_i}}
$$

等价地，一个时间步 $h_i$ 满足这里采用的 CFL 限制，当且仅当：

$$
h_i r_i \leq C.
$$

| 符号 | 含义 |
| --- | --- |
| $D$ | 模拟维数，对应 `hd_core::NDIM` |
| $i$ | active cell 的编号；在 LTS 中对应 active flat 索引 |
| $d$ | 方向编号，范围为 $0 \leq d < D$ |
| $a_{i,d}^{-}$、$a_{i,d}^{+}$ | 第 $d$ 个方向两侧面的波速尺度，单位为长度/时间 |
| $\Delta x_{i,d}$ | 该 cell 沿第 $d$ 个方向的宽度 |
| $r_i$ | 各方向贡献之和，单位为 $1/\text{时间}$，对应代码中的 `rate` 或 `inverse_dt` |
| $C$ | 无量纲 CFL 系数，初始选步时对应 `cfl_` |

若所有方向的波速都为零，则 $r_i=0$，流体波传播不对这个 cell 提供有限的时间步上限，代码将其局部上限设为正无穷。实际步长仍可受到源项、输出时刻或其他时间上限的限制。

## 2. 一维、二维和三维展开

以下将同一方向两个面的波速最大值简写为 $a_x$、$a_y$、$a_z$，省略 cell 下标。

**一维：**

$$
r = \frac{a_x}{\Delta x},
\qquad
\Delta t_{\mathrm{CFL}} = C\frac{\Delta x}{a_x}.
$$

**二维：**

$$
r = \frac{a_x}{\Delta x} + \frac{a_y}{\Delta y},
\qquad
\Delta t_{\mathrm{CFL}} =
\frac{C}{a_x/\Delta x + a_y/\Delta y}.
$$

**三维：**

$$
r = \frac{a_x}{\Delta x} + \frac{a_y}{\Delta y} + \frac{a_z}{\Delta z},
\qquad
\Delta t_{\mathrm{CFL}} =
\frac{C}{a_x/\Delta x + a_y/\Delta y + a_z/\Delta z}.
$$

如果所有方向的单元宽度都为 $\Delta x$，代表波速都为 $a>0$，则：

$$
r = D\frac{a}{\Delta x},
\qquad
\Delta t_{\mathrm{CFL}} = \frac{C\Delta x}{Da}.
$$

这里的维数因子已经包含在各方向求和中，计算后不需要再额外除以 `NDIM`。

## 3. 为什么各方向求和，而不是分别计算后取最小值？

当前 Godunov 空间算子在同一次更新中叠加所有方向的通量散度，所采用的 CFL 约束也累加各方向的传播率。

如果分别计算各方向允许的时间步，再取最小值，得到的是：

$$
\min_d\left(C\frac{\Delta x_{i,d}}{a_{i,d}}\right)
= \frac{C}{\max_d\left(a_{i,d}/\Delta x_{i,d}\right)}.
$$

而当前实现使用的是：

$$
\Delta t_{i,\mathrm{CFL}}
= \frac{C}{\sum_d a_{i,d}/\Delta x_{i,d}}.
$$

对于非零波速，后者不大于前者。只对各方向取最大传播率，会漏掉当前 CFL 定义中其他方向的贡献。

例如，某个二维 cell 的参数为：

| 方向 | 单元宽度 | 两侧面波速最大值 | 波速/宽度 |
| --- | --- | --- | --- |
| x | $2$ | $4$ | $2$ |
| y | $1$ | $3$ | $3$ |

取 $C=0.4$，则：

$$
r=2+3=5,
\qquad
\Delta t_{\mathrm{CFL}}=\frac{0.4}{5}=0.08.
$$

如果只取各方向时间步的最小值，则得到 $\min(0.2,\;0.4/3)\approx0.1333$，它对应的 $h r\approx0.6667$，已经超过这里设定的 $0.4$。

上述公式描述本项目采用的 CFL 规则；允许的 CFL 系数仍取决于空间离散、时间积分方法和实现中的限制。

## 4. 面波速与面速度的区别

LTS 中的 `speed_[face_id]` 是黎曼求解器返回的波速尺度，供 CFL 计算使用。`face_velocity_[face_id]` 是面上的流体速度，供压力功等计算使用。两者的物理含义不同。

对于理想气体 Euler 方程，单个状态沿第 $d$ 个方向的特征速度为 $u_d-c_s$、$u_d$、$u_d+c_s$，其最大绝对值为 $|u_d|+c_s$，其中：

$$
c_s=\sqrt{\frac{\gamma p}{\rho}}.
$$

这说明 CFL 需要计入流体运动和声波传播。实际 Godunov 面波速由重构后的左右状态及具体黎曼求解器确定，计算时直接使用求解器返回的 `speed_`。

## 5. 与 LTS RK2 代码的对应关系

对一个 active cell，先取得它包含 ghost offset 的存储坐标，再通过 `cell_faces` 找到各方向两侧面的独立 `face_id`：

```cpp
const auto index = mesh_.coordinates(storage_index(cell));

double rate = 0.0;
for (int d = 0; d < hd_core::NDIM; ++d)
{
    const auto& pair = mesh_.cell_faces[cell][d];

    const double direction_speed =
        std::max(speed_[pair.lower], speed_[pair.upper]);

    rate += direction_speed / mesh_.cell_width(d, index[d]);
}

local_dt_[cell] = rate > 0.0 ? cfl_ / rate : inf;
```

| 代码 | 对应含义 |
| --- | --- |
| `cell` | active flat 索引 |
| `storage_index(cell)` | 该 cell 在包含 ghost 的 `cells` 数组中的 storage flat 索引 |
| `index[d]` | 第 $d$ 个方向的存储坐标，已经包含 ghost offset |
| `pair.lower`、`pair.upper` | 第 $d$ 个方向两侧面的独立 `face_id` |
| `speed_[pair.lower]`、`speed_[pair.upper]` | $a_{i,d}^{-}$、$a_{i,d}^{+}$ |
| `mesh_.cell_width(d, index[d])` | $\Delta x_{i,d}$ |
| `local_dt_[cell]` | 初始状态给出的局部 CFL 时间步上限 |

`cell_width(d, index[d])` 使用该方向相邻两个面坐标之差：

```cpp
mesh_.face[d][index[d] + 1] - mesh_.face[d][index[d]]
```

因此这里需要传入多维坐标的第 `d` 个分量，不能直接把 active flat 索引当作该方向的坐标。面数组 `speed_` 则始终用独立 `face_id` 索引。

## 6. 全局时间步与 LTS 局部时间步

若所有 active cell 使用同一个时间步，则：

$$
r_{\max}=\max_{i\in\mathrm{active}}r_i,
\qquad
\Delta t_{\mathrm{global}}
=\min_i\Delta t_{i,\mathrm{CFL}}
=\frac{C}{r_{\max}}.
$$

这里先在**同一个 cell 内对方向求和**，再在**所有 active cell 之间取最大值**。`hydro::max_inverse_dt` 保存的就是这个最大传播率。

LTS 则根据每个 cell 的局部上限分配层级。设同步大步为 $H$，cell $i$ 的层级为 $\ell_i$，则它实际使用的局部步长为：

$$
h_i=\frac{H}{2^{\ell_i}},
\qquad
h_i\leq\Delta t_{i,\mathrm{CFL}}.
$$

层级缓冲可能进一步提高 $\ell_i$，从而缩短局部步长。设最大层级为 $L$，现有 LTS 的初始同步大步选择为：

$$
H=\min\left(
\Delta t_{\mathrm{limit}},\;
2^L\min_i\Delta t_{i,\mathrm{CFL}},\;
\Delta t_{\mathrm{source}}
\right).
$$

同步大步可以大于最小局部 CFL 上限，因为较细层会在一个 $H$ 内执行多个子步。

阶段中的波速可能增大，因此初始选步之后还需要阶段校验。多维阶段校验需要遵循相同的方向求和规则，并使用与被检查阶段相符的波速或明确的保守上界。LTS 每次可能只重新求解部分面，不能直接将缓存中不同阶段的面波速当作同一阶段的完整数据。

初始选步使用 `cfl_`，阶段校验阈值使用 `cfg_.cfl_max`；当前 LTS 输入检查要求 $0<C\leq C_{\max}\leq0.5$。

## 7. 对应源码

- `src/hydro/godunov.cpp`：普通 Godunov 中基于各方向面波速的 CFL 求和。
- `src/hydro/wave_timestep.h`：`hydro::max_inverse_dt` 的定义说明。
- `src/hydro/timestep_CFL.cpp`：由全局最大传播率计算时间步；另有基于 cell 状态估计波速的重载。
- `src/time/lts_rk2.cpp`：`run()` 中的初始局部时间步和同步大步选择，以及 `check_wave_speeds()` 中的阶段检查。

撰写本文时，LTS 的初始 `rate` 计算已改为多维求和；阶段校验、邻居分层等多维改造仍在继续，输入检查仍限制为 1D。
