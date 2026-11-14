# 多层守恒 RK2 局部时间步（`lts_rk2`）
这份算法解释主要配合docs文件夹的lts_rk2.cpp使用,这是一个1D的旧版本,与源代码文件夹的3D新版本有一些细节差异，但算法思路一致。
## 类实例化


independent_face_count_: 若第一维是周期边界,则独立面数为n_+1，否则为n_(第0个active_cell)
### 与面有关的hints_和flux_

`hints_` 保存各个单元面上一次求解得到的**星区压力 \(p_*\)**，用于下一次求解该面黎曼问题时的迭代初猜，以加快收敛。缓存值为 `0` 表示没有可用的压力猜测。

它按法向方向分组，通过以下方式访问：

```cpp
hints_[d][face_flat]
```

- `d`：面的法向方向，`0、1、2` 分别对应 x、y、z。
- `face_flat`：该方向面场内部的展平索引，由面的完整多维坐标计算，包含切向位置的信息，**不是仅沿法向方向的坐标索引**。

计算 `face_flat` 时，先从面存储坐标中减去各方向的 ghost 偏移，再按 active 面场的尺寸展平。该面场的各方向尺寸为 `nactive`，其中法向方向多一个面位置。

例如，二维网格包含 `2×2` 个 active cell。对于 x 向面，面场尺寸为：

```text
(nactive[0] + 1) × nactive[1] = 3 × 2
```

令 `(i,j)` 为减去 ghost 偏移后的面坐标，则：

```cpp
face_flat = i + 3 * j;
```

同在 `x = x₁` 这条分界线上的两个面片，对应不同的缓存位置：

| 面片 | 面坐标 `(i,j)` | `face_flat` | 缓存位置 |
| --- | --- | --- | --- |
| 下方一行的面片 | `(1,0)` | `1` | `hints_[0][1]` |
| 上方一行的面片 | `(1,1)` | `4` | `hints_[0][4]` |

虽然两个面片的 x 坐标相同，但它们位于不同的行，连接不同的 cell，因此压力猜测分别存储。

需要区分以下几套索引：

| 数组 | 索引含义 |
| --- | --- |
| `hints_[d][face_flat]` | 方向 `d` 内的 active 面场展平索引 |
| <span style="color: red;">原Godunov</span> 的 `flux[d].faces[face_flat]`(注意不是LTS_RK2后面的flux_) | 与 `hints_` 使用相同的索引布局 |
| LTS 的 `flux_[face_id]`、`face_velocity_[face_id]`、`speed_[face_id]` | 所有方向的独立面统一编号，与 `mesh.face_info[face_id]` 对应 |
| `mesh.face_area[d].faces[flat]` | 包含 ghost 区域的完整面场展平索引 |

此外，`hints_` 的布局保留周期边界的两个端点槽位，而 `face_info` 将周期边界两端合并为一个独立面。因此，访问压力缓存时需要将 `face_id` 转换为对应的方向和 `face_flat`。

flux_的索引同理,其保存的是黎曼问题求解出来的通量。


## 入口
## run函数
### 第一步:初始扫描
```cpp
    solve_faces(start_, all_faces_);
    initial_flux_ = flux_;
    initial_speed_ = speed_;
    initial_velocity_ = face_velocity_;
    initial_hints_ = hints_;
    const double source_limit = evaluate_sources(start_, time_, 0.0);
    initial_source_ = source_;

    double minimum_local_dt = inf;
    hydro::max_inverse_dt = 0.0;
    for (int cell = 0; cell < n_; ++cell)
    {
        const double rate =
            std::max(speed_[cell], speed_[cell + 1])
            / mesh_.cell_width(0, ng_ + cell);
        if (!std::isfinite(rate)) throw std::runtime_error("Non-finite LTS RK2 CFL rate");
        hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
        local_dt_[cell] = rate > 0.0 ? cfl_ / rate : inf;
        minimum_local_dt = std::min(minimum_local_dt, local_dt_[cell]);
    }
```
   这里我们遍历所有活跃区,然后为每一cell通过cfl计算局部dt，找出最小值，为接下来初始大步做准备。注意在这里我们用初始状态(0时刻)解了一次黎曼问题,将所有的面通量(还有面上的波速等等)存在了flux_[]数组里,接下来我们要用.
### 第二步: 选取完整步macro_step
```cpp
    macro_step_ = std::min({
        dt_limit_, std::ldexp(minimum_local_dt, max_depth_), source_limit});
```
完整步选取的逻辑是取下列项的最小值:
1. dt_limit(硬性结束时间限制)
2. 源项限制
3. $\mathrm{dt_{min}\times 2^{d_{max}}}$,对于8层即$256\mathrm{dt_{min}}$

### 第三步: 尝试推进
- 首先检查macro_step合法性
- 接着进行重置(reset_attempt)
- 分配面的层级
- 推进一步




#### 重置(reset_attempt)逻辑
- 首先设置时间步层级,若CFL步长小于二分后的结果,就再次升级,直到大于为止
- 然后我们需要扩展同级缓冲区，将周围的低等级邻居(宽度由buffer参数指定)提升到与中心cell同等级，然后再平滑过渡到更远处的低层级处，使得邻居之间最多只相差1级，不会出现超过1级的跳跃。
- **先扩展同层缓冲区，覆盖重构所需的邻域。**

  根据每个 cell 的 CFL 要求确定其层级 `depth`，然后将自身及沿各坐标方向、距离不超过 `level_buffer_` 的邻居提升到**至少 `depth` 层**。已经处于更高层级的单元保持原样。

  使用 **PCM** 时，`level_buffer_ = 1`，需要提升两侧一格的邻居。虽然 PCM 不计算空间斜率，但一个面的黎曼问题需要左右两个 cell 分别提供状态，因此计算当前 cell 的面通量仍然需要相邻 cell。

  使用 **PLM / PLM_corrected** 时，`level_buffer_ = 2`，需要提升两侧各两格的邻居。除了面两侧的 cell，还需要它们各自的邻居来计算重构斜率。例如，计算 cell `i` 右侧面的通量，需要重构 cell `i` 和 `i+1` 的状态，其中 cell `i+1` 的斜率还会用到 cell `i+2`。

  多维情况下，沿每个坐标方向分别进行上述扩展。

- **再从最高层向下逐层处理，使层级平滑过渡。**

  对每个处于 `depth` 层的 cell，将缓冲范围内的邻居提升到**至少 `depth - 1` 层**。例如，3 级单元周围至少为 2 级，随后处理 2 级单元，使其周围至少为 1 级。这里使用的缓冲宽度同样由 `level_buffer_` 决定，已有的更高层级不会被降低。

  两个步骤都通过面拓扑寻找邻居：周期边界的首尾 active cell 已经连接，会自动绕到另一侧；遇到非周期边界外不存在的 active 邻居，则停止沿该侧扩展。当前多维版本不需要在这个函数中手动取模。

注意，**扩展同层级区域**和**平滑层级过渡**都在 `assign_time_levels()` 中通过调用 `promote_neighbors(cell, depth)` 完成，区别在于传入的目标层级。

- **第一步：扩展同层级区域。**

  根据当前 cell 的 CFL 限制计算所需层级 `depth`，然后调用：

  ```cpp
  promote_neighbors(cell, depth);
  ```

  将当前 cell 及缓冲范围内的邻居提升到至少 `depth` 层。

- **第二步：平滑层级过渡。**

  从最高层向下遍历，对于当前处于 `depth` 层的 cell，调用：

  ```cpp
  promote_neighbors(cell, depth - 1);
  ```

  将其缓冲范围内的邻居提升到至少 `depth - 1` 层。例如，3 级单元周围至少为 2 级，随后处理 2 级单元，使其周围至少为 1 级。

两次调用都通过 `std::max` 更新层级，因此只会提升层级，已有的更高层级不会降低。


- 更新：兼容多维后,promote_neighbor的实现方式修改，因为每一个cell有2个/3个维度与其他cell相邻,每个维度通常有2个方向的邻居(非周期边界只有一个active邻居)：
```cpp
void LTSRK2Step::promote_neighbors(int cell, int depth)
{
    level_[cell] = std::max(level_[cell], depth);

    for (int d = 0; d < hd_core::NDIM; ++d)
    {
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

                level_[neighbor] =
                    std::max(level_[neighbor], depth);
            }
        }
    }
}
```
以三维为例
1. 第一层for循环遍历三个维度
2. 第二层for循环遍历向低侧和向高侧两个方向,ide=-1 表示向坐标减小的低侧移动，side=1 表示向坐标增大的高侧移动。每次开始处理一侧时，都将 neighbor 重置为中心 cell,
3. 第三层循环先沿着低侧方向向外移动neighbor指向的cell,寻找邻居并将其promote,然后再向外移动neighbor,最多走到buffer格外的邻居处,或遇到非周期边界时,结束这一次循环,neighbor退回中间，退回到第二层循环。然后再沿高侧方向向外移动neighbour指向的cell,promote邻居,同理继续循环直到距离达到或遇到非周期边界。然后退回到第一层循环，更换维度继续重复2和3。
4. 第一层所有维度循环结束后遍历了所有邻居。例如，中心 cell 的坐标为 `(i,j,k)`，当 `level_buffer_=2` 时，六个方向的访问顺序分别为：

| 移动方向 | 第一次移动 | 第二次移动 |
| --- | --- | --- |
| x 低侧 | `(i-1,j,k)` | `(i-2,j,k)` |
| x 高侧 | `(i+1,j,k)` | `(i+2,j,k)` |
| y 低侧 | `(i,j-1,k)` | `(i,j-2,k)` |
| y 高侧 | `(i,j+1,k)` | `(i,j+2,k)` |
| z 低侧 | `(i,j,k-1)` | `(i,j,k-2)` |
| z 高侧 | `(i,j,k+1)` | `(i,j,k+2)` |

因此，扩展范围是**穿过中心 cell 的三条坐标轴上的邻域**，不包含斜对角单元。

在远离边界、邻居互不重复的情况下：

- `level_buffer_=1`：中心 cell 加上 6 个邻居，共 7 个 cell。
- `level_buffer_=2`：中心 cell 加上 12 个邻居，共 13 个 cell。

#### 分配面到不同层级逻辑(assign_faces_to_levels)
- 开场的三个二维数组遵循相同结构,第一维是层级,第二维是面或cell的index(注意,层数从0开始,因此第一维是deepest_ + 1个层级)
- cells_by_level 存储每个层级处理哪些cell
- 对于每一个面,取其左右两侧cell中更高的层级,作为其owner的层级.特别地,left=-1或right=-1(小于0)表示这是非周期边界下的最左/最右面,没有左/右侧cell，取另一侧level
- 面的owner层级表示在哪个层级积分这个面的通量
- needed会在所有低于或等于这个面的owner层级的层级里注册这个面,或者说needed[depth]会记录所有层级大于或等于depth的面,方便复用逻辑。
- 重置所有记录

## 推进步逻辑 advance_level(depth,fraction,reuse)
- 这个函数是将层级为depth的所有cell从起点$\mathrm{fraction}$积分推进到终点$\mathrm{fraction+ 2^{-depth}}$(cells_by_level[depth]),同时记录RK2最终计算出来的斜率,作为终点处的斜率临时保存在slope中,便于更低层级处进行预测。
- 整个函数是递归执行，因此会在不同层级不同时刻呼叫这个函数。每一个层级会呼叫两次高一层级。参数 depth 代表本层的层级, fraction代表这次递归调用的起点在长步的几分之几处。而reuse_start_stage则是一个优化开关,表示是否要复用之前的计算结果,这个在后面会解释
- 用到的几个数组: base_ 存储起点处所有cell的状态, view_ 存储所有cell处于同一时刻的状态用于面通量的计算,其中有一些cell的状态是RK2推进计算的准确状态，有一些来自RK2过程中计算的起点/终点斜率外推的预测值。state_里存储的是每一个cell经历过最近一个局部时间步用RK2计算出来的最终的变化率(起点终点斜率的平均值)进行推进得到的终点状态,不同cell其state_中存储的状态对应的时间可能是不一致的(步长不同)。

- 计算本步span以及实际步长macro_step * span, 还有结束时的分数步end_fraction,然后进入RK2的逻辑: 
 推进本层(level=depth)的cell,需要计算起点处cell守恒量状态U的时间变化率(通过解黎曼问题计算面通量),然后再预测终点处的状态U',再次计算终点处的时间变化率,取平均后乘上时间步dt,然后积分推进
### RK2 Phase I——首先计算起点处cell状态的的时间变化率$\frac{dU}{dt}\Big|_{\mathrm{begin}}$
- 首先检测是否为复用步。
```cpp

    if (!reuse_start_stage)
    {
        if (!build_stage_view(depth, fraction, false)) return false;
        solve_faces(view_, faces_needed_by_level_[depth]);
        if (macro_step_ > evaluate_sources(
                view_, time_ + fraction * macro_step_, macro_step_)) return false;
    }
```

在这里复用的意思是,我们是否能够直接使用现在存在flux_数组里的面通量(以及源项贡献的通量,面波速等等,不过主要是通量)?是否能够复用是在后面的代码里决定的,具体来说就是在advance_level的第三个bool reuse参数里


复用逻辑我们稍后讲解.


如果不复用，那么进入一个build_stage_view函数。这个函数的作用在于,让层级比当前更低（步长更长的慢步）的面用他们之前计算的斜率临时预测到当前递归阶段的起点和终点(RK2需要起点和终点的斜率取平均)，参与后面的通量计算。在这里是预测起点。

由于是用斜率进行预测,而不是再次算一次RHS+源项推进，消耗不大。

示例:假如有三个cellA,B,C，分别处于层级0,1,2.当cell C需要从0.75推进到1.00时,此时B已经推进到0.5(层级1,步长为0.5),因此B用0.5处的斜率(之前计算的RHS+源项,临时保存)预测0.75处,而A则用0处的斜率预测0.75处,得到
$$
   U_A^*(0.75) = U(A,0) + 0.75 * \frac{dU_A}{dt}\Big |_{0}
$$
$$
   U_B^*(B,0.75) = U(B,0.5) + 0.25 * \frac{dU_B}{dt}\Big |_{0.5}
$$
带*号表示这是根据斜率计算的中间临时预测值,而且这个值并不会真正更新到网格中。预测的状态会更新到数组view_中

那为什么要计算？例如，如果B与C是相邻的网格，那么B、C之间存在相接面,由于RK2需要起点0.75处cell C左右两侧的通量,因此我们需要计算更低层级的cell在0.75处的值来解黎曼问题(理论上如果只有RHS,只需要算B---如果不是周期边界的话,否则也要算A。但有全局源项,比如自引力的话,那么需要计算所有低层级cell的变化)。而B的一个时间步是0.5,无法在其正常的时间步中计算到0.75处的状态,因此我们用0.5处的时间导数进行大致预测。

每一层advance_level的递归有depth参数和fraction参数,fraction为起点对应的分数步,在上面的例子里就是层级2,起点为0.75,我们计算所有层级0和层级1的cell在0.75处的预测值,而endpoint=false表示这里是起点处预测,我们不需要预测层级2的cell在这里的值
```cpp
        if (!build_stage_view(depth, fraction, false)) return false;
```
层级1(depth=1)的递归中,调用(2,0.5,true)后再调用(2,0.75,false),因此0.75->1.0这一步调用之前,必然已经完成了0.5->0.75,因此0.75处depth=2的所有cell状态已经计算好，不需要预测。

更高层级的cell也不需要预测，因为它们必然在这一层到达fraction之前已经推进到fraction处,与这一层共用起点状态(递归的特性,可以看[下面](#recursion-feature-1))。比如层级0中执行advance_level(1,0),让层级1先从0推进到0.5,在这个过程中层级1会执行两次advance_level(2),把层级2从0推进到0.25,再从0.25推进到0.5。然后层级0再次执行下一行advance_level(1,0.5)(如上,在层级0中这两次advance是紧挨着的),让层级1从0.5推进到1,这个时候层级2以及更高层级必然已经推进到0.5,与层级1共享起点,因此不需要预测。
在build_stage_view函数中,我们可以看到
```cpp
    view_.cells = state_.cells;
    for (int cell_index = 0; cell_index < n_; ++cell_index)
    {
        const bool needs_prediction =
            level_[cell_index] < depth || (endpoint && level_[cell_index] == depth);
        if (!needs_prediction) continue;

        auto& cell = view_.cells[ng_ + cell_index];
        const double elapsed = (fraction - local_step_start_[cell_index]) * macro_step_;
        cell.U = base_[cell_index] + slope_[cell_index] * elapsed;
        if (!is_physical(cell.U, cell.gamma)) return false;
    }
    boundary::FillBoundary(view_, time_ + fraction * macro_step_, hooks_);
    return true;
```
view_首先继承state_中存储的所有cell的状态。state_存储的是所有cell在各自上一个局部时间步推进后的状态,对于不同层级的cell,它们上一次经历的时间步可能不一样,推进的终点也可能不一样。但根据我们刚才的分析,此时所有level>=depth的cell都推进到了共同的终点,也是这一步的起点fraction处,其对应的state_中存储的状态就是fraction处的状态,只有层级更低的cell因为其步长不够精细没有办法推进到fraction处,因此在下面的for循环里我们为更低层级的cell进行外推预测,预测所用的是之前保存的RK2斜率。


做完预测，接下来,用预测状态view_解这一层级的所有面的黎曼问题。

    ```cpp
        solve_faces(view_, faces_needed_by_level_[depth]);
    ```

这个solve_faces函数解出来的面通量，保存在flux_[]数组内. 注意这里我们要计算的是所有needed的face,也就是owner是这一层以及更高层的面(时间步更短)。
首先,为什么更低层级的面不用计算通量?因为我们的每个面的owner取的是这个面的两边cell中更高层级的那个,因此owner_level<escape><</escape>depth 的面必然不是更高level的cell拥有的面(相邻面),而advance(depth,fraction,reuse)这一层递归只会推进层级为depth的cell，所以这样的面没有影响，不用计算其通量。

因此为什么要计算owner_level>=depth的面的通量也很显然了,因为它们有可能是当前level=depth的cell的所属面,而在这个共同起点这些面的flux还都没有计算.不过共同起点处所有面的flux只会计算一次,因为后面高层会由于reuse=true跳过这个if块内的solve_faces,直接用低层计算好的flux.比如在t=0.5处,层级1计算了所有层级1和高于层级1的面的flux,然后通过advance_level(2,0.5,true)和advance_level(2,0.75,false)呼叫层级2,由于0.5处是层级1和层级2的共同起点,而层级1已经把更高层级的面在0.5处的通量计算好了,因此第一个advance_level的reuse为true,而层级2在推进0.5->0.75时，就跳过了这个if块,直接复用在0.5处的flux,避免重复计算owner_depth>=2的面在0.5处的通量。

最后再进行一次源项的检查,确认预测后的状态是否合法:
```cpp
        if (macro_step_ > evaluate_sources(
                view_, time_ + fraction * macro_step_, macro_step_)) return false;
```
    
- 如果可以复用,那么就跳过上面的过程.

举例说明复用的逻辑:

假设我们从depth=0,最低层级的递归出发,最开始fraction=0,也就是所有层级共同的起点,层级0运行完前面的计算逻辑后,让层级1运行fraction=0处的计算走到0.5,起点相同,所以不用重新预测网格,直接沿用之前计算的通量(对于t=0,我们最开始执行run的时候就计算了一次,请见run函数的第一步)。

但接着层级0让层级1从0.5出发,再走一步到1.0达到共同终点,这一步的时候层级1与层级0不在一个起点了,这时候就要重新预测0.5处的所有层级0的cell的状态,然后计算所有owner层级大于等于1的面的通量，并进行接下来的逻辑(计算slope,推进).

同理,层级1会呼叫更高一层的递归2次,第一次是在层级1的0->0.5的歩进中,会呼叫层级2的0->0.25和0.25->0.5,同样的,层级2从0->0.25也会复用0处的flux_,但0.25->0.5不会(需要预测0.25处所有低层级cell的值). 

而在层级1 0.5->1.0的步进中,由于层级1的这次递归(被层级0呼叫的advance_level(1,0.5,false))不复用,因此它预测了一次view_并用其重新计算了一次flux_,所以现在flux_里保存的是0.5时所有面的通量。接下来这次递归会依照上面的代码呼叫层级2的0.5->0.75和0.75->1.0,因此0.5->0.75这步就可以复用,而0.75->1.0这步则不复用,重新用斜率预测.


- 接下来,进行一次面的合法性检查，然后临时保存起点处的面通量。这是由于在RK2最后我们需要起点和终点处的通量取平均，再乘上面积作为最终的状态对时间的导数(斜率),然后乘上时间步得到cell的状态变化量。这里是不包含源项的。对于源项，我们只取终点处的值作为近似，而不用RK2处理。

```cpp
    if (!check_wave_speeds(faces_needed_by_level_[depth])) return false;

    for (int face : faces_owned_by_level_[depth]) first_flux_[face] = flux_[face];
```

- 然后我们计算当前层级所拥有的cell在起点fraction处的斜率slope(即$\frac{dU}{dt}$),并记录其local起点为fraction,用于上文提到的在更高层级处的估算。
```cpp
    for (int cell_index : cells_by_level_[depth])
    {
        const auto& current_cell = state_.cells[ng_ + cell_index];
        base_[cell_index] = current_cell.U;
        local_step_start_[cell_index] = fraction;
        source_start_[cell_index] = source_[ng_ + cell_index];

        pressure_work_start_[cell_index] =
            -(current_cell.gamma - 1.0) * base_[cell_index].internal_energy()
            * (area_[cell_index + 1] * face_velocity_[cell_index + 1]
               - area_[cell_index] * face_velocity_[cell_index])
            / volume_[cell_index];
        slope_[cell_index] =
            (flux_[cell_index] * area_[cell_index]
             - flux_[cell_index + 1] * area_[cell_index + 1])
            / volume_[cell_index] + source_start_[cell_index];
        slope_[cell_index][hd_core::INTERNAL_ENERGY] += pressure_work_start_[cell_index];

        if (!is_physical(base_[cell_index] + local_step * slope_[cell_index],
                         current_cell.gamma)) return false;
    }
```
slope_[cell_index]保存的就是某个cell在当前步的RK2起点斜率,而不是RK2用于推进时间步的最终斜率。这个斜率可以用在后面对cell状态的预测中.
### 向更高层递归
- 接下来是递归的核心,
```cpp

    if (depth < deepest_)
    {
        if (!advance_level(depth + 1, fraction, true)) return false;    //将更高一层从fraction推进到fraction + span/2
        if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;//将更高一层从fraction+span/2推进到fraction +span/2 + span/2,也就是end_fraction
    }
```
因为反复提及多次了所以不再赘述，简而言之就是调用两次下一层级的步进. 如果步进结果为false那么也立刻返回false(最终会让外面的循环continue,再次尝试缩小一倍macro_step).
<a id="recursion-feature-1"></a>
注意,当执行完advance_level(depth + 1, fraction, true)后,所有高于depth的层级都必然推进到fraction+ span /2, 这是递归的特性决定的。

我们举例说明,<span style="color: red;">记depth这一层的步长为span</span>, depth这一层执行depth+1层的第一个半步推进 advance_level(depth+1,fraction)后, 进入depth+1层的递归。

depth+1层在执行完最开始RK2 Phase I的逻辑后会再次执行到这两行(如果depth+1不是最高层的话), 呼叫更高一层 depth+2 的两个1/4步推进 advance_level(depth+2, fraction,true) 和 advance_level(depth+2,fraction + span / 4.0)

同理呼叫的 depth+2 层推进在执行完RK2 Phase I的逻辑后也会继续呼叫更高层,直到递归到最高层后逐层下降返回。为了简化,我们假设 depth+2 就是最高层,那么 depth+2 就不会再呼叫 depth+3 , 接着执行完后面的RK2 Phase II逻辑后, depth+2 层就会返回到 depth+1 层。而我们可以注意到, depth+1 层中执行了两次 depth+2 层的1/4步推进, 第一次将depth+2推进到了fraction + span /4.0处, 第二次将depth+2推进到 fraction + span/4.0 + span/4.0 = fraction + span / 2.0处。 接着 depth+1 层会执行剩下的RK2 Phase II逻辑并返回 depth 层, 此时 depth+1 层也推进到了fraction + span / 2.0处。注意，这里是执行完 depth+1 层第一次半步推进
```cpp
        if (!advance_level(depth + 1, fraction, true)) return false; 
```
之后。当有更多层级时也一样, depth+1层 呼叫的两次 depth+2 层推进会再次呼叫更高层并把它们逐步推进到 fraction + span/2.0处。所以depth + 2 以及更高层必然会比depth +1 先到fraction + span/2.0这第一个半步的终点处。

同理,在执行完下一个半步推进
```cpp
    if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
```
之后,所有比 depth 高的层都推进到了 fraction + span/2.0 + span/2.0 = fraction +span = end_fraction,也即本层( depth 层)即将但还未到达的终点处,然后本层继续执行RK2 Phase II的逻辑,推进到 end_fraction 处。
所以更高的层级一定比本层先推进到这次advance_level的终点,这是递归的特性保证的。
某一层推进到某一个分数步处,意味着state_中属于那个层的cell状态更新到了处于这个分数步时的状态,flux_中属于那个层的面的通量也被更新到了处于这个分数步时的通量。

进行一次跳过检测(本层是否有cell)之后
```cpp
    if (cells_by_level_[depth].empty()) return true;
```
我们进入到RK2的Phase II
### RK2 Phase II——计算终点处cell状态的的时间变化率$\frac{dU}{dt}\Big|_{\mathrm{end}}$

```cpp
    if (!build_stage_view(depth, end_fraction, true)) return false;
```
这里我们预测所有更低层以及当前层的cell在终点end_fraction处的状态,低层为何要预测同理。唯一不同的是这次我们令endpoint=true,在build_stage_view中
```cpp
    view_.cells = state_.cells;
    for (int cell_index = 0; cell_index < n_; ++cell_index)
    {
        const bool needs_prediction =
            level_[cell_index] < depth || (endpoint && level_[cell_index] == depth);
        if (!needs_prediction) continue;

        auto& cell = view_.cells[ng_ + cell_index];
        const double elapsed = (fraction - local_step_start_[cell_index]) * macro_step_;
        cell.U = base_[cell_index] + slope_[cell_index] * elapsed;
        if (!is_physical(cell.U, cell.gamma)) return false;
    }
    boundary::FillBoundary(view_, time_ + fraction * macro_step_, hooks_);
    return true;
```
意味着我们还要为当前层级预测终点状态,这是因为当前层级还没有推进到终点end_fraction,在普通的RK2中这也是需要的——从起点的斜率外推预测一个终点状态，再从这个预测终点状态得到终点斜率,最后与起点斜率取平均作为最终变化率推进。更高层级依然不需要预测终点状态，因为就在之前已经

递归调用了两次下一层级的advance_level,已经将更高层级推进到了终点end_fraction,state_中存储的状态就是end_fraction处的状态,因此直接复制到view_中不用修改。


- 于是我们更新好了view_,其中存储的cell状态为预测的终点处的状态。接下来我们再用终点处的view_状态更新面通量flux_,
```cpp
    solve_faces(view_, faces_owned_by_level_[depth]);
```

- 再经过一个合法性检测后我们到达终点处的斜率计算以及最终推进时间步处
```cpp
    for (int face : faces_owned_by_level_[depth])
    {
        const auto integral =
            (first_flux_[face] + flux_[face]) * (0.5 * local_step * area_[face]);
        const int left = left_cell(face);
        const int right = right_cell(face);
        if (left >= 0) pending_flux_integral_[left] -= integral;
        if (right >= 0) pending_flux_integral_[right] += integral;
    }

    for (int cell_index : cells_by_level_[depth])
    {
        auto& cell = state_.cells[ng_ + cell_index];
        cell.U = base_[cell_index]
            + pending_flux_integral_[cell_index] / volume_[cell_index]
            + (source_start_[cell_index] + source_[ng_ + cell_index])
                * (0.5 * local_step);
        const double pressure_work_end =
            -(cell.gamma - 1.0) * view_.cells[ng_ + cell_index].U.internal_energy()
            * (area_[cell_index + 1] * face_velocity_[cell_index + 1]
               - area_[cell_index] * face_velocity_[cell_index])
            / volume_[cell_index];
        cell.U[hd_core::INTERNAL_ENERGY] +=
            0.5 * local_step * (pressure_work_start_[cell_index] + pressure_work_end);
        if (!is_physical(cell.U, cell.gamma)) return false;
        cell.synchronize_dual_energy();
        pending_flux_integral_[cell_index] = Conserved{};
    }
```
- 第一个循环,我们对owner为当前层级的所有面,计算RK2的最终推进值-通量积分,即面通量的平均值乘上再乘上面积(即RHS，或者说最终的RK2斜率),再乘上当前的时间步长(local_step = span * macrostep)。然后把通量积分的负值添加到面的左侧cell的累积通量积分pending_flux_integral_[left],正值添加到右侧的cell的累积通量pending_flux_integral_[right]中，这些pending_flux_integral_要一直逐层向下传递,直到这个cell的depth推进的那一层,在积分推进中将它加到对应cell的状态中作为积分步里的变化量。

这里的逻辑与前面稍有不同，在前面phase1我们求面通量黎曼问题的逻辑是求owner_level<escape>>=</escape>depth的本层和更高层级的flux,然后在更高层级处通过复用标记reuse flag跳过solve_face避免重复求解,复用低层求出的通量。当然,因为只有owner_level<escape>>=</escape>depth的通量会影响depth层级的cell的状态,我们复用的其实是这一部分。
而在这里,我们前面求解黎曼问题的solve_face只求了owner为本层级的通量作为终点end_fraction处的通量,owner_level<escape><</escape>depth的面的通量不影响depth层级的cell,owner level等于当前层depth的通量由solve_faces求解。而在前面的两步更高一级的递归 advance_level(depth+1,fraction)和advance_level(depth+1,fraction + span/2.0)以及在这两步递归中更高层级的递归里(都推进到了end_fraction),一定把owner_level高于depth的面在终点end_fraction处的通量求了出来。这些面的通量会被加到右侧cell的pending_flux_integral_[right],以及其负值会被加到左侧cell的pending_flux_integral_[left]里,然后pending_flux逐层向下累积传播,到depth这一层再给pending_flux_integral_加上(右侧cell)/减去(左侧cell)刚刚算出的integral,就是所有owner层级大于等于depth的面对于这个cell的通量之和,而更低层owner的面又对这个cell没有影响,因此这就是这个cell的RHS。于是在后面的循环里我们把所有这一层的cell加上pending_flux_integral,就是RK2对这一层的所有cells最终推进结果。然后我们把这一层的cells的pending_flux_integral清空，因为这一depth的cell随这一个时间步的变化已经结算完了,下一步又需要用到新的通量。比如层级0调用两次层级1,第一次是0到0.5，第二次是0.5到1.0，那么层级1会在0到0.5结算一次pending_flux_integral_加到层级1的cell里,然后清空这些pending_flux_integral_(不影响层级0的cell的pending_flux向下传递),然后0.5到1.0又是新的一步,需要重新计算在这新一步里层级1的所有cell的RHS。

- 一切正常进行,返回真,这就是整个LTSRK2的流程
## 流程总结
0. 外部run函数各种各种初始化,包括初始化数组,时间步的分配,cell和face的层级分配,计算出最大步长macro_step(由所有cell中最小的可行dt * $2^{\mathrm{max_depth}}$),然后开始尝试循环。
1. 首先由外部advance_lts_rk2计算一次最初起点处的flux，然后呼叫advance_level(0,0,true).这是起点
2. 每一层的advance_level有层数参数depth,起点参数fraction(分数步长)和是否复用reuse. 计算当前的分数步长span = $2^{-\mathrm{depth}}$和终点end_fraction = fraction + span.先检测是否为复用步,如果不是,复制RK2推进的cell状态结果state_到临时数组view_中,由于步长不一,每个cell在state_里的状态对应的时间是不一致的。通过此前存储的RK2第一部分斜率将所有更低层的cell(因为步长更长,没有推进到这一层的起点)在view_中的状态外推到当前层的起点fraction处。由于递归特性,当前层和更高层级的cell必然已经推进到当前起点fraction处,因此不用修改从state_处复制的值,于是整个view_同步。然后用view_里的cell状态去求一次黎曼问题,求出所有owner为当前层级和更高层级的面上的通量flux_,因为只有当前/更高层级owner的面才可能是当前层级的cell的所属面,通量才会对这些cell有影响。

   如果是复用步,意味着state_中所有的cell处于共同起点,不需要外推,也不需要重新计算面的flux_(上一步已经算过)
3. 保存起点处各个面的通量flux_到first_flux_下。 这个通量乘上时间步和面的面积后,就是起点处不含源项的cell状态U对时间的变化率,也就是RK2第一部分所需要的斜率$k_1$。我们对源项只采用终点处的值进行处理。

 然后遍历本层的所有cell,综合所有的流体通量和源项计算起点处cell的状态对时间的变化率(斜率)$\frac{dU}{dt}\Big |_{begin}$。这个斜率是包含源项的完整的变化率，我们还要存进slope里,为了第二步中提到的外推预测。
4. 两次推进更高一层的递归,先执行第一次advance_level(depth+1,fraction,true),这是让更高一层从fraction走到fraction+span/2,因为共起点所以复用为true。注意这是递归,这里会由下一层advance_level(depth+1)跳转到步骤2继续执行这个流程,直到最高一层计算完毕,或某一层没有cell,才会逐层往上返回到这一层。因此,当执行完这一次advance_level返回后,所有高于depth的层级都推进到了fraction + span/2。

接着这一层继续执行第二次advance_level(depth+1,fraction + span/2,true),这是让更高一层从fraction + span/2 走到end_fraction,由于这次推进起点是本层的起点向后走半个本层的步长,因此本层无法与这一次推进共起点,所以设置reuse=false。同理是递归呼叫,更高一层跳回步骤2继续执行流程，因此只有再所有更高层全部推进到end_fraction后,才会执行本层接下来的所有逻辑。

注意以上过程没有讨论不合法的情形,如果出现合法性检查不通过的层,直接逐级return false,回到run函数里的尝试(attempt)循环中,将最大步长macro_step缩小一半再次尝试。

5.  接下来进入RK2的Phase II。要计算终点处的cell状态对时间的变化率(终点斜率),我们首先通过build_stage_view预测终点end_fraction处的cell状态view_。将state_再次复制到view_中,由于上述提到的递归特性,更高层级的cell已经推进到了终点处,所以更高层级的cell不用修改。而当前层级还没有推进到终点,所以对于当前层级以及更低的层级的cell,都需要用slope外推到终点。

6. 用view_中的cell状态计算solve_face,这一次只计算owner为当前层级的face的flux_。对于owner为更低层级的面,同理,对当前层级的cell的状态没有影响(不可能是其所属相邻面)。owner为更高层的面,其所属的层级已经推进到终点,在这些层级的推进过程中已经计算过终点处的flux_,我们通过一个逐层向下累积传递的pending_flux_integral_[cell_index]给低层级复用。

7. 计算owner为这一层的面产生的通量积分，这里用的是RK2在起点和终点计算的flux取平均值(即平均斜率$\frac{k_1+k_2}{2}$)作为通量，乘上面积和时间步得到通量积分,并把通量积分的负值添加到面的左侧cell的累积通量积分pending_flux_integral_[left],正值添加到右侧的cell的累积通量pending_flux_integral_[right]中，这些pending_flux_integral_要一直逐层向下传递,直到这个cell的depth推进的那一层,在积分推进中将它加到对应cell的状态中作为积分步里的变化量。

8. 积分推进步，遍历属于这一层的cell, 每个cell加上累积通量数组pending_flux_integral_中属于自己的那一个累积通量pending_flux_integral_[cell_index],这是从最高层一直累积到当前层depth的,累积给这个cell的RHS。再加上源项,就完成了这一步推进对这一层所有cell的积分。最后清理掉pending_flux_integral_中属于这一层cell的通量。这是由于后续更低层的积分推进中,已经用不到了(一个cell不可能属于两个层级)，这样所有层级完成后这个数组也全部清零。在中间步的计算中, 如果depth层级的某一次advance_level完成了的话,那么depth层和更高层的cell对应的pending_flux_integral_都会清零，计算更低层(如 depth-1 层)时也不会影响到更低层的cell。
举例:   f1  C1  f2 C2  f3 C3 f4
C1属于第1层,C2属于第0层,C3属于第二层，f1的owner_level维1,f2的owner_level也为1,f3和f4的owner_level为2。在层级2的时候会计算owner level为2的面f3和f4的通量,并把累积通量积分加给相邻cell C2(f3)和C3(f3和f4)的累积通量积分,然后层级2的推进步会令这一层的cell C3用得到的f3和f4给的累积通量积分推进状态，然后清零C3得到的累积通量积分。而C2属于层级0，因此f3给C2的通量不会清零。然后退到层级1，会把面f1和f2的累积通量积分加给C1(f1)和C2(f1,f2)的累积通量积分,然后层级1再令C1用从f1得到的累积通量积分推进状态，然后清零C1得到的累积通量积分，但C2得到的通量仍然不会清零,这时C2得到的累积通量积分已经是面f2和(负的)面f3之和。最后退回层级0,没有属于层级0的面,然后推进结算层级0的cell C2，这时就能结算在高层级逐层传递下来的累积通量积分pending_flux_integral_[C2],把相邻面f2和f3贡献的通量积分之和加载C2的cell状态中，这样尽管C2的相邻面属于更高层级，我们依然正确结算了C2的所有相邻面的贡献,并且每个面的贡献只需计算一次,不重复计算浪费性能。最后全部清零。

9. 返回真值，完成这一步的推进。递归回退到上一级后会逐层完成结算，除非合法性检查不通过，那所有的层级会快速返回false，退回run函数中的attempt循环，尝试0.5倍的宏步。


## 四个 Cell、三个时间层级的 LTS-RK2 执行示例
(这部分因为我太懒不想再写，因此是ChatGPT写的):
假设有四个非周期边界 cell：

```text
边界 | f0 | C0 | f1 | C1 | f2 | C2 | f3 | C3 | f4 | 边界

层级：      0          1          2          1
时间步：    H         H/2        H/4        H/2
```

也就是：

```text
C0：level 0，一个宏步内推进 1 次
C1：level 1，一个宏步内推进 2 次
C2：level 2，一个宏步内推进 4 次
C3：level 1，一个宏步内推进 2 次
```

> 注意：当前实现会根据重构模板建立层级缓冲区，因此真实的 `assign_time_levels()` 可能继续提升相邻 cell。这里把 `[0, 1, 2, 1]` 当作已经确定的层级，用于说明 LTS-RK2 的递归执行过程。

<a id="face-ownership"></a>

### 1. 分配面的负责层级

一个面由左右相邻 cell 中层级更高、时间步更小的一层负责积分。

| 面 | 左侧 | 右侧 | 负责层级 |
|---|---|---|---:|
| `f0` | 非周期边界 | `C0`，level 0 | 0 |
| `f1` | `C0`，level 0 | `C1`，level 1 | 1 |
| `f2` | `C1`，level 1 | `C2`，level 2 | 2 |
| `f3` | `C2`，level 2 | `C3`，level 1 | 2 |
| `f4` | `C3`，level 1 | 非周期边界 | 1 |

因此：

```cpp
faces_owned_by_level_[0] = {f0};
faces_owned_by_level_[1] = {f1, f4};
faces_owned_by_level_[2] = {f2, f3};
```

各层起点阶段需要准备的面为：

```cpp
faces_needed_by_level_[0] = {f0, f1, f2, f3, f4};
faces_needed_by_level_[1] = {f1, f2, f3, f4};
faces_needed_by_level_[2] = {f2, f3};
```

其中：

```text
faces_owned_by_level_：
当前层真正负责进行时间积分的面。

faces_needed_by_level_：
当前层计算 cell 斜率以及为子层准备起点时，需要计算的面。
```

<a id="initial-stage"></a>

### 2. 宏步起点

在宏步起点 $t=0$，所有 cell 都处于相同物理时刻。

程序首先计算：

```text
f0、f1、f2、f3、f4 的初始通量
所有面的波速和面速度
所有 cell 在 t=0 的源项
```

对应：

```cpp
solve_faces(start_, all_faces_);
evaluate_sources(start_, time_, 0.0);
```

然后进入最外层递归：

```cpp
advance_level(0, 0.0, true);
```

<a id="level-0-start"></a>

### 3. Level 0：准备 C0 的局部步

C0 的局部时间步为：

```text
C0：0 → H
```

首先保存 level 0 拥有的面 `f0` 的起点通量：

```cpp
first_flux_[f0] = flux_[f0];
```

然后使用 C0 的左右面 `f0、f1` 和起点源项计算斜率：

$$
\mathrm{slope}_{C0}(0)
=
\frac{
F_{f0}(0)A_{f0}
-
F_{f1}(0)A_{f1}
}{
V_{C0}
}
+
S_{C0}(0)
$$

保存：

```text
base_[C0]  = C0 在 0 时刻的状态
slope_[C0] = C0 在 0 时刻的变化率
```

C0 此时不会立即更新到 $H$。

它必须等待更细层完成，并收集 `f1` 在整个宏步内的通量积分。

<a id="level-1-first-step"></a>

### 4. Level 1 第一步：C1、C3 从 0 推进到 H/2

进入：

```cpp
advance_level(1, 0.0, true);
```

因为 `reuse_start_stage=true`，直接复用宏步起点已经计算好的通量和源项。

保存 level 1 拥有的起点面通量：

```cpp
first_flux_[f1] = flux_[f1];
first_flux_[f4] = flux_[f4];
```

计算 C1 的起点斜率：

$$
\mathrm{slope}_{C1}(0)
=
\frac{
F_{f1}(0)A_{f1}
-
F_{f2}(0)A_{f2}
}{
V_{C1}
}
+
S_{C1}(0)
$$

计算 C3 的起点斜率：

$$
\mathrm{slope}_{C3}(0)
=
\frac{
F_{f3}(0)A_{f3}
-
F_{f4}(0)A_{f4}
}{
V_{C3}
}
+
S_{C3}(0)
$$

C1、C3 暂时不会提交到 $H/2$，因为它们需要等待 C2 完成两个 $H/4$ 小步。

<a id="level-2-first-step"></a>

### 5. Level 2 第一步：C2 从 0 推进到 H/4

进入：

```cpp
advance_level(2, 0.0, true);
```

保存 level 2 拥有的面 `f2、f3` 的起点通量：

```cpp
first_flux_[f2] = flux_[f2];
first_flux_[f3] = flux_[f3];
```

使用 `f2、f3` 和源项计算 C2 的起点斜率：

$$
\mathrm{slope}_{C2}(0)
=
\frac{
F_{f2}(0)A_{f2}
-
F_{f3}(0)A_{f3}
}{
V_{C2}
}
+
S_{C2}(0)
$$

为了计算 RK2 终点通量，程序构造 $H/4$ 时刻的预测状态：

```text
C0：使用 slope_C0(0) 预测到 H/4
C1：使用 slope_C1(0) 预测到 H/4
C2：使用 slope_C2(0) 预测到 H/4
C3：使用 slope_C3(0) 预测到 H/4
```

预测公式为：

$$
U_i^*(H/4)
=
U_i(0)
+
\mathrm{slope}_i(0)\frac{H}{4}
$$

然后使用预测状态求解 `f2、f3` 在 $H/4$ 时刻的 Riemann 问题。

得到终点通量后，对面通量进行 RK2 梯形积分：

$$
I_{f2}^{0\rightarrow H/4}
=
\frac{
F_{f2}(0)+F_{f2}(H/4)
}{2}
A_{f2}\frac{H}{4}
$$

$$
I_{f3}^{0\rightarrow H/4}
=
\frac{
F_{f3}(0)+F_{f3}(H/4)
}{2}
A_{f3}\frac{H}{4}
$$

同一份面通量积分以相反符号写入左右 cell：

```text
f2：写入 C1 和 C2
f3：写入 C2 和 C3
```

C2 的局部步已经结束，因此立即正式提交到 $H/4$。

C1、C3 尚未完成自己的 $H/2$ 局部步，所以它们收到的通量暂存在：

```cpp
pending_flux_integral_[C1]
pending_flux_integral_[C3]
```

<a id="level-2-second-step"></a>

### 6. Level 2 第二步：C2 从 H/4 推进到 H/2

进入：

```cpp
advance_level(2, 0.25, false);
```

因为 `reuse_start_stage=false`，需要重新构造 $H/4$ 时刻的阶段网格：

```text
C0：预测到 H/4
C1：预测到 H/4
C2：使用已经正式到达 H/4 的状态
C3：预测到 H/4
```

然后重新求解 `f2、f3` 在 $H/4$ 的起点通量。

之所以不能直接复用上一个小步的终点通量，是因为：

```text
上一个小步终点使用的是 C2 的 Euler 预测状态；
当前小步起点使用的是 C2 的 RK2 修正状态。
```

随后：

1. 计算 C2 在 $H/4$ 时刻的新斜率；
2. 将相关 cell 预测到 $H/2$；
3. 求解 `f2、f3` 在 $H/2$ 的终点通量；
4. 对 $H/4\rightarrow H/2$ 的面通量做梯形积分；
5. C2 正式提交到 $H/2$；
6. C1、C3 继续累计来自 `f2、f3` 的通量。

此时：

```text
C2：正式位于 H/2
C1、C3：尚未提交，但已经收齐 0 → H/2 的细层界面通量
```

<a id="level-1-first-commit"></a>

### 7. Level 1 第一步结算：C1、C3 到达 H/2

C2 完成两个小步后，构造 level 1 在 $H/2$ 时刻的终点阶段：

```text
C0：使用 slope_C0(0) 预测到 H/2
C1：使用 slope_C1(0) 预测到 H/2
C2：直接使用正式的 H/2 状态
C3：使用 slope_C3(0) 预测到 H/2
```

然后求解 level 1 拥有的面：

```text
f1、f4
```

计算梯形通量积分：

$$
I_{f1}^{0\rightarrow H/2}
=
\frac{
F_{f1}(0)+F_{f1}(H/2)
}{2}
A_{f1}\frac{H}{2}
$$

$$
I_{f4}^{0\rightarrow H/2}
=
\frac{
F_{f4}(0)+F_{f4}(H/2)
}{2}
A_{f4}\frac{H}{2}
$$

此时 C1、C3 已经收齐自己的全部面通量：

```text
C1：左面 f1，加上右面 f2 的两个 H/4 积分
C3：左面 f3 的两个 H/4 积分，加上右面 f4
```

因此 C1、C3 正式提交到 $H/2$。

`f1` 给 C0 的通量贡献继续保存在：

```cpp
pending_flux_integral_[C0]
```

因为 C0 要等到整个宏步结束才正式结算。

<a id="level-1-second-step"></a>

### 8. Level 1 第二步：C1、C3 从 H/2 推进到 H

进入：

```cpp
advance_level(1, 0.5, false);
```

首先构造 $H/2$ 时刻的阶段状态：

```text
C0：从 0 预测到 H/2
C1：使用正式的 H/2 状态
C2：使用正式的 H/2 状态
C3：使用正式的 H/2 状态
```

然后重新求解：

```text
f1、f2、f3、f4
```

并计算 C1、C3 在 $H/2$ 时刻的新斜率。

接着让 C2 再走两个小步：

```text
C2：H/2 → 3H/4
C2：3H/4 → H
```

每个小步都执行：

1. 准备局部步起点状态；
2. 计算起点面通量和源项；
3. 计算 C2 的起点斜率；
4. 将慢层和当前层预测到局部步终点；
5. 求解 `f2、f3` 的终点 Riemann 问题；
6. 对起点、终点通量进行梯形积分；
7. 正式提交 C2；
8. 把界面通量累计给 C1、C3。

完成两个小步后：

```text
C2：正式位于 H
C1、C3：已经收齐 H/2 → H 的细层界面通量
```

随后构造 level 1 在 $H$ 时刻的终点阶段：

```text
C0：预测到 H
C1：预测到 H
C2：直接使用正式的 H 状态
C3：预测到 H
```

求解 `f1、f4` 在 $H$ 时刻的终点通量，并完成它们在 $H/2\rightarrow H$ 内的梯形积分。

最后 C1、C3 正式提交到 $H$。

<a id="level-0-commit"></a>

### 9. Level 0 最终结算：C0 到达 H

此时：

```text
C1、C2、C3：已经正式位于 H
C0：正式状态仍然位于 0
```

但是 C0 已经累计了 `f1` 的两段通量：

```text
f1：0 → H/2
f1：H/2 → H
```

接下来构造 level 0 的 RK2 终点阶段：

```text
C0：使用 slope_C0(0) 预测到 H
C1、C2、C3：使用正式的 H 时刻状态
```

求解 C0 所在层拥有的面：

```text
f0 在 H 时刻的通量
```

C0 最终使用：

```text
f0 在 0 → H 内的梯形通量积分
+
f1 累计的两个 H/2 通量积分
+
C0 起点和终点源项的平均积分
```

正式更新到 $H$。

此时：

```text
C0、C1、C2、C3 全部同步在 H
```

最后填充边界，并将 `state_` 提交给正式网格。

<a id="execution-timeline"></a>

### 10. 整体时间线

```text
C0：|----------------------- H -----------------------|

C1：|---------- H/2 ---------|---------- H/2 ---------|

C2：|-- H/4 --|-- H/4 --|-- H/4 --|-- H/4 --|

C3：|---------- H/2 ---------|---------- H/2 ---------|
```

递归执行顺序为：

```text
准备 C0 在 0 时刻的斜率
│
├─ 准备 C1、C3 在 0 时刻的斜率
│  │
│  ├─ C2：0 → H/4
│  └─ C2：H/4 → H/2
│
├─ 结算 C1、C3 到 H/2
│
├─ 准备 C1、C3 在 H/2 时刻的新斜率
│  │
│  ├─ C2：H/2 → 3H/4
│  └─ C2：3H/4 → H
│
├─ 结算 C1、C3 到 H
│
└─ 结算 C0 到 H
```

<a id="core-principles"></a>

### 11. 核心原则

整个 LTS-RK2 过程遵循以下原则：

1. 每个 cell 根据自己的层级使用不同局部时间步。
2. 每个面只由相邻 cell 中更细的一层负责时间积分。
3. 每个局部步先计算起点斜率，再预测终点状态。
4. 使用预测终点状态求解终点 Riemann 问题和终点源项。
5. 使用起点与终点通量、源项的平均值完成 RK2 更新。
6. 细层先完成多个小步，粗层等待细层界面通量全部到齐。
7. 同一份界面通量以相反符号累计给左右 cell，保证守恒。
8. 所有层最终重新同步在宏步终点 $H$。


## 一些细节和3D扩展
- 在1D架构中, 我们通过cell_index从0~nactive-1遍历活跃区,cell_index为活跃区下标,再加上nghost的offset得到flat_index。拓展到3D情形,我们需要用Index进行变换
- build_stage_view中
```cpp
        auto& cell = view_.cells[ng_ + cell_index];
```
修改为
```cpp

```

### 源项贡献和检查
evaluate_sources函数会逐项将源项贡献加到网格stage上

然后通过source_change_limit计算在给定allowance=0.1 下每个源项限制的最小dt,并且最后再取所有源项之和再计算一次，防止诸如多个加热源项单独计算dt限制那种情形。

allowance的限制是,用源项对cell状态产生的单位时间变化率rate,乘上步长dt,要小于0.1*原状态,即一步之内状态不能因为源项产生超过10%的变化. 这个allowance目前是硬编码在scan_source_changes函数中的。

source_change_limit里,输入的参数rate就是源项对每个cell的单位时间变化率(注意是在evaluate_sources里调用的,每个源项都会调用一次,此时rate是每个源项对cell产生的变化率,最后把源项求和再调用一次,rate是所有源项之和对cell产生的变化率)

if(trial_H>0.0) 的分支就是正常循环里的分支(所有的步长都必然大于0),此时调用scan_source_changes,这个函数计算每个,这里采用一票否决机制,调用中的lambda函数Bound(第三个参数)
```cpp
    scan_source_changes(stage, rates, [&](double allowance, double rate) {
        exceeded |= source_change_exceeds(trial_H, twice_H, allowance, rate);
    });
```
会在scan_source_changes里被执行。scan_source_changes首先遍历每个网格,计算容许值,这里0.1是硬编码的

```cpp
for (int cell_index = 0; cell_index < n_; ++cell_index)
    {
        const auto& rate = rates[ng_ + cell_index];
        for (double value : rate.data)
            if (!std::isfinite(value)) throw std::runtime_error("Non-finite LTS RK2 source");

        if (!source_cache_ready_)
        {
            const auto& cell = stage.cells[ng_ + cell_index];
            const auto primitive = cell.U.to_Primtivies(cell.gamma);
            const double rho = primitive[hd_core::RHO];
            const double pressure = primitive[hd_core::PRE];
            source_rho_[cell_index] = rho;
            source_density_allowance_[cell_index] = 0.1 * rho;
            source_momentum_allowance_[cell_index] = 0.0;
            source_internal_allowance_[cell_index] = 0.0;
            if (rho != 0.0)
            {
                const double velocity = primitive[hd_core::VEL(0)];
                const double sound_speed = std::sqrt(cell.gamma * pressure / rho);
                source_momentum_allowance_[cell_index] =
                    0.1 * (rho * (std::abs(velocity) + sound_speed));
                source_internal_allowance_[cell_index] =
                    0.1 * (pressure / (cell.gamma - 1.0));
            }
        }

        if (source_rho_[cell_index] == 0.0)
        {
            for (double value : rate.data)
                if (value != 0.0)
                    throw std::invalid_argument("LTS RK2 cannot create gas in vacuum");
            continue;
        }

        bound(source_density_allowance_[cell_index], rate[hd_core::RHO]);
        bound(source_momentum_allowance_[cell_index], rate[hd_core::MOMDEN(0)]);
        bound(source_internal_allowance_[cell_index], rate[hd_core::INTERNAL_ENERGY]);
    }
```
计算出的某个cell的三个容许变化量(密度,动量和内能)会传递给bound函数,在这里就是刚才的source_change_exceeds,scan_source_changes会把计算出的容许变化量和函数变化率传递给它,然后它还会使用一倍的步长和两倍的步长作为参数，这个函数的作用就是判断源项在这一格上产生的变化量是否超过了容许值
```cpp
bool source_change_exceeds(double H, double twice_H, double allowance, double rate)
{
    const double magnitude = std::abs(rate);
    if (twice_H * magnitude < allowance) return false;
    return rate != 0.0 && H > allowance / magnitude;
}
```
这里是一个优化浮点除法的小trick,如果两倍的步长产生的变化都没有超过容许值的话,那么必然没有超过,返回false. 否则再回退到原本的除法,即步长是否大于容许值/变化率,如果大于且变化率不等于0就说明超过，返回true,否则返回false表示不超过.

所以前面source_change_limit内
```cpp
    scan_source_changes(stage, rates, [&](double allowance, double rate) {
        exceeded |= source_change_exceeds(trial_H, twice_H, allowance, rate);
    });
```
的含义也很明显,就是只要有任何一个cell的三个状态变化量中的任何一个超出了容许的变化量allowance,那么就会返回true, exceeded就会永远保持true,只有所有cell的所有容许值都没有被超过,那么exceeded才会保持false。
接下来如果exceeded=true,有容许值被超过了，就返回0.0,这时候上级的advance推进函数会返回false,回到attempt循环里,会尝试缩小一半步长。如果exceeded=false,那么就继续执行。

另一个分支,trial_H<=0.0的情形,只在一种情况发生,就是在run()函数做初始化的时候调用了一次evaluate_source,给的步长是0.0,这个时候的目的是试探一下初始条件下源项限制的步长
```cpp
    const double source_limit = evaluate_sources(start_, time_, 0.0);
    initial_source_ = source_;

    double minimum_local_dt = inf;
    hydro::max_inverse_dt = 0.0;
    for (int cell = 0; cell < n_; ++cell)
    {
        const double rate =
            std::max(speed_[cell], speed_[cell + 1])
            / mesh_.cell_width(0, ng_ + cell);
        if (!std::isfinite(rate)) throw std::runtime_error("Non-finite LTS RK2 CFL rate");
        hydro::max_inverse_dt = std::max(hydro::max_inverse_dt, rate);
        local_dt_[cell] = rate > 0.0 ? cfl_ / rate : inf;
        minimum_local_dt = std::min(minimum_local_dt, local_dt_[cell]);
    }

    macro_step_ = std::min({
        dt_limit_, std::ldexp(minimum_local_dt, max_depth_), source_limit});
```
与后面从cell初状态得到的流体允许的最小步长$\times 2^{\mathrm{depth}}$ 以及外层大循环的dt_limit比较,取最小值得到最开始的macro_step尝试值，因此这里scan_source_changes的bound函数里写的是
```cpp
    scan_source_changes(stage, rates, [&](double allowance, double rate) {
        if (rate != 0.0) limit = std::min(limit, allowance / std::abs(rate));
    });
```
它的作用是在遍历每个cell的过程中,对三个状态(密度,动量和内能)分别计算:

$$\frac{\text{状态的容许变化量allowance}}{\text{源项产生的对应状态的变化率rate}}=\text{该状态上的源项步长限制}$$,

并逐步取更小的步长限制，遍历的最终效果就是取最小的源项步长限制并返回给调用者,参与初始步长的猜测和限制。

