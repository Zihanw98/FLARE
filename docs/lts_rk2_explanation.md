# 多层守恒 RK2 局部时间步（`lts_rk2`）

## 类实例化


n_ : mesh.nactive[0] 一维活跃cell数
ng_ : mesh.nghost[0] 一维ghost cell数
independent_face_count_: 若第一维是周期边界,则独立面数为n_+1，否则为n_(第0个active_cell)
## 入口-run函数


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
- 这里有一个提升邻居等级的trick，把高等级周围(宽度由buffer参数指定)提升到同等级，保证数值稳定性。特别地，周期情形还要取模考虑另一侧
- 最后再遍历一遍，让邻居缓冲过渡，例如3级周围要是2级，2级周围要是1级，缓冲宽度取决于buffer参数，特别地，周期情形还要取模考虑另一侧

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

- 计算本步span以及实际步长macro_step * span, 还有结束时的分数步end_fraction,然后进入RK2的逻辑: 推进本层(level=depth)的cell,需要计算起点处cell守恒量状态U的时间变化率(通过解黎曼问题计算面通量),然后再预测终点处的状态U',再次计算终点处的时间变化率,取平均后乘上时间步dt,然后积分推进
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

在这里复用的意思是,我们是否能够直接使用现在存在flux_数组里的面通量(以及源项贡献的通量,面波速等等,不过主要是通量)?是否能够复用是在后面的代码里决定的,具体来说就是在advance_level的第三个bool参数里
```cpp
        if (!advance_level(depth + 1, fraction, true)) return false;
        if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
```
复用逻辑我们稍后讲解.


如果不复用，那么进入一个build_stage_view函数。这个函数的作用在于,让层级比当前更低（步长更长的慢步）的面用他们之前计算的斜率临时预测到当前递归阶段的起点和终点(RK2需要起点和终点的斜率取平均)，参与后面的通量计算。由于是用斜率进行预测,而不是再次算一次RHS+源项推进，消耗不大。
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
更高层级的cell也不需要预测，因为它们必然在这一层到达fraction之前已经推进到fraction处,与这一层共用起点状态。比如层级0中执行advance_level(1,0),让层级1先从0推进到0.5,在这个过程中层级1会执行两次advance_level(2),把层级2从0推进到0.25,再从0.25推进到0.5。然后层级0再次执行下一行advance_level(1,0.5)(如上,在层级0中这两次advance是紧挨着的),让层级1从0.5推进到1,这个时候层级2以及更高层级必然已经推进到0.5,与层级1共享起点,因此不需要预测。
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


- 接下来,进行一次面的合法性检查，然后临时保存起点处的面通量。这是由于在RK2最后我们需要起点和终点处的通量取平均作为最终的时间导数(斜率)。
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

- 接下来是递归的核心,
```cpp

    if (depth < deepest_)
    {
        if (!advance_level(depth + 1, fraction, true)) return false;    //将更高一层从fraction推进到fraction + span/2
        if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;//将更高一层从fraction+span/2推进到fraction +span/2 + span/2,也就是end_fraction
    }
```
因为反复提及多次了所以不再赘述，简而言之就是调用两次下一层级的步进. 如果步进结果为false那么也立刻返回false(最终会让外面的循环continue,再次尝试缩小一倍macro_step).

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
- 第一个循环,我们对owner为当前层级的所有面,计算RK2的最终推进值integral,即面通量的平均值乘上再乘上面积(即RHS),当前的时间步长(local_step = span * macrostep)。除此以外，我们还要加上左右两侧从更高层级累积下来的通量。