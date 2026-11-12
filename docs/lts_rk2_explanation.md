# 多层守恒 RK2 局部时间步（`lts_rk2`）

## 类实例化


n_ : mesh.nactive[0] 一维活跃cell数
ng_ : mesh.nghost[0] 一维ghost cell数
independent_face_count_: 若第一维是周期边界,则独立面数为n_+1，否则为n_(第0个active_cell)
## 入口


### 第一步:初始扫描
```cpp
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
   这里我们遍历所有活跃区,然后为每一cell通过cfl计算局部dt，找出最小值，为接下来初始大步做准备
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
- needed会在所有低于或等于这个面的owner层级的层级里注册这个面,方便计算通量时复用,具体的复用逻辑下面我举个例子:
```
首先，我们的分层是用递归实现的。每一层对下一层通过fraction实现(读完advance_level再写)(结合GPT，讨论了怎么改变fraction以及递归调用子层的半步fraction，以及上面的复用)
```
- 重置所有记录

## 推进步逻辑
- 这个函数是让所有cell从起点$\mathrm{fraction}$到终点$\mathrm{fraction+ 2^{-depth}}$,同时记录RK2最终计算出来的斜率,作为终点处的斜率临时保存在slope中,便于更低层级处进行预测。
- 整个函数是递归执行，因此会在不同层级不同时刻呼叫这个函数。每一个层级会呼叫两次高一层级。参数 depth 代表本层的层级, fraction代表这次递归调用的起点在长步的几分之几处。而reuse_start_stage则是一个优化开关,表示是否要复用之前的计算结果,这个在后面会解释
- 首先计算本步span以及实际步长macro_step * span, 还有结束时的分数步
- 首先检测是否为复用步，如果不复用，那么进入一个build_stage_view函数。这个函数的作用在于,让层级比当前更低（步长更长的慢步）的面用他们之前计算的斜率临时预测到当前递归阶段的起点和终点(RK2需要起点和终点的斜率取平均)，参与后面的通量计算。由于是用斜率进行预测,而不是再次算一次RHS+源项推进，消耗不大。
示例:假如有三个cellA,B,C，分别处于层级0,1,2.当cell C需要从0.75推进到1.00时,此时B已经推进到0.5(层级1,步长为0.5),因此B用0.5处的斜率(之前计算的RHS+源项,临时保存)预测0.75处,而A则用0处的斜率预测0.75处,得到
$$
   U_A^*(0.75) = U(A,0) + 0.75 * \frac{dU_A}{dt}\Big |_{0}
$$
$$
   U_B^*(B,0.75) = U(B,0.5) + 0.25 * \frac{dU_B}{dt}\Big |_{0.5}
$$
带*号表示这是根据斜率计算的中间临时预测值,而且这个值并不会真正更新到网格中。

那为什么要计算？例如，如果B与C是相邻的网格，那么B、C之间存在相接面,由于RK2需要起点0.75处cell C左右两侧的通量,因此我们需要计算更低层级的cell在0.75处的值来解黎曼问题(理论上如果只有RHS,只需要算B---如果不是周期边界的话,否则也要算A。但有全局源项,比如自引力的话,那么需要计算所有低层级cell的变化)。而B的一个时间步是0.5,无法在其正常的时间步中计算到0.75处的状态,因此我们用0.5处的时间导数进行大致预测。

这里的调用,我们设置的fraction是这一层递归起点处的fratcion，在上面的例子里就是层级2,起点为0.75,我们计算所有层级0和层级1的cell在0.75处的预测值,而endpoint=false表示这里是起点处预测,我们不需要预测层级2的cell在这里的值，我们看后面的调用代码:
```cpp
        if (!advance_level(depth + 1, fraction, true)) return false;
        if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
```
层级1(depth=1)的递归中,调用(2,0.5,true)后再调用(2,0.75,false),因此0.75->1.0这一步调用之前,必然已经完成了0.5->0.75,因此0.75处depth=2的所有cell状态已经计算好，不需要预测。

而后面还有一次调用预测终点1.0处的cell状态,用来计算终点处的RHS和源项，进而计算终点处的斜率，用于最后与起点斜率平均作为RK2最终的推进斜率。
```cpp
    if (!build_stage_view(depth, end_fraction, true)) return false;
```
这里endpoint=true,因为要预测1.0处所有cell的状态,而层级2的所有cell才推进到0.75，因此当前层级也需要预测。

注意:更高层级必然会等待更低层级完成更长步之后才会推进，所以更高层级永远不需要预测，也不需要参与当前层级的黎曼问题求解面通量，而是在推进完成后把更高owner层级的面的通量累积到一个pending_flux数组里，之后直接作为RHS的一项加进给当前层的相关cells。比如上面B的例子里,层级1中,B从0.5推进到1.0的过程,一定会等待层级2,以及更高层级完成推进后再推进,例如C从0.5走到0.75，再从0.75走到1.0。在这个过程中,我们会用一个pending_flux承接B与C之间的面在更高层级里收到的通量之和,更新给cell B完成这一步推进。(通量的更新逻辑:如果还有更高层级,比如层级3里的所有面产生的通量,加到pending_flux里,然后给层级2的所有相邻cell更新状态,然后层级2里所有的面更新产生的通量再次累加到pending_flux[faces_owned]里---不重置之前累加的flux，是叠加---再给层级1的cell进行更新)。

总结:
1. 更低层:一定要预测cell
2. 当前层: 起点 不用预测cell| 终点 需要预测cell
3. 更高层: 一定已经推进到后面,cell不需要预测,而推进过程中更高owner层级的面的通量会从高层级向低层级累积,并更新相邻的低层级cell


- 如果可以复用
```cpp
        if (!advance_level(depth + 1, fraction, true)) return false;
        if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
```
假设我们从depth=0,最低层级的递归出发,最开始fraction=0,也就是所有层级共同的起点,层级0运行完前面的计算逻辑后,让层级1运行fraction=0处的计算走到0.5,起点相同,所以不用重新预测网格(可以和层级0这一步用相同的所有cell状态计算所有层级1的面的通量)。但接着层级0让层级1从0 + 1/2.0 =0.5出发,再走一步0.5到1.0,与自己达到共同起点,这一步的时候层级1就是从0.5出发,与层级0不在一个起点了,这时候就要重新预测0.5处的所有层级cell状态(用刚刚层级0在0处计算的RHS+源项结果)。
- 接下来,解这一层级的所有面的黎曼问题。

```cpp
   solve_faces(view_, faces_owned_by_level_[depth]);
```

每一层级递归只用计算owner为当前层级的面的通量。面的owner层级一定取的是左右两边cell更高的层级,因此两侧cell要么与面owner层级，要么更低。相邻同层级的cell会在这一层递归中直接被该面计算出的RHS更新,而相邻低层级的cell会在高层级递归返回到低层级时，被逐层传递下来的累积pending_flux更新(附加在当前层的RHS和源项上)。

