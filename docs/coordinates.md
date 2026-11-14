#Mesh中各个坐标的转换


## Flat和Index的转换
在mesh下,对于(i,j,k)坐标的cell,其对应顺序存储内cells的第flat个cell,由于x,y,z轴坐标的变化速度顺序由快到慢(x最快)---flat = i + Nx * (j + Ny * k);
```cpp
Index index = {i,j,k};
//flat = i + Nx * (j + Ny * k);
stride = 1
flat = 0
for(int d=0;d < NDIM; d++)
{
    flat += index[d] * stride;
    stride *= ntotal(d);
}
```
其中Nx,Ny,Nz分别是ntotal(0),ntotal(1)和ntotal(2),是每一维度上的active cell数量+ghost cell数量
这里的i,j,k和flat都是包含了ghost的,不需要额外的ghost offset,从flat转index,我们用一个逐层取模的方法
```cpp
remaining = flat
Index index = {};
for(int d=0; d < NDIM; d++)
{
    index[d] = remaining % ntotal(d);
    remaining /= ntotal(d);
}
```
## Active Index和Index的转换
active index是减去ghost offset的,相对于active区域最左下角的cell的坐标,显然
```cpp
  for(int i=0;i<NDIM;i++)
  {
    active_index[i] = index[i] - nghost[i];
  }
```
反向转换同理
## Flat和Active Flat的转换
`flat` 是将全域（包含 ghost cell）的多维 cell 数组展平为一维数组后，每个 cell 对应的索引，从 `0` 开始编号(中文我习惯称呼为展平索引,或者顺序索引)。

展平时，**第 0 维（x）变化最快，第 1 维（y）次之，第 2 维（z）最慢**。也就是先沿 x 方向排完一行，再沿 y 方向换行，排完一层后再沿 z 方向换层。

例如，全域尺寸为：

```cpp
Nx = 3;
Ny = 2;
Nz = 2;
```

(注:这里只是为了举例方便,全域实际上不能这么小。全域包含ghost和active,因此Nx,Ny,Nz必然>2)

三维存储坐标 `(i,j,k)` 对应的索引为：

```cpp
flat = i + Nx * (j + Ny * k);
```

各 cell 的 `flat` 编号如下：

```text
第 0 层：k = 0              第 1 层：k = 1

          i                          i
       0  1  2                    0  1  2
j = 1  3  4  5             j = 1  9 10 11
j = 0  0  1  2             j = 0  6  7  8
```

例如，存储坐标 `(i,j,k) = (1,1,1)` 对应：

```cpp
flat = 1 + 3 * (1 + 2 * 1) = 10;
```

因此，这个 cell 存储在 `mesh.cells[10]` 中。这里的 `i,j,k` 是包含 ghost 偏移的存储坐标，`Nx,Ny,Nz` 分别对应 `ntotal(0)、ntotal(1)、ntotal(2)`。

active flat是仅把活跃域的active cell展平到一维数组后每个active cell对应的index,展平逻辑同理,(中文我习惯称呼为活跃域展平索引/顺序索引)。

active flat通常在遍历中出现,先用类似于flat转index的方法,从active flat 得到active index坐标(注意取模和整除的时候要除的是该维度的nactive[d],而不是ntotal(d)),然后再从active index加上ghost offset得到index,再转flat
```cpp
remaining = active_flat
stride = 1
for(int d=0; d < NDIM; d++)
{
    active_index[d] = remaining % nactive(d);
    remaining /= nactive(d);
    index[d] = active_index[d] + nghost[d];
    flat += index[d] * stride;
    stride *= ntotal(d);  
}
```
逆向转换逻辑同理

## 面的转换
在D维下，我们需要D+1个坐标才能确定一个面,以三维为例，对于一维或二维情形只要去掉更高的维度结论同理。

在三维情形下,第{i,j,k}个cell(全域坐标,包含ghost cell)所拥有的面一共有6个,每个方向各有2个,这两个面我们把它叫做d维法向面,因为其法向在这个维度轴上。这两个面中,那个与(d维度)坐标更低的相邻cell相接的那一个叫做下侧面(lower side face),另外一个与坐标更高的相邻face相接的那一个叫做上侧面(upper side face)。例如x法向的下侧面就是左侧面,其与{i-1,j,k}相邻。

我们通过d,i,j,k这四个坐标数可以唯一确定一个面,约定确定的这个面是下侧面。我们以face_area为例，这个mesh成员存储的是某个面的面积,访问方式是:mesh.face_area[d].at({i,j,k}),例如mesh.face_area[0].at({1,1,1})就表示取{1,1,1}这个cell的左侧面的面积。

由于面不像cell状态有完整的存储对象和数组,为了方便管理我们在build_topology中为**所有至少有一侧是active cell的面**生成一个独一无二的face_id,并以此为索引建立数组face_info。我们把这些面叫做**活跃面(active face)**。注意face_info是不记录那些完全在ghost zone里,两侧都是ghost cell的面的,不过这些面的几何信息仍然可以通过face和face_area访问,它们是包含全部的面的。

建立face_id的方法类似于在一维展平的顺序数组中存储多维cell,用flat_id表示一样,我们优先塞入x方向的面,然后y方向,z方向……
例如 **2×2 个 active cell、非周期边界**：

| 面方向 | 面的数量 | `face_id` |
| --- | --- | --- |
| x 向面 | `(2+1)×2 = 6` | `0～5` |
| y 向面 | `2×(2+1) = 6` | `6～11` |

技术细节可以不用在意,因为face_info和另外一个在build_topology中构建的cell_faces为我们提供了相互转换，访问的方法。

```cpp
    struct FaceInfo
    {
        int direction = 0;
        // Includes ghost offsets; usable with face_area[direction].at(...)
        // and reconstruction.face(direction, ...). A periodic seam uses
        // the lower boundary coordinate, including for the last cell's upper face.
        Index storage_coordinates{};
        // Compact active IDs. Physical boundary exteriors have no active cell.
        std::size_t left_active = invalid_active_index;
        std::size_t right_active = invalid_active_index;
    };
    std::vector<FaceInfo> face_info;

    struct CellFacePair
    {
        std::size_t lower = 0; ///< Face ID on the negative side of this cell.
        std::size_t upper = 0; ///< Face ID on the positive side of this cell.
    };
    std::vector<std::array<CellFacePair, NDIM>> cell_faces;
```
`face_info` 通过 `face_id` 索引，记录每个独立单元面的法向方向 `direction`、面存储坐标 `storage_coordinates`，以及两侧相邻 active cell 的编号 `left_active` 和 `right_active`。

`storage_coordinates` 是面在全域存储网格中的多维索引，包含 ghost 偏移。它采用**正侧邻接 cell 的低侧面**这一索引约定：对于内部面，面坐标与其正侧邻接 cell 的存储坐标相同。这里的<span style="color: red;">“正侧”</span>指沿 `direction` 坐标增大的方向；一维中就是右侧。

以**非周期边界、存在 ghost cell 的一维网格**为例：

| 面的位置 | `storage_coordinates[0]` | 正侧（右侧）邻接 cell |
| --- | --- | --- |
| 活跃区的低边界面 | `nghost[0]` | 第一个 active cell |
| 活跃区的高边界面 | `nghost[0] + nactive[0]` | 右侧第一个 ghost cell |

左侧 ghost cell 的存储坐标为 `0, 1, ..., nghost[0]-1`。因此，第一个 active cell 的坐标是 `nghost[0]`，其左侧面就是活跃区的低边界面。同理，右侧第一个 ghost cell 的左侧面就是活跃区的高边界面。

需要注意，`storage_coordinates` 本质上仍然是**面坐标**。即使没有 ghost cell，高边界面也仍然存在，只是其正侧没有对应的存储 cell。

`left_active` 和 `right_active` 分别记录面负侧、正侧邻接 cell 的 **active flat 索引**，即仅在活跃区内展平后的 cell 编号，不包含 ghost。若某一侧没有 active cell，例如非周期边界的外侧，则对应字段为 `invalid_active_index`。

对于周期边界，两端共用一个 `face_id`，`storage_coordinates` 统一采用低边界面的坐标，而 `left_active` 和 `right_active` 分别连接该方向末端和起始端的 active cell。

`cell_faces` 保存了**每个 active cell 在各个方向上的两侧面编号**：

```cpp
std::vector<std::array<CellFacePair, NDIM>> cell_faces;
```

它的大小等于 active cell 总数：

```cpp
cell_faces.size() == active_count;
```

访问时使用 **active flat 索引**，即 active 区域内展平后的 cell 编号，不包含 ghost cell：

```cpp
cell_faces[cell_id]
```

每个 `cell_faces[cell_id]` 包含 `NDIM` 个 `CellFacePair`。其中，第 `d` 个（$0 \le d < \mathrm{NDIM}$）记录该 cell 在方向 `d` 上的两个面编号：

- `lower`：负方向一侧（低侧）的 `face_id`。
- `upper`：正方向一侧（高侧）的 `face_id`。

例如：

```cpp
cell_faces[cell_id][0].lower  // 该 cell 的 x− 面编号
cell_faces[cell_id][0].upper  // 该 cell 的 x+ 面编号
```

因此我们可以利用`cell_faces`取得想要的某个面的`face_id`(需要的信息:包含这个面的cell的active flat,面的法向,以及面在高侧还是低侧)。取得 `face_id` 后，就可以访问对应的面信息：

```cpp
const auto face_id = cell_faces[cell_id][d].lower;
const auto& info = face_info[face_id];
```

相邻 cell 的公共面使用同一个 `face_id`，因此这些面不是某个 cell 独有的。

注：LTS_RK2算法中的hints_和flux_的索引布局与此处的都不一样，请见LTS_RK2的解释文件
## 面的数量，法向轴上的法向面
对于**法向为 `d` 的面场**，其在第 `d` 个方向上的尺寸比 active cell 数多 `1`，其余方向的尺寸与 active cell 数相同。

原因是：沿某个方向排列的 `N` 个 cell，需要 `N+1` 个面位置来划分，包括两端的边界面，以及中间的 `N-1` 个内部面。

例如，一行有 3 个 cell：

```text
面0         面1         面2         面3
 |  cell 0   |  cell 1   |  cell 2   |
```

这 3 个 cell 对应 4 个面位置。


### 展平面场
这里的**展平面场**，指一种**按法向方向分组，再将各方向的面数据分别展平为一维数组**的存储方式，通过两层索引访问：

`flat_facefield[d][index]`

- `d` 表示面法向所在的坐标轴，`0、1、2` 分别对应 x、y、z。
- `flat_facefield[d]` 保存该方向面场中所有面的数据，并按约定的展平顺序排列。
- `index` 是方向 `d` 的面场内部的一维索引；`flat_facefield[d][index]` 访问其中一个具体单元面的数据。

需要注意，`index` 由面的**完整多维坐标**展平得到，包含法向和切向位置的信息。因此，即使两个面的法向坐标相同，只要切向位置不同，它们也对应不同的 `index`。

以 **3×2 个 active cell、非周期边界的二维网格**为例：

```cpp
nactive[0] = 3;  // x 方向有 3 个 cell
nactive[1] = 2;  // y 方向有 2 个 cell
```

对于 **x 向面（`d=0`）**：

- 沿法向 x，有 `3+1=4` 个面位置。
- 沿切向 y，有 2 行 cell，因此每个 x 面位置对应 2 个面片。
- 面场尺寸为 `4×2`，共 8 个面。

对于 **y 向面（`d=1`）**：

- 沿法向 y，有 `2+1=3` 个面位置。
- 沿切向 x，有 3 列 cell，因此每个 y 面位置对应 3 个面片。
- 面场尺寸为 `3×3`，共 9 个面。

| 面场 | 法向 | x 方向尺寸 | y 方向尺寸 | 面的总数 |
| --- | --- | --- | --- | --- |
| x 向面场 | `d=0` | `nactive[0]+1 = 4` | `nactive[1] = 2` | `4×2 = 8` |
| y 向面场 | `d=1` | `nactive[0] = 3` | `nactive[1]+1 = 3` | `3×3 = 9` |

因此，构造方向 `d` 的面场尺寸时，可以先复制 active cell 的尺寸，再仅将法向对应的分量加 `1`：

```cpp
auto shape = mesh.nactive;
++shape[d];
```

这里描述的是 `hints_` 等数组采用的**平展 active 面场布局**。即使采用周期边界，这种布局也保留两端的面位置；而 `face_info` 记录的是独立面，会将周期边界两端合并为同一个 `face_id`。


### face: 存储面的分界坐标的数组

`mesh.face` 保存的是**各个坐标轴上的分界位置，即面所在的物理坐标**：

```cpp
std::array<std::vector<double>, NDIM> face;
```

通过两个下标访问：

```cpp
mesh.face[d][i]
```

- `d`：坐标轴方向，`0、1、2` 分别对应 x、y、z。
- `i`：该方向上的面位置索引，包含 ghost 偏移。
- 返回值：这个面位置的物理坐标，例如 `x_i`、`y_i` 或 `z_i`。

**`i` 是整数索引，`face[d][i]` 是物理坐标，两者含义不同。**

例如，忽略 ghost，考虑一个 **3×2 个 cell 的二维网格**：

```cpp
mesh.face[0] = {0.0, 1.0, 2.0, 3.0};  // x 方向的分界坐标
mesh.face[1] = {0.0, 1.0, 2.0};       // y 方向的分界坐标
```

对应的网格为：

```text
y=2  +-------+-------+-------+
     |       |       |       |
y=1  +-------+-------+-------+
     |       |       |       |
y=0  +-------+-------+-------+
    x=0     x=1     x=2     x=3
```

其中：

```cpp
mesh.face[0][1] == 1.0;
```

表示整条分界线 `x=1` 的位置。这条线经过两行 cell，因此包含上下两个单元面片。两个面片的 x 坐标相同，都使用 `face[0][1]`，但它们连接不同的 cell，可以具有不同的通量。

因此，**`mesh.face` 保存各轴的分界坐标，不为每个单元面片重复保存坐标**：

| 坐标数组 | 坐标数量 | 对应的单元面片数量 |
| --- | --- | --- |
| `face[0]` | 4 个 x 分界位置 | 每个位置分为 2 段，共 `4×2=8` 个 x 向面 |
| `face[1]` | 3 个 y 分界位置 | 每个位置分为 3 段，共 `3×3=9` 个 y 向面 |

三维中也是同样的关系：`face[0][i]` 确定平面 `x=x_i`，不同的 `(j,k)` 则对应这个平面上的不同面片。

包含 ghost 时，每个方向的坐标数组长度为：

```cpp
mesh.face[d].size() == mesh.ntotal(d) + 1;
```

活跃区的两个边界坐标分别是：

```cpp
mesh.face[d][mesh.nghost[d]]                     // 活跃区低边界
mesh.face[d][mesh.nghost[d] + mesh.nactive[d]]    // 活跃区高边界
```

对于法向为 `d`、完整存储坐标为 `I` 的一个具体单元面：

```cpp
mesh.face[d][I[d]]        // 该面在法向上的物理坐标
mesh.face_area[d].at(I)   // 该面片的面积
```

前者只需要法向分量 `I[d]`；后者需要完整坐标 `I`，以确定具体是哪一个面片。

需要注意，`mesh.face` 虽然也采用 `face[d][i]` 的两层索引形式，但它保存的是**各坐标轴上的一维分界坐标，而不是展平后的面场数据**。

由于网格与坐标轴对齐，对于固定的方向 `d` 和法向索引 `i`，位于不同切向位置的单元面片具有相同的 `d` 方向物理坐标。这个坐标只需在 `face[d][i]` 中保存一次，无须为每个面片重复存储。

因此，`face[d]` 只包含 `ntotal(d)+1` 个分界坐标，不会将所有法向为 `d` 的单元面片逐一展平存储。(后续实现非对齐的mesh时,需要将face改为平展面场)