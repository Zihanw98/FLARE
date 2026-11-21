# Conservative Multilevel RK2 Local Time Stepping (`lts_rk2`)

[中文版](lts_rk2_explanation.md)

This document accompanies the [1D reference code](lts_rk2.cpp). It explains how LTS-RK2 assigns time levels, constructs stage states, advances levels recursively, and preserves conservation across level interfaces through shared face-flux integrals. The [current multidimensional implementation](../src/time/lts_rk2.cpp) follows the same design; differences in indexing, neighbor traversal, and CFL checks are identified where relevant.

The main sequence is: **choose a macro step → assign cell levels and face ownership → prepare the start stage → advance finer levels recursively → finalize the current level → synchronize or retry**.

## 1. Time levels and notation

Let a macro step begin at physical time $t_0$ and have duration $H$. The local step size at level $\ell$ is

$$
h_\ell = \frac{H}{2^\ell}.
$$

A higher level uses a shorter step and is called a finer level. A lower level uses a longer step and is called a coarser or slower level. Levels remain fixed within a macro-step attempt and are reassigned on every retry.

| Code variable | Meaning |
| --- | --- |
| `macro_step_` | Macro-step size $H$ |
| `depth` | Current recursion level $\ell$ |
| `fraction` | Local-step start position $f$ within the macro step, corresponding to $t_0+fH$ |
| `span` | Fraction of the macro step covered by this local step, $2^{-\ell}$ |
| `local_step` | Local step size $h_\ell=H\,\mathrm{span}$ |
| `end_fraction` | Local-step endpoint $f+2^{-\ell}$ |
| `max_depth_` | Largest permitted level index, equal to `lts_max_levels - 1` |
| `deepest_` | Highest level actually used in this attempt |

Levels are numbered from 0. Thus, `lts_max_levels = 8` permits levels 0–7 and a maximum step-size ratio of $2^7=128:1$. A ratio of $256:1$ requires `max_depth_ = 8`, or nine levels in total.

Throughout this document, $U_i$ denotes a cell's conserved state per unit volume, and $F_f$ denotes a numerical face flux. Face areas, cell volumes, and time integration connect these quantities: a face flux is not itself a cell's time derivative.

## 2. Working arrays and indexing

### 2.1 Initial state, working state, and stage view

| Array | Stored quantity |
| --- | --- |
| `start_` | Complete grid snapshot at the macro-step start, unchanged across retries |
| `state_` | Each cell's state after its most recently completed local RK2 step; cells may represent different times |
| `view_` | Approximate grid state at one common physical time, constructed for a stage evaluation; copied and predicted in `build_stage_view()` |
| `base_[i]` | Cell $i$'s state $U_i^0$ at the start of its current local step |
| `slope_[i]` | Start-stage time derivative $k_{1,i}$ used for linear prediction |
| `local_step_start_[i]` | Start fraction of that local step |
| `source_start_[i]` | Source term at the local-step start |
| `pressure_work_start_[i]` | Internal-energy pressure work at the local-step start |
| `pending_flux_integral_[i]` | Net face-flux integral accumulated for cell $i$ but not yet applied |

`state_` stores states, not rates. `slope_` stores the start-stage derivative, not the average of the two RK2 derivatives. Linear prediction modifies only `view_`; it does not directly update `state_` or the simulation grid.

At the start of every macro-step attempt, `reset_attempt()` executes `state_.cells = start_.cells`. Phase I copies the current level's local-step starting states from `state_` into `base_` and computes `slope_`, without advancing those cells in `state_`. During recursion, finer levels update their own cells in their own Phase II, so different levels can represent different times.

When Phase II finalizes the current level, `auto& cell = state_.cells[flat]` refers directly to a cell in the trial array; assigning `cell.U` advances that cell to its local-step endpoint. Only after the entire macro step succeeds does `mesh_.cells.swap(state_.cells)` commit the result to the simulation grid. If the attempt fails, the next `reset_attempt()` restarts from `start_`.

Each `build_stage_view()` call first executes `view_.cells = state_.cells`, then extrapolates selected cells only in `view_`: coarser levels for a starting-stage view, and both current and coarser levels for an endpoint view. Finer cells that have already reached the target time retain their copied states.

Face fluxes and sources are evaluated from this `view_` at the common target time. A finer level reusing the start stage reads the existing `flux_`, `speed_`, `face_velocity_`, and `source_` results, which already include the coarser cells' predicted states at that start time. Coarser cells in `state_` may remain at their own local-step starts.

### 2.2 Face quantities and the star-pressure cache

`flux_`, `speed_`, and `face_velocity_` store face fluxes, CFL wave speeds, and face velocities, respectively. `first_flux_` retains each face's start flux for the current local step of its owning level, protecting it from subsequent updates.

`hints_` stores the star-region pressure $p_*$ from the previous Riemann solve as an initial guess for the next solve. A value of `0` indicates that no usable pressure guess is available. In multiple dimensions, it is accessed as

```cpp
hints_[d][face_flat]
```

Here, `d` is the face-normal direction, and `face_flat` is the flattened index within that direction's active face field, including tangential position. The index is formed by subtracting the ghost offset in each direction and flattening with the face-field dimensions: `nactive[d] + 1` in the normal direction and the corresponding `nactive` in every other direction.

For example, a $2\times2$ active grid has a $3\times2$ x-face field. With `(i, j)` denoting face coordinates after removing ghost offsets,

```cpp
face_flat = i + 3 * j;
```

| Face patch | Face coordinates | `face_flat` | Cache entry |
| --- | --- | ---: | --- |
| Patch in the lower row | `(1, 0)` | 1 | `hints_[0][1]` |
| Patch in the upper row | `(1, 1)` | 4 | `hints_[0][4]` |

The two patches have the same x coordinate but connect different cells, so they require separate cache entries.

The following face-index layouts must be distinguished:

| Object | Index layout |
| --- | --- |
| `hints_[d][face_flat]` | Active face field grouped by direction |
| Conventional Godunov `flux[d].faces[face_flat]` | Active face-field layout corresponding to the pressure cache |
| Multidimensional LTS `flux_[face_id]`, `speed_[face_id]`, `face_velocity_[face_id]` | One numbering of independent faces across all directions, matching `mesh.face_info[face_id]` |
| `mesh.face_area[d].faces[flat]` | Full face field including ghost regions |

For $n$ active cells in 1D, nonperiodic boundaries give $n+1$ independent faces, whereas periodic boundaries give $n$. The older 1D face arrays still contain $n+1$ slots, but the two periodic endpoints represent one face, which is solved and accumulated only once. The multidimensional version obtains the number of independent faces directly from `mesh.face_info.size()`.

The pressure cache retains both periodic endpoint slots, while `face_info` merges them into one independent face. Consequently, an LTS `face_id` cannot be used directly as a `hints_` `face_flat`: the solver converts the index and synchronizes the periodic pressure-cache entries.

## 3. `run()`: choosing a macro step and preparing an attempt

### 3.1 Initial scan

First, solve every independent face at the synchronized start and evaluate the source terms:

```cpp
solve_faces(start_, all_faces_);
initial_flux_ = flux_;
initial_speed_ = speed_;
initial_velocity_ = face_velocity_;
initial_hints_ = hints_;
const double source_limit = evaluate_sources(start_, time_, 0.0);
initial_source_ = source_;
```

These results determine the step size and are reused by the initial RK stage. The `initial_*` snapshots support restoration after a failed attempt.

In 1D, cell $i$'s CFL rate and permitted step size are

$$
r_i=\frac{\max(a_{i-1/2},a_{i+1/2})}{\Delta x_i},
\qquad
\Delta t_i^{\mathrm{CFL}}=
\begin{cases}
\mathrm{CFL}/r_i, & r_i>0,\\
\infty, & r_i=0.
\end{cases}
$$

Here, $a_f$ corresponds to `speed_[f]`. The current multidimensional version sums contributions from all coordinate directions:

$$
r_i=\sum_d\frac{\max(a_{i,d,-},a_{i,d,+})}{\Delta x_{i,d}}.
$$

Scanning all active cells produces `local_dt_[i]` and their minimum, $\Delta t_{\min}$.

### 3.2 Choosing the macro step

```cpp
macro_step_ = std::min({
    dt_limit_, std::ldexp(minimum_local_dt, max_depth_), source_limit});
```

Equivalently,

$$
H=\min\left(\mathrm{dt\_limit},\;2^{\ell_{\max}}\Delta t_{\min},\;H_{\mathrm{source}}\right).
$$

The three terms are the externally supplied step-size cap, the hydrodynamic step allowed by the available levels, and the source-term limit. The external `dt_limit_` can prevent stepping beyond a prescribed end time.

### 3.3 Reset, advance, and retry

Each attempt first checks that $H$ can advance time, then calls `reset_attempt()`. This reassigns cell levels and face ownership, restores the working state and stage caches to the synchronized start, and clears pending face-flux integrals. It then calls

```cpp
advance_level(0, 0.0, true);
```

If any stage returns `false`, failure propagates through the recursive calls and invalidates the entire macro step. The outer loop halves $H$ and retries from the synchronized start. Once all levels succeed, it fills the endpoint boundaries and commits `state_` and the pressure cache to the simulation grid and external cache. Exceeding the retry limit, or encountering an exception that is thrown directly, produces an error.

## 4. Assigning cell levels and face ownership

### 4.1 `assign_time_levels()`: satisfy CFL, then build buffers

For each cell, start at level 0 and progressively shorten the step until

$$
H/2^\ell\leq\Delta t_i^{\mathrm{CFL}}.
$$

Two promotion passes then establish the level layout:

1. **Extend a buffer at the required level.** For the cell's CFL-required `depth`, call `promote_neighbors(cell, depth)` to raise the central cell and its buffered neighbors to at least that level.
2. **Smooth level transitions.** Work downward from the highest level. For every cell currently at `depth`, call `promote_neighbors(cell, depth - 1)` so that its neighbors reach at least the next coarser level.

Both operations use `std::max`: they only raise levels and never lower an existing level. Adjacent cells end up differing by at most one level; wider buffers also constrain more distant neighbors.

The buffer width follows the reconstruction stencil:

| Reconstruction | `level_buffer_` | Reason |
| --- | ---: | --- |
| PCM | 1 | A face Riemann problem requires states from both adjacent cells |
| PLM / PLM_corrected | 2 | Reconstructing those states also requires their outer neighbors for spatial slopes |

For example, the right face of cell `i` requires reconstructed states from `i` and `i+1`; the slope in `i+1` also depends on `i+2`.

### 4.2 Traversing neighbors in multiple dimensions

The multidimensional `promote_neighbors()` walks outward along the lower and upper sides of each coordinate axis, up to `level_buffer_` cells:

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

Each new direction starts again from the central cell. Face topology automatically connects the periodic endpoints. At a nonperiodic boundary, a missing active neighbor gives a negative index and stops the walk.

For a 3D central cell `(i,j,k)` with buffer width 2:

| Direction | First neighbor | Second neighbor |
| --- | --- | --- |
| Lower / upper x | `(i±1,j,k)` | `(i±2,j,k)` |
| Lower / upper y | `(i,j±1,k)` | `(i,j±2,k)` |
| Lower / upper z | `(i,j,k±1)` | `(i,j,k±2)` |

A single call directly visits neighbors on the three coordinate axes, excluding diagonal cells. Away from boundaries and repeated neighbors, width 1 covers the center plus six neighbors, and width 2 covers the center plus twelve. Repeated promotions and smoothing can enlarge the overall promoted region.

### 4.3 `assign_faces_to_levels()`: the finer level owns the integral

Each face is assigned to the higher level of its two adjacent active cells:

$$
\ell_f=\max(\ell_{\mathrm{left}},\ell_{\mathrm{right}}).
$$

At a nonperiodic boundary, the missing side has a negative cell index, and the face takes the level of the valid side.

| Level table | Entries at `depth` |
| --- | --- |
| `cells_by_level_[depth]` | Cells whose level equals `depth` |
| `faces_owned_by_level_[depth]` | Faces whose owning level equals `depth` |
| `faces_needed_by_level_[depth]` | Faces whose owning level is greater than or equal to `depth` |

Each table has `deepest_ + 1` level slots. `owned` determines which level integrates a face flux. `needed` supports the current level's start slopes and prepares fluxes for finer levels that share the same start.

`needed` may include faces that do not touch a cell at the current level, because subsequent finer levels also reuse these start-stage results. Conversely, a face owned below `depth` cannot touch a cell at `depth`, so the current level does not need it.

## 5. `build_stage_view()`: placing different levels at one time

When a fine level evaluates a face flux, a neighboring slow cell may not yet have completed its own local step. Its older state in `state_` cannot be used directly. Instead, the method constructs a linear prediction:

$$
\widetilde U_i(t)=U_i^0+(t-t_i)k_{1,i}.
$$

Here, $t_i$ is the start of that cell's current local step. In fractional-time notation, this becomes

```cpp
const double elapsed =
    (fraction - local_step_start_[cell_index]) * macro_step_;
cell.U = base_[cell_index] + slope_[cell_index] * elapsed;
```

For example, let A, B, and C belong to levels 0, 1, and 2, respectively, with C about to advance from $3H/4$ to $H$. A's local step still begins at 0, while B's current local step begins at $H/2$. Their predictions are

$$
\widetilde U_A(3H/4)=U_A(0)+\frac{3H}{4}k_{1,A}(0),
$$

$$
\widetilde U_B(3H/4)=U_B(H/2)+\frac{H}{4}k_{1,B}(H/2).
$$

The common physical start time $t_0$ is omitted here. This level layout is schematic and ignores further promotions that buffers may require. C has completed its previous fine step and uses its RK2-updated state at $3H/4$ directly.

The function first copies `state_` into `view_`, then overwrites cells requiring prediction according to the following rule:

| Stage | Coarser levels, `level < depth` | Current level, `level == depth` | Finer levels, `level > depth` |
| --- | --- | --- | --- |
| Start, `endpoint=false` | Predict to the start time | Use updated state | Use updated state |
| End, `endpoint=true` | Predict to the endpoint | Predict the endpoint using the start slope | Use state already updated to the endpoint |

The corresponding condition is

```cpp
const bool needs_prediction =
    level_[cell_index] < depth || (endpoint && level_[cell_index] == depth);
```

Why can finer levels use `state_` directly? The recursion ensures that the current level and all finer levels share the start when this local step begins. Before the current level constructs its endpoint stage, finer levels have completed every substep covering the interval.

After prediction, the function checks physical validity and fills boundaries at the target physical time. The resulting `view_` is an approximate state at one common time, suitable for both face reconstruction and sources that depend on the full grid, such as self-gravity. PLM must recompute spatial slopes on this stage view instead of combining neighboring states from different times.

## 6. `advance_level()`: one local RK2 step

The call

```cpp
advance_level(depth, fraction, reuse_start_stage);
```

advances current-level cells from $t_0+fH$ to $t_0+(f+2^{-\ell})H$ and recursively advances all finer levels to the same endpoint. It consists of start-stage preparation, finer-level recursion, and current-level finalization.

### 6.1 First stage: prepare start fluxes, sources, and prediction slopes

If the start stage cannot be reused, rebuild the stage view, solve the needed faces, and evaluate sources:

```cpp
if (!reuse_start_stage)
{
    if (!build_stage_view(depth, fraction, false)) return false;
    solve_faces(view_, faces_needed_by_level_[depth]);
    if (macro_step_ > evaluate_sources(
            view_, time_ + fraction * macro_step_, macro_step_)) return false;
}
```

Next, check wave speeds and retain the start fluxes of the faces **owned by this level**:

```cpp
if (!check_wave_speeds(faces_needed_by_level_[depth])) return false;
for (int face : faces_owned_by_level_[depth])
    first_flux_[face] = flux_[face];
```

For each current-level cell, save `base_`, `local_step_start_`, and `source_start_`, then compute the start slope. In 1D,

$$
k_{1,i}=\frac{A_{i-1/2}F_{i-1/2}^{(0)}-A_{i+1/2}F_{i+1/2}^{(0)}}{V_i}
+S_i^{(0)}+Q_i^{(0)}.
$$

$Q_i$ is a pressure-work term whose only nonzero component is the separately evolved internal energy. The code uses

$$
Q_{i,e}=-(\gamma_i-1)e_i\,
\frac{A_{i+1/2}v_{i+1/2}-A_{i-1/2}v_{i-1/2}}{V_i}.
$$

Here, $e_i$ is the separately evolved internal-energy density, and $v_f$ is the face velocity returned by the Riemann solver. The multidimensional version sums flux differences and area-weighted face-velocity differences over all directions.

The resulting $k_{1,i}$ is stored in `slope_[i]` for both the current level's Euler endpoint prediction and the slow-cell predictions needed by finer levels. The code also checks the physical validity of $U_i^0+h_\ell k_{1,i}$. Current-level cells have not yet been updated to the endpoint.

These saves, slope calculations, and wave-speed checks also run when the start stage is reused. The cell loop visits only `cells_by_level_[depth]`, so `base_` and starting pressure come from current-level cells in `state_` that have already reached `fraction`. Neighboring coarser cells contribute through the face quantities and sources evaluated from `view_`. The current level still saves its own `first_flux_` and `source_start_`, and computes its own `pressure_work_start_` and `slope_`.

### 6.2 Start-stage reuse means reusing a stage evaluation

`reuse_start_stage=true` means the required face fluxes, wave speeds, face velocities, and sources have already been evaluated for the same start stage and can be reused.

| Reused array | Use at the current level |
| --- | --- |
| `flux_` | Save owned faces' start fluxes and compute cell flux divergence |
| `speed_` | Check the CFL wave-speed limit on needed faces |
| `face_velocity_` | Compute velocity divergence for starting pressure work |
| `source_` | Save `source_start_` for the starting slope and final two-stage integral |

- The first macro-step call reuses the initial scan in `run()`.
- The first finer half-step shares its parent's start and reuses the results prepared there.
- The second finer half-step begins at a new time and must reconstruct the stage view and evaluate it again.

The caches remain valid even when coarser cells in `state_` have not reached this start time: those cells were first predicted in `view_` when the caches were prepared. Reuse skips rebuilding the starting view and reevaluating face quantities and sources. It does not require every cell in `state_` to be synchronized, and the current level still computes its starting slope.

For example, let the macro step begin at 0 with levels 0, 1, and 2. On entry to `advance_level(1, 0.5, false)`, level 0's `state_` is still at 0, while levels 1 and 2 are at $H/2$. Level 1 first uses `build_stage_view()` to predict level 0 to $H/2$ in `view_`, then evaluates face quantities for `owner >= 1` and sources. When it subsequently calls `advance_level(2, 0.5, true)`, these caches already represent the $H/2$ stage, including the coarser cells' predictions, so level 2 can reuse them. Level 2 takes its own `base_` and pressure from its own `state_`, which has reached $H/2$.

Between preparing the parent's start stage and entering the first finer half-step, the parent only saves local starting states and computes slopes; it neither advances cells nor overwrites these stage caches. The parent's `needed` face set also contains every face needed by the child. The cached time, states, and face coverage therefore satisfy the reuse requirements.

Likewise, an endpoint flux from the previous substep generally cannot serve as the next substep's start flux. The former uses the current level's Euler-predicted state, whereas the latter requires its RK2-corrected state. Equal times do not imply equal states.

<a id="recursion-feature-1"></a>

### 6.3 Recursion: let finer levels complete two half-steps first

```cpp
if (depth < deepest_)
{
    if (!advance_level(depth + 1, fraction, true)) return false;
    if (!advance_level(depth + 1, fraction + span / 2.0, false)) return false;
}
```

The first call advances all finer levels to $f+\mathrm{span}/2$; the second advances them to $f+\mathrm{span}$. Each child call repeats the same two-half-step pattern, so finer levels reach the current call's endpoint before the current level does.

The first call's `true` relies on the start stage already prepared by the parent. Coarser cells' predicted states have been included in that evaluation even if those cells have not advanced. The second call starts at a new time, with finer cells now updated by RK2, so it uses `false` to rebuild the starting view and caches.

Thus, when the current level begins its second stage, all finer levels have reached `end_fraction`, while the current level still retains its local-step start state.

Even if the current level has no cells, finer-level recursion must finish before

```cpp
if (cells_by_level_[depth].empty()) return true;
```

can skip current-level finalization.

### 6.4 Second stage: construct the endpoint view and solve only owned faces

```cpp
if (!build_stage_view(depth, end_fraction, true)) return false;
solve_faces(view_, faces_owned_by_level_[depth]);
if (!check_wave_speeds(faces_owned_by_level_[depth]) ||
    macro_step_ > evaluate_sources(
        view_, time_ + end_fraction * macro_step_, macro_step_)) return false;
```

The endpoint view uses linear predictions for the current and coarser levels, and completed substep states for finer levels. Sources are also reevaluated on this endpoint view.

Here, `build_stage_view()` updates only `view_`; current-level cells in `state_` still hold their local-step starting states. The subsequent face loop also only accumulates flux integrals. The actual write to current-level cells in `state_` occurs in the cell-finalization loop in Section 6.6.

Only `owned` faces are solved here. Faces owned by finer levels have already been integrated during their substeps, with contributions stored in the adjacent cells' accumulation arrays. Integrating them again at the current level would count the same transport twice.

`flux_` is a working cache updated face by face; it does not always represent the entire grid evaluated from one common state. The start-stage `needed` solve and reuse rules supply valid fluxes for the current slopes, while the endpoint solve updates only current-level owned faces.

### 6.5 Shared face-flux integrals

For each face $f$ owned by this level, compute the trapezoidal integral

$$
I_f=A_f\frac{h_\ell}{2}\left(F_f^{(0)}+F_f^{(1)}\right).
$$

Write that same integral with opposite signs to the adjacent active cells:

```cpp
if (left >= 0) pending_flux_integral_[left] -= integral;
if (right >= 0) pending_flux_integral_[right] += integral;
```

A positive flux travels from left to right, decreasing the left cell and increasing the right cell. The two contributions from an internal face cancel in the volume-weighted domain sum. A boundary face contributes only to its existing active neighbor.

When a face connects a coarse cell and a fine cell, the fine cell applies its accumulated contribution after every substep. The coarse cell retains its contribution until its longer local step ends. For example, with a coarse step $H$ and fine steps $H/2$, the coarse cell receives

$$
I_f^{[0,H]}=I_f^{[0,H/2]}+I_f^{[H/2,H]}.
$$

These integrals are stored directly in the shared `pending_flux_integral_` array; they are not copied down the recursion stack on return. Recursion controls when contributions are produced and when each cell applies them.

### 6.6 Finalizing current-level cells

Once every adjacent face contribution is available, update each current-level cell as

$$
U_i^{\mathrm{new}}=U_i^0
+\frac{J_i}{V_i}
+\frac{h_\ell}{2}\left(S_i^{(0)}+S_i^{(1)}\right)
+\frac{h_\ell}{2}\left(Q_i^{(0)}+Q_i^{(1)}\right),
$$

where $J_i=\mathrm{pending\_flux\_integral\_}[i]$ is the signed net face-flux integral. $J_i/V_i$ is already a state increment and must not be multiplied by the time step again. Pressure work contributes only to the separately evolved internal energy.

Both sources and pressure work use the average of the start and endpoint stages. With one level, this gives the usual Heun/RK2 update; across levels, interface transport uses the sum of the finer substep integrals.

This update is written into the current-level cells in `state_`. The corresponding assignments in the current multidimensional implementation are:

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

Before the pressure-work addition, `pressure_work_end` is computed using pressure from the endpoint `view_` and endpoint face velocities. Because `cell` is a reference, these assignments directly modify `state_`; the predicted `view_` is not copied back wholesale. The loop visits only `cells_by_level_[depth]`, advancing the current level to `end_fraction`, retaining completed finer-level updates, and leaving coarser cells for their own Phase II.

After updating, check physical validity, synchronize the dual-energy variables, and clear `pending_flux_integral_` only for **current-level cells**. Coarser cells retain their accumulated contributions until their own local steps end.

Shared face integrals ensure discrete conservation of transport. Sources can change domain-integrated conserved quantities, and dual-energy synchronization can correct total energy. These effects are separate from the cancellation of internal-face contributions.

## Process summary

0. The outer run function performs all the various initialization steps, including initializing arrays, assigning time steps, assigning levels to cells and faces, and calculating the maximum step size macro_step (the smallest feasible dt among all cells * $2^{\mathrm{max_depth}}$). It then begins the attempt loop.
1. First, the outer advance_lts_rk2 computes the flux at the initial starting point once, then calls advance_level(0,0,true). This is the starting point. Before calling `advance_level` on each attempt, `reset_attempt()` restores the macro-step starting state with `state_.cells = start_.cells`; local updates from a failed attempt do not carry into the next attempt.
2. Each level's advance_level has a level parameter depth, a starting-point parameter fraction (the fractional time step), and a reuse flag. It calculates the current fractional step size span = $2^{-\mathrm{depth}}$ and the endpoint end_fraction = fraction + span. It first checks whether this is a reuse step. If not, it calls `build_stage_view(depth, fraction, false)`, which copies the cell states produced by RK2 advancement from state_ into the temporary array view_. Because the step sizes differ, the states of different cells in state_ correspond to different times. Using the previously stored slopes from the first part of RK2, it extrapolates the states in view_ of all lower-level cells (which have longer step sizes and have not reached this level's starting point) to the current level's start fraction. Because of the recursion, cells at the current and higher levels must already have reached the current start fraction, so their copied values from state_ do not need modification. The entire view_ is therefore synchronized. The cell states in view_ are then used to solve the Riemann problems and obtain flux_ for all faces owned by the current or higher levels, because only faces owned by the current or higher levels can belong to cells at the current level and affect those cells through their fluxes.

   A reuse step means the required starting-stage flux_, speed_, face_velocity_, and source_ results are already available, so build_stage_view, solve_faces, and evaluate_sources are skipped. Lower-level cells in state_ may still represent earlier times, but their predictions at the current start time have already contributed through the previously prepared view_. Those predictions and evaluations need not be repeated for this level; its own starting states must still be saved and its slopes computed.
3. Save the starting flux_ of faces owned by this level into first_flux_. The incoming-minus-outgoing sum of area-weighted face fluxes divided by cell volume gives the state derivative excluding sources and pressure work. Adding those terms gives the complete RK2 starting slope $k_1$. The final integral uses the average of the starting and endpoint sources and pressure work.

 Then iterate over current-level cells, save their local starting states from state_ into base_, save sources from source_ into source_start_, and compute pressure_work_start_ and slope_ using the available flux_, face_velocity_, and current-level starting pressure. These operations also run during reuse. The loop visits only cells_by_level_[depth], so the states and pressures read belong to current-level cells that have reached fraction; earlier coarser states are not used as this level's starting state. Phase I only saves starting states and slopes; current-level state_ is updated in its own Phase II.
4. Recursively advance the next higher level twice. First execute advance_level(depth+1,fraction,true), which advances the next higher level from fraction to fraction+span/2. Reuse is true because the parent has already prepared face quantities and sources for the shared start, with coarser predictions included, and its needed face set covers all faces required by the child. Note that this is recursive: the next level's advance_level(depth+1) returns to step 2 and continues this process until the highest level finishes, before returning level by level to the current level; even an empty level must complete finer-level recursion before skipping its own finalization. Thus, after this advance_level call returns, all levels higher than depth have advanced to fraction + span/2.

Next, the current level executes the second call, advance_level(depth+1,fraction + span/2,false), which advances the next higher level from fraction + span/2 to end_fraction. This call begins half of the current level's step after the current level's starting point, so it cannot share that starting point, and reuse=false is set. Even if the previous finer step evaluated endpoint fluxes at this time, those used predicted endpoint states. The finer cells have now completed their RK2 updates, so a new starting view_ must combine those updated states with new predictions for slower cells. This is also recursive: the higher level returns to step 2 and continues the process. Therefore, the remaining logic at the current level is executed only after all higher levels have advanced to end_fraction.

The discussion above does not cover invalid cases. If a validity check fails at any level, return false propagates level by level back to the attempt loop in run, which halves the maximum step size macro_step and retries.

5. Next comes Phase II of RK2. To calculate the time derivative of the cell state at the endpoint (the endpoint slope), we first use build_stage_view to predict the cell states view_ at end_fraction. Copy state_ into view_ again. Because of the recursion described above, cells at higher levels have already reached the endpoint and do not need modification. The current level has not yet reached the endpoint, so cells at the current and lower levels all need to be extrapolated to the endpoint using slope. This updates only `view_`; current-level cells in `state_` retain their local-step starting states until finalization in step 8.

6. Use view_ to call solve_faces, evaluating endpoint flux_ only for faces owned by the current level. Faces with lower owners cannot touch current-level cells. Faces with higher owners have already been integrated during finer substeps, and their contributions remain in the shared pending_flux_integral_ until the corresponding cells are finalized. Applying those accumulated finer-step integrals is a separate operation from reusing starting face quantities and sources through reuse_start_stage in Phase I.

7. For faces owned by this level, average starting first_flux_ and endpoint flux_, multiply by face area and the local step, and add the negative integral to pending_flux_integral_[left] and the positive integral to pending_flux_integral_[right]. These integrals remain in the shared array until the corresponding cell's level finalizes it, dividing by cell volume to obtain a state increment. They are not copied down the recursion stack on return.

8. Iterate over `cells_by_level_[depth]` and reference each current-level cell with `auto& cell = state_.cells[flat]`. Assign `cell.U = ...` starting from `base_[cell_index]`, adding `pending_flux_integral_[cell_index]` divided by volume and the two-stage source integral, then add the two-stage pressure-work integral to the separate internal energy. These assignments directly advance current-level cells in `state_` to `end_fraction`. Higher levels retain the states completed during recursion; lower levels wait for their own Phase II. Finally, clear the entries in pending_flux_integral_ belonging to cells at this level. They will not be needed in the subsequent integration of lower levels (a cell cannot belong to two levels), so the entire array is cleared once all levels finish. During intermediate steps, if an advance_level call at depth has completed, the pending_flux_integral_ entries for cells at depth and all higher levels will have been cleared, without affecting the entries of lower-level cells when those levels (such as depth-1) are processed.
For example:   f1  C1  f2 C2  f3 C3 f4
C1 belongs to level 1, C2 to level 0, and C3 to level 2. The owner_level of f1 is 1, that of f2 is also 1, and those of f3 and f4 are 2. At level 2, the fluxes of faces f3 and f4, whose owner level is 2, are calculated, and their flux integrals are added to the accumulated integrals of the adjacent cells C2 (f3) and C3 (f3 and f4). The level-2 advancement then updates C3 using the accumulated integrals received from f3 and f4 and clears C3's accumulated integral. C2 belongs to level 0, so the contribution from f3 to C2 is not cleared. Execution then returns to level 1, where the accumulated integrals of f1 and f2 are added to the accumulated integrals of C1 (f1) and C2 (f1,f2). Level 1 then advances C1 using the accumulated integral received from f1 and clears C1's accumulated integral. The contributions received by C2 are still not cleared; at this point, C2's accumulated integral is the sum of the integral from f2 and the negative integral from f3. Finally, execution returns to level 0, which owns no faces, and finalizes the level-0 cell C2. The pending_flux_integral_[C2] accumulated and passed down from higher levels can now be applied, adding the sum of the contributions from the adjacent faces f2 and f3 to C2's state. Thus, even though C2's adjacent faces belong to higher levels, all of their contributions are correctly accounted for, and each face contribution needs to be calculated only once, avoiding redundant computation. Finally, everything is cleared.

9. Return true to complete this advancement. As recursion returns to the preceding level, finalization proceeds level by level. If a validity check fails, however, all levels quickly return false, returning to the attempt loop in run to try a macro step that is 0.5 times as large. Once the outermost `advance_level(0, 0, true)` succeeds and all levels have completed the macro step, `run()` fills the endpoint boundaries of `state_` and commits to the simulation grid with `mesh_.cells.swap(state_.cells)`. Otherwise, the next `reset_attempt()` restarts from `start_`.


## 7. Complete example: four cells and three levels

Consider four cells with nonperiodic boundaries and prescribed levels `[0, 1, 2, 1]`:

```text
Boundary | f0 | C0 | f1 | C1 | f2 | C2 | f3 | C3 | f4 | Boundary
Level:          0         1         2         1
Step size:      H        H/2       H/4       H/2
```

This layout illustrates the recursion order. In an actual run, `assign_time_levels()` may promote additional neighbors to satisfy buffer requirements.

<a id="face-ownership"></a>

### 7.1 Face ownership

| Face | Adjacent active cells | Owning level |
| --- | --- | ---: |
| `f0` | `C0` | 0 |
| `f1` | `C0`, `C1` | 1 |
| `f2` | `C1`, `C2` | 2 |
| `f3` | `C2`, `C3` | 2 |
| `f4` | `C3` | 1 |

The face sets are therefore

| Level | `faces_owned_by_level_` | `faces_needed_by_level_` |
| --- | --- | --- |
| 0 | `{f0}` | `{f0, f1, f2, f3, f4}` |
| 1 | `{f1, f4}` | `{f1, f2, f3, f4}` |
| 2 | `{f2, f3}` | `{f2, f3}` |

<a id="initial-stage"></a>
<a id="level-0-start"></a>
<a id="level-1-first-step"></a>

### 7.2 Preparing levels at the common start

For brevity, label the macro-step start as 0. The initial scan evaluates all face fluxes, wave speeds, face velocities, and cell sources. Execution then enters

```text
advance_level(0, 0.0, true)
└─ advance_level(1, 0.0, true)
   └─ advance_level(2, 0.0, true)
```

All three levels reuse the same start-stage results, while each saves its own face fluxes and cell prediction data:

| Level | Start fluxes saved | Cells whose start states and slopes are saved |
| --- | --- | --- |
| 0 | `f0` | `C0` |
| 1 | `f1`, `f4` | `C1`, `C3` |
| 2 | `f2`, `f3` | `C2` |

For example, `C1`'s start slope is the flux balance between its left face `f1` and right face `f2`, divided by its volume, plus sources and internal-energy pressure work. Although level 2 owns the integration of `f2`, level 1 still needs its start flux to construct `C1`'s prediction slope.

<a id="level-2-first-step"></a>
<a id="level-2-second-step"></a>

### 7.3 `C2` completes two fine steps: $0\rightarrow H/4\rightarrow H/2$

The first fine step has no finer level to call, so it constructs the endpoint view at $H/4$ directly. All four cells are predicted to $H/4$ using their start slopes. The code solves `f2` and `f3` and computes

$$
I_f^{[0,H/4]}=A_f\frac{H}{8}\left(F_f^{(0)}+F_f^{(1)}\right),
\qquad f\in\{f2,f3\}.
$$

The stage superscripts denote the start and predicted endpoint of this fine step. Contributions are accumulated with the following signs:

| Face | Left-cell contribution | Right-cell contribution |
| --- | --- | --- |
| `f2` | `C1`: $-I_{f2}$ | `C2`: $+I_{f2}$ |
| `f3` | `C2`: $-I_{f3}$ | `C3`: $+I_{f3}$ |

`C2` applies its net integral, two-stage sources, and pressure work, then clears its own accumulation entry. `C1` and `C3` retain their contributions because their local steps are not yet complete.

The second fine step calls

```cpp
advance_level(2, 0.25, false);
```

In the start view at $H/4$, `C2` uses its newly completed RK2 state, while the other cells still use slow-cell predictions. The start fluxes of `f2` and `f3` must therefore be recomputed; the preceding step's predicted-endpoint fluxes cannot be reused.

Prediction, solution, integration, and updating then repeat to advance `C2` to $H/2$. At this point, `C1` and `C3` have received both $H/4$ integrals from the finer faces but have not yet been updated to $H/2$.

<a id="level-1-first-commit"></a>

### 7.4 Finalizing level 1: `C1` and `C3` reach $H/2$

The level-1 endpoint view predicts `C0`, `C1`, and `C3` to $H/2$ and uses `C2`'s updated state at that time. It then solves the owned faces `f1` and `f4` and computes their trapezoidal integrals over $[0,H/2]$.

Let $I_f^{[a,b]}$ denote the completed face integral over that interval. The pending net integrals are

$$
J_{C1}=I_{f1}^{[0,H/2]}
-I_{f2}^{[0,H/4]}-I_{f2}^{[H/4,H/2]},
$$

$$
J_{C3}=I_{f3}^{[0,H/4]}+I_{f3}^{[H/4,H/2]}
-I_{f4}^{[0,H/2]}.
$$

After adding their sources and pressure work, `C1` and `C3` are updated to $H/2$ and clear their own accumulation entries. The contribution $-I_{f1}^{[0,H/2]}$ to `C0` remains pending.

<a id="level-1-second-step"></a>

### 7.5 The second half of the macro step: $H/2\rightarrow H$

On entry to `advance_level(1, 0.5, false)`, `C1`, `C2`, and `C3` are synchronized at $H/2$; only `C0` needs prediction from 0 to $H/2$. Level 1 resolves the needed faces `f1`–`f4` and saves new start slopes for `C1` and `C3`.

It then executes

```cpp
advance_level(2, 0.5, true);   // H/2 -> 3H/4
advance_level(2, 0.75, false); // 3H/4 -> H
```

The first fine step reuses `flux_`, `speed_`, `face_velocity_`, and `source_` just prepared by level 1. Although `C0` in `state_` is still at 0, these caches were evaluated from a `view_` containing its predicted state at $H/2$. `C2` saves `base_` from its own updated state at $H/2$ and computes its starting slope, without predicting `C0` or evaluating the starting face quantities and sources again.

The second fine step begins at $3H/4$, uses the updated `C2`, and predicts `C0`, `C1`, and `C3` in a new `view_` to reevaluate the start stage. Once `C2` reaches $H$, level 1 constructs its endpoint view, completes the second-half integrals for `f1` and `f4`, and updates `C1` and `C3` to $H$.

<a id="level-0-commit"></a>

### 7.6 Finalizing `C0` last

`C1`, `C2`, and `C3` have now reached $H$. The working state of `C0` is still at 0, but it has received the `f1` contributions from both half-steps.

Level 0 predicts `C0` to $H$ and uses the other cells' updated states directly. It solves the endpoint flux at `f0` and completes that face's macro-step integral. The net integral for `C0` is

$$
J_{C0}=I_{f0}^{[0,H]}
-I_{f1}^{[0,H/2]}-I_{f1}^{[H/2,H]}.
$$

After adding `C0`'s two-stage sources and pressure work, all four cells are synchronized at $H$. The final boundary fill and macro-step commit follow.

<a id="execution-timeline"></a>

### 7.7 Recursive execution order

```text
Prepare C0's slope at 0
│
├─ Prepare C1 and C3 slopes at 0
│  ├─ C2: 0 -> H/4
│  └─ C2: H/4 -> H/2
├─ Finalize C1 and C3 at H/2
│
├─ Prepare new C1 and C3 slopes at H/2
│  ├─ C2: H/2 -> 3H/4
│  └─ C2: 3H/4 -> H
├─ Finalize C1 and C3 at H
│
└─ Finalize C0 at H
```

## 8. Source evaluation and step-size checks

### 8.1 `evaluate_sources()`: check individual terms and their sum

This function reads the stage grid `stage` and accumulates source contributions into `source_`. It does not directly add them to the state in `stage`.

Each call clears the source workspace and invalidates the scale cache. For every source term, it evaluates the contribution, completes the internal-energy source as needed, checks the callback's step-size limit, checks the source-induced state-change limits, and adds the contribution to the total. It checks the summed source once more after all terms have been evaluated.

Individual checks prevent opposing sources from hiding large component contributions. The sum check detects excessive combined changes from sources acting in the same direction.

### 8.2 Allowed changes and caching within a stage

`scan_source_changes()` uses a hard-coded 10% coefficient to compute, for each cell,

$$
\mathcal A_\rho=0.1\rho,\qquad
\mathcal A_{m_d}=0.1\rho(|v_d|+c),\qquad
\mathcal A_e=0.1\frac{p}{\gamma-1},
\qquad c=\sqrt{\frac{\gamma p}{\rho}}.
$$

These allowances constrain source-induced changes in density, momentum density in each direction, and the separately evolved internal-energy density. The momentum scale includes the sound speed, so its allowance is not simply 10% of the current momentum.

For a source-component rate $R$, the condition is

$$
H|R|\leq\mathcal A.
$$

This constrains the entire macro step $H$, even when a fine substep is being evaluated. It is a step-size restriction based on stage rates, not a strict bound on the change over the full nonlinear evolution.

Within one `evaluate_sources()` call, the first source fills the density and allowance caches. Subsequent sources and the sum check reuse those scales. The next stage evaluation refills them; reusing the same `view_` object does not make the previous scales valid.

All source components are checked for finite values. Zero sources are permitted in vacuum cells, but the current implementation does not support creating gas from vacuum through a nonzero source.

### 8.3 Initial step selection and later trial steps return different limits

| Call context | Behavior of `source_change_limit()` |
| --- | --- |
| Initial scan, `trial_H = 0` | Return the actual step-size bound $\min(\mathcal A/\lvert R\rvert)$ over nonzero rates |
| Later stage, `trial_H > 0` | Test the current $H$: return $\infty$ if acceptable, or 0 if exceeded |

The initial bound participates in macro-step selection. Later stages already know the trial $H$, so they only need an acceptance test rather than an exact minimum bound. Source callbacks' own step-size limits still independently constrain the final result of `evaluate_sources()`.

Later stages use

```cpp
scan_source_changes(stage, rates, [&](double allowance, double rate) {
    exceeded |= source_change_exceeds(trial_H, twice_H, allowance, rate);
});
```

Once any component exceeds its allowance, `exceeded` remains true. The scan nevertheless continues to detect invalid values in other cells or sources. A returned limit of 0 causes the caller to reject the macro step and retry with half the value of $H$.

### 8.4 A fast check that avoids unnecessary division

```cpp
bool source_change_exceeds(double H, double twice_H, double allowance, double rate)
{
    const double magnitude = std::abs(rate);
    if (twice_H * magnitude < allowance) return false;
    return rate != 0.0 && H > allowance / magnitude;
}
```

If $2H|R|<\mathcal A$, the current step is safely within the allowance and can be accepted immediately. Otherwise, the function falls back to the original division-based comparison, $H>\mathcal A/|R|$. A zero rate imposes no restriction. The factor of 2 provides margin for the fast acceptance test; it does not change the actual 10% limit.

## 9. Reading the multidimensional implementation from the 1D reference

The original 1D access

```cpp
auto& cell = view_.cells[ng_ + cell_index];
```

corresponds to the following in the current multidimensional code:

```cpp
const auto flat = storage_index(cell_index);
auto& cell = view_.cells[flat];
```

where

```cpp
std::size_t LTSRK2Step::storage_index(int cell_id) const
{
    return mesh_.active_to_storage[cell_id];
}
```

`cell_index` is an active-cell identifier, while `flat` indexes full storage including ghost regions. `base_`, `slope_`, and the flux-accumulation array use active-cell identifiers. `MeshBlock::cells`, `source_`, and cell volumes use storage indices. A single fixed ghost offset cannot replace this mapping in multiple dimensions.

Other principal correspondences are:

| 1D reference | Current multidimensional implementation |
| --- | --- |
| Faces `i` and `i+1` bound cell `i` | `mesh.cell_faces[i][d]` gives the lower and upper faces in direction `d` |
| One pair of face-flux contributions | Sum of area-weighted flux differences over all directions |
| Pressure work from one pair of face velocities | Sum of area-weighted face-velocity differences over all directions |
| CFL rate from left/right wave speeds and cell width | Sum of directional contributions; dynamic checks use face-speed maxima recorded during the current attempt |
| Buffers extended through 1D neighbors | Neighbor traversal along each coordinate direction through face topology |

Wave-speed checks use each adjacent cell's full local step size. Checking a coarse cell using only the face's finer step would be insufficient. Failed dynamic CFL, source-limit, or physical-validity checks follow the same macro-step rollback procedure.

<a id="core-principles"></a>

Three distinct operations connect the full algorithm: **prediction** supplies states at a common time for stage evaluation; **face integration** supplies shared transport contributions to adjacent cells; and **local finalization** updates a level once all its contributions are available. Recursion orders these operations so that every level ultimately synchronizes at the macro-step endpoint.
