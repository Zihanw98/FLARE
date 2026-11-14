# Multilevel explicit LTS Euler

## 实现 pipeline

`advance_lts_euler` 负责校验输入并启动一次内部 `LTSEulerStep`。
网格副本、面量、源项、层级表和递归推进逻辑都由这个类持有，与 `LTSRK2Step` 的组织方式一致。
一次同步大步按以下顺序完成：

1. 复制同步起点状态、填充边界，计算全部独立面的通量和波速。
2. 计算各单元 CFL 局部步，并在同步起点计算一次需要在整个 `H` 内冻结的源项。
3. 由剩余时间、最细层 CFL 容量和源项上限共同选择同步大步 `H`。
4. 给单元分配 `H/2^level`，给面分配负责层；细层重构前建立同一时刻的阶段视图。
5. `advance_level` 递归执行两个细半步；每个物理面只由最细相邻层积分一次。
6. 面积分以相反符号写入左右单元寄存器，每层结束时结算本层单元并同步双能量。
7. 非物理解或动态 CFL 超限会丢弃完整试算、将 `H` 减半并重新分层。
8. 全部层成功后才提交正式网格和星区压力历史。

Select the separate integrator in `src/time/lts_euler.cpp` with:

```ini
time_integrator = lts_euler
reconstructor = PLM
limiter = MC
nxg = 2
cfl = 0.2
lts_max_levels = 8
```

`lts_max_levels` counts all time levels, including the coarsest. It accepts 1–16
and defaults to **2**, preserving the previous two-level configuration. One level
reduces to global Euler. Eight levels permit a maximum coarse/fine step ratio of
128; sixteen permit 32768. These ratios are ceilings, not guaranteed speedups:
source limits, the requested end time and retries can shorten the macro step.

This experimental implementation supports **1D/2D/3D Cartesian Godunov** and
**1D spherical Godunov** on fixed spatial meshes. Select `PCM`, `PLM`, or
`PLM_corrected`; PCM needs at least one ghost cell on every axis, and both PLM
variants need at least two. Periodic boundaries must be paired on each axis.
Existing global Euler and RK2 algorithms are unchanged. **Time accuracy remains
first order**: PLM supplies linear spatial reconstruction, not RK2 time stepping.

## One synchronized macro step

1. Fill start-time boundaries, solve all faces and retain fluxes and wave speeds.
2. Evaluate each source on the synchronized start state. The macro step `H` is
   bounded by `2^(lts_max_levels-1)` times the minimum hydro CFL step, source callback limits, a
   10-percent source-change estimate, and the caller's remaining-time limit.
3. Assign the smallest level `l` for which `H/2^l` fits each cell's CFL bound.
   Add one neighbour along each axis and side for PCM, or two for PLM. Grade
   surrounding levels through the same axial-neighbour traversal. The initial
   promotion follows coordinate axes, not a full diagonal box. Build cell and
   face lists once per trial, not on each finest substep.
4. At each level, predict its cells' Euler endpoints, recursively advance the
   finer levels twice, then commit that level's cells. PCM holds coarser states
   at their **local** step start. PLM interpolates between each coarser cell's
   local start and provisional endpoint to the actual fine-stage time. Fill
   all directional ghosts in the stage view and reset the reconstruction cache
   before each fresh solve. Accepted local states remain in `state_`.
   The first child step reuses its parent's start-time fluxes; the second solves
   only faces touching that level or finer cells.
5. Integrate each physical face only at the finer of its adjacent levels. Add
   the same `area * flux * dt` with opposite signs to the two cells' pending
   flux registers. Each cell consumes and clears its register when its own step
   ends. Its provisional predictor is replaced by this accumulated update.
   This provides flux synchronization without double-counting a coarse flux.
   The periodic seam is one physical face. Validate states and fill final ghosts
   before committing the synchronized macro step.
6. If a candidate is nonphysical or measured waves exceed a cell's own full local
   step bound (including the coarser side), discard the complete trial and halve
   `H`. `integrator_retry_max` bounds retries. External callback/solver errors
   propagate without committing mesh state.

## Source approximation and limits

**All source rates**, including geometric pressure, gravity acceleration and
gravity work, are frozen at the macro-start state for `H`. This is stronger than
freezing only the gravitational field. Sources are not evaluated on asynchronous
states and self-gravity is recomputed at the next synchronization. Each source
is limited independently before summation, using density, acoustic-scaled
momentum and the instantaneous internal-energy rate; the sum is checked too.
These estimates control first-order source changes, not arbitrary time variation
inside a step. Custom source callbacks must report any additional hard timestep
limit, particularly for prescribed time-dependent forcing. Gas creation from
exact vacuum by a source is unsupported.

The source-change checks cache `NDIM + 2` thermodynamic scales per active cell:
density, one acoustic momentum scale for each direction, and internal energy.
The dual-energy version checks the explicit `S[INTERNAL_ENERGY]` rate. The first
independent source fills these scales; later sources and the sum reuse them.
`LTSEulerStep` owns `source_scales_` directly. Every new macro step recomputes
the scales, while retries reuse the same frozen start-state rates and limits.
With no sources, the source evaluation returns an infinite limit immediately.
Finite-rate/vacuum checks and the 10-percent threshold apply to every source
and to their sum.

This optimization also exists in LTS RK2, where the cache must instead refresh
at every source stage evaluation. Global Euler and RK2 share `evaluate_RHS`,
which uses the source callbacks' reported limits and has no generic 10-percent
source-change scan. Adding such a scan there would change timestep semantics,
so those algorithms are not modified by this cache optimization. Optimizations
to shared source or reconstruction routines can benefit all four integrators.

`source_cache_tests` covers consecutive cooling steps with reused pressure hints,
independent limits for cancelling sources, later-source NaNs, changed mesh/state
after an exception, and zero/nonzero sources in vacuum.

Mass and flux contributions cancel across interior faces; boundary fluxes and
sources still change domain totals. Refluxing does not by itself conserve
self-gravitating gas-plus-potential energy, or radial momentum in spherical
geometry. Assess source errors by reducing the macro step/CFL and comparing
against globally stepped Euler with the same spatial discretization.

## Scope and validation

This is a hierarchy of **time steps on a fixed mesh**, not spatial mesh refinement.
The finest cells still take their explicit CFL-limited steps. Source rates remain
frozen over `H`, so the existing source-change and callback limits remain essential;
more time levels do not remove source accuracy requirements. A buffer and
measured local CFL checks are safeguards, not a general proof of nonlinear LTS
stability; strong-shock accuracy needs case-specific refinement studies.

The existing performance profiler records Riemann solves. Each face is solved
at the frequency of its finer adjacent time level; the synchronized initial solve
is reused by the first substeps. Source evaluation remains at macro-step frequency.
Output and the simulation step counter refer to synchronized macro steps.

## Class-owned state and multidimensional indexing

`LTSEulerStep` owns all temporary storage for one synchronized macro step:

- `start_`: the boundary-filled synchronized start, retained for rollback.
- `state_`: each cell's accepted local state, which may be at a different time.
- `view_`: a copy of `state_` with coarser cells predicted to the requested stage
  time, followed by boundary filling on every axis. Reconstruction reads this view.
- `step_start_` and `candidate_`: each active cell's local start and provisional
  Euler endpoint, used to interpolate PLM coarse states.
- Face results, source rates, flux registers, level lists and reconstruction cache.

There is no external `LTSEulerWorkspace`. The optional pointer after `cfl` in
`advance_lts_euler`, or after `dt_next` in `advance_one_step`, is a
`hydro::StarPressureCache*`, just as for LTS RK2. The main driver passes its shared
pressure cache. Omitting it is supported and starts each macro step with fresh
pressure guesses. Temporary arrays and mesh copies are allocated per call.

Cell arrays use active-flat indices and access full storage through
`mesh.active_to_storage`. Face results use global independent `face_id` values;
`mesh.face_info` supplies direction, storage coordinates and active neighbours.
`mesh.cell_faces[cell][d]` supplies the lower and upper faces. Periodic seam faces
have one shared flux register contribution. Pressure hints use the directional
active-face layout, with shape `nactive + e_d`; the duplicate upper periodic hint
is synchronized to the lower hint. Rebuild mesh topology after changing boundary
conditions, and clear pressure hints after remeshing.

For cell `i`, the initial multidimensional CFL rate is

\[
\lambda_i=\sum_{d=0}^{D-1}
\frac{\max(a_{i,d,-},a_{i,d,+})}{\Delta x_{i,d}},
\qquad h_i=\frac{\mathrm{CFL}}{\lambda_i}.
\]

Dynamic checks use the largest speed seen on each face during the current
attempt, sum over directions again, and check the cell's full local step against
`cfl_max`. A retry restores all face-speed records to their initial values.
Flux divergence and pressure work also sum the area-weighted contributions
from both faces on every axis before division by cell volume. See
[multidimensional CFL](multidimensional_cfl.md) for the shared CFL convention.

## Pressure initial guesses

For `riemann_solver=exact`, Euler, RK2 and both LTS methods store star pressures
in the shared `hydro::StarPressureCache` type (`src/hydro/star_pressure_cache.h`),
with one entry per directional face. The main driver
enables reuse automatically. Euler uses the previous step's pressure; RK2 uses
the previous accepted step in stage 0 and stage 0's pressure in stage 1. A fine
LTS substep uses the most recent pressure solved at that face within the macro step. The
pressure equation, bracket construction, periodic bisection safeguard and
relative bracket convergence tolerance are unchanged. Nonpositive, nonfinite
or out-of-bracket guesses fall back to the original acoustic/midpoint guess.
Vacuum clears the entry. The bisection reference solver ignores the input guess.

The per-face cache is named `star_pressure_cache` throughout the time integrators
and spatial operators. `initial_pressure` and `trial_pressure` are its stage-local
copies; `pressure_hint` is a single face's value passed to the Riemann solver.
The cache provides `ensure_size(direction, face_count)` to preserve matching
storage or initialize a changed face count, `reset(direction, face_count)` to
invalidate a direction even at the same size, and `clear()` to discard all
directions. It owns no geometry or scheduling logic. Each integrator retains its
own stage/substep updates, periodic-face synchronization, retries and commits.
Initial and trial caches keep history transactional: rejected steps
and exceptions before state acceptance cannot commit partial history. Every RK2 retry restores the
stage-0 pressure guesses before evaluating stage 1. LTS periodic boundary faces
share the same physical solve.
Direct callers of all four explicit integrators can pass a persistent `hydro::StarPressureCache` through
the optional `star_pressure_cache` argument; omitting it retains the original
cold-start behavior. Both LTS methods keep local trial copies. Call `clear()` on a caller-owned
cache for a new simulation or remeshing. Face-count changes also reset the
affected direction's entries.
The solver owns no cross-call state. Existing `flux()` callers keep the original
initial guess; `flux_with_pressure_hint()` is an optional interface, with a
fallback to `flux()` for custom solvers. Different initial guesses can change
roundoff-level results, so regression checks use physical error scales rather
than requiring byte-identical output.

Run focused checks from the repository root:

```sh
make -C unit_test -j4 TESTS="lts_euler optimization rk2" DIMENSIONS="1 2 3"
```

Tests cover all-slow/global Euler and all-fast/two-Euler equivalence; nonuniform
periodic conservation including thin end cells; reduced face solves; source
limits and synchronized evaluation; spherical pressure balance; wave-triggered
retries and rollback; unsupported input rejection; and a shock crossing the
fast/slow region compared with globally stepped Euler, for all three
reconstructors. A separate PLM stencil test compares every second-substep
interface state against a globally advanced Euler half-step with heating. It
also checks that stale slow-neighbour states would give a different result.
Additional multilevel tests cover 16:1 macro steps, periodic thin cells at either
end or in the interior, conservation and face-solve savings against global Euler,
one-level and all-fine Euler equivalence, fine-stage PLM prediction, source limits,
shock crossing, rollback and pressure-history reuse. Repeated-step regression tests
also change the number of levels between calls. The timing numbers below describe
the original **default two-level** runs, not multilevel performance.

Repeated-step regression tests cover repeated steps, mesh sizes and ghost counts,
geometry rebuilds, failed trials followed by reuse, and a user boundary that
reads a remote slow cell at the midpoint. Pressure-history tests check invalid
and unrelated old guesses, vacuum, exceptions, periodic faces, mesh changes,
rollback, and agreement with the independent bisection solver. Euler/RK2 tests
also check multidimensional face indexing, reuse across accepted steps, source
exceptions, invalid-state failures, and discarding rejected RK2 stage pressures.

The multidimensional tests use unequal axis sizes and ghost counts. They check
single-step agreement with global Euler for all three reconstructors, mixed
periodic/reflecting boundaries, per-direction pressure hints and rollback,
multilevel conservation with heating and a thin cell on each axis or periodic
seam, comparison with global Euler, and directional source limits. A transient
wave-speed test requires rejection based on the **sum** of directional CFL
rates, then verifies that the next attempt clears rejected wave-speed history.

## Historical 1D measurements

The following numbers were measured with the earlier 1D implementation. They
are retained as historical comparisons, not timings of the current class-owned
multidimensional implementation.

The native CCSN MESA case was also checked through 0.1 physical seconds using
PCM and CFL=0.2 for both methods, with temporary output directories. Its existing
`check.py` passed for both runs: LTS took 9,834 macro steps and 26,640,297 exact
Riemann solves versus Euler's 19,667 steps and 53,100,900 solves. The largest
saved relative mass change in LTS was 1.55e-15. At the final time the maximum
pointwise differences relative to Euler were 5.49e-5 in density and 1.08e-3 in
pressure. This is a short-time regression comparison, not a convergence study
or validation through shock breakout.

With PLM/MC and CFL=0.2, the same 0.1-second MESA regression also passed for both
methods: LTS used 9,832 macro steps and 26,674,215 Riemann solves versus Euler's
19,664 steps and 53,092,800 solves. The largest saved relative mass change was
6.66e-16 for LTS. Final maximum pointwise differences relative to Euler were
1.08e-4 in density and 4.64e-4 in pressure. The ready-to-run configuration is
`tests/CCSN_MESA_day13/parameters_lts_plm.txt`.


After removing shared geometry/primitive caching, a three-pair serial timing
comparison isolated the LTS pressure-initial-guess change on the same 0.1-second
MESA PLM/MC case (CFL=0.2, summary profiling, three snapshots): median total time
was 6.5546 s without pressure reuse and 6.4004 s with reuse (2.35% reduction).
Median integration time was 6.1566 s versus 5.9962 s (2.61% reduction). Both used
9,832 macro steps and 26,674,215 Riemann solves; mean pressure iterations fell
from 1.5075 to 1.1854. Final maximum pointwise relative differences were
5.81e-11 in density and 4.24e-11 in pressure. The original MESA check passed.
These short-run timings do not establish performance at later explosion phases.

The same pressure initial guesses are also enabled for global Euler and RK2.
A serial cold/warm MESA regression pair per method (same PLM/MC, CFL=0.2 and
0.1-second endpoint) gave the following counts:

| Integrator | Accepted steps | Mean iterations, cold | Mean iterations, warm |
| --- | ---: | ---: | ---: |
| euler | 19,664 | 1.4981 | 1.1470 |
| rk2 | 19,664 | 1.4981 | 1.0879 |

Each method retained its original Riemann-solve count. Both MESA checks passed;
maximum pointwise density/pressure differences were below 1e-10. The single-pair
wall times were 12.1640/11.5988 s for Euler and 23.9078/23.2419 s for RK2; these
are regression-run observations, not repeated median speedup measurements.

## Multilevel MESA short-run check

The unchanged native MESA case was advanced for 0.1 physical seconds with
PLM/MC, exact Riemann solves, central gravity and self-gravity. Both the two-level
baseline and the eight-level runs passed the existing `check.py`, including grid
mapping, positive states, injected energy, endpoint and mass conservation.
All four runs used the same executable and wrote three snapshots. Maximum saved
relative mass change was at most 6.67e-16.

| Levels | CFL | Macro steps | Riemann solves | Wall time (s) | Max density difference | Max pressure difference |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 0.2 | 9,832 | 26,674,215 | 6.962 | reference | reference |
| 8 | 0.2 | 154 | 867,156 | 0.512 | 1.33% | 5.14% |
| 8 | 0.1 | 308 | 1,732,400 | 0.872 | 0.649% | 2.58% |
| 8 | 0.05 | 615 | 3,462,172 | 1.654 | 0.303% | 1.17% |

Differences are final pointwise relative differences on active cells against
the two-level CFL=0.2 run, not errors against an exact solution. Decreasing CFL
also shortens the macro step and the frozen-source interval; differences shrink
approximately linearly in this check. The wall times are single serial runs,
not repeated timing statistics. More levels improve scheduling efficiency but
do not increase time accuracy or establish validity through shock breakout.


## Dual energy

The existing `Conserved` now includes internal energy density. Face registers also
transport this component. Each local Euler step adds `-h*p*div(v)` using its start
pressure and the sampled Riemann face velocities, both in the coarse predictor
and the final update. Physical source callbacks remain frozen at macro-step start.
Accepted local states synchronize total/internal energies through
`Cell::synchronize_dual_energy`; rejected trials discard the corrections with the
trial mesh. See [the shared dual-energy rules and tests](lts_rk2.md#独立内能与双能量同步2026-11-10).
