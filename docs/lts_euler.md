# Two-level explicit LTS Euler

Select the separate integrator in `src/time/lts_euler.cpp` with:

```ini
time_integrator = lts_euler
reconstructor = PLM
limiter = MC
nxg = 2
cfl = 0.2
```

This experimental implementation supports **1D Godunov only**, Cartesian or
spherical geometry, and fixed spatial meshes. Select `PCM`, `PLM`, or
`PLM_corrected`; both PLM variants need at least two ghost cells. Unknown
reconstruction names and multiple dimensions are rejected explicitly. Existing
Euler and RK2 algorithms are unchanged. **Time accuracy remains
first order**: PLM supplies linear spatial reconstruction, not RK2 time stepping.

## One synchronized macro step

1. Fill start-time boundaries, solve all faces and retain fluxes and wave speeds.
2. Evaluate each source on the synchronized start state. The macro step `H` is
   bounded by twice the minimum hydro CFL step, source callback limits, a
   10-percent source-change estimate, and the caller's remaining-time limit.
3. Mark cells requiring `H/2` as fast and add one neighbour on each side for PCM,
   or two for PLM's wider stencil. Other cells take `H`. Recompute this partition
   at each synchronized macro step.
4. Advance slow cells provisionally with the start-time fluxes. Fast cells take
   two Euler half steps. PCM holds neighbouring slow states at their macro-start
   values. For PLM, predict slow states to `t+H/2` by averaging their start and
   provisional endpoint conserved states, then fill the needed half-time ghosts
   and reset the reconstruction cache. All cells in each PLM stencil then
   represent the same stage time. Predicted states are never committed to slow
   cells. Reuse initial fluxes for the first half step; solve only faces touching
   fast cells on the second half step.
5. Accumulate `area * flux * dt` on fast faces. Replace each mixed interface's
   provisional slow-side integral with the fine integral. Treat a periodic seam
   as one physical face. Validate the corrected states and fill final ghosts
   before committing the synchronized result.
6. If a candidate is nonphysical or measured waves exceed a cell's own step
   bound (including `H` on the slow side), discard the complete trial and halve
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

Mass and flux contributions cancel across interior faces; boundary fluxes and
sources still change domain totals. Refluxing does not by itself conserve
self-gravitating gas-plus-potential energy, or radial momentum in spherical
geometry. Assess source errors by reducing the macro step/CFL and comparing
against globally stepped Euler with the same spatial discretization.

## Scope and validation

This is **two levels**, not an arbitrary hierarchy: the maximum macro/global
step ratio is two. It is a foundation for multilevel subcycling, not an
orders-of-magnitude solution to the native MESA mesh bottleneck. A buffer and
measured local CFL checks are safeguards, not a general proof of nonlinear LTS
stability; strong-shock accuracy needs case-specific refinement studies.

The existing performance profiler records Riemann solves. An accepted macro
step solves all unique faces once plus fast-adjacent faces once, rather than
solving the full domain twice. State/source scans remain at macro-step frequency.
Output and the simulation step counter refer to synchronized macro steps.

## Reused workspace and pressure initial guesses

The main driver keeps one `hd_time::LTSEulerWorkspace` for the simulation and
passes its address as the `lts_workspace` argument to `advance_one_step`. Direct callers
can do the same with `advance_lts_euler`. Omitting this optional argument remains
valid but creates temporary storage on each call and loses cross-step reuse.
The workspace retains mesh snapshots, source/flux arrays, fast-cell/face lists,
and primitive/slope buffers. Cell states are still copied for rollback; fixed
mesh coordinates, face areas and volumes are copied only when geometry changes.
Separate concurrently advanced simulations must use separate workspaces.

Only slow cells read by fine-face PLM stencils are predicted at the midpoint.
Periodic and reflecting ghost donors are included. If a needed user boundary
can read arbitrary active cells, the implementation conservatively predicts all
slow cells before calling it. Unused boundary sides need no midpoint fill.

Mesh geometry keeps its original coordinate, face-area and volume arrays. No
shared geometry owner or separate physical-state cache is used. Reconstruction
retains its own primitive/slope arrays and invalidates them at every stage;
sources use the original callbacks.

For `riemann_solver=exact`, Euler, RK2 and LTS store the last accepted star
pressure in ordinary `std::array<std::vector<double>, NDIM>` storage, one entry
per directional face. No additional history type is introduced. The main driver
enables reuse automatically. Euler uses the previous step's pressure; RK2 uses
the previous accepted step in stage 0 and stage 0's pressure in stage 1. A fine
LTS substep uses the pressure just solved at that face's macro-step start. The
pressure equation, bracket construction, periodic bisection safeguard and
relative bracket convergence tolerance are unchanged. Nonpositive, nonfinite
or out-of-bracket guesses fall back to the original acoustic/midpoint guess.
Vacuum clears the entry. The bisection reference solver ignores the input guess.

Initial and trial pressure arrays keep history transactional: rejected steps
and exceptions before state acceptance cannot commit partial history. Every RK2 retry restores the
stage-0 pressure guesses before evaluating stage 1. LTS periodic boundary faces
share the same physical solve, and LTS geometry changes clear its pressure array.
Euler/RK2 direct callers can pass a persistent pressure array as the final
`advance_euler`, `advance_rk2` or `advance_one_step` argument; omitting it retains
the original cold-start behavior. Clear this array for a new simulation or
remeshing. Face-count changes also reset the affected direction's entries.
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
Workspace regression tests cover repeated steps, mesh sizes and ghost counts,
geometry rebuilds, failed trials followed by reuse, and a user boundary that
reads a remote slow cell at the midpoint. Pressure-history tests check invalid
and unrelated old guesses, vacuum, exceptions, periodic faces, mesh changes,
rollback, and agreement with the independent bisection solver. Euler/RK2 tests
also check multidimensional face indexing, reuse across accepted steps, source
exceptions, invalid-state failures, and discarding rejected RK2 stage pressures.

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
