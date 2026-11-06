# Spherical Sod in 1D

This example evolves a spherical discontinuity at `r = 0.5` on `0 <= r <= 1`.
The inner state is `(rho, u, p) = (1, 0, 1)` and the outer state is
`(0.125, 0, 0.1)`, with `gamma = 1.4`. The origin is reflective; the outer
ghosts copy the last active cell. The calculation ends before the wave reaches
the outer boundary.

From the repository root:

```sh
make PROBLEM_DIR=tests/Spherical_Sod_1D run
```

The supplied configuration uses 128 active cells, two ghost layers per side,
RK2, PLM with the MC limiter, CFL 0.2, and final time 0.15. Change `nxa` to
256 or 512 to compare resolutions.

This is a spherical problem: the planar Sod exact solution is not its exact
reference. A comparison against a finer spherical run checks resolution
dependence, rather than proving an order of convergence.

Validation on 2026-11-05 completed at all three resolutions, with positive
density and pressure in every saved state. Relative to the 512-cell solution
restricted by shell-volume averaging of conserved states, the volume-weighted
density L1 difference decreased from `1.7295e-3` at 128 cells to `6.8326e-4`
at 256 cells. Across saved states and resolutions, the largest relative mass
and total-energy changes were `2.063e-11` and `2.548e-11`.

These conservation values were computed independently from the CSV states
using `V = 4*pi/3 * (r_right^3 - r_left^3)`. The framework's current printed
`mass` and `energy` diagnostics still use Cartesian cell-width products and
must not be interpreted as spherical totals.

A separate uniform static check, including the origin, advanced 128 cells
to time 0.2 with maximum speed `1.439e-11`. It also verified that the registered
geometric source changes only active radial momentum entries and leaves
ghost RHS entries unchanged.

## Longer evolution and animation

The long configuration runs to `t = 0.6` on `0 <= r <= 2`, with 1024 cells.
This keeps `dr` equal to the 512-cell short run and lets the outgoing shock
remain inside the domain. On the original `r <= 1` domain the shock reaches
the outer boundary around `t = 0.316`; the later solution is affected by that
boundary and its mass/energy decrease as fluid leaves.

```sh
make PROBLEM_DIR=tests/Spherical_Sod_1D CONFIG=tests/Spherical_Sod_1D/parameters_long.txt run
python3 animate.py tests/Spherical_Sod_1D/output_long spherical_sod --plot spherical --stride 3
python3 animate.py tests/Spherical_Sod_1D/output_long spherical_sod --plot spherical --plotly --stride 3
```

`--plot spherical` interprets CSV `x` as radius and `vx` as radial velocity.
Each left panel displays a symmetry-plane slice of the radial scalar field;
each right panel displays its radial profile. The projection uses inferred
shell edges, masks any central cavity and the exterior, and excludes ghosts.
The default Cartesian plotting mode is unchanged. GIF/MP4/HTML use the same
radial mapping, and the default spherical output names include `_spherical`.
For older snapshots whose filenames still begin with `mesh_`, use
`mesh_spherical_sod` as the filename prefix instead.

Long-run checks completed at 512, 1024, and 2048 cells. An independent NumPy
HLL + minmod PLM + SSP-RK2 implementation at 4096 cells provided a numerical
cross-check. At 2048 cells, relative volume-weighted L1 differences from that
reference were about 0.036% in density, 0.081% in radial velocity, and 0.017%
in pressure. These are numerical comparisons, not errors against an analytic
exact solution. In the extended domain, the maximum saved-state relative
mass and total-energy changes were below `6.1e-10`.
