# SNEC profile013 thermal-bomb comparison

This case imports SNEC's actual initial spherical shells into the current FLARE
framework and compares early shock propagation through 100 physical seconds.
All adapters and analysis live in `unit_test/`; neither solver's source is edited.
See [VALIDATION.md](VALIDATION.md) for measured results and limits.

## Reproduce

From the FLARE repository root:

```bash
# Rebuild both current sources and run SNEC ideal, SNEC Paczynski+radiation,
# FLARE exact, and FLARE HLL on a common radially resolved initial grid.
python3 unit_test/compare_snec_thermal_bomb.py run \
  --output tests/SNEC_Thermal_Bomb/output_new --radial-cap 0.02

# Resolution check: refine both codes' initial grid, retaining the original
# profile interpolation procedure and the same prescribed explosion.
python3 unit_test/compare_snec_thermal_bomb.py run \
  --output tests/SNEC_Thermal_Bomb/output_new_fine --radial-cap 0.01 --hydro-only

# Regenerate figures and numerical summaries from completed runs.
python3 unit_test/compare_snec_thermal_bomb.py analyze \
  --output tests/SNEC_Thermal_Bomb/output_new

# Run the bundled, already-converted FLARE initial state by itself.
make -C tests/SNEC_Thermal_Bomb -j4 run
```

Use a new output directory for each comparison. The runner refuses to mix runs.
The default SNEC source is the sibling `Research/SNEC/SNEC_NEW`; override with
`--snec-root /path/to/SNEC`. A private build copy is made under `.build/snec`.
The source directory and its historical outputs are only read. Original SNEC
opacity tables are referenced by symlink. Requirements: C++17 compiler,
gfortran, LAPACK, Python, NumPy and Matplotlib.

`make compare` runs the full comparison; `make analyze` regenerates its plots.
For the standalone FLARE run, change `output_directory` in `parameters.txt`
before rerunning if existing snapshots must be preserved.

## Matched setup

| Quantity | Setting |
| --- | --- |
| Progenitor | Local `snec_input13.short`, 4147 rows; matching isotope file |
| Mass cut | 1.502 SNEC solar masses; `M_sun = 1.98892e33 g` |
| Geometry / gravity | Spherical 1D, central mass plus envelope self gravity |
| Prescribed injected energy | `1.6e51 erg`, SNEC `bomb_mode=2` |
| Heating duration | 0–0.1 s |
| Heating mass parameter | 0.1 solar masses, SNEC's discrete shell selection |
| Temporal / mass endpoint ratio | 100 / 100 |
| Matched EoS | Ideal gas, gamma = 1.4, SNEC's existing ideal-EoS setting |
| Matched radiation / nickel | Both off; composition smoothing off |
| FLARE scheme | PLM, minmod, SSP RK2, CFL 0.4; exact and HLL |
| Inner boundary | SNEC fixed wall / FLARE reflective at the same radius |
| Outer boundary | SNEC moving free surface / FLARE fixed-domain outflow |

SNEC's usual `bomb_mode=1` means **target final energy**, with initial binding
energy subtracted to determine the bomb. Mode 2 is used here so the intended
injection is identical even when changing EoS. The extra Paczynski+radiation run
uses the same profile, mass grid and prescribed bomb, but is a physics comparison;
its initial pressure/internal energy and subsequent EoS differ from the ideal runs.
It does not establish agreement between FLARE's and SNEC's radiation physics.

The temporal luminosity is
`L(t) = E_bomb ln(100) exp(-ln(100) t / 0.1) / (0.1 * 0.99)`.
The spatial heating is piecewise constant in each original SNEC mass shell,
proportional to `exp(-A * (m_left - m_cut))`, and normalized by shell mass.
These conventions follow the local `bomb_profile.F90` implementation and the
[SNEC manual](https://stellarcollapse.org/codes/snec_notes-1.00.pdf).

The FLARE hook integrates these weights over each cell's **current enclosed-mass
interval** at every RK stage. Thus heating follows the same mass coordinates
while the Eulerian gas moves. It adds to both total and auxiliary internal energy,
limits heating steps to `1e-4 s` and at most 5% instantaneous internal-energy
change, and is stateless under rejected steps. Runtime source checks verify total
power before/after deliberately redistributing shell masses, and zero power after
the bomb ends.

## Initial-state mapping and grid control

SNEC interpolates rho and temperature at mass centres and integrates radius from
shell mass and density. Simply reading the progenitor's radii into FLARE would
give a different initial model. The runner therefore first executes a tiny SNEC
initialization and exports its full-precision `radius/rho/press/vel.xg` initial
frame. The terminal passive SNEC point is excluded from active cells. Initial
FLARE rho, pressure, shell radii and enclosed mass are checked against SNEC.

The native 1000-point mass grid contains extremely wide radial cells outside
the helium core. `--radial-cap 0.02` adds mass-grid faces until the initial
`delta_r/r_left <= 0.02`. Each refinement re-interpolates the **original** profile
and reintegrates radius, and both codes use this resulting mass grid and state.
Use `--radial-cap 0` to reproduce the native-grid diagnostic. `--refine` additionally
splits only FLARE shells into two equal-volume pieces; this is a separate fixed
initial-state sensitivity test, not a refinement of both codes.

FLARE uses its existing gravity implementation and constants. SNEC uses
`G=6.6742e-8`; FLARE uses `6.67430e-8` (relative difference about `1.5e-5`).
Both have stationary inner boundaries. Outer boundary formulations differ; these
100 s runs measure an early interior shock, not shock breakout or a light curve.

## Outputs and interpretation

- `manifest.json`: FLARE commit, hashes of every solver source and input,
  physical settings, run completion and wall times, source-unchanged check.
- Each run's `parameters*`, `run.log`, initial state, snapshots and performance data.
- `comparison/comparison.png`: shock radius/mass, energy accounting, final velocity.
- `comparison/profiles.png`: final density, velocity and pressure profiles.
- `comparison/summary.json`, `*_history.csv`, `common_times.csv`: numerical results.

The shock estimator uses the largest outward compressive cell-centred velocity
drop with a pressure decrease, identically for both codes. Its adjacent cell
centres give a **grid bracket**, not a statistical confidence interval. Histories
can jump across grid faces. Time interpolation in `common_times.csv` does not
remove this spatial uncertainty. At very early times before a resolved shock is
identified, the history intentionally contains NaN shock positions.

Energy diagnostics integrate internal plus kinetic energy and the exact
gravitational binding of piecewise uniform spherical shells. SNEC's passive
outer point is excluded. The script reconstructs actual injected energy from
accepted-step logs: SNEC left-endpoint heating versus FLARE RK2 trapezoids.
This reconstruction inherits log precision. The reported residual subtracts that
actual heating; it does not add boundary energy fluxes or escaped photons, so
the radiation run's residual is not solely a numerical conservation error.

Prescribed identical heating does not imply bitwise identical numerical deposited
energy. Mesh dependence, discrete heating quadrature, reconstruction/artificial
viscosity, geometry, gravity work and outer boundaries are separate differences.
Use the resolution results before assigning a speed difference to a Riemann solver.
