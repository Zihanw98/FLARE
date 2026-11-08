# Thermal bomb on the radius-fixed profile013 mesh

This example reads
`Research/mesa_snec_profiles/mesa_detailed/profile013_radius_fixed.data` through the
framework MESA initializer. The default path is relative to the problem directory.
The file is also byte-identical to `profile013_radius_fixed.data` currently stored
at the FLARE repository root. The `day13` label refers here to the user's selected
profile013 file: its own header reports `time_seconds=69427.929964800016`, not
thirteen elapsed days, and its 4,147 data rows supersede the stale `num_zones=3563`
header value.

Execution order:

1. The optional unit hook selects the same solar code units as `CCSN_Simple`:
   M=GM_sun/G, L=6.957e10 cm, t=sqrt(L³/GM_sun).
2. The framework reads/reverses MESA rows, converts cgs density/pressure/velocity,
   clips with `mesa_radius_inner` and builds the spherical geometry.
3. The user initializer fills only the empty exterior padding cells, saves the
   pre-bomb state, and adds 1.6e51 erg to the first active shell. It never calls
   `build_mesh` or constructs a polytrope; all active faces remain unchanged.
4. Reflecting inner and zero-gradient, no-inflow outer boundaries are filled,
   then RK2/PLM evolves with central gravity and envelope self-gravity.

Parameters and units:

- `mesa_radius_inner=0.20711213012337312` is a native face at enclosed MESA mass
  4.001260592007234 M_sun, approximately the requested 4 M_sun cut. The same mass
  is used for the fixed central gravitational source.
- The reader's MESA solar-radius convention is 6.96e10 cm; this differs from the
  6.957e10 cm code length unit. The cut is about 1.441e10 cm.
- Retain 2,599 native shells, including their density, total pressure and radial
  velocity. The atmosphere extends from about 1.09146e14 cm to 1.25226e14 cm
  (`mesa_radius_outer=1800`), using the framework's 100 padding cells.
- Exterior density is **1e-15 g/cm³**, converted to code units; initial velocity
  is zero. A positive pressure is set with T=1e4 K and mean molecular weight 0.62.
- The shell-integrated envelope mass is about 8.21586 code masses. It need not
  exactly equal the difference of the MESA mass column: this reader preserves
  supplied radii and densities rather than remapping to the mass coordinate.

Run from this directory:

```sh
make -j4 test
```

The default run ends at **0.1 physical seconds after injection**, an early explosion
test, not shock breakout. The initial cell-centred acoustic CFL minimum is in
the interior near r = 2.46916e11 cm: dr = 962.77 cm and sound speed = 3.33338e7
cm/s give 0.2*dr/(cs+abs(v)) = 5.78 microseconds. The ~70.5 cm surface shell
has a lower sound speed and is not this minimum. Interface wave speeds make
the actual global steps roughly 5 microseconds. At that step size, simulating
days would require billions of steps. No shell coarsening or radial remeshing
is used.

`output/initial_envelope.csv` records the grid and state before injection;
`output/model.json` records unit conversions and bomb parameters. `check.py`
independently checks every retained face and MESA row, the energy increment,
positive states, mass conservation, and reaching the configured end time. It
writes `output/check.json`.

For the extended **1 physical second** run, use `parameters_long.txt`. It keeps
the same initial state, native grid, energy injection and explicit RK2 method,
but writes snapshots every 20,000 steps instead of 2,000, into `output_long/`.
The initial, first and final steps are also saved. The framework currently has
no checkpoint/restart reader, so this reruns from the original initial state;
the original `output/` remains available. To keep the text log equally sparse,
run the following in Bash from this directory:

```sh
set -o pipefail
stdbuf -oL ./flare -par parameters_long.txt 2>&1 |
  awk '/csv =|Loaded parameters|MESA CCSN|Error/ {print; fflush()}' > run_long.log
python3 check.py --parameters parameters_long.txt
```

This filters the saved log only; the numerical code still performs its usual
per-step diagnostics. Check `output_long/check.json` for completion and the
validation results after the run finishes.

The hydrodynamic model remains the Newtonian ideal gas used by `CCSN_Simple`:
MESA total pressure is converted using the configured gamma=5/3; MESA `gamma1`,
specific internal energy, composition and radiation EOS are not evolved. The
single-shell injection is very concentrated (~0.00778 M_sun), so any high speeds
must be interpreted within those model limitations.
