# File-defined 2D Cartesian advection

This example loads both geometry and primitive states through the framework's new
file initialization path. `initialize_units` selects M=2 g, L=3 cm, t=5 s and T=7 K
before construction. The input tables themselves are in **code units**; the generic
file reader does not apply cgs conversions. The later `initialize` callback only
checks the imported state and does not overwrite it.

- `mesh.dat`: whitespace table with `x y` columns, 33 active faces per axis.
  Both axes are nonuniform; the domain is [0,1] × [0,0.75]. Ghost faces are added by
  the framework. The current rectangular table format requires equal numbers of
  faces in all axis columns; it does not represent a 32×48 mesh.
- `profile.dat`: 1,024 active-cell rows, with x varying fastest. Headers may be in
  any order. All six `rho vx vy vz pre gamma` columns are required, even in 2D.
  The density is the cell average of `1 + 0.2 cos(2πx) cos(2πy/0.75)`;
  velocity is (0.2,-0.1), pressure 1 and gamma 1.4.
- Periodic boundaries, PLM/MC and RK2 evolve the wave to t=0.25.

From this directory:

```sh
make -j4 test
```

The checked-in tables are ready to use. `make inputs` regenerates them.
`check.py` verifies initial geometry/state mapping, unit metadata, positive states,
mass/momentum/energy conservation and the analytic advected cell averages. Results
are in `output/check.json`; snapshots are in `output/`.
