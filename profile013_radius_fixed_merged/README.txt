Input: /home/u0190250/Research/FLARE/profile013_radius_fixed.data
Zones: 4147 -> 3892
Minimum CFL dt [s]: 6.633202931e-08 -> 0.01001492259
Target dt [s]: 0.01; CFL = 1; sound-only = False
Solar radius [cm]: 6.96e+10

Geometry: logR gives each outer face; the innermost missing face uses dr.
CFL before/after uses face differences, which can differ from the input dr column.
Merged CFL = group width / maximum original signal speed in that group.
rho and P: volume-weighted; T and velocity: rho*dV-weighted.
This preserves the density-integrated mass and radial momentum.
M: original enclosed mass at each retained outer face, not recomputed from rho.
The dm column is unused. Original and merged tables include source zone ranges.
Thermodynamic profiles use shell midpoints; R-M uses outer face radii.
These are diagnostic averages, not an EOS-consistent MESA restart model.
Internal/total energy and hydrostatic equilibrium are not reconstructed.
