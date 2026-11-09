#!/usr/bin/env python3
"""Plot a CFL=1 timestep estimate for every zone in a MESA profile.

Usage:
    python plot_mesa_cfl.py profile013_radius_fixed.data
    python plot_mesa_cfl.py profile.data --sound-only --show

Requires numpy and matplotlib. Reads the profile directly (no mesa_reader).
Uses native dr [cm], csound [cm/s], and velocity [cm/s] columns. The default
estimate is dr / (csound + abs(velocity)) for a stationary radial mesh;
--sound-only gives dr / csound. These are profile-based estimates, not the
face-wave CFL limits of a numerical Riemann solver. No remeshing is applied.
Column units: https://github.com/MESAHub/mesa/blob/main/star/defaults/profile_columns.list
"""

import argparse
from pathlib import Path

import numpy as np


def read_profile(path):
    """Locate the zone table by column names; support Fortran D exponents."""
    with Path(path).open() as stream:
        for line in stream:
            names = line.split()
            if "zone" in names and "dr" in names and "csound" in names:
                if len(set(names)) != len(names):
                    raise ValueError("Duplicate profile column names")
                rows = (
                    row.replace("D", "E").replace("d", "e")
                    for row in stream
                    if row.strip() and not row.lstrip().startswith("#")
                )
                values = np.loadtxt(rows, ndmin=2)
                if values.size == 0 or values.shape[1] != len(names):
                    raise ValueError("Empty or malformed profile zone table")
                return dict(zip(names, values.T))
    raise ValueError("Profile must contain zone, dr, and csound columns")


def cfl_timestep(profile, sound_only=False):
    """Return one positive finite timestep [s] per zone, with CFL fixed at 1."""
    zone, dr, cs = (profile[name] for name in ("zone", "dr", "csound"))
    if not np.all(np.isfinite(zone)):
        raise ValueError("Non-finite zone numbers")
    if not np.all(np.isfinite(dr) & (dr > 0)):
        raise ValueError("Every dr must be finite and positive")
    if not np.all(np.isfinite(cs) & (cs > 0)):
        raise ValueError("Every csound must be finite and positive")
    if sound_only:
        speed = cs
    else:
        if "velocity" not in profile:
            raise ValueError("Missing velocity column; use --sound-only for dr/csound")
        velocity = profile["velocity"]
        if not np.all(np.isfinite(velocity)):
            raise ValueError("Every velocity must be finite")
        speed = cs + np.abs(velocity)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        dt = dr / speed
    if not np.all(np.isfinite(dt) & (dt > 0)):
        raise ValueError("Computed timesteps must be finite and positive")
    return dt


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("profile", type=Path, help="MESA profile.data file")
    parser.add_argument("-o", "--output", type=Path,
                        help="Output figure (default: <profile_stem>_cfl.png)")
    parser.add_argument("--sound-only", action="store_true",
                        help="Use dr/csound, ignoring bulk velocity")
    parser.add_argument("--show", action="store_true", help="Also display the plot")
    args = parser.parse_args()
    try:
        profile = read_profile(args.profile)
        dt = cfl_timestep(profile, args.sound_only)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    zone = profile["zone"]
    order = np.argsort(zone)
    minimum = int(np.argmin(dt))
    formula = r"$\Delta r/c_s$" if args.sound_only else r"$\Delta r/(c_s+|v|)$"
    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    ax.plot(zone[order], dt[order], lw=1.3, label=formula)
    ax.set_yscale("log")
    ax.axhline(dt[minimum], color="tab:red", ls="--", lw=1,
               label=f"Minimum: {dt[minimum]:.4g} s (zone {zone[minimum]:g})")
    ax.scatter(zone[minimum], dt[minimum], color="tab:red", s=25, zorder=3)
    ax.set_xlabel("MESA zone number (surface to interior)")
    ax.set_ylabel(r"CFL timestep $\Delta t$ [s]")
    ax.set_title(f"{args.profile.name} | CFL = 1.0")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()

    suffix = "_cfl_sound.png" if args.sound_only else "_cfl.png"
    output = args.output or args.profile.with_name(args.profile.stem + suffix)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    print(f"Zones: {len(dt)}")
    print(f"Minimum dt: {dt[minimum]:.10g} s at zone {zone[minimum]:g}")
    print(f"Maximum dt: {dt.max():.10g} s")
    print(f"Saved: {output.resolve()}")
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
