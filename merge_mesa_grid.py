#!/usr/bin/env python3
"""Merge adjacent MESA shells and compare radial profiles (CFL = 1).

Example:
    python merge_mesa_grid.py profile013_radius_fixed.data --min-dt 0.01

Requires numpy and matplotlib, plus plot_mesa_cfl.py in the same directory.
The output is a diagnostic resampling, not a MESA restart model: no EOS or
hydrostatic equilibrium is solved. Density and pressure are volume averages;
temperature and velocity use rho*dV weights. Enclosed mass is retained at
each surviving outer face, independently of the density integral. We do not
use the input dm column (some converted profiles contain q in that column).

Shell faces come from logR, using FLARE's default solar radius (6.96e10 cm).
Only the innermost missing face uses dr. The CFL diagnostic after merging
is total width / max(original csound + abs(velocity)) within each group,
not a sound speed recomputed from an EOS. --sound-only omits velocity.
"""

import argparse
from pathlib import Path

import numpy as np

from plot_mesa_cfl import cfl_timestep, read_profile


def prepare_profile(profile, solar_radius_cm=6.96e10, sound_only=False):
    """Return physical quantities in input (surface-to-interior) order."""
    if not np.isfinite(solar_radius_cm) or solar_radius_cm <= 0:
        raise ValueError("solar-radius-cm must be finite and positive")
    required = ("zone", "logR", "dr", "logT", "mass", "logRho", "logP",
                "velocity", "csound")
    for name in required:
        if name not in profile or not np.all(np.isfinite(profile[name])):
            raise ValueError(f"Missing or non-finite column: {name}")
    cfl_timestep(profile, sound_only)  # Validate positive widths and sound speeds.
    outer = 10.0 ** profile["logR"] * solar_radius_cm
    inner_edge = outer[-1] - profile["dr"][-1]
    tolerance = 16 * np.finfo(float).eps * outer[-1]
    if inner_edge < -tolerance:
        raise ValueError("Innermost dr exceeds its radius; check solar-radius-cm")
    inner = np.r_[outer[1:], max(0.0, inner_edge)]
    if not np.all(np.isfinite(outer) & (outer > inner)):
        raise ValueError("Radii must decrease strictly from surface to interior")
    if not np.all(profile["mass"] > 0) or np.any(np.diff(profile["mass"]) > 0):
        raise ValueError("Enclosed mass must be positive and non-increasing")
    result = {
        "r_inner_cm": inner,
        "r_outer_cm": outer,
        "T_K": 10.0 ** profile["logT"],
        "M_Msun": profile["mass"],
        "rho_g_cm3": 10.0 ** profile["logRho"],
        "P_dyn_cm2": 10.0 ** profile["logP"],
        "v_cm_s": profile["velocity"],
        "signal_max_cm_s": profile["csound"] + (
            0 if sound_only else np.abs(profile["velocity"])),
    }
    if any(not np.all(np.isfinite(a)) for a in result.values()):
        raise ValueError("Non-finite physical quantity after unit conversion")
    for name in ("T_K", "rho_g_cm3", "P_dyn_cm2"):
        if np.any(result[name] <= 0):
            raise ValueError(f"{name} must be positive")
    return result


def shell_volumes(inner, outer):
    """Factored difference avoids cancellation for extremely thin shells."""
    return (4 * np.pi / 3) * (outer - inner) * (outer**2 + outer * inner + inner**2)


def merge_groups(data, min_dt):
    """Greedy contiguous groups, surface inward; absorb a short final group.

    A group is accepted as soon as its width / maximum signal speed meets
    min_dt. A remaining short tail is merged back until it also qualifies.
    This deterministic grouping does not claim a globally optimal mesh.
    """
    if not np.isfinite(min_dt) or min_dt <= 0:
        raise ValueError("min-dt must be finite and positive")
    outer, inner, speed = (data[k] for k in (
        "r_outer_cm", "r_inner_cm", "signal_max_cm_s"))
    groups, start, maximum = [], 0, 0.0
    for stop in range(1, len(outer) + 1):
        maximum = max(maximum, speed[stop - 1])
        if (outer[start] - inner[stop - 1]) / maximum >= min_dt:
            groups.append((start, stop))
            start, maximum = stop, 0.0
    if start < len(outer):
        while groups and (outer[start] - inner[-1]) / maximum < min_dt:
            start, previous_stop = groups.pop()
            maximum = max(maximum, np.max(speed[start:previous_stop]))
        if (outer[start] - inner[-1]) / maximum < min_dt:
            raise ValueError("Target cannot be reached even by merging all shells")
        groups.append((start, len(outer)))
    return groups


def resample(data, zones, groups):
    """Average shell quantities and retain exact group boundary locations."""
    volume = shell_volumes(data["r_inner_cm"], data["r_outer_cm"])
    weights = volume * data["rho_g_cm3"]
    rows = []
    for start, stop in groups:
        s = slice(start, stop)
        inner, outer = data["r_inner_cm"][stop - 1], data["r_outer_cm"][start]
        speed = np.max(data["signal_max_cm_s"][s])
        rows.append({
            "zone_first": zones[start], "zone_last": zones[stop - 1],
            "n_zones": stop - start,
            "r_inner_cm": inner, "r_outer_cm": outer,
            "r_mid_cm": inner + 0.5 * (outer - inner), "dr_cm": outer - inner,
            "T_K": np.average(data["T_K"][s], weights=weights[s]),
            "M_Msun": data["M_Msun"][start],
            "rho_g_cm3": np.average(data["rho_g_cm3"][s], weights=volume[s]),
            "P_dyn_cm2": np.average(data["P_dyn_cm2"][s], weights=volume[s]),
            "v_cm_s": np.average(data["v_cm_s"][s], weights=weights[s]),
            "signal_max_cm_s": speed, "dt_cfl_s": (outer - inner) / speed,
        })
    return {name: np.array([row[name] for row in rows]) for name in rows[0]}


def save_table(path, data):
    np.savetxt(path, np.column_stack(list(data.values())), fmt="%.17e",
               header=" ".join(data))


def plot_profiles(original, merged, output_dir, title, min_dt, show=False):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [
        ("T_K", "Temperature [K]", "R-T", True),
        ("M_Msun", r"Enclosed mass [$M_\odot$]", "R-M", True),
        ("rho_g_cm3", r"Density [g cm$^{-3}$]", "R-rho", True),
        ("P_dyn_cm2", r"Pressure [dyn cm$^{-2}$]", "R-pre", True),
        ("v_cm_s", r"Velocity [cm s$^{-1}$]", "R-velocity", False),
        ("dt_cfl_s", "CFL timestep [s]", "R-cfl", True),
    ]

    def draw(ax, panel):
        key, label, name, log_y = panel
        # Enclosed M belongs to the outer face; thermodynamic fields to a shell.
        radial_key = "r_outer_cm" if key == "M_Msun" else "r_mid_cm"
        for data, legend, style in (
                (original, f"Original ({len(original[key])} zones)", "-"),
                (merged, f"Merged ({len(merged[key])} zones)", "--")):
            ax.plot(data[radial_key][::-1], data[key][::-1], style,
                    lw=1.3, label=legend)
        ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        else:
            ax.axhline(0, color="0.5", lw=0.7)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        if key == "dt_cfl_s":
            ax.axhline(min_dt, color="tab:red", ls=":", label=f"Target: {min_dt:g} s")
        ax.set(xlabel="Radius [cm]", ylabel=label, title=name)
        ax.grid(True, which="both", alpha=0.2)
        ax.legend(fontsize=8)

    overview, axes = plt.subplots(3, 2, figsize=(13, 12), layout="constrained")
    for ax, panel in zip(axes.flat, panels):
        draw(ax, panel)
        fig, separate_ax = plt.subplots(figsize=(8, 5), layout="constrained")
        draw(separate_ax, panel)
        fig.savefig(output_dir / f"{panel[2]}.png", dpi=180)
        plt.close(fig)
    overview.suptitle(f"{title} | CFL = 1 | target dt = {min_dt:g} s")
    overview.savefig(output_dir / "profiles_comparison.png", dpi=180)
    if show:
        plt.show()
    plt.close(overview)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("profile", type=Path)
    parser.add_argument("--min-dt", type=float, default=0.01, help="Target CFL timestep [s]")
    parser.add_argument("--solar-radius-cm", type=float, default=6.96e10)
    parser.add_argument("--sound-only", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    try:
        profile = read_profile(args.profile)
        data = prepare_profile(profile, args.solar_radius_cm, args.sound_only)
        groups = merge_groups(data, args.min_dt)
        original = resample(data, profile["zone"], [(i, i + 1) for i in range(len(profile["zone"]))])
        merged = resample(data, profile["zone"], groups)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))

    output_dir = args.output_dir or args.profile.with_name(args.profile.stem + "_merged")
    output_dir.mkdir(parents=True, exist_ok=True)
    save_table(output_dir / "original_profile.dat", original)
    save_table(output_dir / "merged_profile.dat", merged)
    plot_profiles(original, merged, output_dir, args.profile.name, args.min_dt, args.show)
    summary = (
        f"Input: {args.profile.resolve()}\n"
        f"Zones: {len(original['T_K'])} -> {len(merged['T_K'])}\n"
        f"Minimum CFL dt [s]: {original['dt_cfl_s'].min():.10g} -> {merged['dt_cfl_s'].min():.10g}\n"
        f"Target dt [s]: {args.min_dt:g}; CFL = 1; sound-only = {args.sound_only}\n"
        f"Solar radius [cm]: {args.solar_radius_cm:.10g}\n\n"
        "Geometry: logR gives each outer face; the innermost missing face uses dr.\n"
        "CFL before/after uses face differences, which can differ from the input dr column.\n"
        "Merged CFL = group width / maximum original signal speed in that group.\n"
        "rho and P: volume-weighted; T and velocity: rho*dV-weighted.\n"
        "This preserves the density-integrated mass and radial momentum.\n"
        "M: original enclosed mass at each retained outer face, not recomputed from rho.\n"
        "The dm column is unused. Original and merged tables include source zone ranges.\n"
        "Thermodynamic profiles use shell midpoints; R-M uses outer face radii.\n"
        "These are diagnostic averages, not an EOS-consistent MESA restart model.\n"
        "Internal/total energy and hydrostatic equilibrium are not reconstructed.\n"
    )
    (output_dir / "README.txt").write_text(summary)
    print(summary.split("\n\n")[0])
    print(f"Output: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
