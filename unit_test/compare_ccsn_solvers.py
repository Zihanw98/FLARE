#!/usr/bin/env python3
"""Compare saved CCSN solver runs, without assuming either solver's mass flux."""
import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ccsn_bounce_diagnostics import G, C, gravity_weights, read_snapshot

CASE = Path(__file__).resolve().parents[1] / "tests/CCSN_Hybrid_1D"
spec = importlib.util.spec_from_file_location("ccsn_analysis", CASE / "analyze.py")
case = importlib.util.module_from_spec(spec)
spec.loader.exec_module(case)


def load_run(directory):
    model = json.loads((directory / "model.json").read_text())
    grid = np.genfromtxt(directory / "grid.csv", delimiter=",", names=True)
    volume = grid["volume"]
    ga, gb = gravity_weights(grid)
    paths = sorted(directory.glob("bounce_*.csv"))
    if not paths:
        raise ValueError(f"No bounce snapshots in {directory}")
    frames = [read_snapshot(p) for p in paths]
    records = []
    for time, data in frames:
        rho, vel = data["rho"], data["vx"]
        mass = rho * volume
        interior = np.r_[0, np.cumsum(mass)[:-1]]
        potential = -G * np.sum(rho * interior * ga + rho * rho * gb)
        gas = np.sum(data["total_energy"] * volume)
        pc, uc, gc = case.cold_state(rho, model)
        thermal = data["total_energy"] - 0.5 * rho * vel**2 - uc
        assert all(np.all(np.isfinite(data[k])) for k in
                   ("rho", "vx", "p", "total_energy", "internal_energy"))
        assert np.all(rho > 0)
        assert np.all(thermal >= -1e-10 * uc)
        assert np.allclose(data["p"], pc + (model["gamma_thermal"]-1)*thermal, rtol=1e-9)
        shock = case.shock_measurement(data, model)
        records.append((time, rho[0], np.sum(mass), gas, potential, gas + potential,
                        shock[0], shock[1], shock[2], shock[3], shock[4], shock[6],
                        np.max(np.abs(vel))/C,
                        np.sum(data["dual_energy_correction"]*volume)))
    hist = np.asarray(records)
    peak = int(np.argmax(hist[:, 1]))
    crossed = np.flatnonzero(hist[:, 1] > model["rho_nuclear"])
    bounced = 0 < peak < len(hist)-1 and hist[-1, 1] < 0.9*hist[peak, 1]
    # Same finite-width front diagnostic as the case, only after central bounce.
    shock_mask = ((hist[:, 0] > hist[peak, 0]) & (hist[:, 7] > 0)
                  & (hist[:, 8] < 0) & (hist[:, 9] > 1.5)
                  & (hist[:, 10] > 2) & (hist[:, 11] > 1.1))
    ids = np.flatnonzero(shock_mask)
    drift = hist[:, 5]-hist[0, 5]
    summary = dict(cells=len(grid), snapshots=len(frames),
                   steps=int(paths[-1].stem.rsplit("_", 1)[1]),
                   final_time_s=float(hist[-1, 0]),
                   nuclear_crossing_s=float(hist[crossed[0], 0]) if len(crossed) else None,
                   resolved_bounce=bool(bounced),
                   bounce_time_s=float(hist[peak, 0]) if bounced else None,
                   peak_central_density_g_cm3=float(hist[peak, 1]),
                   final_central_density_g_cm3=float(hist[-1, 1]),
                   initial_mass_Msun=float(hist[0, 2]/1.98847e33),
                   max_relative_mass_error=float(np.max(np.abs(hist[:, 2]/hist[0, 2]-1))),
                   initial_total_energy_erg=float(hist[0, 5]),
                   initial_gravitational_energy_erg=float(hist[0, 4]),
                   final_total_energy_erg=float(hist[-1, 5]),
                   final_energy_drift_erg=float(drift[-1]),
                   final_drift_over_abs_initial_total=float(drift[-1]/abs(hist[0, 5])),
                   final_drift_over_abs_initial_gravity=float(drift[-1]/abs(hist[0, 4])),
                   max_abs_velocity_c=float(np.max(hist[:, 12])),
                   final_max_abs_velocity_c=float(hist[-1, 12]),
                   final_compressive_front_km=float(hist[-1, 6]/1e5),
                   first_outgoing_shock_s=float(hist[ids[0], 0]) if len(ids) else None,
                   last_valid_outgoing_shock_km=float(hist[ids[-1], 6]/1e5) if len(ids) else None,
                   max_abs_dual_energy_correction_erg=float(np.max(np.abs(hist[:, 13]))))
    # Energy leaving the outer boundary: recorded outflow state, trapezoidal
    # quadrature in saved time (not an exact discrete face-flux audit).
    outer_flux = np.array([4*np.pi*grid["r_right"][-1]**2 *
                          (d["total_energy"][-1]+d["p"][-1])*d["vx"][-1]
                          for _, d in frames])
    summary["approx_outer_gas_energy_outflow_erg"] = float(np.trapezoid(outer_flux, hist[:, 0]))
    timing = directory.parent / "timing.txt"
    if timing.exists():
        for line in timing.read_text().splitlines():
            fields = line.split()
            if len(fields) == 2 and fields[0] in ("real", "user", "sys"):
                summary[fields[0]+"_seconds"] = float(fields[1])
    return frames, hist, summary, shock_mask


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact", type=Path, required=True)
    parser.add_argument("--hll", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--end", type=float, default=0.0104)
    args = parser.parse_args()
    runs = {"Exact": load_run(args.exact), "HLL": load_run(args.hll)}
    assert np.array_equal(np.genfromtxt(args.exact/"grid.csv", delimiter=",", names=True),
                          np.genfromtxt(args.hll/"grid.csv", delimiter=",", names=True)), "grid"
    for key in ("rho", "vx", "p", "total_energy"):
        assert np.array_equal(runs["Exact"][0][0][1][key], runs["HLL"][0][0][1][key]), key
    args.output.mkdir(parents=True, exist_ok=True)
    common_end = min(run[1][-1, 0] for run in runs.values())
    report = {name: run[2] for name, run in runs.items()}
    report["requested_end_s"] = args.end
    report["both_completed"] = all(abs(run[1][-1, 0]-args.end) < 1e-12 for run in runs.values())
    report["identical_initial_states"] = True
    report["common_end_s"] = float(common_end)
    report["notes"] = [
        "Snapshots are sampled every 10 steps; bounce times use saved density maxima.",
        "Front radii use finite-width velocity compression, not a fitted discontinuity.",
        "Self-gravity energy is integrated exactly for piecewise constant shell densities.",
        "Energy drift is measured from snapshots; no HLL mass-flux reconstruction is applied to Exact.",
        "Partial runs, if present, are not a completed bounce comparison."]
    (args.output/"comparison.json").write_text(json.dumps(report, indent=2)+"\n")
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .2})
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for (name, (frames, hist, summary, shock_mask)), color in zip(runs.items(), ("#c84b31", "#247ba0")):
        style = "-" if name == "Exact" else "--"
        t = hist[:, 0]*1e3
        axes[0, 0].semilogy(t, hist[:, 1], label=name, color=color, ls=style)
        axes[0, 1].plot(t, (hist[:, 5]-hist[0, 5])/1e51, label=name, color=color, ls=style)
        axes[0, 2].plot(t[shock_mask], hist[shock_mask, 6]/1e5, label=name, color=color, ls=style)
        # Display the last saved profile no later than the common completed time.
        idx = np.flatnonzero(hist[:, 0] <= common_end+1e-14)[-1]
        time, data = frames[idx]
        label = f"{name}, {time*1e3:.4f} ms"
        radius = data["x"]/1e5
        visible = (radius >= 0.1) & (radius <= 1500)
        axes[1, 0].loglog(radius[visible], data["rho"][visible], label=label, color=color, ls=style)
        axes[1, 1].semilogx(radius[visible], data["vx"][visible]/C, label=label, color=color, ls=style)
        axes[1, 2].plot(t, hist[:, 12], label=name, color=color, ls=style)
        np.savetxt(args.output/(name.lower()+"_history.csv"), hist, delimiter=",", comments="",
                   header="time_s,central_density,mass_g,gas_erg,gravity_erg,total_erg,front_cm,"
                   "downstream_v,upstream_v,density_jump,pressure_jump,upstream_Mach,max_abs_v_c,dual_energy_erg")
    axes[0, 0].axhline(2e14, color="gray", ls=":", label="Nuclear density")
    titles = ("Central density", "Total energy drift (gas + self-gravity)", "Outgoing compressive front",
              "Density at common final time", "Velocity at common final time", "Maximum |velocity| / c")
    ylabels = ("Density [g/cm³]", "Energy change [10⁵¹ erg]", "Radius [km]", "Density [g/cm³]", "Velocity / c", "|Velocity| / c")
    for i, ax in enumerate(axes.flat):
        ax.set_title(titles[i]); ax.set_ylabel(ylabels[i]); ax.legend(fontsize=8)
        ax.set_xlabel("Radius [km]" if i in (3, 4) else "Time [ms]")
        if i in (3, 4): ax.set_xlim(0.1, 1500)
    status = "" if report["both_completed"] else " (PARTIAL RUN)"
    fig.suptitle("CCSN toy collapse: Exact versus HLL" + status
                 + " — identical grid, PCM, RK2 and CFL")
    fig.savefig(args.output/"comparison.png", dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
