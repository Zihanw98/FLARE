#!/usr/bin/env python3
"""Run/analyze the two 2D gravity cases using only snapshots listed in run.log.

Requires numpy and matplotlib. Hydrostatic checks are strict by default;
--diagnostic-only reports equilibrium drift without claiming preservation.
"""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parameters(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def snapshot(path, cfg, record):
    with path.open() as stream:
        metadata = dict(next(stream).strip()[2:].split("=", 1) for _ in range(3))
        require(int(metadata["step"]) == record[0], "Snapshot step/log mismatch")
        for key, value in zip(("time", "dt"), record[1:3]):
            require(math.isclose(float(metadata[key]), value, rel_tol=1e-9, abs_tol=1e-14),
                    "Snapshot time/log mismatch")
        data = np.genfromtxt(stream, delimiter=",", names=True)
    nx, ny, gx, gy = (int(cfg[k]) for k in ("nxa", "nya", "nxg", "nyg"))
    require(data.size == (nx+2*gx)*(ny+2*gy), f"Incomplete snapshot: {path}")
    require(all(np.isfinite(data[k]).all() for k in data.dtype.names), "Nonfinite snapshot")
    data = data[np.lexsort((data["x"], data["y"]))].reshape(ny+2*gy, nx+2*gx)
    for axis, size, ghost in (("x", nx, gx), ("y", ny, gy)):
        lo, hi = float(cfg[f"{axis}_lower"]), float(cfg[f"{axis}_upper"])
        centers = lo + (np.arange(-ghost, size+ghost)+.5)*(hi-lo)/size
        expected = centers[None, :] if axis == "x" else centers[:, None]
        require(np.allclose(data[axis], expected, rtol=0, atol=1e-12), "Incorrect coordinates")
    ghosts = np.ones(data.shape, dtype=bool)
    ghosts[gy:gy+ny, gx:gx+nx] = False
    require(np.array_equal(data["is_ghost"], ghosts), "Incorrect ghost flags")
    require(np.all(data["rho"] > 0) and np.all(data["p"] > 0), "Nonpositive density/pressure")
    require(np.all(data["gamma"] == float(cfg["gamma"])), "Incorrect gamma")
    return float(metadata["time"]), float(metadata["dt"]), data[gy:gy+ny, gx:gx+nx]


def plot(initial, final, metrics, directory):
    # Keep matplotlib's cache inside a writable disposable directory.
    with tempfile.TemporaryDirectory(prefix="flare-gravity-plot-") as cache:
        os.environ.setdefault("MPLCONFIGDIR", cache)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
        for ax, key in zip(axes[0], ("rho", "p", "vy")):
            for a, label in ((initial, "Initial"), (final, "Final")):
                ax.plot(a["y"][:, 0], a[key].mean(axis=1), label=label)
            ax.set(xlabel="y", ylabel=key, title=f"Horizontal mean: {key}")
            ax.legend()
        for ax, key in zip(axes[1], ("max_speed", "bottom_minus_top_p", "relative_total_energy_drift")):
            ax.plot([m["time"] for m in metrics], [m[key] for m in metrics])
            ax.set(xlabel="Time", ylabel=key)
        for ax in axes.flat:
            ax.grid(alpha=.2)
        fig.savefig(directory/"diagnostics.png", dpi=150)
        plt.close(fig)


def analyze(args, cfg, root):
    pattern = r"step = (\d+), time = ([^,]+), dt = ([^,]+).*?, csv = (.+)$"
    records = []
    for line in (root/"run.log").read_text().splitlines():
        match = re.match(pattern, line)
        if match:
            records.append((int(match[1]), float(match[2]), float(match[3]), root/match[4]))
    require(len(records) >= 3 and records[0][0:3] == (0, 0., 0.) and records[1][0] == 1,
            "Need initial, first-step and final snapshots in run.log")
    require(all(b[0] > a[0] and b[1] > a[1] for a, b in zip(records, records[1:])),
            "run.log must describe one run with increasing steps/times")
    require(math.isclose(records[-1][1], float(cfg["t_end"]), rel_tol=1e-9),
            "Run did not reach t_end")
    dx = (float(cfg["x_upper"])-float(cfg["x_lower"]))/int(cfg["nxa"])
    dy = (float(cfg["y_upper"])-float(cfg["y_lower"]))/int(cfg["nya"])
    gravity, gamma = float(cfg["gravity_y"]), float(cfg["gamma"])
    _, _, initial = snapshot(records[0][3], cfg, records[0])
    y = initial["y"] - float(cfg["y_lower"])
    z = .5*gravity*dy
    reference = (np.exp(gravity*y)*np.sinh(z)/z
                 if args.case == "hydrostatic" else np.ones(y.shape))
    for key in ("rho", "p"):
        require(np.allclose(initial[key], reference, rtol=1e-12, atol=1e-12),
                f"Wrong {args.case} initial {key}")
    require(np.all(initial["vx"] == 0) and np.all(initial["vy"] == 0), "Initial gas must be at rest")
    metrics = []
    for record in records:
        time, dt, a = snapshot(record[3], cfg, record)
        kinetic = .5*a["rho"]*(a["vx"]**2+a["vy"]**2)
        mass = a["rho"].sum()*dx*dy
        energy = (a["p"]/(gamma-1)+kinetic-gravity*a["rho"]*y).sum()*dx*dy
        if not metrics:
            mass0, energy0 = mass, energy
        require(abs(mass/mass0-1) < 1e-10, "Closed-box mass is not conserved")
        # Exact symmetry should be preserved by this grid-aligned, x-independent problem.
        require(max(np.ptp(a[k], axis=1).max() for k in ("rho", "p", "vx", "vy")) < 1e-10
                and np.max(np.abs(a["vx"])) < 1e-10, "Spurious horizontal structure")
        require(np.all(a["vz"] == 0), "Unexpected out-of-plane motion")
        if args.case == "adjustment" and record[0] == 1:
            require(np.allclose(a["vy"][1:-1], gravity*dt, rtol=1e-9, atol=1e-12),
                    "First-step interior acceleration disagrees with gravity")
        metrics.append({
            "step": record[0], "time": time, "mass": mass,
            "relative_mass_drift": mass/mass0-1,
            "gas_plus_potential_energy": energy,
            "relative_total_energy_drift": energy/energy0-1,
            "max_speed": np.max(np.hypot(a["vx"], a["vy"])),
            "max_mach": np.max(np.hypot(a["vx"], a["vy"])/np.sqrt(gamma*a["p"]/a["rho"])),
            "max_relative_rho_change": np.max(np.abs(a["rho"]/initial["rho"]-1)),
            "max_relative_p_change": np.max(np.abs(a["p"]/initial["p"]-1)),
            "bottom_minus_top_p": a["p"][0].mean()-a["p"][-1].mean(),
            "hydrostatic_gradient_residual": np.max(np.abs(
                np.gradient(a["p"].mean(axis=1), dy, edge_order=2)-gravity*a["rho"].mean(axis=1))),
        })
    output = root/cfg["output_directory"]/"diagnostics"
    output.mkdir(parents=True, exist_ok=True)
    with (output/"metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    plot(initial, a, metrics, output)
    print(f"Completed {args.case}: t={metrics[-1]['time']:.8g}, {records[-1][0]} steps")
    for key in ("max_speed", "max_mach", "max_relative_rho_change", "max_relative_p_change",
                "relative_mass_drift", "relative_total_energy_drift"):
        print(f"  max |{key}| over saved snapshots = {max(abs(m[key]) for m in metrics):.6g}")
    print(f"Diagnostics: {output}")
    # A finite-volume gravity coupling need not conserve gas+potential energy exactly.
    # 0.1% is an explicit coarse regression bound, not an exact-conservation claim.
    require(max(abs(m["relative_total_energy_drift"]) for m in metrics) < 1e-3,
            "Gas+potential energy drift exceeds the 0.1% regression bound")
    if args.case == "adjustment":
        weight = abs(gravity)*mass0/(dx*int(cfg["nxa"]))
        require(max(m["bottom_minus_top_p"] for m in metrics) > .25*weight,
                "No appreciable pressure response to column weight; run longer or check gravity")
        print("PASS: acceleration, pressure response, symmetry, positivity and conservation bounds")
        print("This does not assert relaxation to a static endpoint.")
    else:
        error = max(m[k] for m in metrics for k in
                    ("max_mach", "max_relative_rho_change", "max_relative_p_change"))
        if args.diagnostic_only:
            print(f"DIAGNOSTIC ONLY: equilibrium error={error:.6g}; preservation threshold=1e-10")
        else:
            require(error <= 1e-10, f"Hydrostatic equilibrium NOT preserved: error={error:.6g} > 1e-10. "
                    "The current Godunov/source/reflective-wall combination is not well-balanced.")
            print("PASS: hydrostatic equilibrium preserved")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("adjustment", "hydrostatic"), required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--diagnostic-only", action="store_true")
    args = parser.parse_args()
    try:
        path = args.parameters.resolve()
        cfg, root = parameters(path), path.parent
        require(int(cfg["nxa"]) >= 2 and int(cfg["nya"]) >= 4, "Need at least 2x4 active cells")
        require(cfg["constant_gravity_enabled"].lower() in ("true", "1")
                and float(cfg["gravity_x"]) == 0 and float(cfg["gravity_y"]) < 0,
                "Require enabled, downward vertical gravity")
        require(cfg["time_integrator"].lower() in ("euler", "rk2"),
                "These checks support Euler and RK2")
        for axis, bc in (("x", "periodic"), ("y", "reflective")):
            require(all(cfg[f"boundary_{axis}_{side}"] == bc for side in ("lower", "upper")),
                    f"Require {bc} boundaries in {axis}")
        if args.run:
            require(args.executable is not None, "--run requires --executable")
            print(f"Running {args.case}; progress: {root/'run.log'}", flush=True)
            with (root/"run.log").open("w") as stream:
                result = subprocess.run([str(args.executable.resolve()), "-par", str(path)], cwd=root,
                                        stdout=stream, stderr=subprocess.STDOUT)
            require(result.returncode == 0, f"Simulation failed; see {root/'run.log'}")
        analyze(args, cfg, root)
    except (OSError, ValueError, KeyError, StopIteration) as error:
        parser.exit(1, f"FAIL: {error}\n")


if __name__ == "__main__":
    main()
