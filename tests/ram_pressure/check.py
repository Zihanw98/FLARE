#!/usr/bin/env python3
"""Run and check the default 2D ram pressure setup; plot the saved states."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hydro_suite import parameters, read_frame, require, snapshot_files


def check_boundaries(path, cfg):
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                names = line.strip().split(",")
                break
        keys = ("rho", "vx", "vy", "p", "internal_energy", "total_energy")
        data = np.loadtxt(stream, delimiter=",", usecols=[names.index(k) for k in keys])
    nx, ny = int(cfg["nxa"]), int(cfg["nya"])
    gx, gy = int(cfg["nxg"]), int(cfg["nyg"])
    data = data.reshape(ny + 2*gy, nx + 2*gx, len(keys))
    require(np.isfinite(data).all(), "Nonfinite primitive or energy field")
    require(np.all(data[..., 0] > 0) and np.all(data[..., 3] > 0), "Nonpositive ghost/active state")
    g = float(cfg["gamma"])
    p0 = (1./.8)**2 / g
    wind = np.array([1., 1., 0., p0, p0/(g-1), p0/(g-1)+.5])
    require(np.allclose(data[:, :gx], wind, rtol=0, atol=2e-13), "Left inflow changed")
    for j in range(ny + 2*gy):
        for i in range(gx, nx + 2*gx):
            if gy <= j < gy+ny and i < gx+nx:
                continue
            nearest = data[np.clip(j, gy, gy+ny-1), np.clip(i, gx, gx+nx-1)]
            require(np.array_equal(data[j, i], nearest), "Incorrect zero-gradient outflow")
    return data[gy:gy+ny, gx:gx+nx]


def plot(cfg, paths, first, last, output, animate):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import PillowWriter

    extent = [float(cfg[k]) for k in ("x_lower", "x_upper", "y_lower", "y_upper")]
    figure_width = max(11., 2. * (extent[1]-extent[0]) / (extent[3]-extent[2]))
    diagnostics = output / "diagnostics"
    fig, axes = plt.subplots(3, 1, figsize=(figure_width, 8), constrained_layout=True)
    for ax, key, label in zip(axes, ("rho", "vx", "p"), ("Density", "x velocity", "Pressure")):
        im = ax.imshow(first[key][0], origin="lower", extent=extent, aspect="equal",
                       interpolation="nearest", cmap="viridis",
                       vmin=0 if key != "rho" else float(first["rho"].min()),
                       vmax={"rho": float(first["rho"].max()), "vx": 1,
                             "p": (1./.8)**2/float(cfg["gamma"])}[key])
        ax.set(title=label, ylabel="y")
        fig.colorbar(im, ax=ax, pad=.015, shrink=.85)
    axes[-1].set_xlabel("x")
    fig.suptitle("Initial state: stationary light bubble, ambient Mach = 0.8", fontsize=14)
    fig.savefig(diagnostics / "initial.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(figure_width, 6), constrained_layout=True)
    vmax = max(float(first["rho"].max()), float(last["rho"].max()))
    vmin = min(float(first["rho"].min()), float(last["rho"].min()))
    for ax, frame in zip(axes, (first, last)):
        im = ax.imshow(frame["rho"][0], origin="lower", extent=extent, aspect="equal",
                       interpolation="nearest", cmap="magma", vmin=vmin, vmax=vmax)
        ax.set(title=f"Density, t = {frame['time']:.3f}", xlabel="x", ylabel="y")
        fig.colorbar(im, ax=ax, pad=.015, shrink=.9, label="Density")
    fig.savefig(diagnostics / "comparison.png", dpi=150)
    plt.close(fig)

    if animate:
        # Use the same run manifest as the checker, excluding old CSV files.
        chosen = np.unique(np.linspace(0, len(paths)-1, min(len(paths), 80), dtype=int))
        frames = [read_frame(paths[i], cfg)[1] for i in chosen]
        vmax = max(float(f["rho"].max()) for f in frames)
        vmin = min(float(f["rho"].min()) for f in frames)
        fig, ax = plt.subplots(figsize=(figure_width, 3.5), constrained_layout=True)
        im = ax.imshow(first["rho"][0], origin="lower", extent=extent, aspect="equal",
                       interpolation="nearest", cmap="magma", vmin=vmin, vmax=vmax)
        ax.set(xlabel="x", ylabel="y")
        fig.colorbar(im, ax=ax, pad=.015, shrink=.85, label="Density")
        destination = output / "animations" / "ram_pressure.gif"
        destination.parent.mkdir(parents=True, exist_ok=True)
        writer = PillowWriter(fps=10)
        with writer.saving(fig, str(destination), dpi=100):
            for frame in frames:
                im.set_data(frame["rho"][0])
                ax.set_title(f"Ram pressure | t = {frame['time']:.3f}")
                writer.grab_frame()
        plt.close(fig)
        print(f"Animation: {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=Path(__file__).with_name("parameters.txt"))
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--animate", action="store_true")
    args = parser.parse_args()
    path = args.parameters.resolve()
    cfg = parameters(path)
    output = path.parent / cfg["output_directory"]
    log = output / "run.log"
    if args.run:
        exe = (args.executable or Path(__file__).with_name("flare")).resolve()
        output.mkdir(parents=True, exist_ok=True)
        print(f"Running ram pressure case; progress: {log}", flush=True)
        with log.open("w") as stream:
            result = subprocess.run([str(exe), "-par", str(path)], cwd=path.parent,
                                    stdout=stream, stderr=subprocess.STDOUT)
        require(result.returncode == 0, f"Simulation failed; see {log}")
    paths = snapshot_files(output, cfg["output_suffix"], log)
    axes, first = read_frame(paths[0], cfg)
    require(first["step"] == 0 and first["time"] == 0, "Missing initial snapshot")
    x, y = np.meshgrid(axes[0], axes[1])
    cloud = x*x + y*y < 1.
    require(cloud.any() and (~cloud).any(), "Grid must resolve both cloud and ambient gas")
    p0 = (1./.8)**2 / float(cfg["gamma"])
    expected = {"rho": np.where(cloud, .1, 1.), "vx": np.where(cloud, 0., 1.),
                "vy": 0., "vz": 0., "p": p0}
    for key, value in expected.items():
        require(np.allclose(first[key][0], value, rtol=0, atol=2e-13), f"Incorrect initial {key}")
    wind_sound_speed = np.sqrt(float(cfg["gamma"]) * first["p"][0][~cloud] / first["rho"][0][~cloud])
    wind_mach = np.hypot(first["vx"][0][~cloud], first["vy"][0][~cloud]) / wind_sound_speed
    require(np.allclose(wind_mach, .8, rtol=0, atol=2e-13), "Initial wind must be Mach 0.8 in the ambient gas")
    initial = check_boundaries(paths[0], cfg)
    thermal = p0 / (float(cfg["gamma"])-1)
    require(np.allclose(initial[..., 4], thermal, rtol=0, atol=2e-13), "Incorrect initial internal energy")
    require(np.allclose(initial[..., 5], thermal + np.where(cloud, 0., .5), rtol=0, atol=2e-13),
            "Wind kinetic energy missing from initial total energy")
    min_rho, min_p, previous = float("inf"), float("inf"), -1.
    for source in paths:
        _, last = read_frame(source, cfg)
        require(last["time"] > previous, "Snapshot times must increase")
        previous = last["time"]
        min_rho = min(min_rho, float(last["rho"].min()))
        min_p = min(min_p, float(last["p"].min()))
    require(np.isclose(last["time"], float(cfg["t_end"]), rtol=0, atol=1e-12), "Run did not reach t_end")
    check_boundaries(paths[-1], cfg)
    report = {"snapshots": len(paths), "steps": last["step"], "time": last["time"],
              "initial_pressure_spread": float(np.ptp(first["p"])),
              "initial_cloud_max_speed": float(np.max(np.hypot(first["vx"][0][cloud], first["vy"][0][cloud]))),
              "initial_density_contrast": float(first["rho"][0][cloud].mean() / first["rho"][0][~cloud].mean()),
              "initial_bubble_density": float(first["rho"][0][cloud].mean()),
              "initial_wind_density": float(first["rho"][0][~cloud].mean()),
              "initial_wind_mach": float(wind_mach.mean()),
              "minimum_density": min_rho, "minimum_pressure": min_p,
              "initial_state": "PASS", "initial_energy": "PASS", "boundaries": "PASS"}
    diagnostics = output / "diagnostics"
    diagnostics.mkdir(parents=True, exist_ok=True)
    (diagnostics / "summary.json").write_text(json.dumps(report, indent=2)+"\n")
    plot(cfg, paths, first, last, output, args.animate)
    print("PASS " + json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
