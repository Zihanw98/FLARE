#!/usr/bin/env python3
"""Check and animate a 2D double blast using only CSVs listed by this run."""
import argparse
from contextlib import ExitStack
import csv
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from animate import check_grid, edges, read_snapshot
from animate3d import snapshot_files


def require(condition, message):
    if not condition:
        raise ValueError(message)


def plot(frames, records, output, make_movie):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.animation import PillowWriter, FFMpegWriter, writers
    from matplotlib.colors import LogNorm

    first = frames[0]
    xe, ye = edges(first.x), edges(first.y)
    mid = int(np.argmin(abs(first.y - (ye[0] + ye[-1]) / 2)))
    rho_min = min(frame.values[0].min() for frame in frames)
    rho_max = max(frame.values[0].max() for frame in frames)
    p_min = min(frame.values[3].min() for frame in frames)
    p_max = max(frame.values[3].max() for frame in frames)
    times = np.array([frame.time for frame in frames])
    center_p = np.array([row["center_pressure"] for row in records])
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.4), layout="constrained",
                           gridspec_kw={"height_ratios": (1.65, 1)})
    density = ax[0,0].pcolormesh(xe, ye, first.values[0], cmap="inferno", vmin=rho_min, vmax=rho_max)
    pressure = ax[0,1].pcolormesh(xe, ye, first.values[3], cmap="magma", norm=LogNorm(p_min, p_max))
    fig.colorbar(density, ax=ax[0,0], label="Density")
    fig.colorbar(pressure, ax=ax[0,1], label="Pressure (log scale)")
    for panel, label in zip(ax[0], ("Density", "Pressure")):
        panel.set(xlabel="x", ylabel="y", title=label, aspect="equal")
        panel.axhline(first.y[mid], color="white", alpha=.5, ls="--", lw=.7)
    profile, = ax[1,0].plot(first.x, first.values[3,mid], color="#ab2669", lw=1.7)
    ax[1,0].set(xlim=(xe[0], xe[-1]), ylim=(0, p_max*1.05), xlabel="x", ylabel="p",
                title=f"Horizontal pressure cut, y={first.y[mid]:.5f}")
    ax[1,1].plot(times, center_p, color="#176b87", lw=1.6)
    marker, = ax[1,1].plot([times[0]], [center_p[0]], "o", color="#c9472c")
    cursor = ax[1,1].axvline(times[0], color="#c9472c", alpha=.4)
    ax[1,1].axhline(1., color="gray", ls="--", lw=.7)
    ax[1,1].set(xlim=(times[0], times[-1]), ylim=(0, center_p.max()*1.1),
                xlabel="Time", ylabel="Mean p in central 2 x 2 cells", title="Collision-region pressure")
    for panel in ax[1]:
        panel.grid(alpha=.2)
    title = fig.suptitle("")

    def update(i):
        frame = frames[i]
        density.set_array(frame.values[0].ravel())
        pressure.set_array(frame.values[3].ravel())
        profile.set_ydata(frame.values[3,mid])
        marker.set_data([frame.time], [center_p[i]])
        cursor.set_xdata([frame.time, frame.time])
        title.set_text(f"Two expanding blasts | t = {frame.time:.5f} | step {frame.step}")

    diagnostics = output / "diagnostics"
    update(len(frames)-1)
    fig.savefig(diagnostics / "final.png", dpi=130)
    peak = int(center_p.argmax())
    update(peak)
    fig.savefig(diagnostics / "collision_peak.png", dpi=130)
    if make_movie:
        directory = output / "animations"
        directory.mkdir(exist_ok=True)
        movie_writers = [(PillowWriter(fps=12), directory / "double_blast.gif")]
        if writers.is_available("ffmpeg"):
            movie_writers.append((FFMpegWriter(fps=12, codec="libx264", extra_args=["-pix_fmt", "yuv420p"]),
                                  directory / "double_blast.mp4"))
        with ExitStack() as stack:
            for writer, path in movie_writers:
                stack.enter_context(writer.saving(fig, str(path), dpi=90))
            for i in range(len(frames)):
                update(i)
                for writer, _ in movie_writers:
                    writer.grab_frame()
                if i % 20 == 0:
                    print(f"Rendered {i+1}/{len(frames)} frames", flush=True)
        print(f"Saved animations in {directory}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=Path(__file__).with_name("parameters.txt"))
    parser.add_argument("--executable", type=Path, default=Path(__file__).with_name("flare"))
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--animate", action="store_true")
    args = parser.parse_args()
    parameters = args.parameters.resolve()
    cfg = {}
    for line in parameters.read_text().splitlines():
        line = line.split("#",1)[0].strip()
        if line:
            key, value = line.split("=",1)
            cfg[key.strip()] = value.strip()
    output = parameters.parent / cfg["output_directory"]
    log = output / "run.log"
    if args.run:
        output.mkdir(parents=True, exist_ok=True)
        with log.open("w") as stream:
            result = subprocess.run([str(args.executable.resolve()), "-par", str(parameters)],
                                    cwd=parameters.parent, stdout=stream, stderr=subprocess.STDOUT)
        require(result.returncode == 0, f"Simulation failed; inspect {log}")
    paths = snapshot_files(output, "mesh_" + cfg["output_suffix"], log)
    frames = [read_snapshot(path, int(path.stem.rsplit("_",1)[1])) for path in paths]
    first = frames[0]
    require(first.step == 0 and first.time == 0, "Missing initial snapshot")
    require(tuple(field[0] for field in first.fields) == ("rho", "vx", "vy", "p"), "Expected 2D fields")
    nx, ny = int(cfg["nxa"]), int(cfg["nya"])
    require(first.values.shape == (4, ny, nx) and nx % 2 == ny % 2 == 0,
            "Checker requires an even number of active cells in both directions")
    lower = np.array([float(cfg[f"{d}_lower"]) for d in "xy"])
    upper = np.array([float(cfg[f"{d}_upper"]) for d in "xy"])
    lengths = upper - lower
    for d, (axis, n) in enumerate(((first.x, nx), (first.y, ny))):
        require(np.allclose(axis, lower[d]+(np.arange(n)+.5)*lengths[d]/n, rtol=0, atol=1e-13), "Wrong grid coordinates")
    x, y = np.meshgrid(first.x, first.y)
    radius = .08*min(lengths)
    hot = (np.minimum((x-lower[0]-.375*lengths[0])**2, (x-lower[0]-.625*lengths[0])**2)
           + (y-lower[1]-.5*lengths[1])**2 < radius**2)
    require(np.array_equal(first.values[3], np.where(hot, 10., 1.)), "Incorrect initial pressure disks")
    require(np.all(first.values[0] == 1) and np.all(first.values[1:3] == 0), "Incorrect initial rho/velocity")
    gamma = float(cfg["gamma"])
    area = np.diff(edges(first.y))[:,None] * np.diff(edges(first.x))[None,:]
    records = []
    reference = None
    for frame in frames:
        check_grid(frame, first)
        rho, vx, vy, p = frame.values
        require(np.all(rho > 0) and np.all(p > 0), "Nonpositive density/pressure")
        totals = np.array([np.sum(rho*area), np.sum((p/(gamma-1)+.5*rho*(vx*vx+vy*vy))*area)])
        if reference is None:
            reference = totals
        drift = np.abs(totals/reference-1)
        require(drift.max() < 1e-10, f"Conservation failure at step {frame.step}: {drift}")
        symmetry = 0.
        for f, value in enumerate(frame.values):
            for axis, component in ((1,1), (0,2)):
                sign = -1 if f == component else 1
                error = np.max(abs(value-sign*np.flip(value,axis))) / max(1.,np.max(abs(value)))
                symmetry = max(symmetry,error)
        require(symmetry < 1e-9, f"Reflection symmetry failure at step {frame.step}")
        center = np.s_[ny//2-1:ny//2+1, nx//2-1:nx//2+1]
        wall_p = max(np.max(abs(p[0]-1)),np.max(abs(p[-1]-1)),np.max(abs(p[:,0]-1)),np.max(abs(p[:,-1]-1)))
        records.append(dict(step=frame.step, time=frame.time, mass=totals[0], energy=totals[1],
                            mass_relative_drift=drift[0], energy_relative_drift=drift[1], symmetry_error=symmetry,
                            center_pressure=p[center].mean(), center_density=rho[center].mean(),
                            rho_min=rho.min(), rho_max=rho.max(), p_min=p.min(), p_max=p.max(),
                            wall_pressure_deviation=wall_p))
    require(np.isclose(frames[-1].time,float(cfg["t_end"]),rtol=0,atol=1e-12), "Did not reach t_end")
    # These distinguish an actual collision from two unchanged high-pressure disks.
    require(max(row["center_pressure"] for row in records) > 1.5, "No central pressure rise: run long enough for collision")
    require(max(row["center_density"] for row in records) > 1.2, "No central compression detected")
    diagnostics = output / "diagnostics"
    diagnostics.mkdir(exist_ok=True)
    with (diagnostics / "metrics.csv").open("w",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    peak = max(records,key=lambda row:row["center_pressure"])
    print(f"PASS: {nx}x{ny}, {len(frames)} frames, {frames[-1].step} steps, t={frames[-1].time:g}",flush=True)
    for key in ("mass_relative_drift", "energy_relative_drift", "symmetry_error", "wall_pressure_deviation"):
        print(f"  max {key} = {max(row[key] for row in records):.6e}",flush=True)
    print(f"  central pressure peak = {peak['center_pressure']:.6g} at saved t={peak['time']:.6g}",flush=True)
    plot(frames,records,output,args.animate)


if __name__ == "__main__":
    main()
