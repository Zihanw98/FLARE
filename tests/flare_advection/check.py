#!/usr/bin/env python3
"""Check periodic translation of cold FLARE letters in true vacuum.

    make test      # compile, run, validate, save comparison.png + metrics.csv
    make animate   # read this run's CSVs and also save advection.gif
    python3 check.py --run --animate

Dependencies: numpy, matplotlib, Pillow (requirements.txt).

Two references are used: rigid periodic translation of the initial cell-average
image, and an independently computed spatially first-order upwind solution using
the configured Euler/RK2 method and logged time steps. The first is a
cold/pressureless limit, NOT an exact positive-pressure Euler solution.
Numerical diffusion relative to rigid translation is reported;
it does not by itself fail this translation test. Checks target transport phase,
conservation, vacuum states, periodic ghosts, and agreement with cold upwinding.
"""

import argparse
import csv
import math
from pathlib import Path
import re
import subprocess

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.animation import PillowWriter
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


def read_snapshot(path, cfg):
    """Return metadata and primitive arrays including validated periodic ghosts."""
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                columns = line.strip().split(",")
                break
            key, value = line[1:].strip().split("=", 1)
            metadata[key] = int(value) if key == "step" else float(value)
        else:
            raise ValueError(f"Missing CSV table: {path}")
        required = ("x", "y", "z", "gamma", "rho", "vx", "vy", "vz", "p", "is_ghost")
        data = np.loadtxt(stream, delimiter=",", ndmin=2,
                          usecols=[columns.index(key) for key in required])
    require(set(metadata) == {"step", "time", "dt"}, f"Missing snapshot metadata: {path}")
    require(all(math.isfinite(v) and v >= 0 for v in metadata.values()), "Invalid metadata")
    nx, ny, gx, gy = (int(cfg[key]) for key in ("nxa", "nya", "nxg", "nyg"))
    require(data.shape == ((nx+2*gx)*(ny+2*gy), len(required)), f"Incomplete CSV: {path}")
    require(np.isfinite(data).all(), f"Non-finite snapshot: {path}")
    data = data[np.lexsort((data[:, 0], data[:, 1]))].reshape(ny+2*gy, nx+2*gx, -1)
    for axis, size, ghost, column in (("x", nx, gx, 0), ("y", ny, gy, 1)):
        lower, upper = float(cfg[f"{axis}_lower"]), float(cfg[f"{axis}_upper"])
        centers = lower + (np.arange(-ghost, size+ghost)+.5)*(upper-lower)/size
        target = centers[None, :] if axis == "x" else centers[:, None]
        require(np.allclose(data[:, :, column], target, rtol=0, atol=1e-12), "Wrong grid coordinates")
    mask = np.ones(data.shape[:2], dtype=bool)
    mask[gy:gy+ny, gx:gx+nx] = False
    require(np.array_equal(data[:, :, -1], mask), "Incorrect ghost flags")
    require(np.all(data[:, :, (2, 7)] == 0), "Unexpected z/vz in 2D")
    require(np.all(data[:, :, 3] == float(cfg["gamma"])), "Unexpected gamma")
    active = data[gy:gy+ny, gx:gx+nx]
    wrapped = active[np.arange(-gy, ny+gy) % ny][:, np.arange(-gx, nx+gx) % nx]
    require(np.array_equal(data[:, :, 3:9], wrapped[:, :, 3:9]), "Periodic ghosts are stale or incorrect")
    rho, vx, vy, p = (active[:, :, k] for k in (4, 5, 6, 8))
    require(np.all(rho >= 0), "Negative density")
    require(np.all(p[rho > 0] > 0), "Gas pressure must remain positive")
    require(np.all(active[rho == 0][:, [5, 6, 7, 8]] == 0), "Vacuum must have zero velocity and pressure")
    return metadata, (rho, vx, vy, p)


def totals(values, gamma, cell_area):
    rho, vx, vy, p = values
    return np.array([rho.sum(), (rho*vx).sum(), (rho*vy).sum(),
                     (p/(gamma-1)+.5*rho*(vx*vx+vy*vy)).sum()])*cell_area


def translate(image, sx, sy):
    """Exact cell-average overlap for a translated piecewise-constant image."""
    ix, iy = math.floor(sx), math.floor(sy)
    fx, fy = sx-ix, sy-iy
    shifted = (1-fx)*np.roll(image, ix, axis=1)+fx*np.roll(image, ix+1, axis=1)
    return (1-fy)*np.roll(shifted, iy, axis=0)+fy*np.roll(shifted, iy+1, axis=0)


def upwind_step(image, dt, velocity, spacing):
    """Independent unsplit scalar upwind update, with either velocity sign."""
    cx, cy = np.abs(velocity)*dt/spacing
    require(cx+cy <= 1+1e-9, "Unstable scalar reference CFL")
    return ((1-cx-cy)*image
            + cx*np.roll(image, 1 if velocity[0] >= 0 else -1, axis=1)
            + cy*np.roll(image, 1 if velocity[1] >= 0 else -1, axis=0))


def modes(image):
    """First Fourier coefficient along each axis; its phase measures translation."""
    ny, nx = image.shape
    return np.array([np.sum(image*np.exp(-2j*np.pi*np.arange(nx)/nx)[None, :]),
                     np.sum(image*np.exp(-2j*np.pi*np.arange(ny)/ny)[:, None])])


def plot_results(frames, initial, cfg, output, animate, fps):
    bounds = tuple(float(cfg[key]) for key in ("x_lower", "x_upper", "y_lower", "y_upper"))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), layout="constrained")
    names = ("Numerical density", "Rigid translation (cold limit)", "First-order upwind reference")
    images = []
    for axis, name in zip(axes, names):
        image = axis.imshow(initial, origin="lower", extent=bounds, aspect="equal",
                            vmin=0, vmax=initial.max(), cmap="inferno", interpolation="nearest")
        axis.set(title=name, xlabel="x", ylabel="y")
        images.append(image)
    fig.colorbar(images[0], ax=list(axes), shrink=.65, label="rho")
    title = fig.suptitle("", fontsize=13)

    def update(frame):
        for image, data in zip(images, (frame["rho"], frame["rigid"], frame["upwind"])):
            image.set_data(data)
        title.set_text(f"FLARE in vacuum | step {frame['step']} | t = {frame['time']:.6f}")

    if animate:
        indices = np.linspace(0, len(frames)-1, min(len(frames), 100), dtype=int)
        writer = PillowWriter(fps=fps)
        with writer.saving(fig, str(output/"advection.gif"), dpi=100):
            for index in indices:
                update(frames[index])
                writer.grab_frame()
        print(f"Saved {output/'advection.gif'}")
    update(frames[0])
    fig.savefig(output/"initial.png", dpi=120)
    update(frames[-1])
    fig.savefig(output/"comparison.png", dpi=120)
    plt.close(fig)


def analyze(cfg, problem_dir, log_path, animate, fps):
    method = cfg.get("time_integrator", "euler").lower()
    require(method in ("euler", "rk2"), "Cold upwind reference supports Euler and RK2")
    for axis in ("x", "y"):
        require(all(cfg[f"boundary_{axis}_{side}"] == "periodic" for side in ("lower", "upper")),
                "Translation check requires periodic boundaries on both axes")
    records = []
    pattern = r"step = (\d+), time = ([^,]+), dt = ([^,]+)"
    for line in log_path.read_text().splitlines():
        match = re.match(pattern, line)
        if match:
            filename = line.split(", csv = ", 1)[1] if ", csv = " in line else None
            records.append((int(match[1]), float(match[2]), float(match[3]), filename))
    require(len(records) > 1 and [r[0] for r in records] == list(range(len(records))),
            "run.log must contain one complete run with consecutive steps")
    require(records[0][1:3] == (0., 0.) and records[0][3] and records[-1][3],
            "Initial/final snapshot missing")
    require(math.isclose(records[-1][1], float(cfg["t_end"]), rel_tol=1e-9), "Run did not reach t_end")
    require(all(b[1] > a[1] and b[2] > 0 for a, b in zip(records, records[1:])), "Invalid time sequence")

    header, initial = read_snapshot(problem_dir/records[0][3], cfg)
    require(header["step"] == 0 and header["time"] == 0, "Initial snapshot is stale")
    rho0, vx0, vy0, p0 = initial
    gas = rho0 > 0
    require(gas.any() and (~gas).any(), "Initial state must contain gas and true vacuum")
    length = np.array([float(cfg[f"{axis}_upper"])-float(cfg[f"{axis}_lower"]) for axis in ("x", "y")])
    spacing = length/np.array([int(cfg["nxa"]), int(cfg["nya"])])
    gamma = float(cfg["gamma"])
    conserved0 = totals(initial, gamma, np.prod(spacing))
    velocity = conserved0[1:3]/conserved0[0]
    require(np.allclose(vx0[gas], velocity[0], atol=1e-12, rtol=0)
            and np.allclose(vy0[gas], velocity[1], atol=1e-12, rtol=0), "Initial gas velocity is not uniform")
    expansion = 2*np.sqrt(gamma*np.max(p0[gas]/rho0[gas]))/(gamma-1)*float(cfg["t_end"])
    require(expansion < .25*spacing.min(), "Gas is too warm for this cold-translation check")
    mode0 = modes(rho0)
    require(np.all(np.abs(mode0) > 1e-6*rho0.sum()), "Initial image has no usable first Fourier mode")
    scalar = rho0.copy()
    metrics, frames = [], []
    for step, time, dt, filename in records:
        if step:
            stage = upwind_step(scalar, dt, velocity, spacing)
            scalar = (stage if method == "euler" else
                      .5*scalar + .5*upwind_step(stage, dt, velocity, spacing))
        if filename is None:
            continue
        header, values = read_snapshot(problem_dir/filename, cfg)
        require(header["step"] == step and math.isclose(header["time"], time, rel_tol=1e-9, abs_tol=1e-14)
                and math.isclose(header["dt"], dt, rel_tol=1e-9, abs_tol=1e-14), "Snapshot/log mismatch")
        time = header["time"]
        rho = values[0]
        conserved = totals(values, gamma, np.prod(spacing))
        require(np.all(np.abs(conserved-conserved0) < 5e-10*np.maximum(1., np.abs(conserved0))),
                f"Mass, momentum or total energy is not conserved at step {step}")
        coefficient = modes(rho)
        require(np.all(abs(coefficient) > 1e-6*abs(mode0)), "Transport mode was lost")
        expected_phase = np.exp(-2j*np.pi*velocity*time/length)
        phase_error = -np.angle((coefficient/mode0)/expected_phase)*length/(2*np.pi)
        require(np.max(abs(phase_error)/spacing) < 1., f"Translation position error exceeds one cell at step {step}")
        rigid = translate(rho0, *(velocity*time/spacing))
        rigid_l1 = np.sum(abs(rho-rigid))/rho0.sum()
        upwind_l1 = np.sum(abs(rho-scalar))/rho0.sum()
        require(upwind_l1 < .01, f"Cold gas disagrees with independent upwinding at step {step}: {upwind_l1}")
        metrics.append({"step": step, "time": time, "shift_x": float((velocity[0]*time) % length[0]),
                        "shift_y": float((velocity[1]*time) % length[1]),
                        "position_error_x_cells": float(phase_error[0]/spacing[0]),
                        "position_error_y_cells": float(phase_error[1]/spacing[1]),
                        "mass": float(conserved[0]), "momentum_x": float(conserved[1]),
                        "momentum_y": float(conserved[2]), "energy": float(conserved[3]),
                        "rigid_relative_l1": float(rigid_l1), "upwind_relative_l1": float(upwind_l1)})
        frames.append({"step": step, "time": time, "rho": rho.copy(),
                       "rigid": rigid, "upwind": scalar.copy()})
    output = problem_dir/cfg["output_directory"]/"diagnostics"
    output.mkdir(parents=True, exist_ok=True)
    with (output/"metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    plot_results(frames, rho0, cfg, output, animate, fps)
    final = metrics[-1]
    print(f"PASS: {rho0.shape[1]}x{rho0.shape[0]}, {final['step']} steps, t={final['time']:.8g}, "
          f"gas velocity=({velocity[0]:.6g},{velocity[1]:.6g})")
    print(f"Final wrapped displacement: ({final['shift_x']:.6g},{final['shift_y']:.6g})")
    print(f"Position errors in cells: x={final['position_error_x_cells']:.6g}, y={final['position_error_y_cells']:.6g}")
    print(f"Relative density L1: rigid={final['rigid_relative_l1']:.6g} (shape diffusion), "
          f"independent upwind={final['upwind_relative_l1']:.6g}")
    print(f"Initial vacuum-front expansion scale over this run: {expansion:.6g} "
          f"({expansion/spacing.min():.4g} cells); rigid translation is an approximation.")
    print(f"Saved {output/'comparison.png'} and {output/'metrics.csv'}")


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--parameters", type=Path, default=root/"parameters.txt")
    parser.add_argument("--executable", type=Path, default=root/"flare")
    parser.add_argument("--run", action="store_true", help="Run FLARE first; overwrite this problem's run.log")
    parser.add_argument("--animate", action="store_true", help="Also save the density/reference animation")
    parser.add_argument("--fps", type=float, default=6.)
    args = parser.parse_args()
    if not math.isfinite(args.fps) or not 0 < args.fps <= 50:
        parser.error("Require 0 < fps <= 50")
    try:
        path = args.parameters.resolve()
        cfg = parameters(path)
        log = path.parent/"run.log"
        if args.run:
            print("Running FLARE; progress is recorded in " + str(log), flush=True)
            with log.open("w") as stream:
                result = subprocess.run([str(args.executable.resolve()), "-par", str(path)],
                                        cwd=path.parent, stdout=stream, stderr=subprocess.STDOUT)
            require(result.returncode == 0, f"Simulation failed; see {log}")
        analyze(cfg, path.parent, log, args.animate, args.fps)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Advection check failed: {error}\n")


if __name__ == "__main__":
    main()
