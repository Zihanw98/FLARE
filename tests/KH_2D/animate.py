#!/usr/bin/env python3
"""Plot FLARE's 2D KH snapshots as a four-panel GIF and a final-frame PNG.

Run: python3 animate.py
     python3 animate.py --fps 8 --stride 2

Requires numpy, matplotlib and Pillow. Reads existing CSVs only; it can run while
FLARE is running. Run it again to include newly written snapshots. Ghost cells
are excluded. Each frame is one selected snapshot, not an equal physical-time
interval. Color limits are fixed across the selected frames.
"""

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")  # Save figures on machines without a graphical display.
from matplotlib import pyplot as plt
from matplotlib.animation import PillowWriter
from matplotlib.ticker import ScalarFormatter
import numpy as np


FIELDS = (("rho", "Density", "viridis"), ("p", "Pressure", "magma"),
          ("vx", "x velocity", "RdBu_r"), ("vy", "y velocity", "RdBu_r"))
COLUMNS = ("x", "y", "z", *(field[0] for field in FIELDS), "is_ghost")


@dataclass
class Snapshot:
    step: int
    time: float
    dt: float
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray  # (quantity, y, x); never infer orientation from row order.
    total_cells: int


def signature(path):
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def read_snapshot(path):
    """Read and validate one complete, uniform, rectangular 2D mesh."""
    before = signature(path)
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if not line.strip():
                continue
            if not line.startswith("#"):
                columns = line.strip().split(",")
                break
            key, separator, value = line[1:].strip().partition("=")
            key = key.strip()
            if key in ("step", "time", "dt"):
                if not separator or key in metadata:
                    raise ValueError("Invalid or duplicate snapshot metadata")
                metadata[key] = int(value) if key == "step" else float(value)
        else:
            raise ValueError("Missing CSV table")
        if set(metadata) != {"step", "time", "dt"}:
            raise ValueError("Missing step/time/dt comment header")
        if any(not math.isfinite(v) or v < 0 for v in metadata.values()):
            raise ValueError("Invalid step/time/dt")
        if not set(COLUMNS).issubset(columns):
            raise ValueError("Missing coordinates, primitive quantities or is_ghost")
        data = np.loadtxt(stream, delimiter=",", ndmin=2,
                          usecols=[columns.index(key) for key in COLUMNS])
    if before != signature(path):
        raise ValueError("File changed while being read; still being written")
    if data.size == 0 or not np.isfinite(data).all():
        raise ValueError("Empty table or non-finite values")
    if not np.all(data[:, 2] == 0):
        raise ValueError("Expected a 2D mesh with z = 0")
    if not np.isin(data[:, -1], (0, 1)).all():
        raise ValueError("Invalid ghost flags")

    # Check the full grid too: a truncated CSV may already contain active cells
    # but still be missing its final ghost rows.
    all_x, ix = np.unique(data[:, 0], return_inverse=True)
    all_y, iy = np.unique(data[:, 1], return_inverse=True)
    if len(data) != len(all_x) * len(all_y) or \
            len(np.unique(iy * len(all_x) + ix)) != len(data):
        raise ValueError("Incomplete full grid or duplicate coordinates")
    active = data[data[:, -1] == 0]
    x, ix = np.unique(active[:, 0], return_inverse=True)
    y, iy = np.unique(active[:, 1], return_inverse=True)
    if len(x) < 2 or len(y) < 2 or len(active) != len(x) * len(y):
        raise ValueError("Active cells must form a complete 2D rectangle")
    for axis in (x, y):
        spacing = np.diff(axis)
        if not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-13):
            raise ValueError("This visualizer requires a uniform grid")
    if np.any(active[:, 3:5] <= 0):
        raise ValueError("Nonpositive density or pressure")

    values = np.empty((len(FIELDS), len(y), len(x)))
    values[:, iy, ix] = active[:, 3:7].T
    return Snapshot(**metadata, x=x, y=y, values=values, total_cells=len(data))


def check_grid(snapshot, reference):
    if snapshot.total_cells != reference.total_cells or any(
        a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=1e-13)
        for a, b in ((snapshot.x, reference.x), (snapshot.y, reference.y))
    ):
        raise ValueError("Grid changed or snapshot is incomplete")


def find_snapshots(directory, pattern):
    files = {}
    for path in directory.glob(pattern):
        match = re.search(r"_(\d+)\.csv$", path.name)
        if not match:
            continue
        step = int(match[1])
        if step in files:
            raise ValueError(f"Multiple files for step {step}; narrow --pattern")
        files[step] = path
    if not files:
        raise ValueError(f"No snapshots matching {pattern} in {directory}")
    return sorted(files.items())


def prepare_frames(files, stride, max_frames):
    """Validate the latest file before sampling; omit it if still incomplete."""
    reference = read_snapshot(files[0][1])
    if reference.step != files[0][0]:
        raise ValueError("Initial CSV header step disagrees with its filename")
    if len(files) > 1:
        try:
            latest = read_snapshot(files[-1][1])
            check_grid(latest, reference)
            if latest.step != files[-1][0]:
                raise ValueError("CSV header step disagrees with its filename")
        except (OSError, ValueError) as error:
            print(f"Skipping latest snapshot {files[-1][1].name}: {error}", flush=True)
            files = files[:-1]

    selected = files[::stride]
    if selected[-1] != files[-1]:
        selected.append(files[-1])
    if len(selected) > max_frames:
        indices = np.linspace(0, len(selected) - 1, max_frames, dtype=int)
        selected = [selected[i] for i in indices]

    lower = np.full(len(FIELDS), np.inf)
    upper = np.full(len(FIELDS), -np.inf)
    records = []
    previous_time = -1.0
    for index, (step, path) in enumerate(selected):
        before = signature(path)
        snapshot = read_snapshot(path)
        check_grid(snapshot, reference)
        if before != signature(path):
            raise ValueError(f"Snapshot changed: {path}")
        if snapshot.step != step or snapshot.time <= previous_time:
            raise ValueError(f"Inconsistent steps/times: {path}; check for old run files")
        previous_time = snapshot.time
        lower = np.minimum(lower, snapshot.values.min(axis=(1, 2)))
        upper = np.maximum(upper, snapshot.values.max(axis=(1, 2)))
        records.append((path, before, step, snapshot.time, snapshot.dt))
        if index % 10 == 0 or index == len(selected) - 1:
            print(f"Reading {index + 1}/{len(selected)}: step {step}, "
                  f"t={snapshot.time:.6g}", flush=True)
    return reference, records, lower, upper


def save_animation(reference, records, lower, upper, output_dir, fps, dpi):
    """Read one mesh at a time, keeping a fixed color scale for each quantity."""
    # Diverging velocity colors are centered at zero throughout the movie.
    for k, (name, _, _) in enumerate(FIELDS):
        if name in ("vx", "vy"):
            span = max(abs(lower[k]), abs(upper[k]), 1e-12)
            lower[k], upper[k] = -span, span
        elif upper[k] <= lower[k]:
            pad = max(abs(lower[k]) * 1e-6, 1e-12)
            lower[k], upper[k] = lower[k] - pad, upper[k] + pad

    dx, dy = reference.x[1] - reference.x[0], reference.y[1] - reference.y[0]
    extent = (reference.x[0] - dx / 2, reference.x[-1] + dx / 2,
              reference.y[0] - dy / 2, reference.y[-1] + dy / 2)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.3))
    fig.subplots_adjust(left=0.07, right=0.90, bottom=0.09, top=0.85,
                        wspace=0.36, hspace=0.36)
    title = fig.suptitle("", fontsize=14, y=0.98)
    fig.text(0.5, 0.90, f"{len(reference.x)} x {len(reference.y)} active cells; "
             "fixed color scales; ghosts excluded", ha="center", fontsize=10)
    fig.text(0.5, 0.025, "One selected snapshot per frame; playback is not uniform physical time",
             ha="center", fontsize=9, color="#555555")
    images = []
    for k, (axis, (name, label, cmap)) in enumerate(zip(axes.flat, FIELDS)):
        artist = axis.imshow(reference.values[k], origin="lower", extent=extent,
                             interpolation="nearest", aspect="equal", cmap=cmap,
                             vmin=lower[k], vmax=upper[k])
        axis.set(title=f"{label} ({name})", xlabel="x", ylabel="y")
        formatter = ScalarFormatter(useOffset=False)
        fig.colorbar(artist, ax=axis, pad=0.025, fraction=0.046, format=formatter)
        images.append(artist)

    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / "kh_colormap.gif"
    # Pillow keeps rendered frames in memory; --max-frames bounds their number.
    writer = PillowWriter(fps=fps)
    try:
        with writer.saving(fig, str(destination), dpi=dpi):
            for index, (path, recorded_signature, step, time, dt) in enumerate(records):
                if signature(path) != recorded_signature:
                    raise ValueError(f"Snapshot changed since scanning: {path}; run again")
                snapshot = read_snapshot(path)
                check_grid(snapshot, reference)
                if signature(path) != recorded_signature:
                    raise ValueError(f"Snapshot changed since scanning: {path}; run again")
                title.set_text(f"Kelvin-Helmholtz | step {step} | t = {time:.6f}\n"
                               f"dt = {dt:.4g}")
                for k, artist in enumerate(images):
                    artist.set_data(snapshot.values[k])
                writer.grab_frame()
                if index % 10 == 0 or index == len(records) - 1:
                    print(f"Rendering {index + 1}/{len(records)}", flush=True)
        fig.savefig(output_dir / "kh_colormap_last.png", dpi=dpi)
    finally:
        plt.close(fig)
    print(f"Saved {destination} ({len(records)} frames, {fps:g} fps)")
    print(f"Saved {output_dir / 'kh_colormap_last.png'}")
    print(f"Included through step {records[-1][2]}, t={records[-1][3]:.10g}; "
          "run again to include newer CSVs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parent / "output")
    parser.add_argument("--output", type=Path, help="Default: INPUT/animations")
    parser.add_argument("--pattern", default="mesh_*.csv", help="Includes current mesh_sod_*.csv names")
    parser.add_argument("--fps", type=float, default=6.0)
    parser.add_argument("--stride", type=int, default=1, help="Select every Nth snapshot, plus the last")
    parser.add_argument("--max-frames", type=int, default=120,
                        help="Cap GIF memory use by sampling; keep first/last (default: 120)")
    parser.add_argument("--dpi", type=int, default=100)
    args = parser.parse_args()
    if not math.isfinite(args.fps) or not 0 < args.fps <= 50 or \
            args.stride < 1 or args.max_frames < 2 or args.dpi < 1:
        parser.error("Require 0 < fps <= 50, stride >= 1, max-frames >= 2, dpi >= 1")
    try:
        files = find_snapshots(args.input, args.pattern)
        print(f"Found {len(files)} snapshot files at startup", flush=True)
        reference, records, lower, upper = prepare_frames(files, args.stride, args.max_frames)
        output_dir = args.output if args.output else args.input / "animations"
        save_animation(reference, records, lower, upper, output_dir, args.fps, args.dpi)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Animation failed: {error}\n")


if __name__ == "__main__":
    main()
