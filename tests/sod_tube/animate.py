#!/usr/bin/env python3
"""Read FLARE's 1D Sod CSV snapshots and save colormap/line-plot GIFs.

Ghost cells are excluded. Times come from the CSV comment header, with run.log as a
fallback for older CSVs; time is never estimated from step * dt.
Each saved snapshot is one animation frame; playback is by snapshot, not
uniform physical time. The colormap repeats the 1D data vertically.
"""

import argparse
from dataclasses import dataclass
import fnmatch
import math
from pathlib import Path
import re
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # Also works on servers without a display.
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
from snapshot_io import read_snapshot_csv


FIELDS = (("rho", "Density", "viridis", "#176b87"),
          ("vx", "Velocity", "coolwarm", "#b55132"),
          ("p", "Pressure", "magma", "#7351a1"))


@dataclass
class Snapshot:
    step: int
    time: Optional[float]  # None when neither CSV metadata nor a log is available.
    values: np.ndarray  # shape: (3 physical quantities, active cell count)
    dt: Optional[float] = None


def load_snapshots(input_dir, pattern="mesh_sod_*.csv", log_path=None):
    """Return cell centers and ordered snapshots, selecting this run's files.

    If log_path is supplied, use only snapshots listed in that log. This avoids
    mixing in stale CSVs from previous runs with a different output cadence.
    """
    files = {}
    for path in input_dir.glob(pattern):
        match = re.search(r"_(\d+)\.csv$", path.name)
        if match:
            step = int(match[1])
            if step in files:
                raise ValueError(f"More than one CSV for step {step}; narrow --pattern")
            files[step] = path
    if not files:
        raise ValueError(f"No snapshots matching {pattern} in {input_dir}; run make test first")

    times = {}
    if log_path is not None:
        for line in log_path.read_text().splitlines():
            match = re.match(r"step = (\d+), time = ([^,]+),.*?, csv = (.+)$", line)
            if not match or not fnmatch.fnmatch(Path(match[3]).name, pattern):
                continue
            step, time = int(match[1]), float(match[2])
            if step in times:
                raise ValueError("The log contains repeated steps; use a log from a single run")
            if step not in files or files[step].name != Path(match[3]).name:
                raise ValueError(f"Snapshot listed in the log is missing: {match[3]}")
            if not math.isfinite(time) or time < 0:
                raise ValueError(f"Invalid time for step {step}")
            times[step] = time
        if not times:
            raise ValueError(f"No matching snapshot records in {log_path}; check --log or use --no-log")
        ignored = len(files) - len(times)
        if ignored:
            print(f"Ignoring {ignored} CSV(s) absent from this run's log")
        files = {step: files[step] for step in times}

    x = None
    snapshots = []
    required = {"x", "y", "z", "rho", "vx", "p", "is_ghost"}
    for step, path in sorted(files.items()):
        header, fields, raw_rows = read_snapshot_csv(path)
        if not required.issubset(fields):
            raise ValueError(f"Missing required CSV columns in {path}")
        # Also accept the briefly used per-row metadata format for existing files.
        metadata_columns = {"step", "time", "dt"}
        present = metadata_columns.intersection(fields)
        if present and present != metadata_columns:
            raise ValueError(f"Incomplete step/time/dt metadata in {path}")
        columns = required | present
        rows = [{key: int(row[key]) if key == "step" else float(row[key])
                 for key in columns} for row in raw_rows]
        if not rows or any(not math.isfinite(value) for row in rows for value in row.values()):
            raise ValueError(f"Empty or non-finite CSV: {path}")
        time, dt = times.get(step), None
        if present:
            metadata = tuple(rows[0][key] for key in ("step", "time", "dt"))
            if metadata[0] != step or metadata[1] < 0 or metadata[2] < 0:
                raise ValueError(f"Invalid snapshot step/time/dt in {path}")
            if any(tuple(row[key] for key in ("step", "time", "dt")) != metadata for row in rows):
                raise ValueError(f"Snapshot metadata differs between cells in {path}")
            # Full-precision CSV metadata takes precedence over rounded log times.
            _, time, dt = metadata
        if header:
            if header["step"] != step:
                raise ValueError(f"Header step disagrees with the filename: {path}")
            if present and (header["time"], header["dt"]) != (time, dt):
                raise ValueError(f"Header and row metadata disagree: {path}")
            time, dt = header["time"], header["dt"]
        if any(row["is_ghost"] not in (0, 1) for row in rows):
            raise ValueError(f"Invalid ghost flags in {path}")
        rows = sorted((row for row in rows if row["is_ghost"] == 0), key=lambda row: row["x"])
        if len(rows) < 2 or any(row["y"] != 0 or row["z"] != 0 for row in rows):
            raise ValueError(f"Expected a 1D mesh with at least two active cells: {path}")
        current_x = np.array([row["x"] for row in rows])
        spacing = np.diff(current_x)
        if np.any(spacing <= 0) or not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-14):
            raise ValueError(f"This Sod visualizer expects a uniform 1D mesh: {path}")
        if x is None:
            x = current_x
        elif x.shape != current_x.shape or not np.allclose(x, current_x, rtol=0, atol=1e-13):
            raise ValueError(f"The grid changed between snapshots: {path}")
        values = np.array([[row[key] for row in rows] for key, *_ in FIELDS])
        if np.any(values[[0, 2]] <= 0):
            raise ValueError(f"Nonpositive density or pressure in {path}")
        snapshots.append(Snapshot(step, time, values, dt))
    if len(snapshots) < 2:
        raise ValueError("At least two snapshots are needed for an evolution animation")
    timed = [snapshot for snapshot in snapshots if snapshot.time is not None]
    if any(b.time <= a.time for a, b in zip(timed, timed[1:])):
        raise ValueError("Snapshot times must increase strictly")
    return x, snapshots


def save_animation(x, snapshots, output_dir, kind, fps, dpi):
    """Save one three-panel GIF and its final-frame PNG with fixed scales."""
    values = np.stack([snapshot.values for snapshot in snapshots])
    lower = values.min(axis=(0, 2))
    upper = values.max(axis=(0, 2))
    # Include zero and allow room around the curves; never rescale per frame.
    curve_lower = np.minimum(lower, 0.0)
    spans = np.maximum(upper - curve_lower, 1e-8)
    fig, axes = plt.subplots(3, 1, figsize=(9, 7.2), sharex=True)
    fig.subplots_adjust(left=0.12, right=0.88 if kind == "colormap" else 0.96,
                        bottom=0.12, top=0.88, hspace=0.35)
    title = fig.suptitle("", fontsize=15, fontweight="bold", y=0.97)
    subtitle = "1D data repeated vertically; vertical direction is illustrative" if kind == "colormap" else "Active cells only; fixed axes throughout the evolution"
    fig.text(0.5, 0.92, subtitle, ha="center", color="#555555", fontsize=10)
    fig.text(0.5, 0.025, "One saved snapshot per frame; playback intervals are not physical time intervals",
             ha="center", color="#555555", fontsize=9)
    dx = x[1] - x[0]
    bounds = (x[0] - dx / 2, x[-1] + dx / 2)
    artists = []
    for k, (axis, (_, label, cmap, color)) in enumerate(zip(axes, FIELDS)):
        axis.set_xlim(*bounds)
        axis.set_ylabel(label)
        if kind == "colormap":
            image = axis.imshow(np.tile(snapshots[0].values[k], (24, 1)),
                                origin="lower", aspect="auto", interpolation="nearest",
                                extent=(*bounds, 0, 1), cmap=cmap,
                                vmin=lower[k], vmax=max(upper[k], lower[k] + 1e-8))
            axis.set_yticks([])
            fig.colorbar(image, ax=axis, pad=0.025, fraction=0.035)
            artists.append(image)
        else:
            line, = axis.plot(x, snapshots[0].values[k], color=color, linewidth=2)
            axis.set_ylim(curve_lower[k] - 0.06 * spans[k], upper[k] + 0.08 * spans[k])
            axis.grid(alpha=0.2)
            axis.spines[["top", "right"]].set_visible(False)
            artists.append(line)
    axes[-1].set_xlabel("x")

    def update(frame):
        snapshot = snapshots[frame]
        time_label = f"t = {snapshot.time:.6f}" if snapshot.time is not None else "time unavailable"
        title.set_text(f"Sod shock tube  |  step {snapshot.step}  |  {time_label}")
        for k, artist in enumerate(artists):
            if kind == "colormap":
                artist.set_data(np.tile(snapshot.values[k], (24, 1)))
            else:
                artist.set_ydata(snapshot.values[k])
        return [title, *artists]

    movie = FuncAnimation(fig, update, frames=len(snapshots), interval=1000 / fps, blit=False)
    destination = output_dir / f"sod_{kind}.gif"
    movie.save(destination, writer=PillowWriter(fps=fps), dpi=dpi)
    update(len(snapshots) - 1)
    fig.savefig(output_dir / f"sod_{kind}_last.png", dpi=dpi)
    plt.close(fig)
    print(f"Saved {destination} ({len(snapshots)} frames, {fps:g} fps)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parent / "output")
    parser.add_argument("--output", type=Path, help="Default: INPUT/animations")
    parser.add_argument("--pattern", default="mesh_sod_*.csv")
    logs = parser.add_mutually_exclusive_group()
    logs.add_argument("--log", type=Path, help="Restrict files to this run; provide times for legacy CSVs")
    logs.add_argument("--no-log", action="store_true", help="Read CSVs only; skip automatic legacy log fallback")
    parser.add_argument("--kind", choices=("both", "colormap", "plot"), default="both")
    parser.add_argument("--fps", type=float, default=4.0)
    parser.add_argument("--stride", type=int, default=1, help="Keep every Nth snapshot; always retain the last")
    parser.add_argument("--dpi", type=int, default=110)
    args = parser.parse_args()
    if not math.isfinite(args.fps) or not 0 < args.fps <= 50 or args.stride < 1 or args.dpi < 1:
        parser.error("Require 0 < fps <= 50, stride >= 1, dpi >= 1")
    input_dir = args.input.resolve()
    output_dir = args.output.resolve() if args.output else input_dir / "animations"
    default_log = input_dir.parent / "run.log"
    log_path = args.log
    try:
        x, snapshots = load_snapshots(input_dir, args.pattern, log_path)
        if not args.no_log and log_path is None and default_log.exists() \
                and any(snapshot.time is None for snapshot in snapshots):
            x, snapshots = load_snapshots(input_dir, args.pattern, default_log)
        if any(snapshot.time is None for snapshot in snapshots):
            print("Some snapshots have no CSV time or run log: showing their step numbers only")
        selected = snapshots[::args.stride]
        if selected[-1].step != snapshots[-1].step:
            selected.append(snapshots[-1])
        output_dir.mkdir(parents=True, exist_ok=True)
        for kind in (("colormap", "plot") if args.kind == "both" else (args.kind,)):
            save_animation(x, selected, output_dir, kind, args.fps, args.dpi)
    except (OSError, ValueError) as error:
        parser.exit(1, f"Animation failed: {error}\n")


if __name__ == "__main__":
    main()
