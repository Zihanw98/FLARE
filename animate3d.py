#!/usr/bin/env python3
"""Independent 3D FLARE animator: three orthogonal colored planes + 2D cuts.

python3 animate3d.py tests/Spherical_Blast_3D/output mesh_blast3d --field p
python3 animate3d.py DIRECTORY PREFIX --field rho --stride 2 --output movie.mp4

Requires numpy, matplotlib, Pillow; MP4 additionally requires ffmpeg.
No dependency on animate.py. Ghosts are excluded; arbitrary CSV row order is
accepted. Colors stay fixed across selected frames. Slices use the nearest
cell center (the lower center on ties); no interpolation of the solution.
The 3D panel displays cuts, not a volume rendering or an isosurface.
Playback has equal frame durations, not equal physical-time intervals.
If DIRECTORY/run.log exists, only its CSV entries are used. --log overrides
that manifest; otherwise numbered files are discovered. New files are not
added during playback. Use a separate directory per simulation.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


@dataclass
class Snapshot:
    step: int
    time: float
    axes: tuple
    fields: dict  # arrays indexed [x, y, z]


def read_snapshot(path):
    stamp = (path.stat().st_size, path.stat().st_mtime_ns)
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                names = line.strip().split(",")
                break
            key, value = line[1:].strip().split("=", 1)
            metadata[key] = value
        else:
            raise ValueError(f"Empty snapshot: {path}")
        required = ("x", "y", "z", "rho", "p", "vx", "vy", "vz", "gamma", "is_ghost")
        if len(names) != len(set(names)) or not set(required).issubset(names):
            raise ValueError(f"Missing or duplicate CSV columns: {path}")
        data = np.loadtxt(stream, delimiter=",", ndmin=2,
                          usecols=[names.index(key) for key in required])
    if stamp != (path.stat().st_size, path.stat().st_mtime_ns):
        raise ValueError(f"File changed while reading: {path}")
    if not np.isfinite(data).all() or not np.isin(data[:, -1], [0, 1]).all():
        raise ValueError(f"Invalid values or ghost flags: {path}")
    data = data[data[:, -1] == 0]
    if not data.size:
        raise ValueError(f"No active cells: {path}")
    unique = [np.unique(data[:, d], return_inverse=True) for d in range(3)]
    axes = tuple(item[0] for item in unique)
    shape = tuple(len(axis) for axis in axes)
    if min(shape) < 2 or np.prod(shape) != len(data):
        raise ValueError(f"Expected a complete 3D Cartesian grid: {path}")
    indices = tuple(item[1] for item in unique)
    if len(np.unique(np.ravel_multi_index(indices, shape))) != len(data):
        raise ValueError(f"Duplicate active cells: {path}")
    fields = {}
    for col, key in enumerate(required[3:-1], start=3):
        values = np.empty(shape)
        values[indices] = data[:, col]
        fields[key] = values
    step, time = int(metadata["step"]), float(metadata["time"])
    if step < 0 or not np.isfinite(time) or time < 0:
        raise ValueError(f"Invalid step/time metadata: {path}")
    return Snapshot(step, time, axes, fields)


def snapshot_files(directory, prefix, log=None):
    pattern = re.compile(re.escape(prefix.rstrip("_")) + r"_(\d+)\.csv$")
    manifest = Path(log) if log else directory / "run.log"
    if manifest.is_file():
        paths = []
        for line in manifest.read_text().splitlines():
            match = re.search(r", csv = (.+)$", line)
            if match:
                # Output and its manifest may be moved together.
                path = directory / Path(match[1].strip()).name
                if pattern.fullmatch(path.name):
                    paths.append(path)
    elif log:
        raise ValueError(f"Missing manifest: {manifest}")
    else:
        paths = list(directory.glob(f"{prefix.rstrip('_')}_*.csv"))
    pairs = [(int(match[1]), path) for path in paths
             if (match := pattern.fullmatch(path.name))]
    pairs.sort()
    if not pairs or len({step for step, _ in pairs}) != len(pairs):
        raise ValueError("No matching snapshots or duplicate steps in manifest")
    # Check ALL timestamps, including frames later skipped by --stride.
    previous = -np.inf
    for step, path in pairs:
        with path.open() as stream:
            metadata = {}
            for line in stream:
                if not line.startswith("#"):
                    break
                key, value = line[1:].strip().split("=", 1)
                metadata[key] = value
        time = float(metadata["time"])
        if int(metadata["step"]) != step or not np.isfinite(time) or time < 0 or time <= previous:
            raise ValueError(f"Inconsistent step/time; possibly mixed runs: {path}")
        previous = time
    return [path for _, path in pairs]


def edges(axis):
    return np.r_[axis[0] - (axis[1] - axis[0]) / 2,
                 (axis[:-1] + axis[1:]) / 2,
                 axis[-1] + (axis[-1] - axis[-2]) / 2]


def field_values(snapshot, field):
    if field == "speed":
        return np.sqrt(sum(snapshot.fields[key] ** 2 for key in ("vx", "vy", "vz")))
    return snapshot.fields[field]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", type=Path)
    parser.add_argument("prefix")
    parser.add_argument("--field", choices=("rho", "p", "vx", "vy", "vz", "speed"), default="rho")
    parser.add_argument("--stride", type=int, default=1, help="Sample saved snapshots; always retain first/last")
    parser.add_argument("--fps", type=int, default=10)
    parser.add_argument("--dpi", type=int, default=100)
    parser.add_argument("--output", type=Path, help="GIF or MP4; also saves *_last.png")
    parser.add_argument("--log", type=Path, help="Simulation log listing this run's CSVs")
    parser.add_argument("--slice-x", type=float)
    parser.add_argument("--slice-y", type=float)
    parser.add_argument("--slice-z", type=float)
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--vmax", type=float)
    args = parser.parse_args()
    if min(args.stride, args.fps, args.dpi) < 1:
        parser.error("stride, fps and dpi must be positive")
    output = args.output or args.directory / "animations" / f"{args.prefix}_{args.field}_3d.gif"
    if output.suffix.lower() not in (".gif", ".mp4"):
        parser.error("output must end in .gif or .mp4")
    paths = snapshot_files(args.directory, args.prefix, args.log)
    selected = paths[::args.stride]
    if selected[-1] != paths[-1]:
        selected.append(paths[-1])
    frames = []
    axes = None
    for path in selected:
        snap = read_snapshot(path)
        if axes is None:
            axes = snap.axes
        elif any(a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=1e-13)
                 for a, b in zip(axes, snap.axes)):
            raise ValueError(f"Grid changed between frames: {path}")
        frames.append((snap.step, snap.time, field_values(snap, args.field)))
    bounds = tuple(edges(axis) for axis in axes)
    cuts = []
    for axis, edge, position in zip(axes, bounds, (args.slice_x, args.slice_y, args.slice_z)):
        position = (edge[0] + edge[-1]) / 2 if position is None else position
        if not np.isfinite(position) or not edge[0] <= position <= edge[-1]:
            parser.error("slice positions must be finite and inside the active domain")
        cuts.append(int(np.argmin(np.abs(axis - position))))
    low = min(values.min() for _, _, values in frames) if args.vmin is None else args.vmin
    high = max(values.max() for _, _, values in frames) if args.vmax is None else args.vmax
    if args.field in ("vx", "vy", "vz") and args.vmin is None and args.vmax is None:
        high = max(abs(low), abs(high))
        low = -high
    if not np.isfinite([low, high]).all() or low > high:
        parser.error("invalid color limits")
    if low == high:
        low, high = low - max(1, abs(low)) * 1e-6, high + max(1, abs(high)) * 1e-6

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.animation import FFMpegWriter, PillowWriter, writers
    from matplotlib.colors import Normalize
    if output.suffix.lower() == ".mp4" and not writers.is_available("ffmpeg"):
        parser.error("MP4 requires ffmpeg; use --output movie.gif instead")
    cmap = plt.get_cmap("RdBu_r" if args.field in ("vx", "vy", "vz") else "inferno")
    norm = Normalize(low, high)
    fig = plt.figure(figsize=(11, 7), layout="constrained")
    grid = fig.add_gridspec(3, 2, width_ratios=(2.1, 1))
    spatial = fig.add_subplot(grid[:, 0], projection="3d")
    panels = [fig.add_subplot(grid[k, 1]) for k in range(3)]
    labels = ("x", "y", "z")
    # Each panel uses (horizontal dimension, vertical dimension, fixed dimension).
    planes = ((0, 1, 2), (0, 2, 1), (1, 2, 0))
    meshes, geometry = [], []
    for ax, (a, b, fixed) in zip(panels, planes):
        values = np.take(frames[0][2], cuts[fixed], axis=fixed)
        meshes.append(ax.pcolormesh(bounds[a], bounds[b], values.T, cmap=cmap, norm=norm))
        ax.set(xlabel=labels[a], ylabel=labels[b],
               title=f"{labels[fixed]} = {axes[fixed][cuts[fixed]]:.5g}", aspect="equal")
        first, second = np.meshgrid(bounds[a], bounds[b], indexing="ij")
        xyz = [None] * 3
        xyz[a], xyz[b] = first, second
        xyz[fixed] = np.full_like(first, axes[fixed][cuts[fixed]])
        geometry.append(xyz)
    spatial.set(xlabel="x", ylabel="y", zlabel="z", title="3D orthogonal cuts (not volume rendering)",
                xlim=(bounds[0][0], bounds[0][-1]), ylim=(bounds[1][0], bounds[1][-1]),
                zlim=(bounds[2][0], bounds[2][-1]))
    spatial.set_box_aspect([edge[-1] - edge[0] for edge in bounds])
    spatial.view_init(elev=25, azim=-55)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=[spatial, *panels],
                 shrink=.8, pad=.02, label=args.field)
    title = fig.suptitle("")
    surfaces = []

    def update(frame):
        step, time, values = frame
        for surface in surfaces:
            surface.remove()
        surfaces.clear()
        for image, xyz, (_, _, fixed) in zip(meshes, geometry, planes):
            section = np.take(values, cuts[fixed], axis=fixed)
            image.set_array(section.T.ravel())
            surfaces.append(spatial.plot_surface(*xyz, facecolors=cmap(norm(section)),
                            rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False))
        title.set_text(f"FLARE 3D | {args.field} | step {step} | t = {time:.6f}")

    output.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=args.fps) if output.suffix.lower() == ".gif" else FFMpegWriter(fps=args.fps)
    print(f"Rendering {len(frames)}/{len(paths)} frames to {output}", flush=True)
    with writer.saving(fig, str(output), dpi=args.dpi):
        for i, frame in enumerate(frames):
            update(frame)
            writer.grab_frame()
            if i % 10 == 0:
                print(f"  frame {i + 1}/{len(frames)}", flush=True)
    last = output.with_name(output.stem + "_last.png")
    fig.savefig(last, dpi=args.dpi)
    plt.close(fig)
    print(f"Saved {output}\nSaved {last}")


if __name__ == "__main__":
    main()
