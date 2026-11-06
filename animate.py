#!/usr/bin/env python3
"""Animate FLARE snapshots: colormaps on the left, profiles on the right.

Examples (run from the project root):
  python3 animate.py tests/sod_tube/output mesh_sod --stride 5
  python3 animate.py tests/KH_2D/output mesh_sod --stride 100 --slice-y 0.5
  python3 animate.py tests/KH_2D/output mesh_sod --stride 100 --slice-x 0.25
  python3 animate.py output test --stride 10
  python3 animate.py tests/sod_tube/output mesh_sod --show
  python3 animate.py tests/KH_2D/output mesh_sod --stride 100 --plotly --slice-x 0.25
  python3 animate.py tests/sod_tube/output mesh_sod --output figures/sod.html
  python3 animate.py tests/Spherical_Sod_1D/output spherical_sod --plot spherical

Requires numpy, matplotlib, Pillow (ffmpeg additionally for MP4).
Interactive HTML requires numpy and plotly: python3 -m pip install plotly
Reads PREFIX_<step>.csv (x, optional y/z, rho, vx, vy, p or pressure,
optional is_ghost) or legacy PREFIX_<step>.txt (rho vx p plus FLARE header).
Prefix may include its trailing underscore. Legacy TXT uses cell indices,
because physical coordinates are not stored in those files.

1D: rho/vx/pressure strips repeated vertically, with synchronized x profiles.
--plot spherical: interpret 1D CSV x as radius and vx as radial velocity;
show spherical symmetry slices and radial profiles. A nonzero inner radius
is shown as a masked cavity. GIF, MP4 and HTML use the same radial mapping.
2D: rho/vx/vy/pressure maps, with x profiles at the nearest --slice-y row
(default: domain midpoint), or y profiles at the nearest --slice-x column.
These options are mutually exclusive. A dashed line marks the cut on each map.
--streamlines overlays instantaneous velocity streamlines on 2D rho/pressure
maps, including HTML output. Nonuniform grids are interpolated for tracing.
Scales are fixed over selected frames; ghost cells are excluded. Stride counts
saved files, not solver steps. First/last frames are retained. Playback uses
equal frame durations, not equal physical-time intervals.
--cgs converts CSV coordinates to cm, time to s, density to g/cm^3,
velocity to cm/s and pressure to dyn/cm^2 using the four base-unit headers.
With --cgs, --slice-x/--slice-y are in cm. CSV table values remain unchanged.
Old snapshots without unit metadata can still be plotted without --cgs.
--log uses logarithmic color scales and profile y axes for rho and pressure;
velocity and spatial axes stay linear. All selected rho/pressure values must
be positive. Units and physical values, including HTML hover values, are unchanged.

Default output: DIRECTORY/animations/PREFIX_evolution.gif and *_last.png.
Use --output movie.mp4 for MP4, or --show to play in a GUI without saving.
Use --plotly or --output viewer.html for an offline interactive HTML viewer.
Its six controls are first/previous/play/pause/next/last. It starts paused and
stops at the last frame. HTML embeds all selected data: increase --stride for
large runs. --fps controls playback speed; --dpi applies only to GIF/MP4/PNG.
Files are discovered once at startup; rerun to include new simulation output.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


FIELDS = (("rho", "rho", "viridis"), ("vx", "vx", "RdBu_r"),
          ("vy", "vy", "RdBu_r"), ("p", "pressure", "magma"))


@dataclass
class Snapshot:
    step: int
    time: object
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray  # (field, y, x), with one y row for 1D
    fields: tuple
    xlabel: str = "x"
    units: tuple = ()  # Base-unit values in g, cm, s, K; None for missing entries.
    cgs: bool = False


def signature(path):
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def find_files(directory, prefix, stride):
    if not directory.is_dir():
        raise ValueError(f"Not a directory: {directory}")
    pattern = re.compile(re.escape(prefix.rstrip("_")) + r"_?(\d+)\.(csv|txt)$")
    files = {}
    for path in directory.iterdir():
        match = pattern.fullmatch(path.name)
        if match and path.is_file():
            step = int(match[1])
            if step in files:
                raise ValueError(f"Multiple snapshots for step {step}; use a narrower prefix")
            files[step] = path
    ordered = sorted(files.items())
    if not ordered:
        raise ValueError(f"No {prefix}_<step>.csv/.txt snapshots in {directory}")
    selected = ordered[::stride]
    if selected[-1] != ordered[-1]:
        selected.append(ordered[-1])
    print(f"Found {len(ordered)} snapshots; selected {len(selected)} frames", flush=True)
    return selected


def read_snapshot(path, step, cgs=False):
    before = signature(path)
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if not line.strip():
                continue
            if not line.lstrip().startswith("#"):
                break
            metadata.update(re.findall(r"(\w+)\s*=\s*([^\s]+)", line))
        else:
            raise ValueError(f"Empty snapshot: {path}")

        if path.suffix == ".txt":
            # The first non-comment line is already numeric in legacy TXT.
            from itertools import chain
            data = np.loadtxt(chain((line,), stream), ndmin=2)
            nx, ng, na = (int(metadata[k]) for k in ("nx", "nxg", "nxa"))
            if ng < 0 or na < 2 or nx != na + 2 * ng or data.shape != (nx, 3):
                raise ValueError(f"Invalid legacy TXT grid: {path}")
            data = data[ng:ng + na]
            columns = {"x": np.arange(na, dtype=float), "rho": data[:, 0],
                       "vx": data[:, 1], "p": data[:, 2]}
        else:
            names = [name.strip() for name in line.strip().split(",")]
            if len(names) != len(set(names)):
                raise ValueError(f"Duplicate CSV columns: {path}")
            data = np.loadtxt(stream, delimiter=",", ndmin=2)
            if not data.size or data.shape[1] != len(names):
                raise ValueError(f"Empty or incomplete CSV table: {path}")
            columns = dict(zip(names, data.T))
            # Accept older files that store metadata in each row.
            for key in ("step", "time", "dt"):
                if key in columns:
                    value = columns[key][0]
                    if not np.all(columns[key] == value):
                        raise ValueError(f"Inconsistent {key} column: {path}")
                    if key in metadata and float(metadata[key]) != value:
                        raise ValueError(f"Conflicting {key} metadata: {path}")
                    metadata[key] = str(value)

    if signature(path) != before:
        raise ValueError(f"Snapshot changed while reading: {path}; rerun after output finishes")
    header_step = float(metadata.get("step", metadata.get("timestep", step)))
    if header_step != step:
        raise ValueError(f"Filename and header step disagree: {path}")
    time = float(metadata["time"]) if "time" in metadata else None
    if time is not None and (not np.isfinite(time) or time < 0):
        raise ValueError(f"Invalid time: {path}")
    if "p" not in columns and "pressure" in columns:
        columns["p"] = columns["pressure"]
    if not {"x", "rho", "vx", "p"}.issubset(columns):
        raise ValueError(f"Missing x/rho/vx/p (or pressure) columns: {path}")
    unit_keys = ("unit_mass", "unit_length", "unit_time", "unit_temperature")
    units = tuple(float(metadata[key]) if key in metadata else None for key in unit_keys)
    if any(value is not None and (not np.isfinite(value) or value <= 0) for value in units):
        raise ValueError(f"Base units must be finite and positive: {path}")
    if cgs:
        if path.suffix == ".txt":
            raise ValueError("--cgs requires CSV physical coordinates, not legacy cell indices")
        missing = [key for key, value in zip(unit_keys, units) if value is None]
        if missing:
            raise ValueError(f"--cgs requires unit metadata ({', '.join(missing)}): {path}; regenerate snapshots with the updated solver")
        mass, length, duration, _temperature = units
        try:
            factors = {"x": length, "y": length, "z": length,
                       "rho": mass / length**3, "p": mass / (length * duration**2),
                       "vx": length / duration, "vy": length / duration, "vz": length / duration}
            if any(not np.isfinite(value) or value <= 0 for value in factors.values()):
                raise ValueError("Invalid derived CGS conversion factor")
            with np.errstate(over="raise", invalid="raise"):
                columns = {key: values * factors.get(key, 1.0) for key, values in columns.items()}
            if time is not None:
                time *= duration
                if not np.isfinite(time):
                    raise ValueError("Non-finite CGS time")
        except (OverflowError, ZeroDivisionError, FloatingPointError) as error:
            raise ValueError(f"CGS conversion overflow/underflow: {path}") from error
    if "is_ghost" in columns:
        if not np.isin(columns["is_ghost"], (0, 1)).all():
            raise ValueError(f"Invalid is_ghost flags: {path}")
        active = columns["is_ghost"] == 0
        columns = {key: values[active] for key, values in columns.items()}
    if not len(columns["x"]):
        raise ValueError(f"No active cells: {path}")
    if "z" in columns and np.unique(columns["z"]).size != 1:
        raise ValueError(f"3D snapshots are not supported: {path}")
    x, ix = np.unique(columns["x"], return_inverse=True)
    y, iy = np.unique(columns.get("y", np.zeros_like(columns["x"])), return_inverse=True)
    fields = tuple(field for field in FIELDS if len(y) > 1 or field[0] != "vy")
    required = ["x", *[field[0] for field in fields]]
    required += [key for key in ("y", "z") if key in columns]
    if any(key not in columns for key in required):
        raise ValueError(f"2D snapshots require vy: {path}")
    if any(not np.isfinite(columns[key]).all() for key in required):
        raise ValueError(f"Non-finite coordinates or field values: {path}")
    if len(x) < 2 or len(ix) != len(x) * len(y) or len(np.unique(iy * len(x) + ix)) != len(ix):
        raise ValueError(f"Expected a complete Cartesian grid without duplicate cells: {path}")
    values = np.empty((len(fields), len(y), len(x)))
    for k, (key, _, _) in enumerate(fields):
        values[k, iy, ix] = columns[key]
    return Snapshot(step, time, x, y, values, fields,
                    "active cell index" if path.suffix == ".txt" else "x", units, cgs)


def check_grid(snapshot, reference):
    if snapshot.units != reference.units or snapshot.cgs != reference.cgs:
        raise ValueError("Units changed between snapshots; use files from a single unit system")
    if snapshot.fields != reference.fields or snapshot.xlabel != reference.xlabel or any(
        a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=1e-13)
        for a, b in ((snapshot.x, reference.x), (snapshot.y, reference.y))
    ):
        raise ValueError("Grid changed between snapshots; use files from a single fixed-grid run")


def edges(centers):
    """Infer cell edges from centers, including nonuniform Cartesian grids."""
    return np.r_[centers[0] - (centers[1] - centers[0]) / 2,
                 (centers[:-1] + centers[1:]) / 2,
                 centers[-1] + (centers[-1] - centers[-2]) / 2]


def spherical_projection(snapshot, pixels=241):
    """Map a 1D radial grid onto a symmetry plane without interpolating states.

    CSV stores centers, so shell edges are inferred as in the Cartesian plots.
    The returned mask excludes both the exterior and any inner cavity.
    """
    if len(snapshot.y) != 1:
        raise ValueError("--plot spherical requires a 1D radial snapshot")
    if snapshot.xlabel != "x":
        raise ValueError("--plot spherical requires CSV physical radii, not legacy cell indices")
    radial_edges = edges(snapshot.x)
    tolerance = 1e-12 * max(abs(radial_edges[-1]), 1.0)
    if snapshot.x[0] < 0 or radial_edges[0] < -tolerance:
        raise ValueError("--plot spherical requires nonnegative active radii and shell edges")
    radial_edges[0] = max(0.0, radial_edges[0])
    axis = np.linspace(-radial_edges[-1], radial_edges[-1], pixels)
    radius = np.hypot(axis[:, None], axis[None, :])
    indices = np.searchsorted(radial_edges, radius, side="right") - 1
    outside = (radius < radial_edges[0]) | (radius > radial_edges[-1])
    return axis, np.clip(indices, 0, len(snapshot.x) - 1), outside


def spherical_values(values, projection):
    """Project one field while keeping radial velocity a scalar radial component."""
    _, indices, outside = projection
    return np.ma.array(values[0, indices], mask=outside)


def display_label(name, label, spherical, cgs=False):
    label = "radial velocity" if spherical and name == "vx" else label
    unit = {"rho": "g/cm^3", "vx": "cm/s", "vy": "cm/s", "p": "dyn/cm^2"}
    return f"{label} [{unit[name]}]" if cgs else label


def coordinate_label(label, cgs):
    return f"{label} [cm]" if cgs else label


def logarithmic_limits(lower, upper):
    """Positive limits shared by both backends, including constant fields."""
    if not (np.isfinite(lower) and np.isfinite(upper) and 0 < lower <= upper):
        raise ValueError("--log requires strictly positive finite rho and pressure; zero/negative values cannot be displayed on a log scale")
    if lower == upper:
        return lower / 1.1, upper * 1.1
    return lower, upper


def prepare(files, cgs=False, log=False):
    reference = None
    records = []
    previous_time = -np.inf
    for i, (step, path) in enumerate(files):
        stamp = signature(path)
        snapshot = read_snapshot(path, step, cgs=cgs)
        if log:
            for k, (name, _, _) in enumerate(snapshot.fields):
                if name in ("rho", "p") and np.any(snapshot.values[k] <= 0):
                    raise ValueError(f"--log requires strictly positive {name}: {path}")
        if signature(path) != stamp:
            raise ValueError(f"Snapshot changed while scanning: {path}")
        if reference is None:
            reference = snapshot
            lower = np.full(len(snapshot.fields), np.inf)
            upper = -lower.copy()
        check_grid(snapshot, reference)
        if snapshot.time is not None:
            if snapshot.time <= previous_time:
                raise ValueError(f"Times do not increase: {path}; check for files from older runs")
            previous_time = snapshot.time
        lower = np.minimum(lower, snapshot.values.min(axis=(1, 2)))
        upper = np.maximum(upper, snapshot.values.max(axis=(1, 2)))
        records.append((step, path, stamp))
        if i % 20 == 0 or i == len(files) - 1:
            print(f"Scanning {i + 1}/{len(files)}", flush=True)
    return reference, records, lower, upper


def select_profile(reference, args):
    """Use the same slice and coordinate labels in both rendering backends."""
    is_2d = len(reference.y) > 1
    if args.streamlines and not is_2d:
        raise ValueError("--streamlines requires a 2D snapshot")
    if not is_2d and (args.slice_x is not None or args.slice_y is not None):
        raise ValueError("--slice-x/--slice-y require a 2D snapshot")
    along_y = args.slice_x is not None
    fixed_axis = reference.x if along_y else reference.y
    target = args.slice_x if along_y else args.slice_y
    if target is None:
        target = (fixed_axis[0] + fixed_axis[-1]) / 2
    if is_2d and not edges(fixed_axis)[0] <= target <= edges(fixed_axis)[-1]:
        raise ValueError(f"--slice-{'x' if along_y else 'y'} lies outside the domain")
    cut_index = int(np.argmin(abs(fixed_axis - target)))
    profile_axis = reference.y if along_y else reference.x
    profile_label = "radius r" if args.plot == "spherical" else ("y" if along_y else reference.xlabel)
    profile_label = coordinate_label(profile_label, reference.cgs)
    cut_label = (f" at {'x' if along_y else 'y'} = {fixed_axis[cut_index]:.6g}"
                 + (" cm" if reference.cgs else "")) if is_2d else ""
    return is_2d, along_y, fixed_axis, cut_index, profile_axis, profile_label, cut_label


def streamline_segments(snapshot):
    """Shared line/arrow geometry for both backends, in physical coordinates.

    Matplotlib traces on a uniform lattice. Resample only the velocity used
    for tracing; the scalar maps and profiles retain their original samples.
    """
    from matplotlib.figure import Figure

    x, y = snapshot.x, snapshot.y
    if len(y) < 2:
        raise ValueError("--streamlines requires a 2D snapshot")
    fields = {name: snapshot.values[k] for k, (name, _, _) in enumerate(snapshot.fields)}
    u, v = fields["vx"], fields["vy"]
    if not np.any(u) and not np.any(v):
        return []
    uniform_x = np.linspace(x[0], x[-1], len(x))
    uniform_y = np.linspace(y[0], y[-1], len(y))

    def resample(values):
        rows = np.array([np.interp(uniform_x, x, row) for row in values])
        return np.array([np.interp(uniform_y, y, col) for col in rows.T]).T

    figure = Figure(figsize=(6, 4))
    try:
        ax = figure.subplots()
        ax.set(xlim=(x[0], x[-1]), ylim=(y[0], y[-1]), aspect="equal")
        ax.apply_aspect()
        flow = ax.streamplot(uniform_x, uniform_y, resample(u), resample(v),
                             density=1, color="white", linewidth=0.7,
                             arrowstyle="->", arrowsize=0.9)
        segments = list(flow.lines.get_segments())
        # streamplot adds individual arrow patches to the axes. Their public
        # paths are in data coordinates, independent of the output backend.
        for arrow in ax.patches:
            segments.extend(arrow.get_path().to_polygons(closed_only=False))
        return segments
    finally:
        figure.clear()


def animate(reference, records, lower, upper, args):
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter, writers

    is_2d, along_y, fixed_axis, cut_index, profile_axis, profile_label, cut_label = select_profile(reference, args)
    spherical = args.plot == "spherical"
    projection = spherical_projection(reference) if spherical else None
    profile_edges = edges(profile_axis)

    def profile_values(snapshot, k):
        return snapshot.values[k, :, cut_index] if along_y else snapshot.values[k, cut_index, :]

    nfields = len(reference.fields)
    fig, axes = plt.subplots(nfields, 2, figsize=(12, 3 * nfields), layout="constrained")
    title = fig.suptitle("")
    maps, lines = [], []
    streamlines = []
    xedges = edges(projection[0]) if spherical else edges(reference.x)
    yedges = xedges if spherical else (edges(reference.y) if is_2d else np.linspace(0, 1, 25))

    def map_values(snapshot, k):
        if spherical:
            return spherical_values(snapshot.values[k], projection)
        return snapshot.values[k] if is_2d else np.repeat(snapshot.values[k], 24, axis=0)

    for k, (name, label, cmap) in enumerate(reference.fields):
        label = display_label(name, label, spherical, reference.cgs)
        lo, hi = lower[k], upper[k]
        logarithmic = args.log and name in ("rho", "p")
        if logarithmic:
            lo, hi = logarithmic_limits(lo, hi)
        elif name in ("vx", "vy"):
            hi = max(abs(lo), abs(hi), 1e-12)
            lo = -hi
        elif lo == hi:
            pad = max(abs(lo) * 0.01, 1e-12)
            lo, hi = lo - pad, hi + pad
        left, right = axes[k]
        color_scale = {"norm": LogNorm(vmin=lo, vmax=hi)} if logarithmic else {"vmin": lo, "vmax": hi}
        mesh = left.pcolormesh(xedges, yedges, map_values(reference, k),
                               shading="flat", cmap=cmap, **color_scale)
        fig.colorbar(mesh, ax=left, pad=0.02)
        left.set(title=label + (" (spherical slice)" if spherical else ""),
                 xlabel=coordinate_label("x" if spherical else reference.xlabel, reference.cgs))
        if spherical:
            outer = edges(reference.x)[-1]
            left.set(ylabel=coordinate_label("z", reference.cgs), aspect="equal", xlim=(-outer, outer), ylim=(-outer, outer))
        elif is_2d:
            left.set(ylabel=coordinate_label("y", reference.cgs), aspect="equal")
            mark_cut = left.axvline if along_y else left.axhline
            mark_cut(fixed_axis[cut_index], color="black", linestyle="--", linewidth=0.8)
        else:
            left.set_yticks([])
            left.set_ylabel("repeated 1D")
        line, = right.plot(profile_axis, profile_values(reference, k), linewidth=1.5)
        pad = max((upper[k] - lower[k]) * 0.05, max(abs(lower[k]), abs(upper[k])) * 0.01, 1e-12)
        if logarithmic:
            right.set_yscale("log")
        right.set(xlim=(profile_edges[0], profile_edges[-1]), ylim=(lo, hi) if logarithmic else (lower[k] - pad, upper[k] + pad),
                  xlabel=profile_label, ylabel=label,
                  title=f"{label}({'r' if spherical else 'y' if along_y else 'x'})" + cut_label)
        right.grid(alpha=0.25)
        maps.append(mesh)
        lines.append(line)
        if args.streamlines and name in ("rho", "p"):
            from matplotlib.collections import LineCollection
            import matplotlib.patheffects as effects
            overlay = LineCollection([], colors="white", linewidths=0.7, zorder=3,
                                     path_effects=[effects.Stroke(linewidth=1.5, foreground="#333333"),
                                                   effects.Normal()])
            left.add_collection(overlay)
            streamlines.append(overlay)

    def update(index):
        step, path, stamp = records[index]
        if signature(path) != stamp:
            raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
        snapshot = read_snapshot(path, step, cgs=reference.cgs)
        check_grid(snapshot, reference)
        if signature(path) != stamp:
            raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
        time_label = f"t = {snapshot.time:.7g}" + (" s" if reference.cgs else "") if snapshot.time is not None else "time unavailable"
        title.set_text(f"{args.prefix} | {'1D spherical' if spherical else '2D' if is_2d else '1D'} | step {step} | {time_label}"
                       f" | frame {index + 1}/{len(records)}")
        for k, (mesh, line) in enumerate(zip(maps, lines)):
            mesh.set_array(map_values(snapshot, k).ravel())
            line.set_ydata(profile_values(snapshot, k))
        if streamlines:
            segments = streamline_segments(snapshot)
            for overlay in streamlines:
                overlay.set_segments(segments)
        return [title, *maps, *lines, *streamlines]

    try:
        if args.show:
            update(0)
            movie = FuncAnimation(fig, update, frames=len(records), interval=1000 / args.fps,
                                  blit=False, cache_frame_data=False)
            plt.show()
        else:
            suffix = ("_spherical" if spherical else "") + ("_cgs" if reference.cgs else "") + ("_log" if args.log else "")
            destination = args.output or args.directory / "animations" / f"{args.prefix.rstrip('_')}{suffix}_evolution.gif"
            if destination.suffix.lower() == ".mp4":
                if not writers.is_available("ffmpeg"):
                    raise ValueError("MP4 requires ffmpeg; use --output movie.gif instead")
                writer = FFMpegWriter(fps=args.fps)
            else:
                writer = PillowWriter(fps=args.fps)
            destination.parent.mkdir(parents=True, exist_ok=True)
            # Reload one snapshot at a time. Pillow still buffers rendered GIF frames;
            # use a larger stride or MP4 for long/high-resolution animations.
            with writer.saving(fig, str(destination), dpi=args.dpi):
                for index in range(len(records)):
                    update(index)
                    writer.grab_frame()
                    if index % 20 == 0 or index == len(records) - 1:
                        print(f"Rendering {index + 1}/{len(records)}", flush=True)
            last = destination.with_name(destination.stem + "_last.png")
            fig.savefig(last, dpi=args.dpi)
            print(f"Saved {destination}\nSaved {last}")
    finally:
        plt.close(fig)


PLOTLY_VIEWER = r"""
<style>
body { margin:0; color:#213547; background:#f5f7fa; font:15px system-ui,sans-serif; }
header { position:sticky; top:0; z-index:10; background:#fffffff5; padding:14px 22px;
         border-bottom:1px solid #d8dfe8; }
h1 { font-size:19px; margin:0 0 10px; }
nav { display:flex; gap:8px; flex-wrap:wrap; align-items:center; }
button { font:inherit; padding:7px 14px; border:1px solid #bac9da; border-radius:6px;
         color:#213547; background:white; cursor:pointer; }
button:hover:enabled { background:#e8f0fa; }
button:disabled { opacity:.4; cursor:default; }
#play { background:#176b87; color:white; }
#status { margin-top:10px; font-variant-numeric:tabular-nums; }
#chart { width:100%; background:white; }
.note { margin:6px 0 0; font-size:12px; color:#637083; }
</style></head><body>
<header><h1 id="heading"></h1><nav aria-label="Frame Playing Control">
<button id="first">⏮ First Frame</button><button id="previous">◀ Previous Frame</button>
<button id="play">▶ Play</button><button id="pause">⏸ Pause</button>
<button id="next">Next Frame▶</button><button id="last">Last Frame ⏭</button>
</nav><div id="status" role="status" aria-live="polite">Loading…</div>
<p class="note">Each frame corresponds to a selected output file. Color scales and axis limits remain fixed. Drag to zoom; double-click to reset.</p>
</header><main><div id="chart"></div></main>
<script>
const chart = document.getElementById('chart');
const status = document.getElementById('status');
const buttons = Object.fromEntries(['first','previous','play','pause','next','last']
    .map(id => [id, document.getElementById(id)]));
document.getElementById('heading').textContent = viewer.heading;
const end = frames.length - 1;
let current = -1, target = 0, playing = false, drawing = false, timer = null;
function controls() {
    buttons.first.disabled = buttons.previous.disabled = target <= 0;
    buttons.last.disabled = buttons.next.disabled = target >= end;
    buttons.play.disabled = playing || end === 0;
    buttons.pause.disabled = !playing;
}
function pause() {
    playing = false;
    clearTimeout(timer);
    controls();
}
function schedule() {
    clearTimeout(timer);
    if (!playing || drawing) return;
    if (current >= end) { pause(); return; }
    timer = setTimeout(() => { target = current + 1; draw(); }, 1000 / viewer.fps);
}
async function draw() {
    controls();
    if (drawing) return;
    drawing = true;
    try {
        // Serialize redraws. Rapid navigation records the latest requested frame,
        // so an older asynchronous redraw cannot overwrite a newer selection.
        while (current !== target) {
            const next = target, frame = frames[next];
            const data = viewer.figure.data.map((trace, i) => {
                if (trace.meta && trace.meta.streamlines) {
                    return {...trace, x: frame.streamlines.x, y: frame.streamlines.y};
                }
                const field = frame.values[Math.floor(i / 2)];
                if (i % 2 === 0) {
                    const physical = viewer.spherical
                        ? viewer.radialMap.map(row => row.map(j => j < 0 ? null : field[0][j]))
                        : viewer.is2d ? field : [field[0], field[0]];
                    const z = viewer.logFields[Math.floor(i / 2)]
                        ? physical.map(row => row.map(v => v === null ? null : Math.log10(v)))
                        : physical;
                    return {...trace, z, customdata: physical};
                }
                const y = viewer.alongY ? field.map(row => row[viewer.cut]) : field[viewer.cut];
                return {...trace, y};
            });
            await Plotly.react(chart, data, viewer.figure.layout,
                {responsive:true, displaylogo:false, scrollZoom:true});
            current = next;
            const time = frame.time === null ? 'Time not provided' : `t = ${frame.time.toPrecision(7)}${viewer.cgs ? ' s' : ''}`;
            status.textContent = `Frame ${current + 1} / ${frames.length}  ·  step ${frame.step}  ·  ${time}`;
        }
    } catch (error) {
        pause();
        status.textContent = `Plot Fail：${error.message}`;
        console.error(error);
    } finally {
        drawing = false;
        controls();
    }
    schedule();
}
function jump(index) {
    pause();
    target = Math.max(0, Math.min(end, index));
    draw();
}
buttons.first.onclick = () => jump(0);
buttons.previous.onclick = () => jump(target - 1);
buttons.next.onclick = () => jump(target + 1);
buttons.last.onclick = () => jump(end);
buttons.pause.onclick = pause;
buttons.play.onclick = () => {
    if (playing || end === 0) return;
    playing = true;
    if (target >= end) target = 0;
    controls();
    if (current !== target) draw(); else schedule();
};
draw();
</script></body></html>
"""


def save_plotly(reference, records, lower, upper, args):
    """Write a standalone viewer; stream snapshots rather than duplicating frames."""
    import html
    import json
    try:
        import plotly.graph_objects as go
        from plotly.offline import get_plotlyjs
        from plotly.subplots import make_subplots
    except ImportError as error:
        raise ValueError("Interactive HTML requires Plotly: python3 -m pip install plotly") from error

    is_2d, along_y, fixed_axis, cut, profile_axis, profile_label, cut_label = select_profile(reference, args)
    spherical = args.plot == "spherical"
    projection = spherical_projection(reference) if spherical else None
    titles = []
    for name, label, _ in reference.fields:
        label = display_label(name, label, spherical, reference.cgs)
        titles.extend((label + (" (spherical slice)" if spherical else ""),
                       f"{label}({'r' if spherical else 'y' if along_y else 'x'}){cut_label}"))
    fig = make_subplots(rows=len(reference.fields), cols=2, horizontal_spacing=0.18,
                        vertical_spacing=0.09, subplot_titles=titles)
    xedges = edges(projection[0]).tolist() if spherical else edges(reference.x).tolist()
    yedges = xedges if spherical else (edges(reference.y).tolist() if is_2d else [0, 0.5, 1])
    profile_edges = edges(profile_axis)
    for k, (name, label, cmap) in enumerate(reference.fields):
        label = display_label(name, label, spherical, reference.cgs)
        row = k + 1
        lo, hi = float(lower[k]), float(upper[k])
        logarithmic = args.log and name in ("rho", "p")
        color_ticks = {}
        if logarithmic:
            lo, hi = logarithmic_limits(lo, hi)
            lo, hi = float(np.log10(lo)), float(np.log10(hi))
            tick_step = max(1, int(np.ceil((hi-lo)/6)))
            ticks = np.arange(np.ceil(lo), np.floor(hi)+1, tick_step)
            if len(ticks) < 2:
                ticks = np.linspace(lo, hi, 3)
            color_ticks = dict(tickmode="array", tickvals=ticks.tolist(),
                               ticktext=[f"{10.0**tick:.2g}" for tick in ticks])
        elif name in ("vx", "vy"):
            hi = max(abs(lo), abs(hi), 1e-12)
            lo = -hi
        elif lo == hi:
            pad = max(abs(lo) * 0.01, 1e-12)
            lo, hi = lo - pad, hi + pad
        axis_number = 2 * k + 1
        axis_suffix = str(axis_number) if axis_number > 1 else ""
        domain = fig.layout[f"yaxis{axis_suffix}"].domain
        fig.add_trace(go.Heatmap(
            x=xedges, y=yedges, z=[], name=label,
            colorscale=cmap, zmin=lo, zmax=hi, zauto=False, zsmooth=False,
            colorbar=dict(x=0.435, y=sum(domain) / 2, len=domain[1] - domain[0],
                          thickness=12, outlinewidth=0, tickfont=dict(size=10), **color_ticks),
            hovertemplate=f"x=%{{x:.6g}}" + ("<br>z=%{y:.6g}" if spherical else "<br>y=%{y:.6g}" if is_2d else "")
                          + f"<br>{label}=%{{customdata:.6g}}<extra></extra>"), row=row, col=1)
        fig.add_trace(go.Scatter(x=profile_axis.tolist(), y=[], mode="lines", name=label,
                                line=dict(color="#176b87", width=2),
                                hovertemplate=f"{profile_label}=%{{x:.6g}}<br>{label}=%{{y:.6g}}<extra></extra>"),
                      row=row, col=2)
        fig.update_xaxes(title_text=coordinate_label("x" if spherical else reference.xlabel, reference.cgs), range=[xedges[0], xedges[-1]],
                         autorange=False, row=row, col=1)
        fig.update_yaxes(title_text=coordinate_label("z" if spherical else "y", reference.cgs) if (spherical or is_2d) else "repeated 1D", range=[yedges[0], yedges[-1]],
                         autorange=False, showticklabels=is_2d or spherical, row=row, col=1)
        if is_2d or spherical:
            fig.update_yaxes(scaleanchor=f"x{axis_suffix}", scaleratio=1, constrain="domain", row=row, col=1)
        if is_2d:
            # exclude_empty_subplots=False: data arrays are populated in the browser.
            if along_y:
                fig.add_vline(x=float(fixed_axis[cut]), line_dash="dash", line_width=1,
                              row=row, col=1, exclude_empty_subplots=False)
            else:
                fig.add_hline(y=float(fixed_axis[cut]), line_dash="dash", line_width=1,
                              row=row, col=1, exclude_empty_subplots=False)
        pad = max(float(upper[k] - lower[k]) * 0.05,
                  max(abs(float(lower[k])), abs(float(upper[k]))) * 0.01, 1e-12)
        fig.update_xaxes(title_text=profile_label, range=[float(profile_edges[0]), float(profile_edges[-1])],
                         autorange=False, row=row, col=2)
        fig.update_yaxes(title_text=label, type="log" if logarithmic else "linear",
                         range=[lo, hi] if logarithmic else [float(lower[k]) - pad, float(upper[k]) + pad],
                         autorange=False, row=row, col=2)
    # Append overlays after all heatmap/profile pairs to preserve their indices.
    if args.streamlines:
        for k, (name, _, _) in enumerate(reference.fields):
            if name in ("rho", "p"):
                for color, width in (("#333333", 2.2), ("white", 1)):
                    fig.add_trace(go.Scatter(x=[], y=[], mode="lines", name="Streamlines",
                                            line=dict(color=color, width=width),
                                            meta=dict(streamlines=True), hoverinfo="skip"),
                                  row=k + 1, col=1)
    fig.update_layout(height=310 * len(reference.fields), margin=dict(l=65, r=35, t=45, b=55),
                      template="plotly_white", showlegend=False, uirevision="fixed-grid",
                      font=dict(family="Arial, sans-serif", size=12))
    heading = f"{args.prefix} · {'1D spherical' if spherical else '2D' if is_2d else '1D'}"
    settings = dict(heading=heading, fps=args.fps, is2d=is_2d, alongY=along_y, cut=cut,
                    spherical=spherical, cgs=reference.cgs,
                    logFields=[bool(args.log and name in ("rho", "p")) for name, _, _ in reference.fields],
                    figure=json.loads(fig.to_json()))
    if spherical:
        settings["radialMap"] = np.where(projection[2], -1, projection[1]).tolist()

    def safe_json(value):
        # A filename prefix may contain HTML/script delimiters.
        return json.dumps(value, ensure_ascii=True, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")

    suffix = ("_spherical" if spherical else "") + ("_cgs" if reference.cgs else "") + ("_log" if args.log else "")
    destination = args.output or args.directory / "animations" / f"{args.prefix.rstrip('_')}{suffix}_evolution.html"
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as stream:
        stream.write('<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">'
                     '<meta name="viewport" content="width=device-width,initial-scale=1">'
                     f'<title>{html.escape(heading)}</title><script>')
        stream.write(get_plotlyjs())
        stream.write('</script><script>const viewer=')
        stream.write(safe_json(settings))
        stream.write(';const frames=[')
        for index, (step, path, stamp) in enumerate(records):
            if signature(path) != stamp:
                raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
            snapshot = read_snapshot(path, step, cgs=reference.cgs)
            check_grid(snapshot, reference)
            if signature(path) != stamp:
                raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
            if index:
                stream.write(",")
            frame = dict(step=step, time=snapshot.time, values=snapshot.values.tolist())
            if args.streamlines:
                sx, sy = [], []
                for segment in streamline_segments(snapshot):
                    sx.extend(segment[:, 0].tolist() + [None])
                    sy.extend(segment[:, 1].tolist() + [None])
                frame["streamlines"] = dict(x=sx, y=sy)
            stream.write(safe_json(frame))
            if index % 20 == 0 or index == len(records) - 1:
                print(f"Embedding {index + 1}/{len(records)}", flush=True)
        stream.write('];</script>')
        stream.write(PLOTLY_VIEWER)
    print(f"Saved {destination} ({destination.stat().st_size / 1024**2:.1f} MiB)\n"
          "Open this HTML file in a browser; no server or internet connection required.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", type=Path, help="Directory containing snapshots")
    parser.add_argument("prefix", help="Filename prefix, e.g. mesh_sod or mesh_sod_")
    parser.add_argument("--stride", type=int, default=1, help="One frame per N saved outputs; always keep last")
    parser.add_argument("--fps", type=float, default=8, help="Playback frames per second (default: 8)")
    parser.add_argument("--dpi", type=int, default=100)
    parser.add_argument("--plot", choices=("cartesian", "spherical"), default="cartesian",
                        help="Map style: Cartesian strips/maps (default), or 1D spherical symmetry slices")
    parser.add_argument("--cgs", action="store_true",
                        help="Convert plots to CGS using CSV unit metadata; slice coordinates are then in cm")
    parser.add_argument("--log", action="store_true",
                        help="Logarithmic color scales and profile y axes for positive rho/pressure; velocity stays linear")
    parser.add_argument("--plotly", action="store_true", help="Save an offline interactive HTML viewer")
    parser.add_argument("--streamlines", action="store_true",
                        help="Overlay velocity streamlines on 2D density and pressure maps")
    slices = parser.add_mutually_exclusive_group()
    slices.add_argument("--slice-y", type=float, help="Fix y and plot along x (nearest row; default: middle y)")
    slices.add_argument("--slice-x", type=float, help="Fix x and plot along y (nearest column; 2D only)")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--output", type=Path, help="Destination .gif, .mp4 or interactive .html")
    mode.add_argument("--show", action="store_true", help="Play in a GUI window instead of saving")
    args = parser.parse_args()
    if args.stride < 1 or args.dpi < 1 or not np.isfinite(args.fps) or args.fps <= 0:
        parser.error("Require stride >= 1, dpi >= 1 and finite fps > 0")
    if not args.prefix.rstrip("_") or Path(args.prefix).name != args.prefix:
        parser.error("prefix must be a filename prefix without directory components")
    if any(value is not None and not np.isfinite(value) for value in (args.slice_x, args.slice_y)):
        parser.error("slice-x and slice-y must be finite")
    if args.output and args.output.suffix.lower() not in (".gif", ".mp4", ".html"):
        parser.error("output must end in .gif, .mp4 or .html")
    if args.plotly and (args.show or (args.output and args.output.suffix.lower() != ".html")):
        parser.error("--plotly requires HTML output and cannot be combined with --show")
    try:
        files = find_files(args.directory, args.prefix, args.stride)
        reference, records, lower, upper = prepare(files, cgs=args.cgs, log=args.log)
        if args.plotly or (args.output and args.output.suffix.lower() == ".html"):
            save_plotly(reference, records, lower, upper, args)
        else:
            animate(reference, records, lower, upper, args)
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f"Animation failed: {error}\n")


if __name__ == "__main__":
    main()
