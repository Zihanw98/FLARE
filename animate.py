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

Requires numpy, matplotlib, Pillow (ffmpeg additionally for MP4).
Interactive HTML requires numpy and plotly: python3 -m pip install plotly
Reads PREFIX_<step>.csv (x, optional y/z, rho, vx, vy, p or pressure,
optional is_ghost) or legacy PREFIX_<step>.txt (rho vx p plus FLARE header).
Prefix may include its trailing underscore. Legacy TXT uses cell indices,
because physical coordinates are not stored in those files.

1D: rho/vx/pressure strips repeated vertically, with synchronized x profiles.
2D: rho/vx/vy/pressure maps, with x profiles at the nearest --slice-y row
(default: domain midpoint), or y profiles at the nearest --slice-x column.
These options are mutually exclusive. A dashed line marks the cut on each map.
Scales are fixed over selected frames; ghost cells are excluded. Stride counts
saved files, not solver steps. First/last frames are retained. Playback uses
equal frame durations, not equal physical-time intervals.

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


def read_snapshot(path, step):
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
                    "active cell index" if path.suffix == ".txt" else "x")


def check_grid(snapshot, reference):
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


def prepare(files):
    reference = None
    records = []
    previous_time = -np.inf
    for i, (step, path) in enumerate(files):
        stamp = signature(path)
        snapshot = read_snapshot(path, step)
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
    profile_label = "y" if along_y else reference.xlabel
    cut_label = f" at {'x' if along_y else 'y'} = {fixed_axis[cut_index]:.6g}" if is_2d else ""
    return is_2d, along_y, fixed_axis, cut_index, profile_axis, profile_label, cut_label


def animate(reference, records, lower, upper, args):
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter, writers

    is_2d, along_y, fixed_axis, cut_index, profile_axis, profile_label, cut_label = select_profile(reference, args)
    profile_edges = edges(profile_axis)

    def profile_values(snapshot, k):
        return snapshot.values[k, :, cut_index] if along_y else snapshot.values[k, cut_index, :]

    nfields = len(reference.fields)
    fig, axes = plt.subplots(nfields, 2, figsize=(12, 3 * nfields), layout="constrained")
    title = fig.suptitle("")
    maps, lines = [], []
    xedges = edges(reference.x)
    yedges = edges(reference.y) if is_2d else np.linspace(0, 1, 25)

    def map_values(snapshot, k):
        return snapshot.values[k] if is_2d else np.repeat(snapshot.values[k], 24, axis=0)

    for k, (name, label, cmap) in enumerate(reference.fields):
        lo, hi = lower[k], upper[k]
        if name in ("vx", "vy"):
            hi = max(abs(lo), abs(hi), 1e-12)
            lo = -hi
        elif lo == hi:
            pad = max(abs(lo) * 0.01, 1e-12)
            lo, hi = lo - pad, hi + pad
        left, right = axes[k]
        mesh = left.pcolormesh(xedges, yedges, map_values(reference, k),
                               shading="flat", cmap=cmap, vmin=lo, vmax=hi)
        fig.colorbar(mesh, ax=left, pad=0.02)
        left.set(title=label, xlabel=reference.xlabel)
        if is_2d:
            left.set(ylabel="y", aspect="equal")
            mark_cut = left.axvline if along_y else left.axhline
            mark_cut(fixed_axis[cut_index], color="black", linestyle="--", linewidth=0.8)
        else:
            left.set_yticks([])
            left.set_ylabel("repeated 1D")
        line, = right.plot(profile_axis, profile_values(reference, k), linewidth=1.5)
        pad = max((upper[k] - lower[k]) * 0.05, max(abs(lower[k]), abs(upper[k])) * 0.01, 1e-12)
        right.set(xlim=(profile_edges[0], profile_edges[-1]), ylim=(lower[k] - pad, upper[k] + pad),
                  xlabel=profile_label, ylabel=label,
                  title=f"{label}({'y' if along_y else 'x'})" + cut_label)
        right.grid(alpha=0.25)
        maps.append(mesh)
        lines.append(line)

    def update(index):
        step, path, stamp = records[index]
        if signature(path) != stamp:
            raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
        snapshot = read_snapshot(path, step)
        check_grid(snapshot, reference)
        if signature(path) != stamp:
            raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
        time_label = f"t = {snapshot.time:.7g}" if snapshot.time is not None else "time unavailable"
        title.set_text(f"{args.prefix} | {'2D' if is_2d else '1D'} | step {step} | {time_label}"
                       f" | frame {index + 1}/{len(records)}")
        for k, (mesh, line) in enumerate(zip(maps, lines)):
            mesh.set_array(map_values(snapshot, k).ravel())
            line.set_ydata(profile_values(snapshot, k))
        return [title, *maps, *lines]

    try:
        if args.show:
            update(0)
            movie = FuncAnimation(fig, update, frames=len(records), interval=1000 / args.fps,
                                  blit=False, cache_frame_data=False)
            plt.show()
        else:
            destination = args.output or args.directory / "animations" / f"{args.prefix.rstrip('_')}_evolution.gif"
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
                const field = frame.values[Math.floor(i / 2)];
                if (i % 2 === 0) {
                    return {...trace, z: viewer.is2d ? field : [field[0], field[0]]};
                }
                const y = viewer.alongY ? field.map(row => row[viewer.cut]) : field[viewer.cut];
                return {...trace, y};
            });
            await Plotly.react(chart, data, viewer.figure.layout,
                {responsive:true, displaylogo:false, scrollZoom:true});
            current = next;
            const time = frame.time === null ? 'Time not provided' : `t = ${frame.time.toPrecision(7)}`;
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
    titles = []
    for _, label, _ in reference.fields:
        titles.extend((label, f"{label}({'y' if along_y else 'x'}){cut_label}"))
    fig = make_subplots(rows=len(reference.fields), cols=2, horizontal_spacing=0.18,
                        vertical_spacing=0.09, subplot_titles=titles)
    xedges = edges(reference.x).tolist()
    yedges = edges(reference.y).tolist() if is_2d else [0, 0.5, 1]
    profile_edges = edges(profile_axis)
    for k, (name, label, cmap) in enumerate(reference.fields):
        row = k + 1
        lo, hi = float(lower[k]), float(upper[k])
        if name in ("vx", "vy"):
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
                          thickness=12, outlinewidth=0, tickfont=dict(size=10)),
            hovertemplate=f"x=%{{x:.6g}}" + ("<br>y=%{y:.6g}" if is_2d else "")
                          + f"<br>{label}=%{{z:.6g}}<extra></extra>"), row=row, col=1)
        fig.add_trace(go.Scatter(x=profile_axis.tolist(), y=[], mode="lines", name=label,
                                line=dict(color="#176b87", width=2),
                                hovertemplate=f"{profile_label}=%{{x:.6g}}<br>{label}=%{{y:.6g}}<extra></extra>"),
                      row=row, col=2)
        fig.update_xaxes(title_text=reference.xlabel, range=[xedges[0], xedges[-1]],
                         autorange=False, row=row, col=1)
        fig.update_yaxes(title_text="y" if is_2d else "repeated 1D", range=[yedges[0], yedges[-1]],
                         autorange=False, showticklabels=is_2d, row=row, col=1)
        if is_2d:
            fig.update_yaxes(scaleanchor=f"x{axis_suffix}", scaleratio=1, constrain="domain", row=row, col=1)
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
        fig.update_yaxes(title_text=label, range=[float(lower[k]) - pad, float(upper[k]) + pad],
                         autorange=False, row=row, col=2)
    fig.update_layout(height=310 * len(reference.fields), margin=dict(l=65, r=35, t=45, b=55),
                      template="plotly_white", showlegend=False, uirevision="fixed-grid",
                      font=dict(family="Arial, sans-serif", size=12))
    heading = f"{args.prefix} · {'2D' if is_2d else '1D'}"
    settings = dict(heading=heading, fps=args.fps, is2d=is_2d, alongY=along_y, cut=cut,
                    figure=json.loads(fig.to_json()))

    def safe_json(value):
        # A filename prefix may contain HTML/script delimiters.
        return json.dumps(value, ensure_ascii=True, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")

    destination = args.output or args.directory / "animations" / f"{args.prefix.rstrip('_')}_evolution.html"
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
            snapshot = read_snapshot(path, step)
            check_grid(snapshot, reference)
            if signature(path) != stamp:
                raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
            if index:
                stream.write(",")
            stream.write(safe_json(dict(step=step, time=snapshot.time, values=snapshot.values.tolist())))
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
    parser.add_argument("--plotly", action="store_true", help="Save an offline interactive HTML viewer")
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
        reference, records, lower, upper = prepare(files)
        if args.plotly or (args.output and args.output.suffix.lower() == ".html"):
            save_plotly(reference, records, lower, upper, args)
        else:
            animate(reference, records, lower, upper, args)
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f"Animation failed: {error}\n")


if __name__ == "__main__":
    main()
