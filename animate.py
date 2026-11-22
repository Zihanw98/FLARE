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
  python3 animate.py tests/CCSN_MESA_day13/output ccsn_mesa --plot spherical --cgs --log --rmax 1.4e13
  python3 animate.py tests/Hydrogen_Profile_2D/output hydrogen_profile --composition

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
--radial-grid reads exact active shell edges from a CSV with r_left,r_right
columns in code units, in the same radial order as the snapshot. Use this for
stretched meshes whose edges cannot be inferred accurately from cell centers.
--rmax sets the displayed outer radius for spherical maps and profiles, in cm
with --cgs, otherwise in code length units. Maps are sampled at this radius
range; color scales and profile y limits still use all cells and selected frames.
--profile-range MIN MAX (alias --xlim) overrides the horizontal length range
of the right-hand profiles only. --log-x (alias --log-radius) uses a logarithmic
profile length axis; limits must be positive. Without explicit limits the lower
bound is the first positive cell center when the domain reaches zero.
Left-hand spherical maps keep signed linear coordinates.
2D: rho/vx/vy/pressure maps, with x profiles at the nearest --slice-y row
(default: domain midpoint), or y profiles at the nearest --slice-x column.
These options are mutually exclusive. A dashed line marks the cut on each map.
--composition adds every comp_<index> mass fraction, labeled by the CSV
comp_<index>_name metadata (for example HI, HII, He_I). These dimensionless
fields use a fixed linear [0,1] scale, including all-zero species, with --cgs/--log too.
--thermodynamics adds mean_molecular_weight and temperature. With --cgs,
temperature is converted to kelvin; molecular weight remains dimensionless.
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
--density-threshold RHO (with --log) emphasizes density near RHO using a fixed
piecewise-log blue/yellow/red scale; yellow marks RHO. RHO is in displayed units.
--density-min RHO sets the density color scale and profile lower limit in displayed
units. Lower densities use the bottom color; physical data are unchanged.
--density-max RHO similarly sets the upper display limit. Density limits cannot
be combined with --density-threshold.
--map-pixels controls spherical map sampling (default: 241); increase it together
with --dpi for sharper output. It does not change the simulation grid.

Default output: DIRECTORY/animations/PREFIX_evolution.gif and *_last.png.
Use --output movie.mp4 for MP4, or --show to play in a GUI without saving.
Use --plotly or --output viewer.html for an offline interactive HTML viewer.
Its six controls are first/previous/play/pause/next/last. It starts paused and
stops at the last frame. HTML embeds all selected data: increase --stride for
large runs. --fps controls playback speed; --dpi applies only to GIF/MP4/PNG.
Files are discovered once at startup; rerun to include new simulation output.
"""

import argparse
from dataclasses import dataclass, field as dataclass_field
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
    metadata: dict = dataclass_field(default_factory=dict)
    radial_edges: object = None  # Optional exact shell edges in displayed length units.


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


def read_snapshot(path, step, cgs=False, fields=None, composition=False, thermodynamics=False):
    """Read standard fields, or caller-selected (column, label, colormap) tuples."""
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
                       "vx": length / duration, "vy": length / duration, "vz": length / duration,
                       "temperature": _temperature}
            for key in ("internal_energy", "total_energy", "dual_energy_correction"):
                factors[key] = factors["p"]  # Energy densities, not specific energies.
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
    if fields is None:
        fields = FIELDS
        if composition:
            components = sorted((key for key in columns if re.fullmatch(r"comp_\d+", key)),
                                key=lambda key: int(key.split("_")[1]))
            if not components:
                raise ValueError(f"--composition requires comp_<index> columns: {path}")
            fields += tuple((key, f"{metadata.get(key + '_name', key)} mass fraction", "viridis")
                            for key in components)
        if thermodynamics:
            fields += (("mean_molecular_weight", "mean molecular weight", "viridis"),
                       ("temperature", "temperature", "inferno"))
    fields = tuple(field for field in fields
                   if len(y) > 1 or field[0] != "vy")
    required = ["x", *[field[0] for field in fields]]
    required += [key for key in ("y", "z") if key in columns]
    if any(key not in columns for key in required):
        raise ValueError(f"Missing requested field columns: {path}")
    if any(not np.isfinite(columns[key]).all() for key in required):
        raise ValueError(f"Non-finite coordinates or field values: {path}")
    if len(x) < 2 or len(ix) != len(x) * len(y) or len(np.unique(iy * len(x) + ix)) != len(ix):
        raise ValueError(f"Expected a complete Cartesian grid without duplicate cells: {path}")
    values = np.empty((len(fields), len(y), len(x)))
    for k, (key, _, _) in enumerate(fields):
        values[k, iy, ix] = columns[key]
    return Snapshot(step, time, x, y, values, fields,
                    "active cell index" if path.suffix == ".txt" else "x", units, cgs, metadata)


def check_grid(snapshot, reference):
    if snapshot.units != reference.units or snapshot.cgs != reference.cgs:
        raise ValueError("Units changed between snapshots; use files from a single unit system")
    if snapshot.fields != reference.fields:
        raise ValueError("Fields or species changed between snapshots; check for files from older runs")
    if snapshot.xlabel != reference.xlabel or any(
        a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=1e-13)
        for a, b in ((snapshot.x, reference.x), (snapshot.y, reference.y))
    ):
        raise ValueError("Grid changed between snapshots; use files from a single fixed-grid run")


def edges(centers):
    """Infer cell edges from centers, including nonuniform Cartesian grids."""
    return np.r_[centers[0] - (centers[1] - centers[0]) / 2,
                 (centers[:-1] + centers[1:]) / 2,
                 centers[-1] + (centers[-1] - centers[-2]) / 2]


def read_radial_grid(path, snapshot):
    """Read fixed active shell geometry, checking it against snapshot centers."""
    before = signature(path)
    table = np.genfromtxt(path, names=True, delimiter=",", ndmin=1)
    if signature(path) != before:
        raise ValueError("Radial grid changed while reading; rerun after output finishes")
    if not {"r_left", "r_right"}.issubset(table.dtype.names or ()):
        raise ValueError("--radial-grid requires CSV columns r_left,r_right")
    left, right = table["r_left"], table["r_right"]
    if snapshot.cgs:
        left, right = left * snapshot.units[1], right * snapshot.units[1]
    if (len(snapshot.y) != 1 or len(left) != len(snapshot.x)
        or not np.all(np.isfinite(left)) or not np.all(np.isfinite(right))
        or np.any(left < 0) or np.any(right <= left)
        or not np.allclose(right[:-1], left[1:], rtol=1e-12, atol=0)):
        raise ValueError("--radial-grid requires one contiguous positive-width shell per active cell")
    if not np.allclose(0.5*left+0.5*right, snapshot.x, rtol=1e-12, atol=0):
        raise ValueError("--radial-grid shell centers do not match the snapshot")
    return np.r_[left, right[-1]]


def spherical_projection(snapshot, pixels=241, rmax=None):
    """Map a 1D radial grid onto a symmetry plane without interpolating states.

    Use supplied shell edges; otherwise infer them as in the Cartesian plots.
    The returned mask excludes both the exterior and any inner cavity.
    """
    if not isinstance(pixels, (int, np.integer)) or pixels < 3:
        raise ValueError("--map-pixels must be an integer >= 3")
    if len(snapshot.y) != 1:
        raise ValueError("--plot spherical requires a 1D radial snapshot")
    if snapshot.xlabel != "x":
        raise ValueError("--plot spherical requires CSV physical radii, not legacy cell indices")
    radial_edges = (snapshot.radial_edges.copy() if snapshot.radial_edges is not None
                    else edges(snapshot.x))
    tolerance = 1e-12 * max(abs(radial_edges[-1]), 1.0)
    if snapshot.x[0] < 0 or radial_edges[0] < -tolerance:
        raise ValueError("--plot spherical requires nonnegative active radii and shell edges")
    radial_edges[0] = max(0.0, radial_edges[0])
    outer = radial_edges[-1] if rmax is None else rmax
    if not np.isfinite(outer) or outer <= radial_edges[0]:
        raise ValueError("--rmax must be finite and greater than the inner shell edge")
    axis = np.linspace(-outer, outer, pixels)
    radius = np.hypot(axis[:, None], axis[None, :])
    indices = np.searchsorted(radial_edges, radius, side="right") - 1
    outside = (radius < radial_edges[0]) | (radius > min(outer, radial_edges[-1]))
    return axis, np.clip(indices, 0, len(snapshot.x) - 1), outside


def spherical_values(values, projection):
    """Project one field while keeping radial velocity a scalar radial component."""
    _, indices, outside = projection
    return np.ma.array(values[0, indices], mask=outside)


def display_label(name, label, spherical, cgs=False):
    label = "radial velocity" if spherical and name == "vx" else label
    unit = {"rho": "g/cm^3", "vx": "cm/s", "vy": "cm/s", "p": "dyn/cm^2", "temperature": "K"}
    return f"{label} [{unit[name]}]" if cgs and name in unit else label


def coordinate_label(label, cgs):
    return f"{label} [cm]" if cgs else label


def logarithmic_limits(lower, upper):
    """Positive limits shared by both backends, including constant fields."""
    if not (np.isfinite(lower) and np.isfinite(upper) and 0 < lower <= upper):
        raise ValueError("--log requires strictly positive finite rho and pressure; zero/negative values cannot be displayed on a log scale")
    if lower == upper:
        return lower / 1.1, upper * 1.1
    return lower, upper


DENSITY_COLORS = ((0., "#081d58"), (.25, "#41b6c4"), (.5, "#ffffcc"),
                  (.75, "#fd8d3c"), (1., "#800026"))


def density_display_limits(lower, upper, minimum=None, maximum=None):
    """Apply optional display limits without clipping physical samples."""
    if minimum is None and maximum is None:
        return lower, upper
    lo = lower if minimum is None else minimum
    hi = upper if maximum is None else maximum
    if not np.isfinite(lo) or not np.isfinite(hi) or not 0 < lo < hi:
        raise ValueError("Density display limits must be finite and satisfy 0 < minimum < maximum")
    return lo, hi


def density_color_scale(lower, upper, threshold):
    """Keep all densities, allocating 60% of color space to 0.5–1.5 threshold.

    Bounds include the threshold window even if the selected states do not
    cross it. Knots are fixed across frames; no snapshot values are modified.
    """
    logarithmic_limits(lower, upper)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("--density-threshold must be finite and positive")
    densities = np.array([min(lower, .05*threshold), .1*threshold,
                          .5*threshold, threshold, 1.5*threshold,
                          max(upper, 2*threshold)])
    if not np.isfinite(densities).all() or np.any(densities <= 0):
        raise ValueError("--density-threshold window exceeds floating-point range")
    return np.log10(densities), np.array([0., .08, .25, .5, .85, 1.])


def density_normalization(knots, positions):
    """Invertible physical-density normalization, including projection masks."""
    from matplotlib.colors import FuncNorm

    def forward(values):
        values = np.ma.asarray(values)
        logs = np.ma.log10(values)
        return np.ma.array(np.interp(logs.filled(knots[0]), knots, positions),
                           mask=np.ma.getmaskarray(logs))

    def inverse(values):
        values = np.ma.asarray(values)
        return np.ma.array(10.**np.interp(values.filled(0.), positions, knots),
                           mask=np.ma.getmaskarray(values))

    return FuncNorm((forward, inverse), vmin=10.**knots[0], vmax=10.**knots[-1], clip=True)


def prepare(files, cgs=False, log=False, composition=False, thermodynamics=False):
    reference = None
    records = []
    previous_time = -np.inf
    for i, (step, path) in enumerate(files):
        stamp = signature(path)
        snapshot = read_snapshot(path, step, cgs=cgs, composition=composition,
                                 thermodynamics=thermodynamics)
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


def profile_display_limits(axis, cell_edges, logarithmic=False, bounds=None):
    """Shared physical-coordinate limits for Matplotlib and Plotly profiles."""
    lo, hi = (cell_edges[0], cell_edges[-1]) if bounds is None else bounds
    if logarithmic and bounds is None and lo <= 0:
        positive = axis[axis > 0]
        if not len(positive):
            raise ValueError("--log-x requires positive profile coordinates")
        lo = positive[0]
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        raise ValueError("--profile-range requires finite increasing bounds")
    if logarithmic and lo <= 0:
        raise ValueError("--log-x requires strictly positive --profile-range bounds")
    if hi < cell_edges[0] or lo > cell_edges[-1]:
        raise ValueError("--profile-range does not overlap the displayed domain")
    return float(lo), float(hi)


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
    from matplotlib.colors import LogNorm, LinearSegmentedColormap
    from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter, writers

    is_2d, along_y, fixed_axis, cut_index, profile_axis, profile_label, cut_label = select_profile(reference, args)
    spherical = args.plot == "spherical"
    projection = spherical_projection(reference, pixels=getattr(args, "map_pixels", 241),
                                      rmax=args.rmax) if spherical else None
    profile_edges = (reference.radial_edges.copy() if spherical and reference.radial_edges is not None
                     else edges(profile_axis))
    if spherical and args.rmax is not None and getattr(args, "profile_range", None) is None:
        profile_edges[-1] = args.rmax
    log_x = getattr(args, "log_x", False)
    profile_limits = profile_display_limits(profile_axis, profile_edges, log_x,
                                             getattr(args, "profile_range", None))
    plotted_axis = np.ma.masked_less_equal(profile_axis, 0) if log_x else profile_axis

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
        if name == "rho":
            lo, hi = density_display_limits(lo, hi, getattr(args, "density_min", None),
                                             getattr(args, "density_max", None))
        density_bounds = name == "rho" and any(getattr(args, option, None) is not None
                                               for option in ("density_min", "density_max"))
        logarithmic = args.log and name in ("rho", "p")
        if logarithmic:
            lo, hi = logarithmic_limits(lo, hi)
        elif name.startswith("comp_"):
            lo, hi = 0.0, 1.0
        elif name in ("vx", "vy"):
            hi = max(abs(lo), abs(hi), 1e-12)
            lo = -hi
        elif lo == hi:
            pad = max(abs(lo) * 0.01, 1e-12)
            lo, hi = lo - pad, hi + pad
        left, right = axes[k]
        color_scale = {"norm": LogNorm(vmin=lo, vmax=hi)} if logarithmic else {"vmin": lo, "vmax": hi}
        threshold = getattr(args, "density_threshold", None) if name == "rho" else None
        if threshold is not None:
            knots, positions = density_color_scale(lo, hi, threshold)
            color_scale = {"norm": density_normalization(knots, positions)}
            cmap = LinearSegmentedColormap.from_list("density_threshold", DENSITY_COLORS, N=1024)
        if spherical:
            # Projection pixels are uniform: draw the same shell samples directly
            # instead of rasterizing hundreds of thousands of quadrilaterals.
            mesh = left.imshow(map_values(reference, k), origin="lower", interpolation="nearest",
                               extent=(xedges[0], xedges[-1], yedges[0], yedges[-1]),
                               cmap=cmap, **color_scale)
        else:
            mesh = left.pcolormesh(xedges, yedges, map_values(reference, k),
                                   shading="flat", cmap=cmap, **color_scale)
        colorbar = fig.colorbar(mesh, ax=left, pad=0.02)
        if threshold is not None:
            ticks = 10.**knots
            colorbar.set_ticks(ticks, labels=[f"{value:.2g}" + (" (threshold)" if i == 3 else "")
                                             for i, value in enumerate(ticks)])
            colorbar.ax.minorticks_off()
            colorbar.ax.axhline(threshold, color="#222222", lw=1.2)
            colorbar.set_label("Density: threshold-focused scale")
            right.axhline(threshold, color="#ad2844", ls="--", lw=1.2,
                          label=f"Threshold = {threshold:.2g}")
            right.legend(loc="lower left", fontsize=8)
        left.set(title=label + (" (spherical slice)" if spherical else ""),
                 xlabel=coordinate_label("x" if spherical else reference.xlabel, reference.cgs))
        if spherical:
            outer = projection[0][-1]
            left.set(ylabel=coordinate_label("z", reference.cgs), aspect="equal", xlim=(-outer, outer), ylim=(-outer, outer))
        elif is_2d:
            left.set(ylabel=coordinate_label("y", reference.cgs), aspect="equal")
            mark_cut = left.axvline if along_y else left.axhline
            mark_cut(fixed_axis[cut_index], color="black", linestyle="--", linewidth=0.8)
        else:
            left.set_yticks([])
            left.set_ylabel("repeated 1D")
        line, = right.plot(plotted_axis, profile_values(reference, k), linewidth=1.5)
        if name.startswith("comp_"):
            # Keep X=0 and X=1 curves visible on the fixed range's boundary.
            line.set(zorder=3, clip_on=False)
        pad = max((upper[k] - lower[k]) * 0.05, max(abs(lower[k]), abs(upper[k])) * 0.01, 1e-12)
        if logarithmic:
            right.set_yscale("log")
        if log_x:
            right.set_xscale("log")
        right.set(xlim=profile_limits,
                  ylim=(lo, hi) if logarithmic or name.startswith("comp_") or density_bounds
                  else (lower[k] - pad, upper[k] + pad),
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
        snapshot = read_snapshot(path, step, cgs=reference.cgs,
                                 composition=any(name.startswith("comp_") for name, _, _ in reference.fields),
                                 thermodynamics=any(name == "temperature" for name, _, _ in reference.fields))
        check_grid(snapshot, reference)
        if signature(path) != stamp:
            raise ValueError(f"Snapshot changed since scanning: {path}; rerun")
        time_label = f"t = {snapshot.time:.7g}" + (" s" if reference.cgs else "") if snapshot.time is not None else "time unavailable"
        title.set_text(f"{args.prefix} | {'1D spherical' if spherical else '2D' if is_2d else '1D'} | step {step} | {time_label}"
                       f" | frame {index + 1}/{len(records)}")
        for k, (mesh, line) in enumerate(zip(maps, lines)):
            if spherical:
                mesh.set_data(map_values(snapshot, k))
            else:
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
function densityColor(value) {
    if (value === null) return null;
    const {knots, positions} = viewer.densityScale;
    const x = Math.log10(value);
    if (x <= knots[0]) return 0;
    for (let j = 1; j < knots.length; j++) {
        if (x <= knots[j]) return positions[j-1] + (positions[j]-positions[j-1])
            * (x-knots[j-1]) / (knots[j]-knots[j-1]);
    }
    return 1;
}
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
                    const z = viewer.densityScale && i === viewer.densityScale.field * 2
                        ? physical.map(row => row.map(densityColor))
                        : viewer.logFields[Math.floor(i / 2)]
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
    projection = spherical_projection(reference, pixels=getattr(args, "map_pixels", 241),
                                      rmax=args.rmax) if spherical else None
    titles = []
    for name, label, _ in reference.fields:
        label = display_label(name, label, spherical, reference.cgs)
        titles.extend((label + (" (spherical slice)" if spherical else ""),
                       f"{label}({'r' if spherical else 'y' if along_y else 'x'}){cut_label}"))
    fig = make_subplots(rows=len(reference.fields), cols=2, horizontal_spacing=0.18,
                        vertical_spacing=min(0.09, 0.5 / max(len(reference.fields) - 1, 1)),
                        subplot_titles=titles)
    xedges = edges(projection[0]).tolist() if spherical else edges(reference.x).tolist()
    yedges = xedges if spherical else (edges(reference.y).tolist() if is_2d else [0, 0.5, 1])
    profile_edges = (reference.radial_edges.copy() if spherical and reference.radial_edges is not None
                     else edges(profile_axis))
    if spherical and args.rmax is not None and getattr(args, "profile_range", None) is None:
        profile_edges[-1] = args.rmax
    log_x = getattr(args, "log_x", False)
    profile_limits = profile_display_limits(profile_axis, profile_edges, log_x,
                                             getattr(args, "profile_range", None))
    plotted_axis = [None if log_x and value <= 0 else float(value) for value in profile_axis]
    density_scale = None
    for k, (name, label, cmap) in enumerate(reference.fields):
        label = display_label(name, label, spherical, reference.cgs)
        row = k + 1
        lo, hi = float(lower[k]), float(upper[k])
        if name == "rho":
            lo, hi = density_display_limits(lo, hi, getattr(args, "density_min", None),
                                             getattr(args, "density_max", None))
        density_bounds = name == "rho" and any(getattr(args, option, None) is not None
                                               for option in ("density_min", "density_max"))
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
        elif name.startswith("comp_"):
            lo, hi = 0.0, 1.0
        elif name in ("vx", "vy"):
            hi = max(abs(lo), abs(hi), 1e-12)
            lo = -hi
        elif lo == hi:
            pad = max(abs(lo) * 0.01, 1e-12)
            lo, hi = lo - pad, hi + pad
        profile_lo, profile_hi = lo, hi
        threshold = getattr(args, "density_threshold", None) if name == "rho" else None
        if threshold is not None:
            knots, positions = density_color_scale(float(lower[k]), float(upper[k]), threshold)
            density_scale = dict(field=k, knots=knots.tolist(), positions=positions.tolist())
            cmap = list(DENSITY_COLORS)
            lo, hi = 0., 1.
            color_ticks = dict(tickmode="array", tickvals=positions.tolist(),
                               ticktext=[f"{10.**value:.2g}" + (" (threshold)" if i == 3 else "")
                                         for i, value in enumerate(knots)])
            fig.add_hline(y=threshold, line_dash="dash", line_color="#ad2844",
                          row=row, col=2, exclude_empty_subplots=False)
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
        fig.add_trace(go.Scatter(x=plotted_axis, y=[], mode="lines", name=label,
                                line=dict(color="#176b87", width=2),
                                hovertemplate=f"{profile_label}=%{{x:.6g}}<br>{label}=%{{y:.6g}}<extra></extra>"),
                      row=row, col=2)
        fig.update_xaxes(title_text=coordinate_label("x" if spherical else reference.xlabel, reference.cgs), range=[xedges[0], xedges[-1]],
                         autorange=False, row=row, col=1)
        fig.update_yaxes(title_text=coordinate_label("z" if spherical else "y", reference.cgs) if (spherical or is_2d) else "repeated 1D", range=[yedges[0], yedges[-1]],
                         autorange=False, showticklabels=is_2d or spherical, row=row, col=1)
        if spherical:
            outer = float(projection[0][-1])
            fig.update_xaxes(range=[-outer, outer], row=row, col=1)
            fig.update_yaxes(range=[-outer, outer], row=row, col=1)
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
        fig.update_xaxes(title_text=profile_label, type="log" if log_x else "linear",
                         range=np.log10(profile_limits).tolist() if log_x else list(profile_limits),
                         autorange=False, row=row, col=2)
        fig.update_yaxes(title_text=label, type="log" if logarithmic else "linear",
                         range=[profile_lo, profile_hi] if logarithmic or name.startswith("comp_") or density_bounds
                         else [float(lower[k]) - pad, float(upper[k]) + pad],
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
                    densityScale=density_scale,
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
            snapshot = read_snapshot(path, step, cgs=reference.cgs,
                                     composition=any(name.startswith("comp_") for name, _, _ in reference.fields),
                                     thermodynamics=any(name == "temperature" for name, _, _ in reference.fields))
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
    parser.add_argument("--map-pixels", type=int, default=241,
                        help="Samples per axis in spherical maps (default: 241; e.g. 801 for sharper maps)")
    density_scale = parser.add_mutually_exclusive_group()
    density_scale.add_argument("--density-min", type=float,
                        help="Density color/profile lower limit in displayed units; smaller values use the bottom color")
    density_scale.add_argument("--density-threshold", type=float,
                        help="Emphasize density near this value in displayed units; requires --log")
    parser.add_argument("--density-max", type=float,
                        help="Density color/profile upper limit in displayed units")
    parser.add_argument("--log-x", "--log-radius", action="store_true",
                        help="Logarithmic horizontal length axis on right-hand profiles; spherical maps stay linear")
    parser.add_argument("--profile-range", "--xlim", type=float, nargs=2, metavar=("MIN", "MAX"),
                        help="Right-hand profile length range (cm with --cgs, otherwise code units)")
    parser.add_argument("--plot", choices=("cartesian", "spherical"), default="cartesian",
                        help="Map style: Cartesian strips/maps (default), or 1D spherical symmetry slices")
    parser.add_argument("--rmax", type=float,
                        help="Spherical map/profile outer radius (cm with --cgs, otherwise code units)")
    parser.add_argument("--radial-grid", type=Path,
                        help="Exact active shell CSV (r_left,r_right in code units); requires --plot spherical")
    parser.add_argument("--cgs", action="store_true",
                        help="Convert plots to CGS using CSV unit metadata; slice coordinates are then in cm")
    parser.add_argument("--composition", action="store_true",
                        help="Add all species mass fractions, labeled from CSV metadata, on a fixed [0,1] scale")
    parser.add_argument("--thermodynamics", action="store_true",
                        help="Add mean molecular weight and temperature (kelvin with --cgs)")
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
    if args.map_pixels < 3:
        parser.error("--map-pixels must be >= 3")
    if args.density_min is not None and (not np.isfinite(args.density_min) or args.density_min <= 0):
        parser.error("--density-min must be finite and positive")
    if args.density_max is not None:
        if not np.isfinite(args.density_max) or args.density_max <= 0:
            parser.error("--density-max must be finite and positive")
        if args.density_threshold is not None:
            parser.error("--density-max cannot be combined with --density-threshold")
    if args.profile_range is not None:
        lo, hi = args.profile_range
        if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
            parser.error("--profile-range requires finite increasing bounds")
        if args.log_x and lo <= 0:
            parser.error("--log-x requires strictly positive --profile-range bounds")
    if args.density_threshold is not None:
        if not args.log or not np.isfinite(args.density_threshold) or args.density_threshold <= 0:
            parser.error("--density-threshold requires --log and a finite positive value")
    if args.rmax is not None:
        if args.plot != "spherical":
            parser.error("--rmax requires --plot spherical")
        if not np.isfinite(args.rmax) or args.rmax <= 0:
            parser.error("--rmax must be finite and positive")
    if args.radial_grid is not None and args.plot != "spherical":
        parser.error("--radial-grid requires --plot spherical")
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
        reference, records, lower, upper = prepare(files, cgs=args.cgs, log=args.log,
                                                  composition=args.composition,
                                                  thermodynamics=args.thermodynamics)
        if args.radial_grid is not None:
            reference.radial_edges = read_radial_grid(args.radial_grid, reference)
        if args.plotly or (args.output and args.output.suffix.lower() == ".html"):
            save_plotly(reference, records, lower, upper, args)
        else:
            animate(reference, records, lower, upper, args)
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f"Animation failed: {error}\n")


if __name__ == "__main__":
    main()
