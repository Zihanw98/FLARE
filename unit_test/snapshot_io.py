"""Shared manifest validation and 3D CSV reading for simulation checks.

These data readers are retained from the removed 3D animation module; rendering
is intentionally separate. Only snapshots listed by a supplied manifest are used.
"""
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
