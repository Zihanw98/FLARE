#!/usr/bin/env python3
"""Run/check the cubic 3D blast: positivity, conservation, symmetry, evolution."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from animate3d import edges, read_snapshot, snapshot_files


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=Path(__file__).with_name("parameters.txt"))
    parser.add_argument("--executable", type=Path, default=Path(__file__).with_name("flare"))
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    parameters = args.parameters.resolve()
    cfg = {}
    for line in parameters.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    output = parameters.parent / cfg["output_directory"]
    log = output / "run.log"
    if args.run:
        output.mkdir(parents=True, exist_ok=True)
        with log.open("w") as stream:
            result = subprocess.run([str(args.executable.resolve()), "-par", str(parameters)],
                                    cwd=parameters.parent, stdout=stream, stderr=subprocess.STDOUT)
        require(result.returncode == 0, f"Simulation failed; inspect {log}")
    files = snapshot_files(output, "mesh_" + cfg["output_suffix"], log)
    initial = read_snapshot(files[0])
    require(initial.step == 0 and initial.time == 0, "Missing initial snapshot")
    shape = tuple(int(cfg[f"n{axis}a"]) for axis in "xyz")
    require(len(set(shape)) == 1, "Symmetry checker requires equal x/y/z resolution")
    lengths = [float(cfg[f"{axis}_upper"]) - float(cfg[f"{axis}_lower"]) for axis in "xyz"]
    require(np.allclose(lengths, lengths[0], rtol=0, atol=1e-13), "Checker requires a cubic domain")
    for d, axis in enumerate("xyz"):
        lower = float(cfg[f"{axis}_lower"])
        centers = lower + (np.arange(shape[d]) + .5) * lengths[d] / shape[d]
        require(initial.axes[d].shape == centers.shape and
                np.allclose(initial.axes[d], centers, rtol=0, atol=1e-13), "Incorrect grid")
    positions = np.meshgrid(*initial.axes, indexing="ij")
    offsets = [positions[d] - .5 * (float(cfg[f"{axis}_lower"]) + float(cfg[f"{axis}_upper"]))
               for d, axis in enumerate("xyz")]
    radius2 = sum(offset ** 2 for offset in offsets)
    expected_pressure = np.where(radius2 < (.12 * min(lengths)) ** 2, 5., 1.)
    require(np.array_equal(initial.fields["p"], expected_pressure), "Incorrect initial high-pressure sphere")
    require(np.all(initial.fields["rho"] == 1), "Initial density must be one")
    velocity_keys = ("vx", "vy", "vz")
    require(all(np.all(initial.fields[key] == 0) for key in velocity_keys), "Initial gas must be at rest")
    volumes = (np.diff(edges(initial.axes[0]))[:, None, None]
               * np.diff(edges(initial.axes[1]))[None, :, None]
               * np.diff(edges(initial.axes[2]))[None, None, :])
    reference = None
    records = []
    for path in files:
        snap = read_snapshot(path)
        require(all(np.array_equal(a, b) for a, b in zip(initial.axes, snap.axes)), "Grid changed")
        fields = snap.fields
        rho, pressure = fields["rho"], fields["p"]
        gamma = fields["gamma"]
        require(np.all(gamma == float(cfg["gamma"])), "Unexpected gamma")
        require(np.all(rho > 0) and np.all(pressure > 0), "Nonpositive density/pressure")
        velocity2 = sum(fields[key] ** 2 for key in velocity_keys)
        mass = np.sum(rho * volumes)
        energy = np.sum((pressure / (gamma - 1) + .5 * rho * velocity2) * volumes)
        if reference is None:
            reference = np.array([mass, energy])
        drift = np.abs(np.array([mass, energy]) / reference - 1)
        require(np.max(drift) < 1e-10, f"Mass/energy conservation failed at step {snap.step}: {drift}")
        symmetry = 0.
        for key in ("rho", "p", *velocity_keys):
            field = fields[key]
            scale = max(1., np.max(np.abs(field)))
            for axis in range(3):
                sign = -1 if key == velocity_keys[axis] else 1
                symmetry = max(symmetry, np.max(np.abs(field - sign * np.flip(field, axis))) / scale)
            for a, b in ((0, 1), (0, 2), (1, 2)):
                other = velocity_keys[b] if key == velocity_keys[a] else (
                    velocity_keys[a] if key == velocity_keys[b] else key)
                symmetry = max(symmetry, np.max(np.abs(field - fields[other].swapaxes(a, b))) / scale)
        require(symmetry < 1e-9, f"Reflection/axis permutation symmetry failed at step {snap.step}")
        records.append(dict(step=snap.step, time=snap.time, mass=mass, energy=energy,
                            mass_relative_drift=drift[0], energy_relative_drift=drift[1],
                            symmetry_error=symmetry, rho_min=rho.min(), rho_max=rho.max(),
                            p_min=pressure.min(), p_max=pressure.max(),
                            speed_max=np.sqrt(velocity2.max())))
    require(np.isclose(snap.time, float(cfg["t_end"]), rtol=0, atol=1e-12), "Run did not reach t_end")
    # Prevent a no-op solver from passing conservation and symmetry alone.
    radial_flux = sum(fields[key] * offsets[d] for d, key in enumerate(velocity_keys))
    require(np.sum(rho * radial_flux * volumes) > 1e-6, "No net outward expansion")
    require(rho.max() > 1.05 and pressure.max() < 5, "Missing compression/pressure relaxation")
    diagnostics = output / "diagnostics"
    diagnostics.mkdir(exist_ok=True)
    with (diagnostics / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    print(f"PASS: {shape[0]}^3 cells, {len(records)} snapshots, final t={snap.time:g}, step={snap.step}")
    for key in ("mass_relative_drift", "energy_relative_drift", "symmetry_error"):
        print(f"  max {key} = {max(record[key] for record in records):.6e}")
    print(f"  final rho=[{rho.min():.6g}, {rho.max():.6g}], p=[{pressure.min():.6g}, {pressure.max():.6g}]")
    print(f"Saved {diagnostics / 'metrics.csv'}")


if __name__ == "__main__":
    main()
