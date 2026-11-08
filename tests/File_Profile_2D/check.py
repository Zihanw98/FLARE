#!/usr/bin/env python3
"""Check file mapping, unit hook, conservation, and the advected cell averages."""
import csv
import json
import math
from pathlib import Path
from generate_inputs import cosine_average

ROOT = Path(__file__).resolve().parent


def snapshot(path):
    lines = path.read_text().splitlines()
    meta = dict(line[2:].split("=", 1) for line in lines if line.startswith("# "))
    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(line for line in lines if not line.startswith("#"))]
    assert all(math.isfinite(v) for row in rows for v in row.values())
    assert all(row["rho"] > 0 and row["p"] > 0 for row in rows)
    return meta, [row for row in rows if row["is_ghost"] == 0]


def main():
    grid = [list(map(float, line.split())) for line in (ROOT / "mesh.dat").read_text().splitlines()[1:]]
    x, y = zip(*grid)
    text = (ROOT / "profile.dat").read_text().splitlines()
    profile = [dict(zip(text[0].split(), map(float, line.split()))) for line in text[1:]]
    paths = sorted((ROOT / "output").glob("file_profile_*.csv"))
    meta0, initial = snapshot(paths[0])
    assert [float(meta0[f"unit_{key}"]) for key in ["mass", "length", "time", "temperature"]] == [2, 3, 5, 7]
    assert len(initial) == len(profile) == 1024
    volumes = []
    for index, (row, expected) in enumerate(zip(initial, profile)):
        i, j = index % 32, index // 32
        volumes.append((x[i + 1] - x[i]) * (y[j + 1] - y[j]))
        assert math.isclose(row["x"], (x[i] + x[i + 1]) / 2, abs_tol=1e-14)
        assert math.isclose(row["y"], (y[j] + y[j + 1]) / 2, abs_tol=1e-14)
        for key in ["rho", "vx", "vy", "vz", "gamma"]:
            assert math.isclose(row[key], expected[key], rel_tol=1e-14, abs_tol=1e-14)
        assert math.isclose(row["p"], expected["pre"], rel_tol=1e-14)
    def totals(rows):
        return [sum(volume * value for volume, value in zip(volumes, component)) for component in [
            [r["rho"] for r in rows],
            [r["rho"] * r["vx"] for r in rows],
            [r["rho"] * r["vy"] for r in rows],
            [r["p"] / (r["gamma"] - 1) + 0.5 * r["rho"] * (r["vx"]**2 + r["vy"]**2) for r in rows],
        ]]
    before = totals(initial)
    for path in paths:
        _, rows = snapshot(path)
        assert len(rows) == 1024
        assert all(math.isclose(a, b, rel_tol=2e-11, abs_tol=2e-12) for a, b in zip(before, totals(rows)))
    meta, final = snapshot(paths[-1])
    time = float(meta["time"])
    assert time == 0.25
    error = 0.0
    for index, row in enumerate(final):
        i, j = index % 32, index // 32
        exact = 1 + 0.2 * cosine_average(x[i], x[i + 1], 1, 0.2 * time) * cosine_average(y[j], y[j + 1], 0.75, -0.1 * time)
        error += volumes[index] * abs(row["rho"] - exact)
    error /= sum(volumes)
    assert error < 0.003, error
    report = {"snapshots": len(paths), "steps": int(meta["step"]), "time": time,
              "density_L1": error, "initial_mass": before[0], "final_mass": totals(final)[0],
              "initial_mapping_and_units": "PASS", "conservation": "PASS"}
    (ROOT / "output/check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
