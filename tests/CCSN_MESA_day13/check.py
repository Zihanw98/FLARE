#!/usr/bin/env python3
"""Verify the native MESA mesh, cgs conversion, thermal bomb, and completed run."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_csv(path):
    lines = path.read_text().splitlines()
    meta = dict(line[2:].split("=", 1) for line in lines if line.startswith("# "))
    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(line for line in lines if not line.startswith("#"))]
    return meta, rows


def parameters(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=ROOT / "parameters.txt")
    args = parser.parse_args()
    cfg = parameters(args.parameters)
    output = ROOT / cfg["output_directory"]
    info = json.loads((output / "model.json").read_text())
    source = (ROOT / cfg["init_mesa_profile_name"]).resolve()
    lines = source.read_text().splitlines()
    header = dict(zip(lines[1].split(), lines[2].split()))
    columns = lines[5].split()
    native = [dict(zip(columns, map(float, line.split()))) for line in lines[6:] if line.strip()][::-1]
    # Independently reconstruct only native outer faces. The selected inner cut
    # lies on a native face, so there is no interpolation or remapping in this case.
    length = info["length_unit_cm"]
    radii = [(10 ** row["logR"] * 6.96e10) / length for row in native]
    cut = float(cfg["mesa_radius_inner"])
    cut_index = min(range(len(radii)), key=lambda i: abs(radii[i] - cut))
    assert math.isclose(radii[cut_index], cut, rel_tol=2e-14)
    retained = native[cut_index + 1:]
    assert len(retained) == info["native_cells"] == 2599
    assert info["atmosphere_cells"] == 100
    assert math.isclose(native[cut_index]["mass"], info["core_mass_code"], rel_tol=1e-14)
    assert abs(info["core_mass_code"] - 4) < 0.002
    _, before = read_csv(output / "initial_envelope.csv")
    assert len(before) == len(retained) + 100
    assert before[0]["r_left"] == cut
    assert before[-1]["r_right"] == float(cfg["mesa_radius_outer"])
    for i, row in enumerate(before):
        a, b = row["r_left"], row["r_right"]
        assert b > a
        if i:
            assert a == before[i - 1]["r_right"]
        expected_volume = 4 * math.pi / 3 * (b - a) * (a * a + a * b + b * b)
        assert math.isclose(row["volume"], expected_volume, rel_tol=2e-14)
        if i < len(retained):
            original = retained[i]
            assert math.isclose(b, radii[cut_index + 1 + i], rel_tol=2e-14)
            assert math.isclose(row["rho"] * info["density_unit_g_cm3"], 10 ** original["logRho"], rel_tol=2e-13)
            assert math.isclose(row["p"] * info["pressure_unit_dyn_cm2"], 10 ** original["logP"], rel_tol=2e-13)
            assert math.isclose(row["vx"] * info["velocity_unit_cm_s"], original["velocity"], rel_tol=2e-13, abs_tol=1e-10)
            assert row["is_atmosphere"] == 0
        else:
            assert row["is_atmosphere"] == 1
            assert math.isclose(row["rho"] * info["density_unit_g_cm3"], 1e-15, rel_tol=2e-14)
    paths = sorted(output.glob("ccsn_mesa_*.csv"))
    _, initial_all = read_csv(paths[0])
    initial = [r for r in initial_all if r["is_ghost"] == 0]
    assert len(initial) == len(before)
    added_energy = 0.0
    for i, (row, old) in enumerate(zip(initial, before)):
        assert row["rho"] == old["rho"]
        assert row["vx"] == old["vx"]
        energy = row["p"] / (row["gamma"] - 1) + 0.5 * row["rho"] * row["vx"]**2
        difference = energy - old["energy_density_before"]
        if i == 0:
            added_energy = difference * old["volume"] * info["energy_unit_erg"]
        else:
            assert math.isclose(energy, old["energy_density_before"], rel_tol=2e-13)
    assert math.isclose(added_energy, 1.6e51, rel_tol=2e-13)
    initial_mass = sum(r["rho"] * old["volume"] for r, old in zip(initial, before))
    maximum_mass_error = 0.0
    maximum_speed = 0.0
    for path in paths:
        meta, all_rows = read_csv(path)
        assert all(math.isfinite(v) for row in all_rows for v in row.values()), path
        assert all(row["rho"] > 0 and row["p"] > 0 for row in all_rows), path
        rows = [r for r in all_rows if r["is_ghost"] == 0]
        assert len(rows) == len(before)
        for row, old in zip(rows, before):
            assert math.isclose(row["x"], (old["r_left"] + old["r_right"]) / 2, rel_tol=2e-14)
        mass = sum(r["rho"] * old["volume"] for r, old in zip(rows, before))
        maximum_mass_error = max(maximum_mass_error, abs(mass / initial_mass - 1))
        maximum_speed = max(maximum_speed, max(abs(r["vx"]) for r in rows) * info["velocity_unit_cm_s"])
    assert maximum_mass_error < 1e-9, maximum_mass_error
    assert math.isclose(float(meta["time"]), float(cfg["t_end"]), rel_tol=2e-14)
    assert max(r["vx"] for r in rows[:30]) > 0.1  # Outgoing motion from the inner thermal bomb.
    report = {
        "profile": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "profile_rows": len(native), "header_num_zones": int(header["num_zones"]),
        "profile_header_time_seconds": float(header["time_seconds"]),
        "core_mass_msun": info["core_mass_code"], "cut_radius_cm": cut * length,
        "retained_native_shells": len(retained), "atmosphere_shells": 100,
        "outer_radius_cm": before[-1]["r_right"] * length,
        "native_envelope_mass_code": info["envelope_mass_code"],
        "mesa_mass_column_envelope_msun": native[-1]["mass"] - native[cut_index]["mass"],
        "injected_energy_erg": added_energy, "grid_and_profile_mapping": "PASS",
        "last_step": int(meta["step"]), "elapsed_physical_time_s": float(meta["time"]) * info["time_unit_s"],
        "snapshots_checked": len(paths), "maximum_relative_mass_change": maximum_mass_error,
        "maximum_saved_speed_cm_s": maximum_speed, "maximum_saved_speed_over_c": maximum_speed / 2.99792458e10,
        "scope": "Early thermal-bomb evolution on the unchanged native MESA grid; not shock breakout.",
    }
    (output / "check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
