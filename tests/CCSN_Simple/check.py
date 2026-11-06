"""Check the mass cut, instantaneous thermal injection, and saved evolution."""
import argparse
import csv
import json
import math
from pathlib import Path


def rows(path):
    with path.open() as stream:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(line for line in stream if not line.startswith("#"))
        ]


def check(output):
    meta = json.loads((output / "model.json").read_text())
    profile = rows(output / "initial_envelope.csv")
    snapshots = sorted(output.glob("ccsn_*.csv"))
    assert snapshots and snapshots[0].name == "ccsn_00000000.csv"
    initial = [row for row in rows(snapshots[0]) if row["is_ghost"] == 0]
    assert len(initial) == len(profile)
    core = meta["core_mass_code"]
    assert math.isclose(core, 4.0, abs_tol=1e-12)
    assert math.isclose(meta["enclosed_mass_at_cut_code"], core, abs_tol=1e-12)
    assert math.isclose(profile[0]["r_left"], meta["cut_radius_code"], abs_tol=1e-12)
    initial_mass = 0.0
    injection = 0.0
    for index, (row, before) in enumerate(zip(initial, profile)):
        a, b = before["r_left"], before["r_right"]
        volume = 4 * math.pi / 3 * (b - a) * (a*a + a*b + b*b)
        assert math.isclose(volume, before["volume"], rel_tol=1e-12)
        assert math.isclose(row["rho"], before["rho"], rel_tol=1e-14)
        assert row["vx"] == 0.0
        assert math.isclose(row["x"], (a+b)/2, rel_tol=1e-14)
        initial_mass += row["rho"] * volume
        added = (row["p"] / (row["gamma"]-1) - before["energy_density_before"]) * volume
        if index == 0:
            injection = added
        else:
            assert math.isclose(row["p"], before["p"], rel_tol=1e-13)
    # Includes a tiny atmosphere, hence the tolerance on the six-solar-mass envelope.
    assert abs(initial_mass - 6.0) < 1e-5, initial_mass
    injected_erg = injection * meta["energy_unit_erg"]
    assert math.isclose(injected_erg, 1.6e51, rel_tol=1e-12), injected_erg
    maximum_speed = 0.0
    maximum_fast_mass = 0.0
    for path in snapshots:
        saved = rows(path)
        assert len(saved) >= len(profile)
        for row in saved:
            assert all(math.isfinite(value) for value in row.values()), path
            assert row["rho"] > 0 and row["p"] > 0, path
        active = [row for row in saved if row["is_ghost"] == 0]
        assert len(active) == len(profile)
        maximum_speed = max(maximum_speed, max(row["vx"] for row in active))
        fast_mass = sum(row["rho"] * cell["volume"] for row, cell in zip(active, profile)
                        if abs(row["vx"]) * meta["velocity_unit_cm_s"] > 0.1 * 29979245800.0)
        maximum_fast_mass = max(maximum_fast_mass, fast_mass)
    final = [row for row in rows(snapshots[-1]) if row["is_ghost"] == 0]
    final_mass = sum(row["rho"] * cell["volume"] for row, cell in zip(final, profile))
    last_time = next(float(line.split("=", 1)[1]) for line in snapshots[-1].read_text().splitlines()
                     if line.startswith("# time="))
    report = {
        "cut_radius_solar": meta["cut_radius_code"],
        "core_mass_solar": core,
        "initial_envelope_mass_solar": initial_mass,
        "injected_energy_erg": injected_erg,
        "injection_mass_solar": meta["injection_mass_code"],
        "snapshots_checked": len(snapshots),
        "last_time_code": last_time,
        "final_domain_mass_solar": final_mass,
        "maximum_saved_outward_speed_km_s": maximum_speed * meta["velocity_unit_cm_s"] / 1e5,
        "maximum_saved_outward_speed_over_c": maximum_speed * meta["velocity_unit_cm_s"] / 29979245800.0,
        "maximum_saved_mass_above_0_1c_solar": maximum_fast_mass,
        "high_speed_limit_note": "Newtonian ideal-gas calculation; high-speed outer layers require relativistic treatment for quantitative predictions.",
    }
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=Path(__file__).parent / "output")
    parser.add_argument("--expect-time", type=float, default=100.0)
    args = parser.parse_args()
    report = check(args.output)
    assert math.isclose(report["last_time_code"], args.expect_time, rel_tol=1e-12, abs_tol=1e-14), report
    if args.expect_time > 0:
        assert report["maximum_saved_outward_speed_km_s"] > 0, report
    (args.output / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
