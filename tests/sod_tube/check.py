#!/usr/bin/env python3
"""Run the standard Sod case and check CSV output against its exact solution.

Uses only Python's standard library. The reference is evaluated independently
of FLARE's C++ Riemann solver using shock/rarefaction relations; see
https://www.clawpack.org/riemann_book/html/Euler.html .
This check targets the supplied [0,1], gamma=1.4, t_end=0.2 configuration.
"""

import argparse
import csv
import math
from pathlib import Path
import re
import subprocess
from snapshot_io import read_snapshot_csv


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_parameters(path):
    parameters = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            parameters[key.strip()] = value.strip()
    return parameters


def exact_sod(x, time):
    """Return rho, u, p at cell centers (not exact cell averages)."""
    gamma = 1.4
    c_left = math.sqrt(gamma)
    c_right = math.sqrt(gamma * 0.1 / 0.125)

    def left_wave(p):
        return 2.0 * c_left / (gamma - 1.0) * (p ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)

    def right_wave(p):
        return (p - 0.1) * math.sqrt(
            (2.0 / ((gamma + 1.0) * 0.125))
            / (p + (gamma - 1.0) / (gamma + 1.0) * 0.1)
        )

    lower, upper = 0.1, 1.0
    for _ in range(80):
        mid = (lower + upper) / 2.0
        if left_wave(mid) + right_wave(mid) > 0.0:
            upper = mid
        else:
            lower = mid
    p_star = (lower + upper) / 2.0
    u_star = (right_wave(p_star) - left_wave(p_star)) / 2.0
    c_star = c_left * p_star ** ((gamma - 1.0) / (2.0 * gamma))
    shock_speed = c_right * math.sqrt(
        (gamma + 1.0) / (2.0 * gamma) * (p_star / 0.1)
        + (gamma - 1.0) / (2.0 * gamma)
    )
    b = (gamma - 1.0) / (gamma + 1.0)
    rho_star_right = 0.125 * ((p_star / 0.1) + b) / (b * (p_star / 0.1) + 1.0)
    result = []
    for position in x:
        speed = (position - 0.5) / time
        if speed <= -c_left:
            state = (1.0, 0.0, 1.0)
        elif speed < u_star - c_star:
            u = 2.0 / (gamma + 1.0) * (c_left + speed)
            c = 2.0 / (gamma + 1.0) * (c_left - (gamma - 1.0) / 2.0 * speed)
            state = ((c / c_left) ** (2.0 / (gamma - 1.0)), u,
                     (c / c_left) ** (2.0 * gamma / (gamma - 1.0)))
        elif speed <= u_star:
            state = (p_star ** (1.0 / gamma), u_star, p_star)
        elif speed < shock_speed:
            state = (rho_star_right, u_star, p_star)
        else:
            state = (0.125, 0.0, 0.1)
        result.append(state)
    return result


def read_snapshot(path, count, ghosts, step, time, dt):
    metadata, _, raw_rows = read_snapshot_csv(path)
    rows = [{key: float(value) for key, value in row.items()} for row in raw_rows]
    require(len(rows) == count + 2 * ghosts, "CSV cell count is incorrect")
    require({"step", "time", "dt"}.issubset(metadata), "CSV snapshot header missing")
    require(metadata["step"] == step, "CSV step does not match the snapshot")
    # Console diagnostics have ten significant digits; CSV preserves full precision.
    require(math.isclose(metadata["time"], time, rel_tol=1e-9, abs_tol=1e-14), "CSV time disagrees with log")
    require(math.isclose(metadata["dt"], dt, rel_tol=1e-9, abs_tol=1e-14), "CSV dt disagrees with log")
    for row in rows:
        require(all(math.isfinite(value) for value in row.values()), "CSV contains NaN/inf")
        require(row["rho"] > 0 and row["p"] > 0, "Density/pressure lost positivity")
        require(row["gamma"] == 1.4, "Unexpected gamma")
        require(all(row[key] == 0 for key in ("y", "z", "vy", "vz")), "Nonzero unused direction")
    active = [row for row in rows if row["is_ghost"] == 0]
    require(len(active) == count, "CSV active cell count is incorrect")
    for i, row in enumerate(active):
        require(abs(row["x"] - (i + 0.5) / count) < 1e-12, "Incorrect cell center")
    for row in rows[:ghosts]:
        require(row["is_ghost"] == 1, "Lower ghost flag missing")
        require(all(row[key] == active[0][key] for key in ("rho", "vx", "p", "gamma")),
                "Lower outflow ghost is stale")
    for row in rows[-ghosts:]:
        require(row["is_ghost"] == 1, "Upper ghost flag missing")
        require(all(row[key] == active[-1][key] for key in ("rho", "vx", "p", "gamma")),
                "Upper outflow ghost is stale")
    return active


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=Path(__file__).parent / "flare")
    parser.add_argument("--parameters", type=Path, default=Path(__file__).parent / "parameters.txt")
    args = parser.parse_args()
    parameter_file = args.parameters.resolve()
    problem_dir = parameter_file.parent
    cfg = read_parameters(parameter_file)
    count, ghosts = int(cfg["nxa"]), int(cfg["nxg"])
    end_time = float(cfg["t_end"])
    require(float(cfg["gamma"]) == 1.4 and end_time == 0.2
            and float(cfg["x_lower"]) == 0.0 and float(cfg["x_upper"]) == 1.0
            and count >= 200 and ghosts >= 1,
            "make test expects the standard Sod domain, gamma=1.4, t_end=0.2 and >=200 cells")
    run = subprocess.run([str(args.executable.resolve()), "-par", str(parameter_file)],
                         cwd=problem_dir, capture_output=True, text=True)
    (problem_dir / "run.log").write_text(run.stdout + run.stderr)
    require(run.returncode == 0, "Simulation failed; see " + str(problem_dir / "run.log"))

    records = []
    snapshots = {}
    for line in run.stdout.splitlines():
        match = re.match(r"step = (\d+), time = ([^,]+), dt = ([^,]+)", line)
        if not match:
            continue
        step, time, dt = int(match[1]), float(match[2]), float(match[3])
        records.append((step, time, dt))
        if ", csv = " in line:
            snapshots[step] = problem_dir / line.split(", csv = ", 1)[1]
    require(records and records[0] == (0, 0.0, 0.0), "Initial diagnostics missing")
    require(records[-1][1] == end_time, "Simulation did not reach t_end")
    require([record[0] for record in records] == list(range(len(records))), "Missing step diagnostics")
    require(all(record[2] > 0 for record in records[1:]), "Nonpositive time step")
    last_step = records[-1][0]
    interval = int(cfg["output_step_every"])
    expected_steps = {0, 1, last_step} | set(range(interval, last_step + 1, interval))
    require(set(snapshots) == expected_steps, "Snapshot cadence is incorrect")
    states = {step: read_snapshot(path, count, ghosts, *records[step])
              for step, path in snapshots.items()}
    for row in states[0]:
        expected = (1.0, 1.0) if row["x"] < 0.5 else (0.125, 0.1)
        require(row["rho"] == expected[0] and abs(row["p"] - expected[1]) < 1e-14
                and row["vx"] == 0.0, "Incorrect initial state")

    final = states[last_step]
    reference = exact_sod([row["x"] for row in final], end_time)
    errors = {key: sum(abs(row[key] - exact[k]) for row, exact in zip(final, reference)) / count
              for k, key in enumerate(("rho", "vx", "p"))}
    # First-order Godunov has expected diffusion near the contact and shock.
    for key, limit in (("rho", 0.02), ("vx", 0.04), ("p", 0.015)):
        require(errors[key] < limit, f"Sod {key} L1 error {errors[key]} exceeds {limit}")
    dx = 1.0 / count
    mass = sum(row["rho"] for row in final) * dx
    momentum = sum(row["rho"] * row["vx"] for row in final) * dx
    energy = sum(row["p"] / 0.4 + 0.5 * row["rho"] * row["vx"] ** 2 for row in final) * dx
    require(abs(mass - 0.5625) < 1e-9, "Mass balance failed")
    require(abs(energy - 1.375) < 1e-9, "Energy balance failed")
    # Fixed end pressures exert a net force; total momentum is not constant.
    require(abs(momentum - 0.9 * end_time) < 1e-9, "Momentum balance failed")

    reference_path = snapshots[last_step].parent / "sod_reference.csv"
    with reference_path.open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("x", "rho", "vx", "p"))
        writer.writerows((row["x"], *exact) for row, exact in zip(final, reference))
    print(f"Sod passed: {count} active cells, {last_step} steps, t={end_time}, {len(snapshots)} snapshots")
    print("L1 errors: " + ", ".join(f"{key}={value:.6g}" for key, value in errors.items()))
    print(f"Mass={mass:.12g}, momentum={momentum:.12g}, energy={energy:.12g}")
    print(f"Final CSV: {snapshots[last_step]}")
    print(f"Exact reference: {reference_path}")


if __name__ == "__main__":
    main()
