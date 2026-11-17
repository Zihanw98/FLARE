#!/usr/bin/env python3
"""Run/check all saved H/He mixing states, including species conservation and motion."""
import argparse
import csv
import json
from pathlib import Path
import subprocess

import numpy as np

from mixing import (check_grid, check_initial, diagnostics, load_frame, parameters,
                    require, run_paths)


def analyze(parameter_file):
    cfg, output, paths = run_paths(parameter_file)
    reference, initial = load_frame(paths[0])
    require(reference.step == 0 and reference.time == 0, "Missing initial state")
    check_initial(reference, initial, cfg)
    records = []
    species_metadata = {k: v for k, v in reference.metadata.items() if k.startswith("comp_")}
    for path in paths:
        frame, fields = load_frame(path)
        check_grid(frame, reference)
        require({k: v for k, v in frame.metadata.items() if k.startswith("comp_")} == species_metadata,
                "Species metadata changed between snapshots")
        for axis in (frame.x, frame.y):
            require(np.allclose(np.diff(axis), axis[1]-axis[0], rtol=1e-10, atol=1e-13),
                    "Requires uniform grid")
        record = diagnostics(frame, fields)
        if records:
            for key in ("mass", "hydrogen_mass", "helium_mass", "energy"):
                drift = abs(record[key]/records[0][key]-1)
                require(drift < 2e-9, f"{key} conservation failed at step {frame.step}: {drift}")
        records.append(record)
    require(len(records) >= 3 and np.isclose(frame.time, float(cfg["t_end"]), rtol=0, atol=1e-12),
            "Incomplete run: did not reach t_end")
    motion = float(np.mean(abs(fields["H"]-initial["H"])))
    require(motion > 1e-6, "Composition did not advect")
    growth = max(r["transverse_kinetic_energy"] for r in records)/records[0]["transverse_kinetic_energy"]
    if float(cfg["t_end"]) >= 2:
        require(growth > 2, "Kelvin-Helmholtz perturbation did not grow")
        require(records[-1]["mass_weighted_mixing"] > records[0]["mass_weighted_mixing"],
                "No increase in cell-scale mixing")
    summary = dict(status="PASS", snapshots=len(records), final_step=frame.step, final_time=frame.time,
                   grid=[len(frame.x), len(frame.y)], composition_mean_absolute_change=motion,
                   peak_transverse_energy_growth=growth,
                   max_fraction_sum_error=max(r["fraction_sum_error"] for r in records),
                   fraction_min=min(r["fraction_min"] for r in records),
                   fraction_max=max(r["fraction_max"] for r in records),
                   initial_mixing=records[0]["mass_weighted_mixing"],
                   final_mixing=records[-1]["mass_weighted_mixing"],
                   max_relative_drift={key: max(abs(r[key]/records[0][key]-1) for r in records)
                                       for key in ("mass", "hydrogen_mass", "helium_mass", "energy")},
                   temperature_model="derived neutral monatomic ideal gas; mu=1/(X_H+X_He/4)")
    destination = output / "diagnostics"
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (destination / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def run(parameter_file, executable):
    cfg = parameters(parameter_file)
    output = (parameter_file.parent / cfg["output_directory"]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    log = output / "run.log"
    print(f"Running {executable}; log: {log}", flush=True)
    with log.open("w") as stream:
        result = subprocess.run([str(executable.resolve()), "-par", str(parameter_file)],
                                cwd=parameter_file.parent, stdout=stream, stderr=subprocess.STDOUT)
    require(result.returncode == 0, f"Simulation failed ({result.returncode}); see {log}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    case = Path(__file__).resolve().parent
    parser.add_argument("--parameters", type=Path, default=case / "parameters.txt")
    parser.add_argument("--executable", type=Path, default=case / "flare")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--smoke", action="store_true", help="32x32, t=0.15 in output/smoke")
    args = parser.parse_args()
    parameter_file = args.parameters.resolve()
    try:
        if args.smoke:
            require(args.run, "--smoke requires --run")
            cfg = parameters(parameter_file)
            output = (parameter_file.parent/cfg["output_directory"] / "smoke").resolve()
            cfg.update(nxa="32", nya="32", t_end="0.15", step_every="5", output_directory=str(output),
                       composition_activate_file=str((parameter_file.parent/cfg["composition_activate_file"]).resolve()))
            output.mkdir(parents=True, exist_ok=True)
            parameter_file = output / "parameters.txt"
            parameter_file.write_text("".join(f"{key} = {value}\n" for key, value in cfg.items()))
        if args.run:
            run(parameter_file, args.executable)
        analyze(parameter_file)
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f"H/He check failed: {error}\n")


if __name__ == "__main__":
    main()
