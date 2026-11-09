#!/usr/bin/env python3
"""Supervise this isolated run and persist progress without retaining every step log."""
import csv
import datetime
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parent
META = json.loads((ROOT / "manifest.json").read_text())
UNIT = META["time_unit_seconds"]
TARGET = META["target_physical_seconds"]
START = time.monotonic()
status = {"status": "starting", "supervisor_pid": os.getpid(),
          "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "target_physical_seconds": TARGET, "time_unit_seconds": UNIT,
          "integrator": "lts_rk2", "max_levels": 8, "cfl": 0.2,
          "physical_time_seconds": 0.0, "step": 0}

def save():
    status["updated_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    status["elapsed_wall_seconds"] = time.monotonic() - START
    tmp = ROOT / "progress.json.tmp"
    tmp.write_text(json.dumps(status, indent=2) + "\n")
    tmp.replace(ROOT / "progress.json")

def verify_final():
    paths = list((ROOT / "output").glob("ccsn_mesa_*.csv"))
    path = max(paths, key=lambda p: int(p.stem.rsplit("_", 1)[1]))
    lines = path.read_text().splitlines()
    meta = dict(line[2:].split("=", 1) for line in lines if line.startswith("# "))
    rows = list(csv.DictReader(line for line in lines if not line.startswith("#")))
    active = [row for row in rows if int(row["is_ghost"]) == 0]
    finite = all(math.isfinite(float(v)) for row in rows for v in row.values())
    positive = all(float(row["rho"]) > 0 and float(row["p"]) > 0 for row in rows)
    physical_time = float(meta["time"]) * UNIT
    reached = math.isclose(physical_time, TARGET, rel_tol=2e-14)
    report = {"last_snapshot": str(path), "last_step": int(meta["step"]),
              "physical_time_seconds": physical_time, "target_reached": reached,
              "all_states_finite": finite, "density_and_pressure_positive": positive,
              "active_cells": len(active)}
    (ROOT / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    status["physical_time_seconds"] = physical_time
    return reached and finite and positive and len(active) == META["expected_active_cells"]

save()
pattern = re.compile(r"^step = (\d+), time = ([^,]+), dt = ([^,]+),")
last_save = last_logged = START
try:
    with (ROOT / "run.log").open("a", buffering=1) as log:
        process = subprocess.Popen(
            ["stdbuf", "-oL", "-eL", str(ROOT / "flare"), "-par", str(ROOT / "parameters.txt"),
             "--performance-file", str(ROOT / "performance.json")],
            cwd=ROOT, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1)
        status.update(status="running", solver_pid=process.pid)
        save()
        for line in process.stdout:
            now = time.monotonic()
            match = pattern.match(line)
            if match:
                step = int(match[1])
                physical = float(match[2]) * UNIT
                status.update(step=step, physical_time_seconds=physical,
                              last_macro_step_seconds=float(match[3]) * UNIT,
                              last_diagnostic=line.strip())
                elapsed = now - START
                if physical > 0 and elapsed > 0:
                    status["mean_physical_seconds_per_wall_second"] = physical / elapsed
                    status["estimated_remaining_wall_seconds_at_mean_rate"] = max(0, TARGET-physical) * elapsed / physical
                if step <= 5 or step % 1000 == 0 or "csv =" in line or now-last_logged >= 60:
                    log.write(line)
                    last_logged = now
            else:
                log.write(line)
            if now-last_save >= 10:
                save()
                last_save = now
        status["exit_code"] = process.wait()
        status["status"] = "failed"
        if status["exit_code"] == 0:
            status["status"] = "completed" if verify_final() else "verification_failed"
        save()
except Exception as error:
    status.update(status="supervisor_error", error=repr(error))
    save()
    raise
