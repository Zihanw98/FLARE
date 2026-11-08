#!/usr/bin/env python3
"""Build and byte-compare unchanged legacy inputs against the current framework.

Run from any directory. Requires git, make, g++ and Python's standard library.
Only the selected output directory is written; no checkout or existing simulation
output is modified. The default baseline is the last buildable pre-initialization
commit: 78d6452 itself fails because read_table repeated default arguments.
"""
import argparse
import concurrent.futures
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[2]
CASES = ["sod_tube", "Periodic_Acoustic_2D", "CCSN_Simple"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default="ee4aa5b")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "output/compatibility")
    args = parser.parse_args()
    work = args.output.resolve()
    # Keep earlier results intact; choose a new output directory for a repeat.
    work.mkdir(parents=True, exist_ok=False)
    revision = subprocess.check_output(["git", "rev-parse", args.baseline], cwd=ROOT, text=True).strip()
    archive = subprocess.check_output(["git", "archive", revision], cwd=ROOT)
    old = work / "baseline"
    old.mkdir()
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
        stream.extractall(old, filter="data")
    legacy = work / "legacy_inputs"
    shutil.copytree(old / "tests", legacy)

    def build_run(framework, problem):
        with (problem / "build.log").open("w") as log:
            subprocess.run(["make", "-f", str(framework / "Makefile"), f"PROBLEM_DIR={problem}", "-j3"],
                           cwd=framework, stdout=log, stderr=subprocess.STDOUT, check=True)
        with (problem / "run.log").open("w") as log:
            subprocess.run([str(problem / "flare"), "-par", "parameters.txt"], cwd=problem,
                           stdout=log, stderr=subprocess.STDOUT, check=True)

    def compare(case):
        before, after = old / "tests" / case, legacy / case
        build_run(old, before)
        build_run(ROOT, after)
        before_files = sorted(p.relative_to(before / "output") for p in (before / "output").rglob("*") if p.is_file())
        after_files = sorted(p.relative_to(after / "output") for p in (after / "output").rglob("*") if p.is_file())
        assert before_files == after_files, f"{case}: output file lists differ"
        manifest = {}
        for name in before_files:
            a, b = (before / "output" / name).read_bytes(), (after / "output" / name).read_bytes()
            assert a == b, f"{case}: byte mismatch in {name}"
            manifest[str(name)] = {"bytes": len(a), "sha256": hashlib.sha256(a).hexdigest()}
        for name in ["parameters.txt", "user_problem.cpp", "user_config.h"]:
            assert (before / name).read_bytes() == (after / name).read_bytes()
        result = {"case": case, "file_count": len(manifest), "byte_identical": True,
                  "legacy_inputs_unchanged": True, "files": manifest}
        print(f"{case}: {len(manifest)} files byte-identical", flush=True)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(compare, CASES))
    report = {"baseline": revision, "current_head": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "current_worktree_diff": subprocess.check_output(["git", "diff"], cwd=ROOT, text=True),
        "cases": results}
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
