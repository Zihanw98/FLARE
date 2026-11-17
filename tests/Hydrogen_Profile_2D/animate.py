#!/usr/bin/env python3
"""Plot this case's fluid fields and species; extra options go to the root animate.py."""
import argparse
from pathlib import Path
import subprocess
import sys

CASE_DIR = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("--parameters", type=Path, default=CASE_DIR / "parameters.txt")
    args, options = parser.parse_known_args()
    cfg = {}
    for line in args.parameters.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    command = [sys.executable, str(CASE_DIR.parents[1] / "animate.py"),
               str(CASE_DIR / cfg.get("output_directory", "output")),
               cfg.get("output_suffix", "hydrogen_profile"), "--composition", *options]
    return subprocess.run(command).returncode


if __name__ == "__main__":
    sys.exit(main())
