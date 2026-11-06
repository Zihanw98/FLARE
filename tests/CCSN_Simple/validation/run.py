#!/usr/bin/env python3
"""Isolated numerical validation of the current CCSN_Simple experiment.

Run from anywhere: python3 tests/CCSN_Simple/validation/run.py
Creates a new output_validation/<timestamp> directory; never overwrites output/.
Requires numpy and matplotlib. The production source and parameters are preserved.
"""
import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"Source changed; inspect validation adapter: {old!r}")
    return source.replace(old, new)


def adapted_source(source):
    source = replace_once(source, "constexpr double atmosphere_density =", "double atmosphere_density =")
    source = replace_once(source, "constexpr double injection_energy_erg = 1.6e51;",
                          "double injection_energy_erg = 1.6e51;\ndouble validation_heating_mass = 0;")
    source = replace_once(source, """    mesh.cells[injection_flat].U[hd_core::ENERGY]+=
        injected_energy/mesh.volume_cell[injection_flat];""", """    const double first_mass=mesh.cells[injection_flat].U[hd_core::RHO]*mesh.volume_cell[injection_flat];
    const double heating_mass=validation_heating_mass>0 ? validation_heating_mass : first_mass;
    double remaining=heating_mass, heating_outer=0;
    int heating_cells=0;
    for(int i=0; i<mesh.nactive[0] && remaining>1e-14*heating_mass; ++i)
    {
        const auto flat=mesh.flat_index({i},true);
        const double dm=mesh.cells[flat].U[hd_core::RHO]*mesh.volume_cell[flat];
        const double weight=std::min(dm,remaining);
        mesh.cells[flat].U[hd_core::ENERGY]+=injected_energy*(weight/heating_mass)/mesh.volume_cell[flat];
        remaining-=weight;
        heating_outer=mesh.cell_face_pos(0,i+mesh.nghost[0]).second;
        ++heating_cells;
    }
    if(remaining>1e-12*heating_mass) throw std::runtime_error("Heating mass exceeds envelope");""")
    # Preserve all initializer behavior except the explicit test controls above.
    source = replace_once(source, r'<<",\n  \"injection_cell_count\": 1"',
                          r'<<",\n  \"injection_cell_count\": "<<heating_cells')
    source = replace_once(source, '<<mesh.cell_face_pos(0,mesh.nghost[0]).second', '<<heating_outer')
    source = replace_once(source, '<<mesh.cells[injection_flat].U[hd_core::RHO]*mesh.volume_cell[injection_flat]',
                          '<<heating_mass')
    return source


def parameters(source, updates):
    for key, value in updates.items():
        source, count = re.subn(rf"(?m)^{re.escape(key)}\s*=.*$", f"{key} = {value}", source)
        if count != 1:
            raise ValueError(f"Expected one parameter {key}, got {count}")
    return source


def experiments():
    cases = []
    def add(name, n=512, heating=0.1, atmosphere=1, cfl=0.2, outer=1100, energy=1.6e51):
        cases.append(dict(name=name, n=n, heating=heating, atmosphere=atmosphere,
                          cfl=cfl, outer=outer, energy=energy))
    for n in (256, 512, 1024):
        add(f"single_{n}", n=n, heating=0)
        add(f"mass01_{n}", n=n)
    add("mass01_2048", n=2048)
    add("mass01_512_atm001", atmosphere=0.01)
    add("mass01_512_atm100", atmosphere=100)
    add("mass01_512_cfl01", cfl=0.1)
    add("mass01_wide", n=752, outer=1500)
    add("mass003_512", heating=0.03)
    add("mass03_512", heating=0.3)
    add("no_bomb_512", heating=0, energy=0)
    return cases


def run_case(case, dest, binary, original_parameters):
    work = dest / case["name"]
    work.mkdir()
    par = parameters(original_parameters, dict(nxa=case["n"], x_upper=case["outer"],
        cfl=case["cfl"], maxsteps=200000, t_end=100, step_every=1,
        output_directory="output", output_suffix="validation"))
    (work / "parameters.txt").write_text(par)
    (work / "experiment.json").write_text(json.dumps(case, indent=2) + "\n")
    start = time.monotonic()
    with (work / "run.log").open("w") as log:
        subprocess.run([str(binary), "parameters.txt", str(case["atmosphere"]),
                        str(case["heating"]), str(case["energy"])], cwd=work,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f"Completed {case['name']} ({time.monotonic()-start:.1f} s)", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--cases", nargs="+", help="Names of selected experiments")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    cases = experiments()
    if args.cases:
        unknown = set(args.cases)-{c["name"] for c in cases}
        if unknown: parser.error(f"Unknown cases: {sorted(unknown)}")
        cases = [c for c in cases if c["name"] in args.cases]
    dest = (args.output or CASE / "output_validation" / dt.datetime.now().strftime("%Y%m%d_%H%M%S")).resolve()
    dest.mkdir(parents=True, exist_ok=False)
    build = dest / "build"
    build.mkdir()
    original = (CASE / "user_problem.cpp").read_text()
    original_parameters = (CASE / "parameters.txt").read_text()
    (build / "original_user_problem.cpp.txt").write_text(original)
    (build / "original_parameters.txt").write_text(original_parameters)
    (build / "user_problem.cpp").write_text(adapted_source(original))
    shutil.copy2(CASE / "user_config.h", build)
    shutil.copy2(HERE / "driver.cpp", build)
    shutil.copy2(HERE / "run.py", build / "run_validation.py")
    shutil.copy2(HERE / "analyze.py", build / "analyze_validation.py")
    sources = sorted(p for p in (ROOT / "src").rglob("*.cpp") if p.name != "main.cpp")
    provenance_paths = sorted(set(list((ROOT/"src").rglob("*.cpp")) + list((ROOT/"src").rglob("*.h"))
        + [CASE/"user_problem.cpp", CASE/"parameters.txt", CASE/"user_config.h", HERE/"run.py", HERE/"driver.cpp", HERE/"analyze.py", ROOT/"Makefile"]))
    (dest / "provenance.json").write_text(json.dumps({str(p.relative_to(ROOT)):
        hashlib.sha256(p.read_bytes()).hexdigest() for p in provenance_paths}, indent=2) + "\n")
    binary = build / "flare"
    print(f"Validation directory: {dest}", flush=True)
    with (dest / "build.log").open("w") as log:
        subprocess.run(["make", "-f", str(ROOT / "Makefile"), f"-j{min(args.jobs,4)}",
                        f"PROBLEM_DIR={build}", "USER_SOURCES=driver.cpp",
                        f"CPPFLAGS=-I{ROOT / 'src'}",
                        "FRAMEWORK_SOURCES="+" ".join(map(str,sources))],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_case,c,dest,binary,original_parameters) for c in cases]
        for f in concurrent.futures.as_completed(futures): f.result()
    # Fail if production inputs changed during the suite, instead of silently
    # attributing results to different source files.
    if (CASE/"user_problem.cpp").read_text()!=original or (CASE/"parameters.txt").read_text()!=original_parameters:
        raise RuntimeError("Production inputs changed during validation; results use the saved copies")
    print(f"Finished {len(cases)} experiments in {dest}", flush=True)


if __name__ == "__main__":
    main()
