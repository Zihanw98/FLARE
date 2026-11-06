#!/usr/bin/env python3
"""Check solar-scale hydrostatic initial data and measured evolution drift."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

CASE = Path(__file__).resolve().parent


def config(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0]
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


def snapshot(path):
    with path.open() as f:
        meta = dict(f.readline()[2:].strip().split("=", 1) for _ in range(3))
        data = np.genfromtxt((line for line in f if not line.lstrip().startswith("#")), delimiter=",", names=True)
    return float(meta["time"]), data[data["is_ghost"] == 0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--parameters", type=Path, default=CASE / "parameters.txt")
    parser.add_argument("--executable", type=Path, default=CASE / "flare")
    args = parser.parse_args()
    cfg = config(args.parameters)
    directory = CASE / cfg["output_directory"]
    prefix = cfg["output_suffix"]
    if args.run:
        directory.mkdir(parents=True, exist_ok=True)
        # Only replace this diagnostic's generated snapshots, never source files.
        for path in directory.glob(f"{prefix}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].csv"):
            path.unlink()
        with (directory / "run.log").open("w") as log:
            subprocess.run([str(args.executable.resolve()), "-par", str(args.parameters.resolve())],
                           cwd=CASE, stdout=log, stderr=subprocess.STDOUT, check=True)
    paths = sorted(directory.glob(f"{prefix}_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].csv"))
    if len(paths) < 2:
        raise RuntimeError("Need initial and evolved snapshots; run make test first")
    model = json.loads((directory / "model.json").read_text())
    states = [snapshot(path) for path in paths]
    t0, initial = states[0]
    n = int(cfg["nxa"])
    faces = np.linspace(float(cfg["x_lower"]), float(cfg["x_upper"]), n + 1)
    a, b = faces[:-1], faces[1:]
    volume = 4*np.pi/3*(b**3-a**3)
    rc = 0.75*(b**4-a**4)/(b**3-a**3)
    dm0 = initial["rho"]*volume
    mass0 = dm0.sum()
    G = model["G_code"]

    def radius99(dm):
        target = 0.99*dm.sum()
        cumulative = np.cumsum(dm)
        i = int(np.searchsorted(cumulative, target))
        before = cumulative[i-1] if i else 0.0
        fraction = (target-before)/dm[i]
        return (a[i]**3+fraction*(b[i]**3-a[i]**3))**(1/3)

    def total_energy(data):
        rho, p, v = data["rho"], data["p"], data["vx"]
        dm = rho*volume
        inner = np.cumsum(dm)-dm
        # Exact potential energy for piecewise uniform-density spherical shells.
        shell_factor = dm/(b**3-a**3)
        W = -4*np.pi*G*np.sum(rho*(inner*(b*b-a*a)/2 + shell_factor*
            ((b**5-a**5)/5-a**3*(b*b-a*a)/2)))
        return np.sum((p/(data["gamma"]-1)+0.5*rho*v*v)*volume)+W

    initial_energy = total_energy(initial)
    initial_r99 = radius99(dm0)
    np.savetxt(directory / "initial_profile_cgs.csv",
               np.column_stack((a*model["length_unit_cm"], b*model["length_unit_cm"],
                                rc*model["length_unit_cm"],
                                initial["rho"]*model["density_unit_g_cm3"],
                                initial["p"]*model["pressure_unit_dyn_cm2"],
                                np.cumsum(dm0)*model["mass_unit_g"])),
               delimiter=",", comments="",
               header="r_inner_cm,r_outer_cm,r_volume_center_cm,rho_g_cm3,p_dyn_cm2,m_enclosed_outer_g")
    rows = []
    for t, data in states:
        if len(data) != n or not all(np.all(np.isfinite(data[k])) for k in data.dtype.names):
            raise AssertionError("Invalid snapshot shape or nonfinite values")
        if np.any(data["rho"] <= 0) or np.any(data["p"] <= 0):
            raise AssertionError("Nonpositive density or pressure")
        dm = data["rho"]*volume
        mach = data["vx"]/np.sqrt(data["gamma"]*data["p"]/data["rho"])
        rows.append([t, abs(dm.sum()/mass0-1),
                     np.sum(abs(data["rho"]-initial["rho"])*volume)/mass0,
                     np.sqrt(np.sum(dm*mach*mach)/dm.sum()),
                     abs(radius99(dm)/initial_r99-1),
                     abs(total_energy(data)/initial_energy-1),
                     data["rho"][0]/initial["rho"][0]-1])
    rows = np.asarray(rows)

    # Independent continuum-HSE diagnostic; this is NOT the numerical flux RHS.
    enclosed = np.cumsum(dm0)-dm0 + dm0*(rc**3-a**3)/(b**3-a**3)
    force = initial["rho"]*G*enclosed/rc**2
    dpdr = np.gradient(initial["p"], rc, edge_order=2)
    interior = (rc < 0.9) & (np.arange(n)>1)
    hse = np.sqrt(np.sum((dpdr[interior]+force[interior])**2*volume[interior]) /
                  np.sum(force[interior]**2*volume[interior]))
    metrics = {
        "cells": n, "end_time_code": float(rows[-1,0]),
        "end_time_seconds": float(rows[-1,0]*model["time_unit_s"]),
        "initial_mass_solar": float(mass0), "initial_r99_solar": float(initial_r99),
        "central_density_g_cm3": model["central_density_code"]*model["density_unit_g_cm3"],
        "central_pressure_dyn_cm2": model["central_pressure_code"]*model["pressure_unit_dyn_cm2"],
        "initial_continuum_hse_relative_rms": float(hse),
        "max_mass_relative_drift": float(rows[:,1].max()),
        "max_density_relative_L1_drift": float(rows[:,2].max()),
        "max_mass_weighted_rms_mach": float(rows[:,3].max()),
        "max_r99_relative_drift": float(rows[:,4].max()),
        "max_total_energy_relative_drift": float(rows[:,5].max()),
        "max_central_density_relative_drift": float(abs(rows[:,6]).max()),
        "final_max_atmosphere_speed_km_s": float(np.max(abs(states[-1][1]["vx"][rc>=1]))*
                                                model["velocity_unit_cm_s"]/1e5),
    }
    # Finite-time near-equilibrium acceptance, not exact well-balanced preservation.
    checks = {
        "starts_at_rest": t0 == 0 and np.max(abs(initial["vx"])) == 0,
        "completed": abs(rows[-1,0]-float(cfg["t_end"])) < 1e-10,
        "solar_mass": abs(mass0-1) < 1e-5,
        "emden_surface": abs(model["xi1"]-6.89684862) < 1e-6,
        "emden_mass": abs(model["mu1"]-2.01823595) < 1e-6,
        "initial_hse": hse < 1e-3,
        "mass_drift": rows[:,1].max() < 1e-5,
        "density_drift": rows[:,2].max() < 0.01,
        "rms_mach": rows[:,3].max() < 0.01,
        "r99_drift": rows[:,4].max() < 0.01,
        "energy_drift": rows[:,5].max() < 0.01,
    }
    metrics["checks"] = {k: bool(v) for k,v in checks.items()}
    (directory / "metrics.json").write_text(json.dumps(metrics, indent=2)+"\n")
    np.savetxt(directory / "metrics.csv", rows, delimiter=",",
               header="time,mass_drift,density_L1,rms_mach,r99_drift,energy_drift,central_density_change")

    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "flare-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2,2, figsize=(11,8), constrained_layout=True)
    final = states[-1][1]
    for data,label,style in [(initial,"Initial","--"),(final,"Evolved","-")]:
        axes[0,0].semilogy(rc,data["rho"]*model["density_unit_g_cm3"],style,label=label)
        axes[0,1].semilogy(rc,data["p"]*model["pressure_unit_dyn_cm2"],style,label=label)
    axes[0,0].set(ylabel="Density [g cm$^{-3}$]",xlabel="Radius / solar radius")
    axes[0,1].set(ylabel="Pressure [dyn cm$^{-2}$]",xlabel="Radius / solar radius")
    axes[1,0].plot(rc,final["vx"]*model["velocity_unit_cm_s"]/1e5)
    axes[1,0].set(xlabel="Radius / solar radius",ylabel="Final radial velocity [km/s]")
    axes[1,1].plot(rows[:,0],rows[:,2]*100,label="Density L1 drift")
    axes[1,1].plot(rows[:,0],rows[:,4]*100,label="99% mass-radius drift")
    axes[1,1].set(xlabel="Time / dynamical time",ylabel="Drift [%]")
    for ax in axes.flat:
        ax.grid(alpha=0.2)
    axes[0,0].legend(); axes[0,1].legend(); axes[1,1].legend()
    fig.suptitle("Solar-scale n=3 polytrope — hydrostatic preservation")
    fig.savefig(directory / "equilibrium.png", dpi=160)
    plt.close(fig)
    print(json.dumps(metrics, indent=2))
    if not all(checks.values()):
        raise SystemExit("FAILED: "+", ".join(k for k,v in checks.items() if not v))
    print("PASS: finite-time hydrostatic preservation within stated tolerances")


if __name__ == "__main__":
    main()
