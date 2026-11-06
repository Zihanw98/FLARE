#!/usr/bin/env python3
"""Analyze a completed run.py suite; all diagnostics are in solar code units."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np

C = 2.99792458e10


def snapshot(path):
    metadata = {}
    with path.open() as stream:
        while True:
            line = stream.readline()
            if not line.startswith("#"):
                break
            key, value = line[1:].strip().split("=", 1)
            metadata[key] = float(value)
        data = np.genfromtxt(stream, delimiter=",", names=line.strip().split(","))
    return metadata, data[data["is_ghost"] == 0]


def crossing(d, field, radius):
    reached = np.flatnonzero(d[field] >= radius)
    if not len(reached):
        return None
    i = int(reached[0])
    if i == 0:
        return float(d["time"][0])
    return float(np.interp(radius, d[field][i-1:i+1], d["time"][i-1:i+1]))


def analyze_case(work):
    output = work / "output"
    settings = json.loads((work / "experiment.json").read_text())
    model = json.loads((output / "model.json").read_text())
    d = np.genfromtxt(output / "diagnostics.csv", delimiter=",", names=True)
    assert d["time"][-1] == 100, f"Incomplete run: {work}"
    assert np.all(np.diff(d["time"]) > 0)
    assert np.all(d["min_rho"] > 0) and np.all(d["min_p"] > 0)
    for name in d.dtype.names:
        if not name.startswith("v_"):
            assert np.isfinite(d[name]).all(), (work, name)
    profiles = {}
    for p in sorted(output.glob("validation_*.csv")):
        meta, data = snapshot(p)
        profiles[meta["time"]] = data
        for name in data.dtype.names:
            assert np.isfinite(data[name]).all(), (p, name)
        assert np.all(data["rho"] > 0) and np.all(data["p"] > 0)
    assert set(profiles) == {0,10,20,30,40,50,75,100}, work
    initial = np.genfromtxt(output / "initial_envelope.csv", delimiter=",", names=True)
    volumes = initial["volume"]
    eu, vu = model["energy_unit_erg"], model["velocity_unit_cm_s"]
    before = np.sum(initial["energy_density_before"] * volumes)
    after = np.sum(profiles[0]["p"]/(profiles[0]["gamma"]-1)*volumes)
    actual = (after-before)*eu
    np.testing.assert_allclose(actual, settings["energy"], rtol=2e-13, atol=1e35)
    np.testing.assert_allclose(model["enclosed_mass_at_cut_code"],4,rtol=1e-12)
    np.testing.assert_allclose(d["mass"][0],np.sum(profiles[0]["rho"]*volumes),rtol=2e-13)
    if settings["heating"]:
        np.testing.assert_allclose(model["injection_mass_code"],settings["heating"],rtol=1e-13)
    for time, profile in profiles.items():
        record=d[np.flatnonzero(d["time"]==time)[0]]
        mass=np.sum(profile["rho"]*volumes)
        kinetic=np.sum(0.5*profile["rho"]*profile["vx"]**2*volumes)
        internal=np.sum(profile["p"]/(profile["gamma"]-1)*volumes)
        np.testing.assert_allclose([record["mass"],record["kinetic"],record["internal"]],
                                   [mass,kinetic,internal],rtol=2e-12,atol=1e-24)
    # Screening tolerances for numerical bookkeeping, not physical accuracy.
    assert np.max(abs(d["mass_residual"]))/d["mass"][0]<1e-5, work
    normalization=max(settings["energy"]/eu,before)
    assert np.max(abs(d["energy_residual"]))/normalization<1e-4, work
    peak = int(np.argmax(d["vmax"]))
    before_surface = d["front10"] < 1000
    row = dict(name=settings["name"], n=settings["n"], atmosphere=settings["atmosphere"],
        heating_mass=model["injection_mass_code"], heating_cells=model["injection_cell_count"],
        heated_cell_support_mass=float(np.sum(initial["rho"][:model["injection_cell_count"]]*volumes[:model["injection_cell_count"]])) if settings["energy"] else 0.0,
        injected_erg=actual, steps=int(d["step"][-1]),
        peak_v_c=float(d["vmax"][peak]*vu/C), peak_time_code=float(d["time"][peak]),
        peak_radius=float(d["vmax_r"][peak]), peak_cell_mass=float(d["vmax_cell_mass"][peak]),
        max_mass_gt_01c=float(np.max(d["mass_gt_01c"])),
        max_kinetic_gt_01c_frac=float(np.max(d["kinetic_gt_01c"])*eu/1.6e51),
        mass_residual_relative=float(np.max(abs(d["mass_residual"]))/d["mass"][0]),
        energy_residual_over_16e51=float(np.max(abs(d["energy_residual"]))*eu/1.6e51),
        presurface_energy_residual_over_16e51=float(np.max(abs(d["energy_residual"][before_surface]))*eu/1.6e51),
        final_mass=float(d["mass"][-1]), final_escaped_energy_erg=float(d["escaped_energy"][-1]*eu))
    for key in ("v_0001", "v_001", "v_01", "v_1"):
        row[f"max_{key}_km_s"] = float(np.nanmax(d[key])*vu/1e5)
        row[f"post_t10_max_{key}_km_s"] = float(np.nanmax(d[key][d["time"]>=10])*vu/1e5)
        row[f"t40_{key}_km_s"] = float(d[key][np.flatnonzero(d["time"]==40)[0]]*vu/1e5)
    for threshold in (10,100):
        times = {r:crossing(d,f"front{threshold}",r) for r in (600,800,950,1000,1050)}
        for a,b in ((600,800),(800,950),(950,1000),(1000,1050)):
            row[f"front{threshold}_{a}_{b}_km_s"] = (
                (b-a)/(times[b]-times[a])*vu/1e5
                if times[a] is not None and times[b] is not None else None)
        row[f"surface_time{threshold}_code"] = times[1000]
    if settings["energy"] == 0:
        row["max_abs_v_snapshot_km_s"] = max(float(np.max(abs(p["vx"]))) for p in profiles.values())*vu/1e5
        row["final_density_l1_change"] = float(np.sum(abs(profiles[100]["rho"]-profiles[0]["rho"])*volumes)/d["mass"][0])
    return row, d, model, profiles, volumes


def convergence(cases):
    rows = []
    for lo, hi in ((256,512),(512,1024),(1024,2048)):
        low, high = f"mass01_{lo}", f"mass01_{hi}"
        if low not in cases or high not in cases:
            continue
        _, _, _, lp, lv = cases[low]
        _, _, _, hp, hv = cases[high]
        np.testing.assert_allclose(hv.reshape(-1,2).sum(axis=1),lv,rtol=1e-10)
        for time in sorted((set(lp)&set(hp))-{0}):
            a,b = lp[time],hp[time]
            fields = lambda p: {
                "rho":p["rho"], "momentum":p["rho"]*p["vx"],
                "gas_energy":p["p"]/(p["gamma"]-1)+0.5*p["rho"]*p["vx"]**2}
            for field, fine in fields(b).items():
                restricted = (fine*hv).reshape(-1,2).sum(axis=1)/lv
                err = np.sum(abs(fields(a)[field]-restricted)*lv)/np.sum(abs(fine)*hv)
                rows.append(dict(coarse=lo,fine=hi,time=time,field=field,relative_l1=float(err)))
    return rows


def plots(dest, cases):
    os.environ.setdefault("MPLCONFIGDIR",str(dest/"matplotlib_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":10,"axes.grid":True,"grid.alpha":0.2})
    fig, axes = plt.subplots(2,2,figsize=(12,8),layout="constrained")
    for n in (256,512,1024,2048):
        for mode, style in (("single","--"),("mass01","-")):
            name=f"{mode}_{n}"
            if name not in cases: continue
            _, d, m, _, _=cases[name]
            axes[0,0].plot(d["time"]*m["time_unit_s"]/86400,d["vmax"]*m["velocity_unit_cm_s"]/C,
                           style,label=f"{mode}, N={n}",lw=1.1)
    axes[0,0].set(xlabel="Time [day]",ylabel="Maximum gas velocity / c",title="Resolution and heated mass")
    axes[0,0].legend(fontsize=8,ncol=2)
    for name, label in (("mass01_512_atm001","atmosphere x0.01"),("mass01_512","atmosphere x1"),
                        ("mass01_512_atm100","atmosphere x100"),("mass01_wide","outer radius 1500 Rsun")):
        if name not in cases: continue
        _,d,m,_,_=cases[name]
        axes[0,1].plot(d["time"]*m["time_unit_s"]/86400,d["vmax"]*m["velocity_unit_cm_s"]/C,label=label)
    axes[0,1].set(xlabel="Time [day]",ylabel="Maximum gas velocity / c",title="Sensitivity of the nearly empty exterior")
    axes[0,1].legend(fontsize=8)
    if "single_512" in cases:
        _,d,m,_,_=cases["single_512"]
        for key,mass in (("v_0001","0.001"),("v_001","0.01"),("v_01","0.1"),("v_1","1")):
            axes[1,0].plot(d["time"]*m["time_unit_s"]/86400,d[key]*m["velocity_unit_cm_s"]/1e5,label=f">= {mass} Msun faster")
    axes[1,0].set(xlabel="Time [day]",ylabel="Velocity threshold [km/s]",title="Original single-cell heating, N=512")
    axes[1,0].legend(fontsize=8)
    for name in ("single_512","mass01_512","mass01_1024","mass01_512_cfl01"):
        if name not in cases: continue
        _,d,m,_,_=cases[name]
        axes[1,1].plot(d["time"]*m["time_unit_s"]/86400,
                       d["energy_residual"]*m["energy_unit_erg"]/1.6e51,label=name)
    axes[1,1].set(xlabel="Time [day]",ylabel="(Egas + W + Eescaped - Einitial) / Ebomb",title="Energy accounting (boundary flux: trapezoids)")
    axes[1,1].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    axes[1,1].legend(fontsize=8)
    fig.savefig(dest/"validation.png",dpi=180)
    fig.savefig(dest/"validation.pdf")
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    args=parser.parse_args()
    dest=args.directory.resolve()
    cases={p.parent.name:analyze_case(p.parent) for p in sorted(dest.glob("*/experiment.json"))}
    if not cases: raise ValueError("No experiments found")
    rows=[v[0] for v in cases.values()]
    conv=convergence(cases)
    (dest/"summary.json").write_text(json.dumps({"cases":rows,"convergence":conv},indent=2,allow_nan=False)+"\n")
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with (dest/"summary.csv").open("w") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    plots(dest,cases)
    print("case                         vmax/c     mass>0.1c      E residual/Ebomb   v(>=0.01Msun) km/s")
    for row in rows:
        print(f"{row['name']:28} {row['peak_v_c']:8.4f} {row['max_mass_gt_01c']:13.4g} "
              f"{row['energy_residual_over_16e51']:18.3g} {row['max_v_001_km_s']:18.1f}")
    print(f"All {len(rows)} cases passed completion, positive-state, mass-cut, injection-budget, snapshot and conservation-screening checks.")
    print(f"Summary and plots: {dest}")
    (dest/"analysis_sha256.txt").write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()+"\n")


if __name__=="__main__":
    main()
