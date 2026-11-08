#!/usr/bin/env python3
"""Render the two initialization examples from their saved output (matplotlib)."""
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/flare-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TESTS = Path(__file__).resolve().parents[1]
FIGURES = Path(__file__).resolve().parent / "figures"


def snapshot(path):
    lines = path.read_text().splitlines()
    meta = dict(line[2:].split("=", 1) for line in lines if line.startswith("# "))
    data = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(line for line in lines if not line.startswith("#"))]
    data = [row for row in data if not row.get("is_ghost", 0)]
    return meta, {key: np.array([row[key] for row in data]) for key in data[0]}


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    problem = TESTS / "File_Profile_2D"
    paths = sorted((problem / "output").glob("file_profile_*.csv"))
    _, initial = snapshot(paths[0])
    meta, final = snapshot(paths[-1])
    grid = np.loadtxt(problem / "mesh.dat", skiprows=1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, data, title in zip(axes, [initial, final], ["File-loaded initial state", f"Advected state: t = {float(meta['time']):g}"]):
        picture = ax.pcolormesh(grid[:, 0], grid[:, 1], data["rho"].reshape(32, 32),
                                cmap="viridis", vmin=0.8, vmax=1.2, shading="flat", rasterized=True)
        ax.set(title=title, xlabel="x [code units]", ylabel="y [code units]", aspect="equal")
    fig.colorbar(picture, ax=axes, label="Density [code units]", shrink=0.85)
    fig.suptitle("Nonuniform 32 x 32 Cartesian mesh + profile input")
    fig.savefig(FIGURES / "file_profile_2d.png", dpi=180)
    plt.close(fig)

    problem = TESTS / "CCSN_MESA_day13"
    info = json.loads((problem / "output/model.json").read_text())
    paths = sorted((problem / "output").glob("ccsn_mesa_*.csv"))
    _, initial = snapshot(paths[0])
    meta, final = snapshot(paths[-1])
    _, before = snapshot(problem / "output/initial_envelope.csv")
    radius = initial["x"]
    seconds = float(meta["time"]) * info["time_unit_s"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for data, label, style in [(initial, "Immediately after injection", "--"), (final, f"t = {seconds:g} s", "-")]:
        axes[0, 0].loglog(radius * info["length_unit_cm"], data["rho"] * info["density_unit_g_cm3"], style, label=label)
        inner = radius < 0.4
        axes[0, 1].semilogy(radius[inner], data["p"][inner] * info["pressure_unit_dyn_cm2"], style, label=label)
        axes[1, 0].plot(radius[inner], data["vx"][inner] * info["velocity_unit_cm_s"] / 1e5, style, label=label)
    axes[0, 0].axvline(info["surface_radius_code"] * info["length_unit_cm"], color="0.6", lw=1, label="Native surface")
    axes[0, 0].set(xlabel="Radius [cm]", ylabel="Density [g/cm³]", title="Native envelope + external atmosphere")
    axes[0, 1].set(xlabel="Radius [code length]", ylabel="Pressure [dyn/cm²]", title="Thermal bomb and outgoing wave")
    axes[1, 0].set(xlabel="Radius [code length]", ylabel="Radial velocity [km/s]", title="Inner envelope response")
    volumes = before["volume"]
    mass0 = np.sum(initial["rho"] * volumes)
    times, errors = [], []
    for path in paths:
        metadata, data = snapshot(path)
        times.append(float(metadata["time"]) * info["time_unit_s"])
        errors.append((np.sum(data["rho"] * volumes) / mass0 - 1) * 1e12)
    axes[1, 1].plot(times, errors, "o-", color="#087e8b")
    axes[1, 1].set(xlabel="Time since injection [s]", ylabel="Relative mass change × 10¹²", title="Mass conservation after wall-geometry fix")
    for ax in axes.flat:
        ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].legend(fontsize=8)
    fig.suptitle("MESA profile013: 4 M_sun cut, 1.6e51 erg thermal bomb\nUnchanged native active mesh; early evolution only")
    fig.savefig(FIGURES / "ccsn_mesa_early.png", dpi=180)
    plt.close(fig)
    print(FIGURES)


if __name__ == "__main__":
    main()
