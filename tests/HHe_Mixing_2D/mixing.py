"""Shared reader and diagnostics for the neutral H/He example (no EOS side effects)."""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from animate import read_snapshot, check_grid, edges, signature
from animate3d import snapshot_files
sys.path.insert(0, str(ROOT / "tests"))
from hydro_suite import parameters, require

# CODATA constants match src/physcon/constants.h. Species masses use A*m_u.
M_U = 1.66053906892e-24
K_B = 1.380649e-16
RAW_FIELDS = tuple((key, key, "viridis") for key in (
    "rho", "p", "vx", "vy", "gamma", "comp_0", "comp_1",
    "internal_energy", "total_energy", "dual_energy_correction"))
COMPOSITION = (("H", "H mass fraction", "Blues"),
               ("He", "He mass fraction", "Oranges"),
               ("mixing", "Mixing: 4 X(H) X(He)", "magma"))
HYDRO = (("rho", "Density [code]", "viridis"),
         ("p", "Pressure [code]", "magma"),
         ("temperature", "Neutral ideal-gas T [K]", "inferno"),
         ("speed", "Speed [code]", "viridis"),
         ("vx", "x velocity [code]", "RdBu_r"),
         ("vy", "y velocity [code]", "RdBu_r"),
         ("internal_energy", "Internal energy density [code]", "magma"),
         ("total_energy", "Total energy density [code]", "magma"))


def load_frame(path):
    step = int(path.stem.rsplit("_", 1)[1])
    frame = read_snapshot(path, step, fields=RAW_FIELDS)
    require(len(frame.y) > 1, "H/He mixing requires a 2D grid")
    fields = {key: frame.values[k] for k, (key, _, _) in enumerate(frame.fields)}
    require(np.all(fields["rho"] > 0) and np.all(fields["p"] > 0), "Nonpositive rho/p")
    require(np.allclose(fields["gamma"], 5/3, rtol=0, atol=1e-14), "Requires gamma=5/3")
    species = {}
    for k in range(2):
        name = frame.metadata[f"comp_{k}_name"]
        require(name not in species, "Duplicate species metadata")
        species[name] = k
    require(set(species) == {"HI", "He_I"}, "Requires neutral HI and He_I metadata")
    for label, name, mass_number in (("H", "HI", 1), ("He", "He_I", 4)):
        k = species[name]
        require(int(frame.metadata[f"comp_{k}_A"]) == mass_number and
                int(frame.metadata[f"comp_{k}_charge"]) == 0, "Invalid species mass/charge")
        fields[label] = fields[f"comp_{k}"]
        require(np.min(fields[label]) >= -5e-10 and np.max(fields[label]) <= 1+5e-10,
                f"Unbounded {label} mass fraction at step {step}")
    require(np.max(abs(fields["H"] + fields["He"] - 1)) < 5e-10,
            f"Mass fractions do not sum to one at step {step}")
    require(all(value is not None for value in frame.units), "Missing physical base units")
    _, length, duration, _ = frame.units
    # This is a derived neutral ideal-gas diagnostic, not Cell.temperature or an ionization EOS.
    fields["temperature"] = (fields["p"] / fields["rho"] * (length/duration)**2 * M_U/K_B
                              / (fields["H"] + fields["He"]/4))
    fields["mixing"] = 4 * fields["H"] * fields["He"]
    fields["speed"] = np.hypot(fields["vx"], fields["vy"])
    require(all(np.isfinite(value).all() for value in fields.values()), "Nonfinite derived field")
    require(np.all(fields["temperature"] > 0), "Nonpositive temperature")
    return frame, fields


def run_paths(parameter_file):
    parameter_file = parameter_file.resolve()
    cfg = parameters(parameter_file)
    output = (parameter_file.parent / cfg["output_directory"]).resolve()
    # Use this run's completed-output manifest, never stale files from another run.
    paths = snapshot_files(output, cfg["output_suffix"], output / "run.log")
    return cfg, output, paths


def check_initial(frame, fields, cfg):
    expected_axes = []
    for axis in "xy":
        lo, hi = float(cfg[f"{axis}_lower"]), float(cfg[f"{axis}_upper"])
        n = int(cfg[f"n{axis}a"])
        expected_axes.append(lo + (np.arange(n)+.5)*(hi-lo)/n)
    for actual, expected in zip((frame.x, frame.y), expected_axes):
        require(actual.shape == expected.shape and np.allclose(actual, expected, rtol=0, atol=1e-13),
                "Unexpected active grid")
    x = (frame.x[None, :] - float(cfg["x_lower"])) / (float(cfg["x_upper"])-float(cfg["x_lower"]))
    y = (frame.y[:, None] - float(cfg["y_lower"])) / (float(cfg["y_upper"])-float(cfg["y_lower"]))
    layer = .5*(1-np.tanh(np.cos(2*np.pi*y)/.15))
    expected = dict(H=1-layer, He=layer, rho=1+layer, p=2.5, vx=layer-.5,
                    vy=.02*np.sin(4*np.pi*x)*np.exp(-(np.cos(2*np.pi*y)/.3)**2))
    for key, values in expected.items():
        require(np.allclose(fields[key], values, rtol=0, atol=2e-13), f"Incorrect initial {key}")


def diagnostics(frame, fields):
    # Volume per unit depth, using the same inferred Cartesian edges as animate.py.
    volume = np.diff(edges(frame.y))[:, None]*np.diff(edges(frame.x))[None, :]
    mass = fields["rho"] * volume
    return dict(step=frame.step, time=frame.time, mass=float(mass.sum()),
                hydrogen_mass=float((mass*fields["H"]).sum()),
                helium_mass=float((mass*fields["He"]).sum()),
                energy=float((volume*fields["total_energy"]).sum()),
                correction=float((volume*fields["dual_energy_correction"]).sum()),
                fraction_sum_error=float(abs(fields["H"]+fields["He"]-1).max()),
                fraction_min=float(min(fields["H"].min(), fields["He"].min())),
                fraction_max=float(max(fields["H"].max(), fields["He"].max())),
                mass_weighted_mixing=float((mass*fields["mixing"]).sum()/mass.sum()),
                transverse_kinetic_energy=float((.5*mass*fields["vy"]**2).sum()))
