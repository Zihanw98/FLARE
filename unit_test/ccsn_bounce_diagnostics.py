#!/usr/bin/env python3
"""Measure fast ejecta, swept atmosphere, and Newtonian energy in CCSN snapshots.

This is a read-only diagnostic of the case's spherical PCM/HLL simulation.
Relativistic energies below
are comparisons for the saved Newtonian cell masses, not a relativistic rerun.
"""
import argparse
from functools import lru_cache
import importlib.util
import json
from pathlib import Path

import numpy as np

G = 6.67430e-8
C = 2.99792458e10
MSUN = 1.98847e33
CASE = Path(__file__).resolve().parents[1] / "tests/CCSN_Hybrid_1D"


@lru_cache(maxsize=1)
def atmosphere_quadrature():
    return np.polynomial.legendre.leggauss(128)


def read_snapshot(path):
    with path.open() as stream:
        meta = {}
        while True:
            line = stream.readline()
            if not line.startswith("#"):
                names = line.strip().split(",")
                break
            key, value = line[1:].strip().split("=", 1)
            meta[key] = float(value)
        data = np.loadtxt(stream, delimiter=",", dtype=[(n, float) for n in names])
    for key in ("unit_mass", "unit_length", "unit_time"):
        if meta[key] != 1:
            raise ValueError("This case diagnostic expects CGS snapshots")
    return meta["time"], data[data["is_ghost"] == 0]


def gravity_weights(grid):
    """Exact -integral G M(<r)/r dm for piecewise constant shell densities.

    Three-point Gauss quadrature integrates the fourth-degree polynomial exactly.
    The factored cube difference avoids cancellation for thin shells.
    """
    a, b = grid["r_left"], grid["r_right"]
    x, w = np.polynomial.legendre.leggauss(3)
    r = (a + b)[:, None] / 2 + (b - a)[:, None] / 2 * x
    dr = (b - a)[:, None] / 2 * w
    cube_difference = (r - a[:, None]) * (r*r + r*a[:, None] + a[:, None]**2)
    return (4*np.pi*np.sum(r*dr, axis=1),
            16*np.pi**2/3*np.sum(cube_difference*r*dr, axis=1))


def atmosphere_mass(inner, outer):
    """Initial artificial-atmosphere mass between two radii in cm."""
    if outer <= inner:
        return 0.0
    x, w = atmosphere_quadrature()
    r = (inner + outer)/2 + (outer - inner)/2*x
    return float(4*np.pi*(outer-inner)/2*np.sum(w*1e3*r*r/(1+(r/1e8)**4)))


def run(directory, destination):
    spec = importlib.util.spec_from_file_location("ccsn_case_analysis", CASE / "analyze.py")
    case_analysis = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(case_analysis)
    model = json.loads((directory / "model.json").read_text())
    grid = np.genfromtxt(directory / "grid.csv", names=True, delimiter=",")
    volume = grid["volume"]
    ga, gb = gravity_weights(grid)
    # Independent analytic check of the potential-energy integral.
    uniform_mass = np.sum(volume)
    uniform_w = -G*np.sum(np.r_[0, np.cumsum(volume)[:-1]]*ga + gb)
    assert np.isclose(uniform_w, -3*G*uniform_mass**2/(5*grid["r_right"][-1]), rtol=1e-12)
    files = sorted(directory.glob("bounce_*.csv"))
    records = []
    first = final = None
    for i, path in enumerate(files):
        time, data = read_snapshot(path)
        assert len(data) == len(volume)
        rho, vel = data["rho"], data["vx"]
        mass = rho*volume
        interior = np.r_[0, np.cumsum(mass)[:-1]]
        potential = -G*np.sum(rho*interior*ga + rho*rho*gb)
        kinetic = np.sum(0.5*mass*vel**2)
        internal = np.sum(data["internal_energy"]*volume)
        gas_energy = np.sum(data["total_energy"]*volume)
        assert np.isclose(gas_energy, internal+kinetic, rtol=1e-12)
        pc, uc, gc = case_analysis.cold_state(rho, model)
        cs = np.sqrt((gc*pc+model["gamma_thermal"]*(data["p"]-pc))/rho)
        a, b = grid["r_left"], grid["r_right"]
        geom = a*a+a*b+b*b
        alpha = (b*b+2*a*b+3*a*a)/(4*geom)
        acceleration = -3*G*(interior+alpha*mass)/geom
        gravity_work = np.sum(mass*vel*acceleration)
        # Reconstruct the current PCM/HLL mass flux to check the gravity-work
        # coupling, independently of differences between saved energy snapshots.
        sl = np.minimum(vel[:-1]-cs[:-1], vel[1:]-cs[1:])
        sr = np.maximum(vel[:-1]+cs[:-1], vel[1:]+cs[1:])
        fl, fr = rho[:-1]*vel[:-1], rho[1:]*vel[1:]
        face_flux = np.where(sl >= 0, fl, np.where(sr <= 0, fr,
            (sr*fl-sl*fr+sl*sr*(rho[1:]-rho[:-1]))/(sr-sl)))
        hll_diffusion_flux = np.where((sl < 0) & (sr > 0),
            sl*sr*(rho[1:]-rho[:-1])/(sr-sl), 0.0)
        mass_flux = 4*np.pi*np.r_[0, b]**2*np.r_[0, face_flux, rho[-1]*vel[-1]]
        mass_rate = -np.diff(mass_flux)
        outer_potential = np.cumsum((mass*ga/volume)[::-1])[::-1]-mass*ga/volume
        potential_per_mass = -G*(interior*ga/volume+2*mass*gb/volume**2+outer_potential)
        potential_rate = np.sum(potential_per_mass*mass_rate)
        diffusion_work_faces = 4*np.pi*b[:-1]**2*hll_diffusion_flux*np.diff(potential_per_mass)
        diffusion_rate = np.sum(diffusion_work_faces)
        shock = case_analysis.shock_measurement(data, model)
        j = int(np.argmax(vel))
        records.append(dict(time_s=time, shock_radius_cm=shock[0],
                            shock_speed_RH_cm_s=shock[5],
                            downstream_velocity_cm_s=shock[1], upstream_velocity_cm_s=shock[2],
                            compression=shock[3], vmax_cm_s=vel[j], vmax_radius_cm=data["x"][j],
                            vmax_density_g_cm3=rho[j], max_sound_speed_cm_s=np.max(cs),
                            central_density_g_cm3=rho[0], mass_g=np.sum(mass),
                            kinetic_erg=kinetic, internal_erg=internal,
                            cold_internal_erg=np.sum(uc*volume), potential_erg=potential,
                            total_energy_erg=gas_energy+potential,
                            gravity_work_erg_s=gravity_work,
                            potential_rate_HLL_erg_s=potential_rate,
                            semidiscrete_energy_defect_erg_s=gravity_work+potential_rate,
                            HLL_density_diffusion_potential_rate_erg_s=diffusion_rate,
                            advection_gravity_mismatch_erg_s=gravity_work+potential_rate-diffusion_rate,
                            diffusion_inside_50km_erg_s=np.sum(diffusion_work_faces[b[:-1] < 5e6]),
                            diffusion_inside_300km_erg_s=np.sum(diffusion_work_faces[b[:-1] < 3e7]),
                            dual_energy_correction_erg=np.sum(data["dual_energy_correction"]*volume),
                            outer_velocity_cm_s=vel[-1], outer_pressure=data["p"][-1],
                            outer_energy_density=data["total_energy"][-1],
                            swept_initial_atmosphere_g=atmosphere_mass(model["initial_core_radius"], shock[0])))
        if i == 0:
            first = data
        final = data
    history = {k: np.array([row[k] for row in records]) for k in records[0]}
    times = history["time_s"]
    bounce = int(np.argmax(history["central_density_g_cm3"]))
    # A finite time fit suppresses grid quantization of the measured shock radius.
    fits = np.full(len(times), np.nan)
    for i, time in enumerate(times):
        if time <= times[bounce] + 0.00025:
            continue
        take = (abs(times-time) <= 0.00025) & (times > times[bounce])
        if np.count_nonzero(take) >= 3:
            fits[i] = np.polyfit(times[take]-time, history["shock_radius_cm"][take], 1)[0]
    history["shock_speed_position_fit_cm_s"] = fits
    initial_energy = history["total_energy_erg"][0]
    drift = history["total_energy_erg"] - initial_energy
    fast = []
    mass = final["rho"]*volume
    beta = final["vx"]/C
    for threshold in (0.1, 0.3, 0.5, 0.7, 0.8):
        take = beta > threshold
        fast.append(dict(velocity_threshold_c=threshold, mass_g=np.sum(mass[take]),
                         mass_Msun=np.sum(mass[take])/MSUN,
                         kinetic_newtonian_erg=np.sum(0.5*mass[take]*final["vx"][take]**2),
                         kinetic_SR_same_mass_and_velocity_erg=(
                             np.sum(mass[take]*C*C*(1/np.sqrt(1-beta[take]**2)-1))
                             if np.all(beta[take] < 1) else None)))
    shock_radius = records[-1]["shock_radius_cm"]
    a, b = grid["r_left"], grid["r_right"]
    lo, hi = np.maximum(a, model["initial_core_radius"]), np.minimum(b, shock_radius)
    partial_volume = 4*np.pi/3*np.maximum(0, hi**3-lo**3)
    initial_shell_mass = np.sum(first["rho"]*partial_volume)
    initial_outer_volume = 4*np.pi/3*np.maximum(0, b**3-np.maximum(a, model["initial_core_radius"])**3)
    final_outer_volume = 4*np.pi/3*np.maximum(0, b**3-np.maximum(a, shock_radius)**3)
    # In spherical ordering the initial core is behind the final shock. Subtract
    # small exterior masses rather than nearly equal ~1.5 Msun enclosed masses.
    net_atmosphere = np.sum(first["rho"]*initial_outer_volume)-np.sum(final["rho"]*final_outer_volume)
    # The outer boundary is far from the shock; quantify the residual flux scale.
    boundary_luminosity = (4*np.pi*model["outer_radius"]**2*history["outer_velocity_cm_s"]
                          *(history["outer_energy_density"]+history["outer_pressure"]))
    surface_flux = np.trapezoid(boundary_luminosity, times)
    selected = []
    for target in (0, .008, times[bounce], .0095, .010, .011, .012, .013, .014, .015):
        i = int(np.argmin(abs(times-target)))
        selected.append({k: float(v[i]) if np.isfinite(v[i]) else None for k, v in history.items()})
    gravity_integral = np.trapezoid(history["gravity_work_erg_s"], times)
    defect_integral = np.trapezoid(history["semidiscrete_energy_defect_erg_s"], times)
    gas_change = (history["kinetic_erg"]+history["internal_erg"])[-1] - (history["kinetic_erg"]+history["internal_erg"])[0]
    report = dict(
        input_directory=str(directory.resolve()), frames=len(files), cells=len(volume),
        initial_core_radius_km=model["initial_core_radius"]/1e5,
        bounce_time_s=times[bounce], final=records[-1], selected_times=selected,
        final_shock_speed_position_fit_c=fits[-1]/C,
        final_initial_atmosphere_integral_g=atmosphere_mass(model["initial_core_radius"], shock_radius),
        final_initial_shell_mass_from_grid_g=initial_shell_mass,
        final_swept_atmosphere_mass_balance_estimate_g=net_atmosphere,
        initial_atmosphere_to_1500km_g=atmosphere_mass(model["initial_core_radius"], 1.5e8),
        final_fast_ejecta=fast,
        global_vmax_c=float(np.max(history["vmax_cm_s"])/C),
        global_max_sound_speed_c=float(np.max(history["max_sound_speed_cm_s"])/C),
        energy_initial_erg=initial_energy, energy_final_erg=history["total_energy_erg"][-1],
        energy_final_drift_erg=drift[-1], energy_max_abs_drift_erg=np.max(abs(drift)),
        energy_max_abs_drift_over_initial_binding=np.max(abs(drift))/abs(history["potential_erg"][0]),
        energy_max_abs_drift_over_initial_total=np.max(abs(drift))/abs(initial_energy),
        integrated_cell_center_gravity_work_erg=gravity_integral,
        gas_energy_change_minus_gravity_work_erg=gas_change-gravity_integral+surface_flux,
        integrated_HLL_gravity_energy_defect_erg=defect_integral,
        integrated_HLL_density_diffusion_potential_change_erg=np.trapezoid(history["HLL_density_diffusion_potential_rate_erg_s"], times),
        integrated_advection_gravity_mismatch_erg=np.trapezoid(history["advection_gravity_mismatch_erg_s"], times),
        integrated_density_diffusion_inside_50km_erg=np.trapezoid(history["diffusion_inside_50km_erg_s"], times),
        integrated_density_diffusion_inside_300km_erg=np.trapezoid(history["diffusion_inside_300km_erg_s"], times),
        energy_drift_minus_integrated_defect_erg=drift[-1]-defect_integral,
        approximate_outer_gas_energy_flux_integral_erg=surface_flux,
        notes=["Shock RH speed uses finite-width downstream/upstream samples; position speed uses a 0.5 ms full window (one-sided at endpoints).",
               "Swept initial atmosphere is an initial Eulerian shell integral; infall means it is not an exact Lagrangian swept-mass tracer.",
               "The exterior-mass balance estimate includes infall but assumes the original core lies behind the identified front; it remains sensitive to numerical shock width.",
               "Potential energy is exact for piecewise constant saved cell densities and Newtonian spherical gravity.",
               "The energy-defect reconstruction assumes this case's PCM/HLL flux, reflecting origin, outflow outer boundary, and no central point mass.",
               "Relativistic kinetic energies are post-processing comparisons, not self-consistent relativistic solutions."])
    destination.mkdir(parents=True, exist_ok=True)
    (destination/"relativistic_check.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    np.savetxt(destination/"relativistic_history.csv", np.column_stack(list(history.values())),
               header=",".join(history), delimiter=",", comments="")
    plot(history, first, final, volume, model, destination)
    print(json.dumps(report, indent=2, allow_nan=False))


def plot(history, first, final, volume, model, destination):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .2})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    t = history["time_s"]*1e3
    ax = axes[0, 0]
    ax.plot(t, history["vmax_cm_s"]/C, label="Maximum gas velocity")
    ax.plot(t, history["shock_speed_position_fit_cm_s"]/C, label="Shock: radius/time fit")
    valid = t > history["time_s"][np.argmax(history["central_density_g_cm3"])]*1e3+.25
    ax.plot(t[valid], history["shock_speed_RH_cm_s"][valid]/C, alpha=.6, label="Shock: mass jump")
    ax.axhline(1, color="black", ls=":", label="Speed of light")
    bounce_ms = t[np.argmax(history["central_density_g_cm3"])]
    ax.set(xlim=(max(0, bounce_ms-.1), t[-1]), ylim=(0, 1.3), xlabel="Time since start [ms]", ylabel="Velocity / c",
           title="The fast front enters the relativistic regime")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    ax.loglog(first["x"]/1e5, first["rho"], label="Initial density")
    ax.loglog(final["x"]/1e5, final["rho"], label=f"Final density ({t[-1]:g} ms)")
    ax.axhline(1e8, color="gray", ls=":", label="Animation density cutoff")
    ax.axvline(model["initial_core_radius"]/1e5, color="black", ls="--", label="Initial core surface")
    ax.set(xlim=(10, 1500), ylim=(1e2, 1e13), xlabel="Radius [km]", ylabel=r"Density [g cm$^{-3}$]",
           title="A compact core surrounded by artificial atmosphere")
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    beta = np.linspace(0, .9, 181)
    masses = np.array([np.sum(final["rho"][final["vx"] > b*C]*volume[final["vx"] > b*C])/MSUN for b in beta])
    ax.semilogy(beta, np.where(masses > 0, masses, np.nan))
    ax.set(xlabel="Outward velocity threshold / c", ylabel=r"Mass above threshold [$M_\odot$]",
           title="Only a small mass reaches the largest velocities", xlim=(0, .9))
    ax = axes[1, 1]
    for key, label in (("kinetic_erg", "Kinetic"), ("internal_erg", "Internal"),
                       ("potential_erg", "Gravitational"), ("total_energy_erg", "Total")):
        ax.plot(t, history[key]/1e52, label=label)
    ax.set(xlabel="Time since start [ms]", ylabel=r"Energy [$10^{52}$ erg]",
           title="Newtonian energy budget, including self-gravity", xlim=(0, t[-1]))
    ax.legend(fontsize=8)
    delta = (history["total_energy_erg"][-1]-history["total_energy_erg"][0])/1e52
    ax.text(.04, .08, rf"Total-energy drift: +{delta:.2f} $\times 10^{{52}}$ erg",
            transform=ax.transAxes, fontsize=10, color="firebrick")
    fig.savefig(destination/"relativistic_check.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=CASE/"output_bounce_hires")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run(args.input, args.output or args.input/"diagnostics")
