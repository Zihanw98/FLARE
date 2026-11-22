#!/usr/bin/env python3
"""Check nuclear crossing, core bounce, and an outward shock in saved toy data."""
import argparse
import json
from pathlib import Path
import numpy as np


def read_snapshot(path):
    with path.open() as stream:
        lines = stream.readlines()
    meta = dict(line[2:].strip().split('=', 1) for line in lines if line.startswith('# '))
    data = np.genfromtxt([line for line in lines if not line.startswith('#')], names=True, delimiter=',')
    return float(meta['time']), data[data['is_ghost'] == 0]


def cold_state(rho, model):
    rn, g1, g2 = model['rho_nuclear'], model['gamma_low'], model['gamma_stiff']
    pn = model['K_cold'] * rn**g1
    g = np.where(rho > rn, g2, g1)
    pc = pn * (rho/rn)**g
    uc = pc/(g-1) + np.where(rho > rn, rho*pn/rn*(1/(g1-1)-1/(g2-1)), 0)
    return pc, uc, g


def shock_measurement(data, model):
    # Locate the steep compressive velocity drop, then sample three/four cells
    # either side. This is a finite-resolution diagnostic, not an exact jump fit.
    r, v, rho, p = (data[k] for k in ('x', 'vx', 'rho', 'p'))
    region = (r[:-1] > 5e5) & (r[:-1] < 1e8)
    score = np.where(region, v[:-1]-v[1:], -np.inf)
    j = int(np.argmax(score))
    left, right = max(0, j-3), min(len(r)-1, j+4)
    compression = rho[left]/rho[right]
    pressure_jump = p[left]/p[right]
    shock_speed = ((rho[left]*v[left]-rho[right]*v[right])/(rho[left]-rho[right])
                   if compression > 1.01 else 0.0)
    pc, _, gc = cold_state(rho[right], model)
    cs = np.sqrt((gc*pc+model['gamma_thermal']*(p[right]-pc))/rho[right])
    return ((r[j]+r[j+1])/2, v[left], v[right], compression, pressure_jump,
            shock_speed, (shock_speed-v[right])/cs)


def analyze(directory, expected_end=None, write=True):
    model = json.loads((directory/'model.json').read_text())
    grid = np.genfromtxt(directory/'grid.csv', names=True, delimiter=',')
    files = sorted(directory.glob('bounce_*.csv'))
    if len(files) < 3:
        raise RuntimeError('Need initial and time-resolved bounce snapshots')
    frames = [read_snapshot(f) for f in files]
    history = []
    for time, data in frames:
        for name in ('rho', 'vx', 'p', 'internal_energy', 'total_energy'):
            assert np.all(np.isfinite(data[name])), (time, name)
        assert len(data) == len(grid) and np.all(data['rho'] > 0)
        pc, uc, _ = cold_state(data['rho'], model)
        assert np.all(data['p'] >= pc*(1-1e-10)), time
        thermal = data['total_energy']-0.5*data['rho']*data['vx']**2-uc
        assert np.all(thermal >= -1e-10*uc), time
        assert np.allclose(data['p'], pc+(model['gamma_thermal']-1)*thermal, rtol=1e-9)
        mass = np.sum(data['rho']*grid['volume'])
        history.append((time, data['rho'][0], mass, *shock_measurement(data, model)))
    history = np.asarray(history)
    times, central = history[:, 0], history[:, 1]
    peak = int(np.argmax(central))
    crossed = np.flatnonzero(central > model['rho_nuclear'])
    assert len(crossed), 'Center never entered stiff nuclear branch'
    assert 0 < peak < len(central)-1 and central[-1] < 0.9*central[peak], 'No resolved central bounce'
    assert times[peak] >= times[crossed[0]]
    mass_error = np.max(np.abs(history[:, 2]/history[0, 2]-1))
    assert mass_error < 1e-7, ('mass error', mass_error)
    # After the density maximum: downstream outward, upstream inward,
    # compressive pressure/density jumps, and supersonic upstream in shock frame.
    outgoing = ((times > times[peak]) & (history[:, 4] > 0) & (history[:, 5] < 0)
                & (history[:, 6] > 1.5) & (history[:, 7] > 2)
                & (history[:, 8] > 0) & (history[:, 9] > 1.1))
    shock_ids = np.flatnonzero(outgoing)
    assert len(shock_ids) >= 3, 'No sustained outgoing bounce shock'
    first, last = shock_ids[0], shock_ids[-1]
    assert history[last, 3]-history[first, 3] > 5e5, 'Shock has not propagated outward by 5 km'
    if expected_end is not None:
        assert abs(times[-1]-expected_end) < 1e-10*expected_end, 'Run did not reach requested end'
    report = {
        'frames': len(frames), 'cells': len(grid), 't_end_s': float(times[-1]),
        'initial_mass_Msun': float(history[0, 2]/1.98847e33),
        'nuclear_crossing_s': float(times[crossed[0]]), 'bounce_s': float(times[peak]),
        'peak_central_density_g_cm3': float(central[peak]),
        'max_relative_mass_error': float(mass_error),
        'first_outgoing_shock_s': float(times[first]),
        'first_outgoing_shock_radius_km': float(history[first, 3]/1e5),
        'last_shock_radius_km': float(history[last, 3]/1e5),
        'last_pressure_jump': float(history[last, 7]),
        'last_density_jump': float(history[last, 6]),
        'last_upstream_shock_frame_Mach': float(history[last, 9]),
    }
    if write:
        (directory/'validation.json').write_text(json.dumps(report, indent=2)+'\n')
        np.savetxt(directory/'history.csv', history, delimiter=',', comments='',
                   header='time,central_density,mass,shock_radius,downstream_velocity,upstream_velocity,density_jump,pressure_jump,shock_speed_RH,upstream_Mach')
    return model, frames, history, report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).parent/'output_bounce')
    parser.add_argument('--t-end', type=float, default=0.0104)
    args = parser.parse_args()
    print(json.dumps(analyze(args.output, args.t_end)[3], indent=2))
