#!/usr/bin/env python3
"""Check all saved snapshots against periodic translation and EoS identities."""
import argparse
import json
from pathlib import Path

import numpy as np

from generate import exact_fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def snapshot(path):
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith('#'):
                header = line.strip()
                break
            key, value = line[1:].strip().split('=', 1)
            metadata[key] = value
        data = np.genfromtxt(stream, names=header, delimiter=',', ndmin=1)
    return metadata, data[data['is_ghost'] == 0]


def check(directory):
    model = json.loads((directory/'model.json').read_text())
    paths = sorted((directory/'output').glob('blob_*.csv'))
    require(len(paths) >= 3, f'{directory}: need initial, intermediate and final snapshots')
    n = model['n']
    dx = model['length']/n
    volume = dx*dx  # Integrated masses/energies are per unit length in the third direction.
    history = []
    initial_totals = None
    for path in paths:
        meta, data = snapshot(path)
        require(len(data) == n*n, f'{path}: incomplete grid')
        require(all(np.isfinite(data[name]).all() for name in data.dtype.names), f'{path}: nonfinite output')
        time = float(meta['time'])
        require(not history or time > history[-1]['time'], f'{path}: non-increasing time')
        x = data['comp_0']
        exact = exact_fraction(data['x'], data['y'], dx, dx, model, time)
        rho, pressure = data['rho'], data['p']
        totals = np.array([np.sum(rho), np.sum(rho*x), np.sum(rho*data['comp_1']),
                           np.sum(data['total_energy'])]) * volume
        if initial_totals is None:
            require(time == 0, 'Missing initial snapshot')
            initial_totals = totals
            require(np.max(abs(x-exact)) < 1e-11, 'Profile differs from analytic cell averages')
            # Check the analytic integral independently of the disk primitive.
            mean = 0.1+0.8*(35/128)**2 if model['shape'] == 'smooth' else np.pi*model['radius']**2
            require(abs(np.mean(x)-mean) < 1e-11, 'Incorrect initial species mass')
        conservation = abs(totals/initial_totals-1)
        mu = 1/(x+2*(1-x)) if model['eos'] == 'Ionized' else 1/(x+(1-x)/4)
        temperature = pressure*mu/(rho*(1.380649e-16/1.66053906892e-24))
        energy = pressure/(data['gamma']-1)
        if model['eos'] == 'Ionized':
            energy += rho*data['comp_1']/1.66053906892e-24*13.598435*1.602176634e-12
        row = dict(time=time, composition_L1=float(np.mean(abs(x-exact))),
                   composition_Linf=float(np.max(abs(x-exact))),
                   pressure_relative_error=float(np.max(abs(pressure/model['pressure']-1))),
                   density_relative_error=float(np.max(abs(rho/model['density']-1))),
                   velocity_error_scaled=float(max(np.max(abs(data['vx']-model['velocity'])),
                                                   np.max(abs(data['vy']-model['velocity'])))/1e6),
                   composition_sum_error=float(np.max(abs(x+data['comp_1']-1))),
                   composition_min=float(min(x.min(), data['comp_1'].min())),
                   composition_max=float(max(x.max(), data['comp_1'].max())),
                   mu_relative_error=float(np.max(abs(data['mean_molecular_weight']/mu-1))),
                   temperature_relative_error=float(np.max(abs(data['temperature']/temperature-1))),
                   energy_eos_relative_error=float(np.max(abs(data['internal_energy']/energy-1))),
                   mass_error=float(conservation[0]), species_mass_error=float(max(conservation[1:3])),
                   total_energy_error=float(conservation[3]),
                   peak_X0=float(x.max()), mixed_area_fraction=float(np.mean((x>0.1)&(x<0.9))),
                   mixing_measure=float(np.mean(x*(1-x))))
        for key in ('pressure_relative_error', 'density_relative_error', 'velocity_error_scaled',
                    'composition_sum_error', 'mu_relative_error', 'temperature_relative_error',
                    'energy_eos_relative_error', 'mass_error', 'species_mass_error', 'total_energy_error'):
            require(row[key] < 1e-8, f'{path.name}: {key}={row[key]}')
        require(row['composition_min'] >= -1e-12 and row['composition_max'] <= 1+1e-12,
                f'{path.name}: composition outside [0,1]')
        # Tolerances are for these supplied moderate-resolution transport cases,
        # not a claim of pointwise convergence across a discontinuity.
        limit = 1e-10 if model['case'] == 'stationary' and model['solver'] == 'exact' else 0.06
        require(row['composition_L1'] < limit, f'{path.name}: transport L1={row["composition_L1"]}')
        history.append(row)
    require(abs(history[-1]['time']-model['end_time']) < 1e-10, 'Simulation did not reach t_end')
    require(any(0.2*model['end_time'] < h['time'] < 0.8*model['end_time'] for h in history),
            'Need intermediate snapshots: a frozen profile could otherwise pass the full-period test')
    maxima = {key: max(row[key] for row in history) for key in history[0]
              if 'error' in key or key in ('composition_L1', 'composition_Linf')}
    result = dict(case=model['case'], n=n, solver=model['solver'], reconstructor=model['reconstructor'],
                  frames=len(paths), final=history[-1], maxima=maxima)
    np.savetxt(directory/'history.csv', [[row[key] for key in history[0]] for row in history],
               delimiter=',', header=','.join(history[0]), comments='')
    (directory/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directories', type=Path, nargs='+')
    args = parser.parse_args()
    for directory in args.directories:
        check(directory)


if __name__ == '__main__':
    main()
