#!/usr/bin/env python3
"""Check conservation, thermodynamic output, and visible H/He entrainment."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_snapshot(path):
    metadata = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith('#'):
                header = line.strip()
                break
            key, value = line[1:].strip().split('=', 1)
            metadata[key] = value
        data = np.genfromtxt(stream, names=header, delimiter=',', ndmin=1)
    data = data[data['is_ghost'] == 0]
    return metadata, data[np.lexsort((data['x'], data['y']))]


def check(directory, parameters):
    cfg = {}
    for line in parameters.read_text().splitlines():
        line = line.split('#')[0].strip()
        if line:
            key, value = line.split('=', 1)
            cfg[key.strip()] = value.strip()
    paths = sorted(directory.glob(cfg['output_suffix'] + '_*.csv'))
    require(len(paths) >= 3, 'Need initial, intermediate and final snapshots')
    nx, ny = int(cfg['nxa']), int(cfg['nya'])
    lx = float(cfg['x_upper']) - float(cfg['x_lower'])
    ly = float(cfg['y_upper']) - float(cfg['y_lower'])
    volume = lx*ly/(nx*ny)  # Per unit depth in 2D.
    history = []
    initial_totals = None
    for path in paths:
        meta, data = read_snapshot(path)
        require(len(data) == nx*ny, f'{path}: incomplete grid')
        require(len(np.unique(data['x'])) == nx and len(np.unique(data['y'])) == ny,
                f'{path}: unexpected mesh')
        require(all(np.isfinite(data[name]).all() for name in data.dtype.names), f'{path}: nonfinite state')
        require(meta.get('comp_0_name') == 'HI' and meta.get('comp_1_name') == 'He_I',
                f'{path}: expected neutral H/He species')
        time = float(meta['time'])
        require(not history or time > history[-1]['time'], f'{path}: times must increase')
        xh, xhe = data['comp_0'], data['comp_1']
        rho, p = data['rho'], data['p']
        require(np.all(rho > 0) and np.all(p > 0), f'{path}: nonpositive rho or P')
        mass = rho*volume
        totals = np.array([mass.sum(), (mass*xh).sum(), (mass*xhe).sum(),
                           data['total_energy'].sum()*volume])
        momentum = np.array([(mass*data['vx']).sum(), (mass*data['vy']).sum()])
        if initial_totals is None:
            require(time == 0, 'Missing initial snapshot')
            initial_totals = totals
            initial_momentum = momentum
            require(np.max(abs(p/p[0]-1)) < 1e-12, 'Initial pressure is not uniform')
            require(np.max(abs(rho/rho[0]-1)) < 1e-12, 'Initial density is not uniform')
            require(abs(totals[1]/totals[0]-0.5) < 1e-12, 'Initial H band should contain half the mass')
        relative_error = abs(totals/initial_totals-1)
        mu = 1/(xh+xhe/4)
        temperature = p*mu/(rho*(1.380649e-16/1.66053906892e-24))
        yn = (data['y']-float(cfg['y_lower']))/ly
        outside = (yn < 0.25) | (yn > 0.75)
        fraction_grid = xh.reshape(ny, nx)
        row = dict(time=time, mass_error=float(relative_error[0]),
                   hydrogen_mass_error=float(relative_error[1]), helium_mass_error=float(relative_error[2]),
                   energy_error=float(relative_error[3]),
                   momentum_error_scaled=float(np.max(abs(momentum-initial_momentum))/(totals[0]*5e5)),
                   composition_sum_error=float(np.max(abs(xh+xhe-1))),
                   composition_min=float(min(xh.min(), xhe.min())),
                   composition_max=float(max(xh.max(), xhe.max())),
                   mu_error=float(np.max(abs(data['mean_molecular_weight']/mu-1))),
                   temperature_error=float(np.max(abs(data['temperature']/temperature-1))),
                   rho_min=float(rho.min()), p_min=float(p.min()),
                   vy_rms=float(np.sqrt((mass*data['vy']**2).sum()/mass.sum())),
                   hydrogen_outside_initial_band=float((mass[outside]*xh[outside]).sum()/totals[1]),
                   helium_inside_initial_band=float((mass[~outside]*xhe[~outside]).sum()/totals[2]),
                   composition_x_rms=float(np.sqrt(np.mean((fraction_grid-fraction_grid.mean(axis=1)[:, None])**2))),
                   mixed_mass_fraction=float(mass[(xh>0.1)&(xh<0.9)].sum()/mass.sum()),
                   mixing_measure=float((mass*4*xh*xhe).sum()/mass.sum()),
                   dual_energy_correction=float(data['dual_energy_correction'].sum()*volume))
        for key in ('mass_error', 'hydrogen_mass_error', 'helium_mass_error', 'energy_error',
                    'momentum_error_scaled', 'composition_sum_error', 'mu_error', 'temperature_error'):
            require(row[key] < 1e-8, f'{path.name}: {key}={row[key]}')
        require(row['composition_min'] >= -1e-12 and row['composition_max'] <= 1+1e-12,
                f'{path.name}: composition outside [0,1]')
        history.append(row)
    require(abs(history[-1]['time']-float(cfg['t_end'])) < 1e-9, 'Did not reach t_end')
    require(max(r['vy_rms'] for r in history) > 3*history[0]['vy_rms'], 'No clear growth of transverse flow')
    require(max(r['hydrogen_outside_initial_band'] for r in history) > 0.02, 'No clear H entrainment')
    require(max(r['composition_x_rms'] for r in history) > 0.02, 'No clear deformation of the H/He layers')
    # These last checks establish a visible KH demonstration, not convergence
    # or agreement with an analytic KH growth rate or physical diffusivity.
    result = dict(frames=len(paths), grid=[nx, ny], final=history[-1],
                  maxima={key: max(row[key] for row in history) for key in history[0] if 'error' in key},
                  max_vy_rms=max(r['vy_rms'] for r in history), initial_vy_rms=history[0]['vy_rms'])
    np.savetxt(directory/'history.csv', [[row[key] for key in history[0]] for row in history],
               delimiter=',', header=','.join(history[0]), comments='')
    (directory/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?', type=Path, default=HERE/'output')
    parser.add_argument('--parameters', type=Path, default=HERE/'parameters.txt')
    args = parser.parse_args()
    check(args.directory, args.parameters)


if __name__ == '__main__':
    main()
