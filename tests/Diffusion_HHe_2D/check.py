#!/usr/bin/env python3
"""Check conservative H/He transport; all plotting remains in animate.py."""
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


def check(directory, parameters, control):
    cfg = {}
    for line in parameters.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if line:
            key, value = line.split('=', 1)
            cfg[key.strip()] = value.strip()
    paths = sorted(directory.glob(cfg['output_suffix'] + '_*.csv'))
    require(len(paths) >= 3, 'Need initial, intermediate and final snapshots')
    nx, ny = int(cfg['nxa']), int(cfg['nya'])
    lx = float(cfg['x_upper'])-float(cfg['x_lower'])
    ly = float(cfg['y_upper'])-float(cfg['y_lower'])
    volume = lx*ly/(nx*ny)  # cm^2: masses and energies are per unit depth.
    history = []
    for path in paths:
        meta, data = read_snapshot(path)
        require(len(data) == nx*ny, f'{path}: incomplete grid')
        require(len(np.unique(data['x'])) == nx and len(np.unique(data['y'])) == ny,
                f'{path}: unexpected mesh')
        require(all(np.isfinite(data[k]).all() for k in data.dtype.names), f'{path}: nonfinite state')
        require(meta.get('comp_0_name') == 'HI' and meta.get('comp_1_name') == 'He_I',
                f'{path}: expected H I / He I')
        require(all(float(meta[f'unit_{k}']) == 1 for k in ('mass', 'length', 'time', 'temperature')),
                f'{path}: this checker expects CGS')
        time = float(meta['time'])
        require(not history or time > history[-1]['time'], f'{path}: times must increase')
        rho, p = data['rho'], data['p']
        xh, xhe = data['comp_0'], data['comp_1']
        mass = rho*volume
        totals = np.array([mass.sum(), (mass*xh).sum(), (mass*xhe).sum(),
                           data['total_energy'].sum()*volume])
        momentum = np.array([(mass*data['vx']).sum(), (mass*data['vy']).sum()])
        require(np.all(rho > 0) and np.all(p > 0), f'{path}: nonpositive rho or P')
        if not history:
            require(time == 0, 'Missing initial snapshot')
            initial = data.copy()
            initial_totals, initial_momentum = totals.copy(), momentum.copy()
            speed_scale = np.sqrt(float(cfg['gamma'])*p[0]/rho[0])
            require(np.max(abs(p/p[0]-1)) < 1e-12, 'Initial P must be uniform')
            require(np.max(abs(rho/rho[0]-1)) < 1e-12, 'Initial rho must be uniform')
            require(np.max(np.hypot(data['vx'], data['vy'])) == 0, 'Initial velocity must vanish')
            x = (data['x']-float(cfg['x_lower']))/lx-0.5
            y = (data['y']-float(cfg['y_lower']))/ly-0.5
            outside = np.hypot(x, y) > 0.20
        require(np.array_equal(data['x'], initial['x']) and np.array_equal(data['y'], initial['y']),
                f'{path}: grid changed')
        errors = abs(totals/initial_totals-1)
        mu = 1/(xh+xhe/4)
        temperature = p*mu/(rho*(1.380649e-16/1.66053906892e-24))
        row = dict(time=time, mass=float(totals[0]), hydrogen_mass=float(totals[1]),
                   helium_mass=float(totals[2]), total_energy=float(totals[3]),
                   mass_error=float(errors[0]), hydrogen_mass_error=float(errors[1]),
                   helium_mass_error=float(errors[2]), energy_error=float(errors[3]),
                   momentum_error_scaled=float(np.max(abs(momentum-initial_momentum))
                                               /(initial_totals[0]*speed_scale)),
                   composition_sum_error=float(np.max(abs(xh+xhe-1))),
                   mu_error=float(np.max(abs(data['mean_molecular_weight']/mu-1))),
                   temperature_error=float(np.max(abs(data['temperature']/temperature-1))),
                   composition_min=float(min(xh.min(), xhe.min())),
                   composition_max=float(max(xh.max(), xhe.max())),
                   hydrogen_peak=float(xh.max()),
                   hydrogen_outside_initial_circle=float((mass[outside]*xh[outside]).sum()/totals[1]),
                   mixing_measure=float((mass*4*xh*xhe).sum()/totals[0]),
                   max_hydrogen_change=float(np.max(abs(xh-initial['comp_0']))),
                   max_speed=float(np.max(np.hypot(data['vx'], data['vy']))),
                   dual_energy_correction_scaled=float(np.sum(abs(data['dual_energy_correction']))
                                                       *volume/initial_totals[3]))
        for key, value in row.items():
            if 'error' in key or key == 'dual_energy_correction_scaled':
                require(value < 1e-8, f'{path.name}: {key}={value}')
        require(row['composition_min'] >= -1e-12 and row['composition_max'] <= 1+1e-12,
                f'{path.name}: composition outside [0,1]')
        history.append(row)
    require(abs(history[-1]['time']-float(cfg['t_end'])) < 1e-9, 'Did not reach t_end')
    if control:
        require(max(r['max_hydrogen_change'] for r in history) < 1e-10, 'D=0 control must stay stationary')
        require(max(r['max_speed'] for r in history)/speed_scale < 1e-10, 'D=0 control developed flow')
    else:
        require(history[-1]['mixing_measure'] > history[0]['mixing_measure']+0.05, 'No clear mixing')
        require(history[-1]['hydrogen_outside_initial_circle'] >
                history[0]['hydrogen_outside_initial_circle']+0.10, 'No clear outward H transport')
        require(history[-1]['hydrogen_peak'] < history[0]['hydrogen_peak']-0.10, 'H peak did not diffuse')
    result = dict(control=control, integrator=cfg['time_integrator'], frames=len(paths),
                  grid=[nx, ny], initial=history[0], final=history[-1],
                  maxima={k: max(r[k] for r in history) for k in history[0]
                          if 'error' in k or k == 'dual_energy_correction_scaled'})
    np.savetxt(directory/'history.csv', [[r[k] for k in history[0]] for r in history],
               delimiter=',', header=','.join(history[0]), comments='')
    (directory/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?', type=Path, default=HERE/'output')
    parser.add_argument('--parameters', type=Path, default=HERE/'parameters.txt')
    parser.add_argument('--control', action='store_true', help='Expect stationary D=0 output')
    args = parser.parse_args()
    check(args.directory, args.parameters, args.control)
