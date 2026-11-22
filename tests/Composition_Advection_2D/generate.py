#!/usr/bin/env python3
"""Generate editable cell-average profiles and parameters for periodic advection."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def circle_primitive(x, y, radius):
    """Signed disk area between (0,0) and (x,y), evaluated analytically."""
    a, b = np.minimum(abs(x), radius), np.minimum(abs(y), radius)
    cut = np.sqrt(np.maximum(radius**2 - b**2, 0))

    def integral(s):
        return 0.5 * (s * np.sqrt(np.maximum(radius**2 - s**2, 0))
                      + radius**2 * np.arcsin(np.clip(s / radius, 0, 1)))

    area = np.where(a <= cut, a * b, cut * b + integral(a) - integral(cut))
    return np.sign(x) * np.sign(y) * area


def exact_fraction(x, y, dx, dy, model, time):
    """Exact cell average after periodic translation; x and y are center arrays."""
    length = model['length']
    cx = (model['center'][0] * length + model['velocity'] * time) % length
    cy = (model['center'][1] * length + model['velocity'] * time) % length
    if model['shape'] == 'smooth':
        # cos^8(pi*s) has this finite Fourier series. Integrating each mode
        # avoids comparing a finite-volume result with point-sampled data.
        def average(s, width):
            result = np.full_like(s, 35.0 / 128)
            for k, coefficient in enumerate((56, 28, 8, 1), 1):
                result += coefficient / 128 * np.cos(2*np.pi*k*s/length) * np.sinc(k*width/length)
            return result
        return 0.1 + 0.8 * average(x-cx, dx) * average(y-cy, dy)

    area = np.zeros_like(x + y)
    radius = model['radius'] * length
    for ix in (-1, 0, 1):
        for iy in (-1, 0, 1):
            left, right = x-dx/2-cx-ix*length, x+dx/2-cx-ix*length
            bottom, top = y-dy/2-cy-iy*length, y+dy/2-cy-iy*length
            area += (circle_primitive(right, top, radius) - circle_primitive(left, top, radius)
                     - circle_primitive(right, bottom, radius) + circle_primitive(left, bottom, radius))
    return np.clip(area / (dx*dy), 0, 1)


def generate(directory, case, n, reconstructor='PLM', solver='exact'):
    directory.mkdir(parents=True, exist_ok=True)
    if any((directory / name).exists() for name in ('profile.dat', 'parameters.txt', 'model.json', 'composition.txt')):
        raise FileExistsError(f'{directory} already has inputs; choose a new directory to preserve edits')
    ionized = case == 'ionized'
    model = dict(case=case, shape='smooth' if case == 'smooth' else 'circle', n=n,
                 length=1e8, density=1e-10, pressure=100.0, gamma=5/3,
                 velocity=0.0 if case == 'stationary' else 1e6,
                 center=[0.35, 0.5], radius=0.15, end_time=100.0,
                 eos='Ionized' if ionized else 'IdealGas_Mixed',
                 species=['HI', 'HII' if ionized else 'He_I'],
                 reconstructor=reconstructor, solver=solver)
    dx = model['length'] / n
    axis = (np.arange(n) + 0.5) * dx
    x, y = np.meshgrid(axis, axis)
    fraction = exact_fraction(x, y, dx, dx, model, 0).ravel()
    table = np.column_stack((np.full(n*n, model['density']),
                             np.full(n*n, model['velocity']), np.full(n*n, model['velocity']),
                             np.zeros(n*n), np.full(n*n, model['pressure']),
                             np.full(n*n, model['gamma']), fraction, 1-fraction))
    np.savetxt(directory/'profile.dat', table, fmt='%.17g', comments='',
               header='rho vx vy vz pre gamma ' + ' '.join(model['species']))
    (directory/'composition.txt').write_text('h1 0 1\n' if ionized else 'h1 0 0\nhe4 0 0\n')
    settings = {}
    for line in (HERE/'parameters.txt').read_text().splitlines():
        line = line.split('#')[0].strip()
        if line:
            key, value = line.split('=', 1)
            settings[key.strip()] = value.strip()
    settings.update(nxa=str(n), nya=str(n), equation_of_state=model['eos'],
                    reconstructor=reconstructor, riemann_solver=solver)
    # Species filenames resolve relative to parameters.txt; profile/output paths
    # resolve relative to the process working directory, so make those absolute.
    settings['init_profile_file_name'] = str((directory/'profile.dat').resolve())
    settings['output_directory'] = str((directory/'output').resolve())
    (directory/'parameters.txt').write_text(''.join(f'{k} = {v}\n' for k, v in settings.items()))
    (directory/'model.json').write_text(json.dumps(model, indent=2)+'\n')
    print(f'Generated {case}: {n} x {n}, {directory.resolve()}', flush=True)
    return model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('stationary', 'smooth', 'sharp', 'ionized'), default='sharp')
    parser.add_argument('--n', type=int, default=128)
    parser.add_argument('--directory', type=Path)
    parser.add_argument('--reconstructor', choices=('PCM', 'PLM'), default='PLM')
    parser.add_argument('--solver', choices=('exact', 'hll'), default='exact')
    args = parser.parse_args()
    if args.n < 8:
        parser.error('--n must be at least 8')
    generate(args.directory or HERE/f'output_{args.case}_{args.n}', args.case, args.n,
             args.reconstructor, args.solver)


if __name__ == '__main__':
    main()
