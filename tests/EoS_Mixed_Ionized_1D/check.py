#!/usr/bin/env python3
"""Run the supplied cases in a disposable directory and check their invariants."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def snapshot(path):
    lines = path.read_text().splitlines()
    metadata = dict(line[2:].split('=', 1) for line in lines if line.startswith('# '))
    rows = [{k: float(v) for k, v in row.items()}
            for row in csv.DictReader(line for line in lines if not line.startswith('#'))]
    return metadata, [row for row in rows if row['is_ghost'] == 0]


def check(executable, root, name, method):
    directory = root / (name + '_' + method)
    directory.mkdir(parents=True, exist_ok=True)
    text = (HERE / ('parameters_' + name + '.txt')).read_text()
    comp = 'mixed' if name.startswith('mixed') else 'ionized'
    text += f'\ncomposition_activate_file = {HERE / ("composition_" + comp + ".txt")}\n'
    text += f'time_integrator = {method}\noutput_directory = {directory / "output"}\n'
    parameters = directory / 'parameters.txt'
    parameters.write_text(text)
    with (directory / 'run.log').open('w') as log:
        subprocess.run([str(executable), '-par', str(parameters)], cwd=directory,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    paths = sorted((directory / 'output').glob(name + '_*.csv'))
    meta0, initial = snapshot(paths[0])
    meta, final = snapshot(paths[-1])
    heating = 'heating' in name
    time = float(meta['time'])
    assert abs(time - (0.2 if heating else 200)) < 1e-10
    quantities = {
        'mass': lambda row: row['rho'],
        'H_or_HI': lambda row: row['rho'] * row['comp_0'],
        'He_or_HII': lambda row: row['rho'] * row['comp_1'],
        'energy': lambda row: row['total_energy'],
    }
    errors = {}
    for label, value in quantities.items():
        before = math.fsum(map(value, initial))
        expected = before
        if heating and label == 'energy':
            expected += math.fsum(row['rho'] for row in initial) * 1e13 * time
        after = math.fsum(map(value, final))
        errors[label] = abs(after/expected-1)
        assert errors[label] < 2e-8, (name, method, label, errors[label])
    for row in final:
        assert all(math.isfinite(v) for v in row.values())
        assert abs(row['comp_0']+row['comp_1']-1) < 1e-9
        assert 0 <= row['comp_0'] <= 1 and 0 <= row['comp_1'] <= 1
        if not heating:
            assert abs(row['p']/100-1) < 2e-8
            assert abs(row['vx']/1e5-1) < 2e-8
    if heating:
        expected_p = initial[0]['p'] + (initial[0]['gamma']-1)*1e-10*1e13*time
        assert max(abs(row['p']/expected_p-1) for row in final) < 2e-8
        if method != 'implicit':
            first_step, _ = snapshot(paths[1])
            cap = 0.1*initial[0]['internal_energy']/(initial[0]['rho']*1e13)
            assert abs(float(first_step['dt'])/cap-1) < 1e-9
    else:
        # Smooth periodic contact: allow PCM diffusion, but detect a frozen or
        # incorrectly transported composition pattern after a fifth of a period.
        errors['composition_L1'] = math.fsum(abs(row['comp_0']-(0.5+0.3*math.sin(
            2*math.pi*(row['x']-1e5*time)/1e8))) for row in final)/len(final)
        assert errors['composition_L1'] < 0.025
    thermo = list(csv.DictReader((directory/'output/thermodynamics_last.csv').open()))
    assert len(thermo) == len(final)
    for row, cache in zip(final, thermo):
        x = row['comp_0']
        mu = 1/(x+(1-x)/4) if comp == 'mixed' else 1/(x+2*(1-x))
        expected_t = row['p']*mu/(row['rho']*(1.380649e-16/1.66053906892e-24))
        assert abs(float(cache['temperature'])/expected_t-1) < 2e-8
        assert abs(float(cache['mu'])/mu-1) < 1e-9
        if comp == 'mixed':
            assert float(cache['n_e']) == 0
        else:
            expected_ne = row['rho']*row['comp_1']/1.66053906892e-24
            assert abs(float(cache['n_e'])/expected_ne-1) < 2e-8
    # Ionization is frozen, so its specific contribution must follow HII.
    ion_ev = 13.598435  # tables/species.json, HII relative to neutral HI.
    for row in initial + final:
        expected = row['p']/(row['gamma']-1)
        if comp == 'ionized':
            expected += row['rho']*row['comp_1']/1.66053906892e-24*ion_ev*1.602176634e-12
        assert abs(row['internal_energy']/expected-1) < 2e-8
    result = {'case': name, 'integrator': method, 'steps': int(meta['step']), **errors}
    print(json.dumps(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, default=HERE/'flare')
    parser.add_argument('--output-root', type=Path)
    args = parser.parse_args()
    root = args.output_root or Path(tempfile.mkdtemp(prefix='flare-eos-cases-'))
    root = root.resolve()
    results = []
    for name in ['mixed_contact', 'ionized_contact']:
        for method in ['rk2', 'lts_rk2']:
            results.append(check(args.exe.resolve(), root, name, method))
    for method in ['euler', 'rk2', 'lts_euler', 'lts_rk2', 'implicit']:
        results.append(check(args.exe.resolve(), root, 'ionized_heating', method))
    (root/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    print('All EoS cases passed:', root)


if __name__ == '__main__':
    main()
