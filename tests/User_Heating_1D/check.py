#!/usr/bin/env python3
"""Check uniform user heating against its analytic solution; standard library only."""
import argparse
import csv
import math
from pathlib import Path
import re
import subprocess
import tempfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_parameters(path):
    values = {}
    for line in path.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if line:
            key, value = line.split('=', 1)
            values[key.strip()] = value.strip()
    return values


def close(actual, expected, message, tolerance=2e-11):
    require(math.isfinite(actual) and math.isclose(actual, expected,
        rel_tol=tolerance, abs_tol=tolerance), f'{message}: {actual} != {expected}')


def run(executable, parameters):
    with (parameters.parent / 'run.log').open('w') as log:
        subprocess.run([str(executable), '-par', parameters.name],
                       cwd=parameters.parent, stdout=log, stderr=subprocess.STDOUT, check=True)


def analyze(parameters, expect_source_limited=False):
    root = parameters.parent
    cfg = read_parameters(parameters)
    custom = read_parameters(root / 'heating.txt')
    user_rate = float(custom.get('specific_heating_rate', 0.5))
    fraction = float(custom.get('max_fractional_energy_change', 0.1))
    builtin_rate = (float(cfg.get('heating_rate', 0.0))
                    if cfg.get('constant_heating_enabled', 'false').lower() in ('true', '1') else 0.0)
    rate = user_rate + builtin_rate
    # Independent sources each impose their own cap on the same internal energy.
    timescale_factor = min((f / q for q, f in ((user_rate, fraction), (builtin_rate, 0.1))
                            if q > 0.0), default=math.inf)
    gamma, end = float(cfg['gamma']), float(cfg['t_end'])
    length = float(cfg['x_upper']) - float(cfg['x_lower'])
    count, ghost = int(cfg['nxa']), int(cfg['nxg'])
    initial_energy = 1.0 / (gamma - 1.0)
    pattern = re.compile(r'^step = (\d+), time = ([^,]+), dt = ([^,]+)')
    records = []
    metrics = []
    for line in (root / 'run.log').read_text().splitlines():
        match = pattern.match(line)
        if not match:
            continue
        step, time, dt = int(match[1]), float(match[2]), float(match[3])
        if records:
            previous = records[-1]
            require(step == previous[0] + 1 and time > previous[1], 'Missing/unordered step')
            close(time, previous[1] + dt, 'Time increment', tolerance=2e-9)
            if rate > 0:
                cap = timescale_factor * (initial_energy + rate * previous[1])
                require(dt <= cap * (1 + 2e-9), 'Source timestep bound exceeded')
        else:
            require((step, time, dt) == (0, 0.0, 0.0), 'Missing initial record')
        records.append((step, time, dt))
        if ', csv = ' not in line:
            continue
        path = root / line.split(', csv = ', 1)[1]
        with path.open() as stream:
            metadata = dict(next(stream).strip()[2:].split('=', 1) for _ in range(3))
            rows = list(csv.DictReader(line for line in stream if not line.lstrip().startswith('#')))
        require(int(metadata['step']) == step, 'Snapshot step mismatch')
        snapshot_time = float(metadata['time'])
        close(snapshot_time, time, 'Snapshot time mismatch', tolerance=2e-9)
        close(float(metadata['dt']), dt, 'Snapshot dt mismatch', tolerance=2e-9)
        require(len(rows) == count + 2 * ghost, 'Incomplete snapshot')
        expected_pressure = 1.0 + (gamma - 1.0) * rate * snapshot_time
        max_error = 0.0
        total_energy = 0.0
        mass = 0.0
        for i, row in enumerate(rows):
            values = {key: float(value) for key, value in row.items()}
            require(all(math.isfinite(value) for value in values.values()), 'Nonfinite state')
            close(values['x'], float(cfg['x_lower']) + (i - ghost + 0.5) * length / count,
                  'Cell center')
            close(values['gamma'], gamma, 'Gamma')
            close(values['rho'], 1.0, 'Density changed')
            for key in ('vx', 'vy', 'vz'):
                close(values[key], 0.0, 'Spurious velocity')
            close(values['p'], expected_pressure, 'Pressure differs from analytic heating')
            max_error = max(max_error, abs(values['p'] - expected_pressure))
            active = ghost <= i < ghost + count
            require(values['is_ghost'] == int(not active), 'Incorrect ghost flag')
            if active:
                volume = length / count
                mass += values['rho'] * volume
                total_energy += (values['p'] / (gamma - 1.0)
                    + 0.5 * values['rho'] * values['vx'] ** 2) * volume
        close(mass, length, 'Mass conservation')
        close(total_energy, length * (initial_energy + rate * snapshot_time), 'Energy budget')
        metrics.append({'step': step, 'time': snapshot_time, 'analytic_pressure': expected_pressure,
                        'max_pressure_error': max_error, 'mass': mass, 'total_energy': total_energy})
    require(len(records) >= 2 and len(metrics) >= 3, 'Insufficient output')
    close(records[-1][1], end, 'Simulation did not reach t_end')
    close(metrics[-1]['time'], end, 'Missing final snapshot')
    if expect_source_limited:
        close(records[1][2], timescale_factor * initial_energy,
              'First step should be limited by heating', tolerance=2e-9)
    diagnostics = root / cfg.get('output_directory', 'output') / 'diagnostics'
    diagnostics.mkdir(parents=True, exist_ok=True)
    with (diagnostics / 'metrics.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    error = max(item['max_pressure_error'] for item in metrics)
    print(f"PASS {cfg['time_integrator']}: q={rate:g}, {records[-1][0]} steps, "
          f"{len(metrics)} snapshots, final p={metrics[-1]['analytic_pressure']:.12g}, "
          f"max pressure error={error:.3e}")


def suite(executable, parameters):
    baseline = read_parameters(parameters)
    with tempfile.TemporaryDirectory(prefix='flare-user-heating-') as directory:
        root = Path(directory)
        for integrator in ('euler', 'rk2'):
            for rate in (0.0, 0.5, 1.5):
                case = root / f'{integrator}-{rate}'
                case.mkdir()
                cfg = dict(baseline, time_integrator=integrator, output_directory='output',
                           constant_heating_enabled='false')
                path = case / 'parameters.txt'
                path.write_text(''.join(f'{key} = {value}\n' for key, value in cfg.items()))
                fraction = 0.001 if rate == 1.5 else 0.1
                (case / 'heating.txt').write_text(f'specific_heating_rate = {rate}\n'
                    f'max_fractional_energy_change = {fraction}\n')
                run(executable, path)
                analyze(path, expect_source_limited=rate == 1.5)
        # Same physics and timestep safety factor: built-in and user implementations
        # must produce identical snapshots/logs, including source-limited heating.
        for integrator in ('euler', 'rk2'):
            for rate in (0.5, 500.0):
                pair = []
                for builtin in (False, True):
                    case = root / f'compare-{integrator}-{rate}-{builtin}'
                    case.mkdir()
                    cfg = dict(baseline, time_integrator=integrator, t_end='0.02',
                               output_directory='output', output_step_every='1',
                               constant_heating_enabled=str(builtin).lower(),
                               heating_rate=str(rate))
                    path = case / 'parameters.txt'
                    path.write_text(''.join(f'{key} = {value}\n' for key, value in cfg.items()))
                    (case / 'heating.txt').write_text(
                        f'specific_heating_rate = {0.0 if builtin else rate}\n'
                        'max_fractional_energy_change = 0.1\n')
                    run(executable, path)
                    analyze(path, expect_source_limited=rate == 500.0)
                    pair.append(case)
                files = sorted(p.relative_to(pair[0]) for p in (pair[0] / 'output').glob('*.csv'))
                require(files == sorted(p.relative_to(pair[1])
                    for p in (pair[1] / 'output').glob('*.csv')), 'Built-in snapshot list differs')
                for name in files + [Path('run.log')]:
                    require((pair[0] / name).read_bytes() == (pair[1] / name).read_bytes(),
                            f'Built-in/user heating output differs: {name}')
                print(f'PASS built-in/user byte comparison: {integrator}, q={rate:g}, '
                      f'{len(files)} snapshots and run.log identical')
        for index, (settings, message) in enumerate([
            ('specific_heating_rate = -1\n', 'specific_heating_rate must be nonnegative'),
            ('max_fractional_energy_change = 0\n', 'max_fractional_energy_change must be'),
            ('heating_typo = 1\n', 'Unknown parameter'),
        ]):
            case = root / f'invalid-{index}'
            case.mkdir()
            (case / 'parameters.txt').write_text(parameters.read_text())
            (case / 'heating.txt').write_text(settings)
            result = subprocess.run([str(executable), '-par', 'parameters.txt'], cwd=case,
                                    capture_output=True, text=True)
            require(result.returncode != 0 and message in result.stderr,
                    f'Invalid user parameters were not rejected: {settings}')
        print('PASS user parameter validation: negative rate, zero fraction, unknown key')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parameters', type=Path, default=Path(__file__).with_name('parameters.txt'))
    parser.add_argument('--executable', type=Path, default=Path(__file__).with_name('flare'))
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--suite', action='store_true')
    args = parser.parse_args()
    parameters, executable = args.parameters.resolve(), args.executable.resolve()
    if args.run:
        run(executable, parameters)
    analyze(parameters)
    if args.suite:
        suite(executable, parameters)


if __name__ == '__main__':
    main()
