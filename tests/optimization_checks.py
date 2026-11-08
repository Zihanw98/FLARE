#!/usr/bin/env python3
"""Compare legacy and optimized exact solves on three complete short runs.

Build sod_tube, File_Profile_2D, CCSN_MESA_day13 first. --baseline-dir optionally
contains CASE.baseline binaries saved before these optimizations. All run outputs
are temporary; only timing/error reports are saved. Run with no competing jobs.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import tempfile
import time

TESTS = Path(__file__).resolve().parent
CASES = ('sod_tube', 'File_Profile_2D', 'CCSN_MESA_day13')


def snapshot(data):
    lines = data.decode().splitlines()
    meta = dict(line[2:].split('=', 1) for line in lines if line.startswith('# '))
    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(
        line for line in lines if not line.startswith('#'))]
    return meta, [r for r in rows if not r['is_ghost']]


def final_snapshot(files):
    snapshots = [snapshot(data) for name, data in files.items()
                 if re.search(r'_\d{8}\.csv$', name)]
    return max(snapshots, key=lambda item: float(item[0]['time']))


def compare(reference, current):
    oldmeta, old = final_snapshot(reference)
    meta, new = final_snapshot(current)
    assert oldmeta['time'] == meta['time'] and len(old) == len(new)
    result = {}
    for key in ('rho', 'p', 'vx', 'vy', 'vz'):
        scale = max(abs(r[key]) for r in old)
        # Velocity errors use at least the maximum sound speed (avoids dividing by zero).
        if key.startswith('v'):
            scale = max(scale, max(math.sqrt(r['gamma'] * r['p'] / r['rho']) for r in old if r['rho'] > 0))
        diff = max(abs(a[key] - b[key]) for a, b in zip(old, new))
        normalized = diff / scale if scale else diff
        assert normalized < 2e-6, (key, normalized)
        result[key] = {'max_absolute': diff, 'max_scaled': normalized}
        if key in ('rho', 'p'):
            relative = max(abs(a[key] - b[key]) / a[key] for a, b in zip(old, new) if a[key] > 0)
            assert relative < 2e-5, (key, relative)
            result[key]['max_pointwise_relative'] = relative
    for a, b in zip(old, new):
        assert all(a[k] == b[k] for k in ('x', 'y', 'z', 'gamma'))
        assert b['rho'] >= 0 and b['p'] >= 0
        assert all(math.isfinite(v) for v in b.values())
    result['final_step'] = int(meta['step'])
    result['final_time'] = float(meta['time'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--report-dir', type=Path, default=TESTS / 'initialization_validation/output/optimization')
    args = parser.parse_args()
    assert args.repeats > 0
    report_dir = args.report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    with tempfile.TemporaryDirectory(prefix='flare-optimization-') as temp:
        scratch = Path(temp)
        for name in CASES:
            case = TESTS / name
            binary = case / 'flare'
            output = scratch / name
            par = scratch / f'{name}.txt'
            original = (case / 'parameters.txt').read_text()
            if name == 'CCSN_MESA_day13':
                original = re.sub(r'^t_end\s*=.*$', 't_end = 0.0000062780216922807006', original, flags=re.M)
                original = re.sub(r'^output_step_every\s*=.*$', 'output_step_every = 20000', original, flags=re.M)
            modes = [('bisection', binary, 'exact_bisection', 1), ('newton', binary, 'exact', args.repeats)]
            if args.baseline_dir:
                modes.insert(0, ('baseline', args.baseline_dir.resolve() / f'{name}.baseline', 'exact', args.repeats))
            reference = None
            records[name] = {}
            for mode, executable, method, repeats in modes:
                par.write_text(original + f'\nriemann_solver = {method}\noutput_directory = {output}\n')
                timings = []
                for repeat in range(repeats):
                    if output.exists(): shutil.rmtree(output)
                    report_path = report_dir / f'{name}.{mode}.json'
                    start = time.perf_counter()
                    process = subprocess.run([str(executable), '-par', str(par), '--performance-file', str(report_path)],
                                             cwd=case, capture_output=True)
                    timings.append(time.perf_counter() - start)
                    assert process.returncode == 0, process.stderr.decode()
                    files = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob('*') if p.is_file()}
                    if reference is None: reference = files
                    if mode in ('baseline', 'bisection'):
                        assert files == reference, f'{name}: bisection output lost byte compatibility'
                    if repeat == 0:
                        first = files
                        records[name][mode] = compare(reference, files)
                    else:
                        assert files == first, 'non-deterministic physics output'
                    report = json.loads(report_path.read_text())
                    assert report['status'] == 'completed'
                    if mode != 'baseline':
                        roots = report['exact_roots']
                        assert roots['solves'] > 0 and roots['failures'] == 0
                        assert sum(roots['iteration_histogram']) == roots['solves']
                        assert roots['initial_guesses'] + roots['newton_steps'] + roots['bisection_steps'] == roots['iterations']
                        records[name][mode]['roots'] = roots
                    (report_dir / f'{name}.{mode}.txt').write_bytes(process.stderr)
                records[name][mode]['wall_seconds'] = timings
                records[name][mode]['median_seconds'] = statistics.median(timings)
                print(name, mode, records[name][mode], flush=True)
            old = records[name].get('baseline', records[name]['bisection'])
            records[name]['speedup'] = old['median_seconds'] / records[name]['newton']['median_seconds']
            print(f"{name}: speedup {records[name]['speedup']:.3f}x", flush=True)
    (report_dir / 'validation.json').write_text(json.dumps(records, indent=2) + '\n')
    print('PASS: legacy byte equality, final-state tolerances, positive states, root accounting and reproducibility')


if __name__ == '__main__':
    main()
