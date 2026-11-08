#!/usr/bin/env python3
"""Check profiler CLI, byte-identical physics/output, and accounting on three cases.

Build sod_tube, File_Profile_2D and CCSN_MESA_day13 first. Optionally pass a directory
containing CASE.baseline executables saved before instrumentation. No existing case
outputs are overwritten. MESA uses the same mesh but a short 0.01-second evolution.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

TESTS = Path(__file__).resolve().parent
CASES = ('sod_tube', 'File_Profile_2D', 'CCSN_MESA_day13')


def run(binary, case, parameters, flags):
    start = time.perf_counter()
    result = subprocess.run([str(binary), '-par', str(parameters), *map(str, flags)],
                            cwd=case, capture_output=True)
    return result, time.perf_counter() - start


def validate_report(report, stdout, detail):
    assert report['status'] == 'completed'
    assert report['mode'] == ('detail' if detail else 'summary')
    sections = {s['name']: s for s in report['sections']}
    wall = report['wall_seconds']
    assert wall > 0
    assert abs(sum(s['self_seconds'] for s in sections.values()) - wall) < 1e-7 * max(1, wall)
    for s in sections.values():
        assert s['calls'] > 0 and 0 <= s['self_seconds'] <= s['inclusive_seconds'] + 1e-12
    steps = int(re.findall(rb'^step = (\d+)', stdout, re.M)[-1])
    counters = report['counters']
    assert counters['accepted_steps'] == steps
    assert sections['advance']['calls'] == steps
    assert sections['output']['calls'] == steps + 1
    if counters['rk2_attempts']:
        assert counters['rk2_attempts'] - counters['rk2_rejected_attempts'] == steps
    if detail:
        assert sections['riemann']['calls'] == sections['reconstruction']['calls']
        assert sections['riemann']['calls'] > steps
    else:
        assert 'riemann' not in sections and 'reconstruction' not in sections


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path)
    parser.add_argument('--report-dir', type=Path,
                        default=TESTS / 'initialization_validation/output/performance')
    args = parser.parse_args()
    report_dir = args.report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    records = {}
    with tempfile.TemporaryDirectory(prefix='flare-performance-check-') as temporary:
        scratch = Path(temporary)
        for name in CASES:
            case = TESTS / name
            binary = case / 'flare'
            assert binary.is_file(), f'Build {case} first'
            output = scratch / name
            parameters = scratch / f'{name}.txt'
            text = (case / 'parameters.txt').read_text()
            if name == 'CCSN_MESA_day13':
                text = re.sub(r'^t_end\s*=.*$', 't_end = 0.0000062780216922807006', text, flags=re.M)
                text = re.sub(r'^output_step_every\s*=.*$', 'output_step_every = 20000', text, flags=re.M)
            parameters.write_text(text + f'\noutput_directory = {output}\n')
            modes = [('off', binary, [])]
            if args.baseline_dir:
                modes.insert(0, ('baseline', args.baseline_dir.resolve() / f'{name}.baseline', []))
            modes += [('summary', binary, ['--performance-file', report_dir / f'{name}.summary.json']),
                      ('detail', binary, ['--performance-detail', '--performance-file', report_dir / f'{name}.detail.json'])]
            reference = None
            reference_stdout = None
            records[name] = {}
            for mode, executable, flags in modes:
                if output.exists():
                    shutil.rmtree(output)
                result, elapsed = run(executable, case, parameters, flags)
                assert result.returncode == 0, result.stderr.decode()
                files = {str(f.relative_to(output)): f.read_bytes() for f in output.rglob('*') if f.is_file()}
                assert files
                if reference is None:
                    reference, reference_stdout = files, result.stdout
                assert files == reference, f'{name}: physics files differ in {mode} mode'
                assert result.stdout == reference_stdout, f'{name}: stdout differs in {mode} mode'
                if mode in ('off', 'baseline'):
                    assert not result.stderr, result.stderr.decode()
                else:
                    report = json.loads((report_dir / f'{name}.{mode}.json').read_text())
                    validate_report(report, result.stdout, mode == 'detail')
                    assert b'FLARE performance' in result.stderr
                    (report_dir / f'{name}.{mode}.txt').write_bytes(result.stderr)
                records[name][mode] = {'external_wall_seconds': elapsed,
                                       'identical_files': len(files), 'identical_stdout': True}
                print(f'{name} {mode}: {elapsed:.3f}s, {len(files)} identical files', flush=True)
            records[name]['sha256'] = {f: hashlib.sha256(data).hexdigest() for f, data in reference.items()}

        # A simulation error must unwind timers and still emit a partial report.
        case = TESTS / 'sod_tube'
        failure = scratch / 'failure.txt'
        failure.write_text((case / 'parameters.txt').read_text() +
                           f'\nmaxsteps = 1\noutput_directory = {scratch / "failure"}\n')
        failed_json = report_dir / 'failed.json'
        result, _ = run(case / 'flare', case, failure, ['--performance', '--performance-file', failed_json])
        assert result.returncode != 0 and b'Maximum number of steps' in result.stderr
        failed = json.loads(failed_json.read_text())
        assert failed['status'] == 'failed' and failed['counters']['accepted_steps'] == 1
        assert failed['wall_seconds'] > 0
        for flags in (['--performance-file'], ['--unknown-performance-flag']):
            bad = subprocess.run([str(case / 'flare'), *flags], capture_output=True)
            assert bad.returncode != 0 and b'Error:' in bad.stderr
        help_result = subprocess.run([str(case / 'flare'), '--help'], capture_output=True)
        assert help_result.returncode == 0 and b'--performance-detail' in help_result.stdout
        assert not help_result.stderr

    (report_dir / 'validation.json').write_text(json.dumps(records, indent=2) + '\n')
    print('PASS: byte equality, accounting, failure reporting, and CLI validation')


if __name__ == '__main__':
    main()
