#!/usr/bin/env python3
"""Animate saved heating snapshots with the analytic pressure solution."""
import argparse
import csv
import io
import os
from pathlib import Path
import tempfile

from check import analyze, read_parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parameters', type=Path, default=Path(__file__).with_name('parameters.txt'))
    args = parser.parse_args()
    parameters = args.parameters.resolve()
    analyze(parameters)  # Reject stale or incomplete outputs before plotting.
    root = parameters.parent
    cfg = read_parameters(parameters)
    custom = read_parameters(root / 'heating.txt')
    rate = float(custom.get('specific_heating_rate', 0.5))
    if cfg.get('constant_heating_enabled', 'false').lower() in ('true', '1'):
        rate += float(cfg.get('heating_rate', 0.0))
    gamma = float(cfg['gamma'])
    snapshots = []
    for line in (root / 'run.log').read_text().splitlines():
        if ', csv = ' not in line:
            continue
        with (root / line.split(', csv = ', 1)[1]).open() as stream:
            metadata = dict(next(stream).strip()[2:].split('=', 1) for _ in range(3))
            rows = [row for row in csv.DictReader(line for line in stream if not line.lstrip().startswith('#'))
                    if row['is_ghost'] == '0']
        snapshots.append((float(metadata['time']), int(metadata['step']),
                          [float(row['x']) for row in rows], [float(row['p']) for row in rows]))
    output = root / cfg.get('output_directory', 'output') / 'animations'
    output.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix='flare-heating-plot-') as cache:
        os.environ.setdefault('MPLCONFIGDIR', cache)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from PIL import Image

        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                             'axes.labelcolor': '#314358', 'text.color': '#162c44',
                             'xtick.color': '#53657a', 'ytick.color': '#53657a'})
        fig, (space, history) = plt.subplots(1, 2, figsize=(11, 5.4), dpi=110)
        fig.patch.set_facecolor('#f5f7fa')
        fig.subplots_adjust(left=0.08, right=0.97, bottom=0.19, top=0.72, wspace=0.24)
        fig.text(0.08, 0.92, 'Uniform heating | 1D', fontsize=23, weight='bold')
        fig.text(0.08, 0.855, f'q = {rate:g}    |    gamma = {gamma:g}    |    '
                 f'{cfg["time_integrator"].upper()}    |    density = 1, velocity = 0', color='#53657a')
        status = fig.text(0.08, 0.79, '', fontsize=12, weight='bold')
        fig.text(0.08, 0.07, 'Saved simulation snapshots; no interpolation.  '
                 'Analytic solution: p(t) = 1 + (gamma - 1) q t', fontsize=10, color='#53657a')
        end = snapshots[-1][0]
        final_pressure = 1 + (gamma - 1) * rate * end
        pad = max(0.05, (final_pressure - 1) * 0.15)
        for axis in (space, history):
            axis.set_facecolor('white')
            axis.grid(color='#e2e8ef', linewidth=0.8)
            axis.set_axisbelow(True)
            axis.set_ylim(1 - pad, final_pressure + pad)
            axis.set_ylabel('Pressure')
            for spine in axis.spines.values():
                spine.set_color('#d3dce6')
        space.set(xlabel='Position x', xlim=(float(cfg['x_lower']), float(cfg['x_upper'])),
                  title='Spatial pressure profile')
        history.set(xlabel='Time', xlim=(-0.02 * end, 1.02 * end), title='Pressure rise over time')
        analytic_space, = space.plot([], [], '--', color='#de8735', lw=2, label='Analytic')
        numeric_space, = space.plot([], [], 'o', color='#087f8c', ms=4, label='Simulation')
        history.plot([0, end], [1, final_pressure], '--', color='#de8735', lw=2, label='Analytic')
        trajectory, = history.plot([], [], 'o-', color='#087f8c', ms=4, lw=1.5, label='Simulation')
        current, = history.plot([], [], 'o', color='#087f8c', ms=9)
        cursor = history.axvline(0, color='#a4b4c5', lw=1)
        for axis in (space, history):
            axis.legend(loc='upper left', frameon=False, fontsize=9, ncol=2)
        frames, durations = [], []
        for i, (time, step, x, pressure) in enumerate(snapshots):
            expected = 1 + (gamma - 1) * rate * time
            analytic_space.set_data([float(cfg['x_lower']), float(cfg['x_upper'])], [expected, expected])
            numeric_space.set_data(x, pressure)
            trajectory.set_data([s[0] for s in snapshots[:i + 1]],
                                [sum(s[3]) / len(s[3]) for s in snapshots[:i + 1]])
            current.set_data([time], [sum(pressure) / len(pressure)])
            cursor.set_xdata([time, time])
            status.set_text(f't = {time:.3f}     |     step {step}     |     p = {expected:.4f}')
            buffer = io.BytesIO()
            fig.savefig(buffer, format='png', facecolor=fig.get_facecolor())
            buffer.seek(0)
            with Image.open(buffer) as frame:
                frames.append(frame.convert('RGB'))
            durations.append(max(60, round(3500 * (snapshots[i + 1][0] - time)))
                             if i + 1 < len(snapshots) else 1400)
        target = output / 'heating_analytic.gif'
        frames[0].save(target, save_all=True, append_images=frames[1:], duration=durations,
                       loop=0, disposal=2)
        fig.savefig(output / 'heating_analytic_last.png', facecolor=fig.get_facecolor())
        plt.close(fig)
        print(f'Saved {len(frames)} frames: {target}')


if __name__ == '__main__':
    main()
