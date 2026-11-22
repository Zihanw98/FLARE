#!/usr/bin/env python3
"""Animate validated spherical collapse snapshots; no interpolation of states."""
import argparse
import os
import tempfile
from pathlib import Path
import numpy as np
from analyze import analyze, cold_state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).parent/'output_bounce')
    parser.add_argument('--frames', type=int, default=160)
    args = parser.parse_args()
    model, frames, history, report = analyze(args.output)
    destination = args.output/'animations'
    destination.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='flare-bounce-mpl-') as cache:
        os.environ.setdefault('MPLCONFIGDIR', cache)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from PIL import Image
        plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                             'axes.spines.top': False, 'axes.spines.right': False})
        fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), dpi=100)
        fig.patch.set_facecolor('#f5f7fb')
        fig.subplots_adjust(left=.085, right=.97, bottom=.11, top=.82, wspace=.26, hspace=.37)
        fig.text(.085, .94, 'Core collapse → nuclear bounce → outgoing shock', fontsize=20, weight='bold')
        fig.text(.085, .897, '1D spherical · Newtonian self-gravity · no composition · hybrid EoS', color='#52617a')
        clock = fig.text(.085, .853, '', color='#1a3760', weight='bold')
        fig.text(.085, .025, 'Saved simulation states; playback slows near bounce.  '
                 'Low-density gamma = 4/3; stiff gamma = 2.5; thermal gamma = 1.5.', fontsize=9, color='#52617a')
        density, velocity, fraction, evolution = axes.flat
        for ax in axes.flat:
            ax.grid(alpha=.18)
        density.set(xscale='log', yscale='log', xlim=(.12, 500), ylim=(1e5, 4e14),
                    xlabel='Radius [km]', ylabel='Density [g cm$^{-3}$]')
        density.axhline(model['rho_nuclear'], color='#c05e33', ls='--', lw=1.2, label='Nuclear threshold')
        velocity.set(xlim=(0, 180), ylim=(-9, 7), xlabel='Radius [km]',
                     ylabel='Radial velocity [10$^9$ cm s$^{-1}$]')
        velocity.axhline(0, color='#777', lw=.7)
        fraction.set(xlim=(0, 180), ylim=(0, 1.02), xlabel='Radius [km]',
                     ylabel='Thermal pressure / total pressure')
        evolution.set(yscale='log', xlim=(0,history[-1, 0]*1e3), ylim=(8e11,4e14),
                      xlabel='Time [ms]', ylabel='Central density [g cm$^{-3}$]')
        evolution.plot(history[:,0]*1e3, history[:,1], color='#c5cddd', lw=1.5)
        evolution.axhline(model['rho_nuclear'], color='#c05e33', ls='--', lw=1)
        evolution.axvline(report['bounce_s']*1e3, color='#9162a5', ls=':', label='Density maximum')
        initial = frames[0][1]
        density.plot(initial['x']/1e5, initial['rho'], color='#b4bfd0', ls=':', label='Initial')
        rho_line, = density.plot([],[], color='#126c99', lw=2, label='Current')
        v_line, = velocity.plot([],[], color='#126c99', lw=2)
        f_line, = fraction.plot([],[], color='#d57b2c', lw=2)
        track, = evolution.plot([],[], color='#126c99', lw=2)
        current, = evolution.plot([],[], 'o', color='#c05e33', ms=5)
        marker_v = velocity.axvline(0, color='#c05e33', ls='--', label='Compressive front')
        marker_f = fraction.axvline(0, color='#c05e33', ls='--')
        density.legend(fontsize=8, loc='lower left')
        velocity.legend(fontsize=8, loc='upper right')
        evolution.legend(fontsize=8, loc='lower right')
        # Half the frames cover the final 1.4 ms, so the nascent shock is visible.
        split = report['bounce_s']-6e-4
        targets = np.r_[np.linspace(0,split,args.frames//2,endpoint=False),
                        np.linspace(split,history[-1,0],args.frames-args.frames//2)]
        ids = np.unique(np.clip(np.searchsorted(history[:,0], targets),0,len(frames)-1))
        images = []
        for j in ids:
            time, data = frames[j]
            pc, _, _ = cold_state(data['rho'], model)
            rho_line.set_data(data['x']/1e5,data['rho'])
            v_line.set_data(data['x']/1e5,data['vx']/1e9)
            f_line.set_data(data['x']/1e5,(data['p']-pc)/data['p'])
            track.set_data(history[:j+1,0]*1e3,history[:j+1,1])
            current.set_data([time*1e3],[history[j,1]])
            phase = ('COLLAPSE' if time < report['nuclear_crossing_s'] else
                     'NUCLEAR STIFFENING' if time < report['bounce_s'] else 'BOUNCE / OUTGOING SHOCK')
            clock.set_text(f't = {time*1e3:6.3f} ms   |   {phase}   |   central density = {history[j,1]:.2e} g/cm³')
            show_front = time >= report['bounce_s']
            for line in (marker_v, marker_f):
                line.set_visible(show_front)
                line.set_xdata([history[j,3]/1e5]*2)
            fig.canvas.draw()
            images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
        durations = [75]*len(images)
        durations[0], durations[-1] = 650, 1800
        target = destination/'nuclear_bounce.gif'
        images[0].save(target, save_all=True, append_images=images[1:], duration=durations, loop=0)
        fig.savefig(destination/'nuclear_bounce_last.png', dpi=140)
        plt.close(fig)
        print(f'{len(images)} frames: {target.resolve()}')


if __name__ == '__main__':
    main()
