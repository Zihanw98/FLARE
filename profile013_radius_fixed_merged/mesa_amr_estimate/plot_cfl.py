"""Compare fixed-mesh CFL=1 diagnostics on the approximate AMR meshes."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from estimate import run, p, r0, rc0, m0

OUT = Path(__file__).resolve().parent
records = []
datasets = {}

def summarize(name, mesh, config=None):
    ro, ri = mesh['r_outer_cm'], mesh['r_inner_cm']
    dr = ro-ri
    cs, v = mesh['csound_cm_s'], mesh['velocity_cm_s']
    dt = dr/(cs+np.abs(v))
    sound_dt = dr/cs
    assert np.all(np.isfinite(dt)) and np.all(dt>0)
    assert np.all(ro>ri)
    assert np.allclose(ri[:-1], ro[1:], rtol=0, atol=0)
    center = ri+dr/2
    k = int(np.argmin(dt))
    row = dict(name=name,zones=len(dt),min_dt_s=float(dt[k]),
               min_sound_only_dt_s=float(sound_dt.min()),
               limiting_zone_surface_inward=k+1,
               limiting_r_inner_cm=float(ri[k]),limiting_r_outer_cm=float(ro[k]),
               limiting_dr_cm=float(dr[k]),limiting_csound_cm_s=float(cs[k]),
               limiting_velocity_cm_s=float(v[k]),config=config)
    records.append(row)
    columns = dict(zone=np.arange(1,len(dt)+1),**mesh,
                   r_mid_cm=center,dr_cm=dr,dt_cfl1_s=dt,dt_sound_only_s=sound_dt)
    np.savetxt(OUT/(name+'.csv'),np.column_stack(list(columns.values())),
               delimiter=',',header=','.join(columns),comments='',fmt='%.16e')
    datasets[name]=(center,dt,row)
    return row

for cut in (False,True):
    scope = 'cut4msun' if cut else 'whole'
    n = int(np.sum(r0>0.20711213012337312*6.96e10)) if cut else len(r0)
    rc = r0[n] if cut else rc0
    original = dict(r_outer_cm=r0[:n],r_inner_cm=np.r_[r0[1:n],rc],
                    mass_msun=m0[:n],velocity_cm_s=p['velocity'][:n],
                    csound_cm_s=p['csound'][:n])
    summarize(scope+'_original',original)
    for base, log, full, center, tag in (
        (1200,True,False,True,'log_merge'),
        (1200,True,True,True,'log_amr'),
        (1200,True,True,False,'log_amr_no_center_split'),
        (1600,False,False,False,'linear_merge'),
        (1600,False,True,False,'linear_amr')):
        config, mesh = run(base,log,full,cut,'retain',split_center=center,
                           return_mesh=True)
        assert config['converged']
        assert np.isclose(mesh['mass_msun'].sum(),m0[:n].sum(),rtol=2e-14)
        assert np.isclose(np.sum(mesh['mass_msun']*mesh['velocity_cm_s']**2),
                          np.sum(m0[:n]*p['velocity'][:n]**2),rtol=2e-14)
        summarize(scope+'_'+tag,mesh,config)

plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes = plt.subplots(2,2,figsize=(15,10),sharey=True)
for row,scope in enumerate(('whole','cut4msun')):
    for col,mode in enumerate(('log','linear')):
        ax = axes[row,col]
        curves=[('original','Original','#8b95a5','-',1.0),
                (mode+'_merge','Merge only','#e17c05','-',1.5),
                (mode+'_amr','Split + merge','#147d92','-',1.5)]
        if mode=='log' and scope=='whole':
            curves.append(('log_amr_no_center_split','Split + merge; no center split',
                           '#9758a6','--',1.4))
        for tag,label,color,style,width in curves:
            radius,dt,info=datasets[scope+'_'+tag]
            label=f"{label}: N={info['zones']:,}, min={info['min_dt_s']:.2e} s"
            ax.plot(radius[::-1],dt[::-1],color=color,ls=style,lw=width,
                    alpha=.85,label=label)
            k=int(np.argmin(dt))
            ax.scatter(radius[k],dt[k],s=30,color=color,zorder=6)
        title=('Whole profile (4,147 input cells)' if row==0 else
               'After inner ~4 Msun cut (2,599 input cells)')
        ax.set_title(title+'\n'+('Log-radius target: 1,200' if col==0 else
                                   'Linear-radius target: 1,600'),fontsize=12,pad=10)
        ax.set_xscale('log');ax.set_yscale('log')
        ax.set_xlabel('Shell midpoint radius [cm]')
        if col==0: ax.set_ylabel('Local CFL=1 timestep [s]')
        ax.grid(True,which='major',alpha=.2)
        ax.legend(loc='upper left',fontsize=8.5,framealpha=.96)
fig.suptitle(r'Approximate MESA remeshing: $\Delta t=\Delta r/(c_s+|v|)$',
             fontsize=18,y=.98)
fig.text(.5,.025,'Dots mark global minima. Widths use logR shell faces for every curve. '
         'Sound speed is inherited from the retained cell; no EOS reconstruction.\n'
         'Meshes are iterated to a fixed point; these are profile diagnostics, '
         'not timestep measurements from a MESA or FLARE run.',
         ha='center',fontsize=10,color='#444444')
fig.tight_layout(rect=(0,.075,1,.94),h_pad=2.5,w_pad=2)
fig.savefig(OUT/'cfl1_comparison.png',dpi=180)
fig.savefig(OUT/'cfl1_comparison.pdf')
plt.close(fig)
(OUT/'cfl1_summary.json').write_text(json.dumps(dict(
    formula='dt = (r_outer-r_inner)/(csound+abs(velocity)); CFL=1',
    assumptions='Same approximate AMR algorithm as estimate.py, inherited sound speed; '
                'no EOS, no hydrodynamic evolution, no gravitational radius correction.',
    results=records),indent=2)+'\n')
for info in records:
    print(info['name'],info['zones'],f"min_dt={info['min_dt_s']:.10g} s",
          f"limiting_zone={info['limiting_zone_surface_inward']}")
