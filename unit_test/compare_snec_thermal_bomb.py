#!/usr/bin/env python3
"""Rebuild unmodified local SNEC, import its initial state, run and compare.

All generated files stay in tests/SNEC_Thermal_Bomb. SNEC's original directory
is read-only. Requires numpy, matplotlib, g++, gfortran and LAPACK.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "tests/SNEC_Thermal_Bomb"
MSUN = 1.98892e33
ENERGY = 1.6e51
G_SNEC = 6.6742e-8
G_FLARE = 6.67430e-8


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes(root):
    return {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob("*"))
            if p.is_file() and p.suffix in (".cpp", ".h", ".F90", ".dek")}


def launch(cmd, cwd, log):
    stream = log.open("w")
    proc = subprocess.Popen([str(x) for x in cmd], cwd=cwd, stdout=stream, stderr=subprocess.STDOUT)
    stream.close()
    return proc


def run_checked(cmd, cwd, log):
    proc = launch(cmd, cwd, log)
    if proc.wait():
        raise RuntimeError(f"Command failed; see {log}")


def xg(path):
    times, frames = [], []
    rows = []
    for line in path.open():
        if 'Time =' in line:
            if rows:
                frames.append(np.array(rows)); rows = []
            times.append(float(line.split('=')[1].replace('D', 'E')))
        elif line.strip():
            rows.append([float(v.replace('D', 'E')) for v in line.split()])
    if rows:
        frames.append(np.array(rows))
    return np.array(times), np.array(frames)


def snec_parameters(tend, radiation=0, eos=1, n=1000):
    return f'''# Shared injected energy, not an asymptotic final energy.
outdir = "data"
profile_name = "profiles/snec_input13.short"
comp_profile_name = "profiles/snec_input13.iso.dat"
initial_data = "Thermal_Bomb"
bomb_mode = 2
final_energy = {ENERGY:.10e}
bomb_tstart = 0.0
bomb_tend = 0.1
bomb_mass_spread = 0.1
bomb_start_point = 1
imax = {n}
gridding = "from_file_by_mass"
mass_excision = 1
mass_excised = 1.502
radiation = {radiation}
eoskey = {eos}
Ni_switch = 0
Ni_mass = 0.0
Ni_boundary_mass = 11.0
Ni_period = 1.0e4
saha_ncomps = 3
boxcar_smoothing = 0
opacity_floor_envelope = 0.01
opacity_floor_core = 0.24
ntmax = 10000000
tend = {tend:.10e}
dtout = {max(tend/100, 1e-5):.10e}
dtout_scalar = {max(tend/100, 1e-5):.10e}
dtout_check = 1.0e4
ntout = -1
ntout_scalar = -1
ntout_check = -1
ntinfo = 1
dtmin = 1.0e-12
dtmax = 1.0e-2
sedov = 0
'''


def grid_pattern(radial_cap=0):
    grid=np.loadtxt(CASE/'inputs/GridPattern.dat')
    if not radial_cap:
        return grid
    profile=np.loadtxt(CASE/'inputs/snec_input13.short',skiprows=1)
    m0=1.502*MSUN; mt=profile[-1,1]
    m=m0+grid*(mt-m0)
    r0=np.interp(m0,profile[:,1],profile[:,2])
    for _ in range(30):
        rho=np.interp((m[:-1]+m[1:])/2,profile[:,1],profile[:,4])
        r=np.cbrt(np.r_[r0**3,3/(4*np.pi)*np.diff(m)/rho].cumsum())
        splits=np.maximum(1,np.ceil(np.log(r[1:]/r[:-1])/np.log1p(radial_cap)).astype(int))
        if np.max(splits)==1:
            return (m-m0)/(mt-m0)
        segments=[]
        for i,n in enumerate(splits):
            # Insert mass faces at geometrically spaced radii in the current
            # shell; re-evaluate the original profile on the new mass grid.
            f=np.geomspace(r[i],r[i+1],n+1)
            fraction=(f**3-r[i]**3)/(r[i+1]**3-r[i]**3)
            part=m[i]+fraction*(m[i+1]-m[i]);part[0]=m[i]
            segments.append(part[:-1])
        m=np.r_[np.concatenate(segments),mt]
    raise RuntimeError('Radial mesh refinement did not converge')


def setup_snec(work, source, tend, radiation=0, eos=1, radial_cap=0):
    work.mkdir(parents=True)
    (work / "data").mkdir()
    (work / "profiles").symlink_to(CASE / "inputs", target_is_directory=True)
    (work / "tables").mkdir()
    for p in (source / "tables").iterdir():
        if p.is_file() and p.name != "GridPattern.dat":
            (work / "tables" / p.name).symlink_to(p)
    grid = grid_pattern(radial_cap)
    np.savetxt(work / "tables/GridPattern.dat", grid, fmt="%.17e")
    (work / "parameters").write_text(snec_parameters(tend, radiation, eos, len(grid)))


def export_initial(snec, flare, tend, solver="exact", refine=1):
    """SNEC already interpolated rho/T in mass and reintegrated radius.

    Import full precision xg data instead of rounded *_initial.dat. This avoids
    a different radial remap masquerading as a hydrodynamic solver difference.
    """
    flare.mkdir(parents=True)
    data = snec / "data"
    _, rr = xg(data / "radius.xg")
    _, rh = xg(data / "rho.xg")
    _, pp = xg(data / "press.xg")
    _, vv = xg(data / "vel.xg")
    _, ee = xg(data / "eps.xg")
    m, r = rr[0,:,0], rr[0,:,1]
    rho, pressure = rh[0,:-1,1], pp[0,:-1,1]
    np.testing.assert_allclose(pressure/(rho*ee[0,:-1,1]),0.4,rtol=1e-12,
                               err_msg='The selected SNEC ideal EoS is not gamma=1.4')
    match=re.search(r'Total energy of the bomb\s*=\s*([\d.EeDd+\-]+)',(snec/'run.log').read_text())
    if not match or abs(float(match[1].replace('D','E'))/ENERGY-1)>1e-12:
        raise RuntimeError('SNEC did not apply the requested injected bomb energy (bomb_mode=2)')
    vel = (vv[0,:-1,1]+vv[0,1:,1])/2
    # Optional splitting retains exactly the coarse initial piecewise state.
    faces = np.concatenate([np.cbrt(np.linspace(a**3,b**3,refine+1))[:-1] for a,b in zip(r[:-1],r[1:])] + [r[-1:]])
    masses = np.concatenate([np.linspace(a,b,refine+1)[:-1] for a,b in zip(m[:-1],m[1:])] + [m[-1:]])
    shells = np.column_stack((faces[:-1], faces[1:], np.repeat(rho,refine), np.repeat(pressure,refine),
                             np.repeat(vel,refine), masses[:-1],masses[1:]))
    n = len(shells)
    np.savetxt(flare / "initial_shells.dat", shells, header=str(n), comments='', fmt='%.17e')
    nb = int(np.searchsorted(m, m[0]+0.1*MSUN))
    alpha = np.log(100)/(m[nb-1]-m[0]) if nb > 1 else 0
    weights = np.exp(-alpha*(m[:nb]-m[0]))
    weights /= np.dot(weights,np.diff(m)[:nb])
    np.savetxt(flare / "bomb.dat", np.column_stack((m[:nb],m[1:nb+1],weights)),
               header=f'{ENERGY:.17e} 0.1 1e-4 {nb}', comments='', fmt='%.17e')
    params = f'''# CGS. Identical SNEC initial mass/rho/P and gamma=1.4.
nxa = {n}
nxg = 2
x_lower = {faces[0]:.17e}
x_upper = {faces[-1]:.17e}
geometry = spherical
gamma = 1.4
equation_of_state = IdealGas_Default
spherical_gravity_enabled = true
spherical_self_gravity_enabled = true
spherical_gravity_central_mass = {m[0]:.17e}
constant_gravity_enabled = false
constant_heating_enabled = false
boundary_x_lower = reflective
boundary_x_upper = user
t_end = {tend:.17e}
cfl = 0.4
cfl_max = 0.5
maxsteps = 3000000
riemann_solver = {solver}
riemann_tolerance = 1e-10
riemann_max_iterations = 200
reconstructor = PLM
limiter = minmod
time_integrator = rk2
output_step_every = 100
output_directory = output
output_suffix = shock
'''
    (flare / "parameters.txt").write_text(params)
    return dict(active_cells=n, inner_radius_cm=float(r[0]), outer_radius_cm=float(r[-1]),
                excised_mass_g=float(m[0]), envelope_mass_g=float(m[-1]-m[0]),
                bomb_cells_snec=nb, bomb_support_mass_msun=float((m[nb]-m[0])/MSUN),
                bomb_energy_erg=ENERGY, bomb_duration_s=.1, gamma=1.4)


def run(args):
    source = args.snec_root.resolve()
    output = args.output.resolve()
    if output.exists():
        raise RuntimeError(f"Refusing to mix results in existing {output}; choose --output")
    output.mkdir(parents=True)
    before = source_hashes(ROOT / "src")
    source_before = source_hashes(source / "src")
    build = CASE / ".build/snec"
    build.mkdir(parents=True, exist_ok=True)
    (build / "src").mkdir(exist_ok=True)
    for p in (source / "src").iterdir():
        if p.suffix in ('.F90','.dek') or p.name == 'Makefile':
            shutil.copy2(p, build / "src" / p.name)
    for name in ('Makefile','make.inc'):
        shutil.copy2(source / name, build / name)
    run_checked(['make','clean'], build, output / "snec_clean.log")
    run_checked(['make'], build, output / "snec_build.log")
    run_checked(['make','-j4'], CASE, output / "flare_build.log")
    probe = output / "snec_initial"
    setup_snec(probe, source, 1e-8,radial_cap=args.radial_cap)
    run_checked([build / "snec"], probe, probe / "run.log")
    if 'tend reached' not in (probe / 'run.log').read_text():
        raise RuntimeError(f"SNEC initialization failed: {probe / 'run.log'}")
    model = export_initial(probe, output / "flare_exact", args.t_end)
    initial_shells=np.loadtxt(output/'flare_exact/initial_shells.dat',skiprows=1)
    model['max_initial_dr_over_r']=float(np.max(initial_shells[:,1]/initial_shells[:,0]-1))
    if args.radial_cap:
        assert model['max_initial_dr_over_r'] <= args.radial_cap*(1+1e-7)
    export_initial(probe, output / "flare_hll", args.t_end, 'hll')
    if args.refine:
        export_initial(probe, output / "flare_exact_2x", args.t_end, refine=2)
    variants = [('snec_ideal',0,1),('snec_paczynski_rad',1,2)]
    if args.hydro_only:
        variants = variants[:1]
    jobs = []
    for name,rad,eos in variants:
        work = output / name
        setup_snec(work,source,args.t_end,rad,eos,radial_cap=args.radial_cap)
        jobs.append((name,launch([build/'snec'],work,work/'run.log'),time.monotonic()))
    for name in ['flare_exact','flare_hll']+(['flare_exact_2x'] if args.refine else []):
        work = output/name
        jobs.append((name,launch([CASE/'flare','-par','parameters.txt','--performance-file','performance.json'],work,work/'run.log'),time.monotonic()))
    manifest = dict(flare_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                    flare_sources=before, snec_root=str(source), snec_sources=source_before,
                    inputs={name:digest(CASE/'inputs'/name) for name in
                            ('snec_input13.short','snec_input13.iso.dat','GridPattern.dat')},
                    adapter_sha256=digest(ROOT/'unit_test/snec_thermal_bomb.cpp'),
                    runner_sha256=digest(Path(__file__)),
                    t_end=args.t_end, radial_cap=args.radial_cap, model=model, results={})
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    while jobs:
        for name,proc,start in jobs[:]:
            if proc.poll() is not None:
                log = (output/name/'run.log').read_text()
                ok = proc.returncode == 0 and ('tend reached' in log if name.startswith('snec') else json.loads((output/name/'performance.json').read_text())['status']=='completed')
                manifest['results'][name] = dict(success=ok, exit_code=proc.returncode, wall_s=time.monotonic()-start)
                print(name, manifest['results'][name], flush=True)
                jobs.remove((name,proc,start))
                (output/'manifest.json').write_text(json.dumps(manifest,indent=2))
        if jobs:
            time.sleep(2)
    assert before == source_hashes(ROOT/'src'), 'FLARE framework source changed during run'
    assert source_before == source_hashes(source/'src'), 'Original SNEC source changed during run'
    manifest['source_unchanged'] = True
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2))
    analyze(args)


def flare_frames(work):
    grid = np.genfromtxt(work/'output/grid.csv', delimiter=',',names=True)
    frames=[]
    for file in sorted((work/'output').glob('shock_*.csv')):
        with file.open() as f:
            f.readline(); t=float(f.readline().split('=')[1])
        a=np.genfromtxt(file,delimiter=',',names=True,skip_header=7)
        a=a[a['is_ghost']==0]
        dm=a['rho']*grid['volume']
        m=np.r_[grid['mass_left'][0], grid['mass_left'][0]+np.cumsum(dm)]
        frames.append(dict(t=t,r=np.r_[grid['r_left'],grid['r_right'][-1]],m=m,
                           rho=a['rho'],p=a['p'],v=a['vx'],e=a['total_energy']*grid['volume'],
                           u=a['internal_energy']*grid['volume']))
    return frames


def snec_frames(work):
    arrays={k:xg(work/'data'/f'{k}.xg') for k in ('radius','rho','press','vel','eps')}
    ts=arrays['radius'][0]; frames=[]
    for i,t in enumerate(ts):
        r=arrays['radius'][1][i,:,1]; m=arrays['radius'][1][i,:,0]
        rho=arrays['rho'][1][i,:-1,1]; p=arrays['press'][1][i,:-1,1]
        v=(arrays['vel'][1][i,:-1,1]+arrays['vel'][1][i,1:,1])/2
        u=arrays['eps'][1][i,:-1,1]*np.diff(m)
        frames.append(dict(t=t,r=r,m=m,rho=rho,p=p,v=v,u=u,e=u+.5*v*v*np.diff(m)))
    return frames


def shock(frame):
    # Maximum outward compressive velocity drop, using cell-centred velocities
    # for both codes. Retain adjacent volume-centre radii as a grid bracket.
    r=np.cbrt((frame['r'][:-1]**3+frame['r'][1:]**3)/2)
    v=frame['v']; p=frame['p']
    strength=np.maximum(v[:-1]-v[1:],0)
    eligible=(v[:-1]>1e7)&(p[:-1]>1.1*p[1:])&(r[:-1]>frame['r'][0]*1.06)
    if not np.any(eligible):
        return [np.nan]*5
    score=np.where(eligible,strength,0)
    # Largest compressive drop, with its bracket retained as spatial uncertainty.
    j=int(np.argmax(score))
    return [float(frame['r'][j+1]),float(frame['m'][j+1]/MSUN),float(r[j]),float(r[j+1]),float(strength[j])]


def injected_history(work, snec, tend):
    """Reconstruct accepted time quadrature from the unmodified runtime logs.

    SNEC evaluates heating at the beginning of its step, FLARE uses RK2.
    SNEC log times have six significant decimals, so this is an audit estimate.
    """
    log=(work/'run.log').read_text()
    c=np.log(100)/.1
    def power(t):
        return c*np.exp(-c*np.minimum(t,1))/.99*(t<=.1)
    if snec:
        starts=np.array([float(a[1]) for line in log.splitlines()
                         if len(a:=line.split())==5 and a[0].isdigit()])
        ends=np.r_[starts[1:],tend]
        increments=power(starts)*(ends-starts)
    else:
        a=np.array([[float(t),float(dt)] for t,dt in re.findall(
            r'time = ([\deE+.\-]+), dt = ([\deE+.\-]+)',log)])
        ends,dt=a.T
        increments=.5*(power(ends)+power(ends-dt))*dt
    return ends,np.cumsum(increments)


def analyze(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=args.output.resolve(); manifest=json.loads((out/'manifest.json').read_text())
    report=out/'comparison'; report.mkdir(exist_ok=True)
    allframes={}
    for name,result in manifest['results'].items():
        if result['success']:
            allframes[name]=snec_frames(out/name) if name.startswith('snec') else flare_frames(out/name)
    if not {'snec_ideal','flare_exact','flare_hll'} <= allframes.keys():
        raise RuntimeError('A required run failed; see manifest and logs')
    initial=allframes['snec_ideal'][0]; summary={}
    fig,ax=plt.subplots(2,2,figsize=(12,9))
    for name,frames in allframes.items():
        hist=[]
        grav=G_SNEC if name.startswith('snec') else G_FLARE
        ti,ei=injected_history(out/name,name.startswith('snec'),manifest['t_end'])
        def total(f):
            # Exact gravitational integral of each uniform-density spherical
            # shell, including its self gravity and all interior mass. Algebra
            # in shell width avoids subtracting nearly equal high powers.
            a,b=f['r'][:-1],f['r'][1:]; h=b-a; d=a*a+a*b+b*b
            dm=np.diff(f['m'])
            center=1.5*(a+b)/d
            self_term=3*(1.5*a**3+2*a*a*h+a*h*h+.2*h**3)/(d*d)
            return float(np.sum(f['e'])-np.sum(grav*dm*(f['m'][:-1]*center+dm*self_term)))
        e0=total(frames[0]); mass0=frames[0]['m'][-1]-frames[0]['m'][0]
        for f in frames:
            hist.append([f['t'],*shock(f),total(f),(f['m'][-1]-f['m'][0])/mass0-1])
        h=np.array(hist)
        np.savetxt(report/f'{name}_history.csv',h,delimiter=',',
                   header='time_s,shock_radius_cm,shock_mass_msun,shock_bracket_lo_cm,shock_bracket_hi_cm,velocity_drop_cms,total_energy_erg,relative_mass_change',comments='')
        last=frames[-1]
        summary[name]=dict(final_time_s=last['t'],shock_radius_cm=h[-1,1],shock_mass_msun=h[-1,2],
            shock_bracket_cm=h[-1,3:5].tolist(),energy_residual_over_bomb=(total(last)-e0-ENERGY)/ENERGY,
            mass_change=h[-1,7],max_speed_cms=float(np.max(np.abs(last['v']))),
            finite_positive=bool(all(np.all(np.isfinite(f['rho'])) and np.all(f['rho']>0) and np.all(np.isfinite(f['p'])) and np.all(f['p']>0) for f in frames)))
        summary[name]['injected_energy_over_nominal']=float(ei[-1])
        summary[name]['energy_residual_after_actual_heating_over_bomb']=(total(last)-e0)/ENERGY-float(ei[-1])
        assert summary[name]['finite_positive']
        assert abs(last['t']-manifest['t_end']) < 1e-8*max(1,manifest['t_end'])
        if name.startswith('flare'):
            f=frames[0]
            assert json.loads((out/name/'output/source_checks.json').read_text()) == {
                'initial_and_redistributed_mass_power':True,'after_bomb_zero':True}
            if len(f['rho'])==len(initial['rho']):
                errors={key:float(np.max(np.abs(f[key]/initial[key]-1))) for key in ('rho','p','r','m')}
                summary[name]['initial_max_relative_errors']=errors
                assert max(errors.values())<2e-9
        ax[0,0].plot(h[:,0],h[:,1]/1e10,label=name)
        ax[0,1].plot(h[:,0],h[:,2],label=name)
        ax[1,0].plot(h[:,0],(h[:,6]-e0)/ENERGY-np.interp(h[:,0],ti,ei,left=0),label=name)
        r=np.cbrt((last['r'][:-1]**3+last['r'][1:]**3)/2)
        ax[1,1].plot(r/1e10,last['v']/1e8,label=name)
    ref=summary['snec_ideal']['shock_radius_cm']
    for name in summary:
        summary[name]['shock_radius_vs_snec_ideal_percent']=100*(summary[name]['shock_radius_cm']/ref-1)
    (report/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    ax[0,0].set(xlabel='Time [s]',ylabel='Shock radius [1e10 cm]')
    ax[0,1].set(xlabel='Time [s]',ylabel='Enclosed shock mass [SNEC solar mass]')
    ax[1,0].set(xlabel='Time [s]',ylabel='(E - E0 - actual heating) / nominal bomb')
    ax[1,1].set(xlabel='Radius [1e10 cm]',ylabel='Velocity [1000 km/s]',xlim=(0,ref/1e10*1.5))
    for a in ax.flat:
        a.grid(alpha=.2); a.legend(fontsize=8)
    fig.suptitle(f'profile013 | mass cut 1.502 Msun | thermal bomb 1.6e51 erg, 0.1 s',fontsize=13)
    fig.tight_layout();fig.savefig(report/'comparison.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,10),sharex=True)
    for name,frames in allframes.items():
        f=frames[-1]; r=np.cbrt((f['r'][:-1]**3+f['r'][1:]**3)/2)
        axes[0].semilogy(r/1e10,f['rho'],label=name)
        axes[1].plot(r/1e10,f['v']/1e8,label=name)
        axes[2].semilogy(r/1e10,f['p'],label=name)
    for a,y in zip(axes,('Density [g/cm3]','Velocity [1000 km/s]','Pressure [dyn/cm2]')):
        a.set(ylabel=y,xlim=(0,ref/1e10*1.5));a.grid(alpha=.2);a.legend(fontsize=8)
    axes[-1].set_xlabel('Radius [1e10 cm]')
    fig.suptitle(f'Same initial profile and prescribed bomb | t = {manifest["t_end"]:g} s')
    fig.tight_layout();fig.savefig(report/'profiles.png',dpi=180);plt.close(fig)
    # Compare shock trajectories at common physical times, preserving the raw
    # histories above. Temporal interpolation does not improve spatial accuracy.
    rows=[]
    for t in (1.,3.,10.,30.,100.):
        if t>manifest['t_end']: continue
        for name,frames in allframes.items():
            history=np.loadtxt(report/f'{name}_history.csv',delimiter=',',skiprows=1)
            valid=np.isfinite(history[:,1]); h=history[valid]
            if t<h[0,0]: continue
            rows.append(f'{t:g},{name},{np.interp(t,h[:,0],h[:,1]):.10e},{np.interp(t,h[:,0],h[:,2]):.10e}')
    (report/'common_times.csv').write_text('time_s,run,shock_radius_cm,shock_mass_msun\n'+'\n'.join(rows)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['run','analyze'])
    p.add_argument('--snec-root',type=Path,default=ROOT.parent/'SNEC/SNEC_NEW')
    p.add_argument('--output',type=Path,default=CASE/'output_radial02')
    p.add_argument('--t-end',type=float,default=100.)
    p.add_argument('--radial-cap',type=float,default=0.02,help='maximum initial dr/r in both codes; 0 uses native SNEC mass grid')
    p.add_argument('--refine',action='store_true',help='also split each FLARE shell into two equal-volume shells')
    p.add_argument('--hydro-only',action='store_true')
    a=p.parse_args()
    if a.t_end<=.1: p.error('t-end must exceed the 0.1 s bomb duration')
    if a.radial_cap<0: p.error('radial-cap must be nonnegative')
    run(a) if a.command=='run' else analyze(a)
