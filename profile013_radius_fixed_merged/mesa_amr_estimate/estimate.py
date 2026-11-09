"""Diagnostic MESA AMR count estimate, not an EOS-consistent MESA remesh.

Replicates geometric candidate selection, default velocity limit, neighbor
selection, KE-based merged velocity and mass-weighted abundances. Sound speed
is either inherited from the retained cell or approximated by mass-weighted
cs**2. Splits inherit fluid state (no reconstruction); global radius rescaling
is omitted (uniform scaling preserves the geometric size ratios). Mass comes
from enclosed mass differences: this converted profile's dm column is q.
"""
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from plot_mesa_cfl import read_profile

p = read_profile(ROOT / 'profile013_radius_fixed.data')
isotopes = ['neut','h1','prot','he3','he4','c12','n14','o16','ne20','mg24',
            'si28','s32','ar36','ca40','ti44','cr48','cr60','fe52','fe54','fe56',
            'co56','ni56']
r0 = 10**p['logR']*6.96e10
rc0 = r0[-1]-p['dr'][-1]
m0 = p['mass']-np.r_[p['mass'][1:],0.0]
x0 = np.array([p[k] for k in isotopes]).T

def run(nbase, logarithmic, full_amr, cut, sound_mode, velocity_limit=0.1,
        split_center=True, return_mesh=False):
    n = int(np.sum(r0 > 0.20711213012337312*6.96e10)) if cut else len(r0)
    rc = r0[n] if n < len(r0) else rc0
    r, m, u, cs, xa = (a[:n].copy() for a in
                       (r0,m0,p['velocity'],p['csound'],x0))
    xm = m.sum()
    merges = splits = 0
    checkpoints = {}
    for step in range(20000):
        faces = np.r_[r,rc]
        coord = np.log(np.maximum(faces,1.0)) if logarithmic else faces
        dx = coord[:-1]-coord[1:]
        base = (coord[0]-coord[-1])/nbase
        dq = m/xm
        score = np.maximum(base/dx,1e-14/dq)
        du = np.abs(np.diff(u))
        jump = np.maximum(np.r_[0.,du],np.r_[du,0.])/cs
        allowed = (np.arange(len(r))>=2)&(dq<0.2)&(jump<=velocity_limit)
        score[~allowed] = 0.
        # MESA scans from the center; reverse argmax preserves its tie rule.
        im = len(r)-1-int(np.argmax(score[::-1]))
        big = dx/base
        big[dq<=5e-14] = 0.
        if not split_center:
            big[-1] = 0.
        isp = len(r)-1-int(np.argmax(big[::-1]))
        can_m = score[im]>1.5
        can_s = full_amr and big[isp]>1.5
        if not can_m and not can_s:
            break
        # Fortran iter starts at 1; odd prefers split, even prefers merge.
        do_m = can_m and ((step+1)%2==0 or not can_s)
        if do_m:
            i = im
            if 0 < i < len(r)-1:
                if np.argmax(xa[i]) == np.argmax(xa[i-1]):
                    left = r[-1]-rc if i+1 == len(r)-1 else r[i]-r[i+1]
                    right = r[i-1]-r[i]
                    if right < left:
                        i -= 1
            if i == len(r)-1:
                i -= 1
            j = i+1
            total = m[i]+m[j]
            unew = np.sqrt((m[i]*u[i]**2+m[j]*u[j]**2)/total)
            if u[i]<0:
                unew = -unew
            xa[i] = (m[i]*xa[i]+m[j]*xa[j])/total
            if sound_mode == 'mass_weighted_cs2':
                cs[i] = np.sqrt((m[i]*cs[i]**2+m[j]*cs[j]**2)/total)
            m[i],u[i] = total,unew
            r,m,u,cs,xa = (np.delete(a,j,axis=0) for a in (r,m,u,cs,xa))
            merges += 1
        else:
            i = isp
            inner = r[i+1] if i+1<len(r) else rc
            # MESA bisects radius. Here density is held uniform within parent.
            mid = 0.5*(r[i]+inner)
            frac_outer = (r[i]-mid)*(r[i]**2+r[i]*mid+mid**2)/(
                (r[i]-inner)*(r[i]**2+r[i]*inner+inner**2))
            parent_mass = m[i]
            m[i] *= frac_outer
            r = np.insert(r,i+1,mid)
            m,u,cs,xa = (np.insert(a,i+1,a[i],axis=0) for a in (m,u,cs,xa))
            m[i+1] = parent_mass*(1-frac_outer)
            splits += 1
        if step+1 in (1000,2000,5000,10000):
            checkpoints[str(step+1)] = dict(zones=len(r),merges=merges,splits=splits)
    result = dict(initial_zones=n,nbase=nbase,logarithmic=logarithmic,
                full_amr=full_amr,cut=cut,sound_mode=sound_mode,
                velocity_limit=velocity_limit,split_center=split_center,
                final_zones=len(r),
                merges=merges,splits=splits,converged=step<19999,
                checkpoints=checkpoints)
    if return_mesh:
        return result, dict(r_outer_cm=r, r_inner_cm=np.r_[r[1:],rc],
                            mass_msun=m, velocity_cm_s=u, csound_cm_s=cs)
    return result

if __name__ == '__main__':
    results=[]
    for cut in (False,True):
        for nbase,log in ((1200,True),(1600,False)):
            for full in (False,True):
                for sound in ('retain','mass_weighted_cs2'):
                    result=run(nbase,log,full,cut,sound,split_center=log)
                    results.append(result)
                    print(json.dumps(result),flush=True)
    for cut in (False,True):
        result=run(1200,True,True,cut,'retain',split_center=False)
        results.append(result)
        print(json.dumps(result),flush=True)
    out=Path(__file__).with_name('results.json')
    out.write_text(json.dumps({'input':str(ROOT/'profile013_radius_fixed.data'),
                              'results':results},indent=2)+'\n')
