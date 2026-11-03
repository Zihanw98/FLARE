#!/usr/bin/env python3
"""Run/check/animate acoustic, oblique Sod, RM and 3D shock-cloud benchmarks.

References: linearized Euler acoustics and exact Sod Riemann relations.
RM/cloud diagnostics are qualitative regression checks, not analytic solutions.
Only snapshots listed by output/run.log are analyzed. No passive tracer is assumed.
"""
import argparse
from contextlib import ExitStack
import csv
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from animate3d import snapshot_files


def require(condition, message):
    if not condition:
        raise ValueError(message)


def parameters(path):
    cfg = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    return cfg


def read_frame(path, cfg):
    stamp = path.stat().st_size, path.stat().st_mtime_ns
    meta = {}
    with path.open() as stream:
        for line in stream:
            if not line.startswith("#"):
                names = line.strip().split(",")
                break
            key, value = line[1:].strip().split("=", 1)
            meta[key] = float(value)
        keys = ("x", "y", "z", "rho", "vx", "vy", "vz", "p", "gamma", "is_ghost")
        data = np.loadtxt(stream, delimiter=",", usecols=[names.index(k) for k in keys], ndmin=2)
    require(stamp == (path.stat().st_size, path.stat().st_mtime_ns), f"File changed: {path}")
    require(np.isfinite(data).all() and np.isin(data[:,-1], (0,1)).all(), f"Invalid CSV: {path}")
    require(np.all(data[:,8] == float(cfg["gamma"])), "Gamma mismatch")
    dimension = 3 if "nza" in cfg else 2
    count = np.prod([int(cfg[f"n{a}a"])+2*int(cfg[f"n{a}g"]) for a in "xyz"[:dimension]])
    require(len(data) == count, "Incomplete snapshot")
    actual_ghost = np.zeros(len(data), dtype=bool)
    axes = []
    for d, a in enumerate("xyz"):
        if d < dimension:
            lo, hi, n = float(cfg[f"{a}_lower"]), float(cfg[f"{a}_upper"]), int(cfg[f"n{a}a"])
            axes.append(lo+(np.arange(n)+.5)*(hi-lo)/n)
            actual_ghost |= (data[:,d]<lo) | (data[:,d]>hi)
        else:
            axes.append(np.array([0.]))
    require(np.array_equal(actual_ghost, data[:,-1]), "Incorrect ghost flags")
    data = data[~actual_ghost]
    nx, ny, nz = map(len, axes)
    require(len(data) == nx*ny*nz, "Wrong active cell count")
    data = data[np.lexsort((data[:,0],data[:,1],data[:,2]))].reshape(nz,ny,nx,-1)
    z, y, x = np.meshgrid(axes[2],axes[1],axes[0],indexing="ij")
    for d, expected in enumerate((x,y,z)):
        require(np.allclose(data[:,:,:,d],expected,rtol=0,atol=1e-13), "Wrong/duplicate coordinates")
    frame = {key: data[:,:,:,i].copy() for i,key in enumerate(keys) if key in ("rho","vx","vy","vz","p")}
    require(np.all(frame["rho"]>0) and np.all(frame["p"]>0), "Nonpositive density or pressure")
    frame.update(time=meta["time"],step=int(meta["step"]))
    return axes, frame


def sod_points(s, time):
    """Independent scalar pressure solve; vectorized exact planar Sod reference."""
    g = 1.4
    cl, cr = np.sqrt(g), np.sqrt(g*.1/.125)
    left = lambda p: 2*cl/(g-1)*(p**((g-1)/(2*g))-1)
    right = lambda p: (p-.1)*np.sqrt(2/((g+1)*.125)/(p+(g-1)/(g+1)*.1))
    lo, hi = .1, 1.
    for _ in range(70):
        mid = (lo+hi)/2
        if left(mid)+right(mid)>0: hi=mid
        else: lo=mid
    ps=(lo+hi)/2
    us=(right(ps)-left(ps))/2
    if time == 0:
        return np.where(s<0,1.,.125), np.zeros_like(s), np.where(s<0,1.,.1)
    speed=s/time
    cs=cl*ps**((g-1)/(2*g))
    shock=cr*np.sqrt((g+1)/(2*g)*ps/.1+(g-1)/(2*g))
    b=(g-1)/(g+1)
    rho=np.full_like(s,.125);u=np.zeros_like(s);p=np.full_like(s,.1)
    middle=speed<shock
    rho[middle]=.125*(ps/.1+b)/(b*ps/.1+1);u[middle]=us;p[middle]=ps
    rho[speed<=us]=ps**(1/g)
    fan=(speed>=-cl)&(speed<us-cs)
    c=2/(g+1)*(cl-(g-1)*speed[fan]/2)
    rho[fan]=(c/cl)**(2/(g-1));p[fan]=(c/cl)**(2*g/(g-1))
    u[fan]=2/(g+1)*(cl+speed[fan])
    rho[speed<-cl]=1;u[speed<-cl]=0;p[speed<-cl]=1
    return rho,u,p


def sod_average(axes, cfg, time):
    x,y=np.meshgrid(axes[0],axes[1]);nx,ny=2/np.sqrt(5),1/np.sqrt(5)
    cx=(float(cfg["x_lower"])+float(cfg["x_upper"]))/2
    cy=(float(cfg["y_lower"])+float(cfg["y_upper"]))/2
    s=nx*(x-cx)+ny*(y-cy)
    rho=np.zeros_like(s);momentum=rho.copy();energy=rho.copy()
    dx,dy=axes[0][1]-axes[0][0],axes[1][1]-axes[1][0]
    for b in range(8):
        for a in range(8):
            r,u,p=sod_points(s+nx*((a+.5)/8-.5)*dx+ny*((b+.5)/8-.5)*dy,time)
            rho+=r;momentum+=r*u;energy+=p/.4+.5*r*u*u
    rho/=64;momentum/=64;energy/=64
    u=momentum/rho;p=.4*(energy-.5*momentum*u)
    return s,rho,u,p


def rm_interface(axes, frame):
    # Entropy proxy, NOT a material tracer: light gas q~1, heavy gas q~0.215.
    q=frame["p"][0]/frame["rho"][0]**1.4
    interface=[]
    for row in q:
        crossing=np.flatnonzero((row[:-1]>.6)&(row[1:]<=.6))
        require(len(crossing)>0, "Entropy-proxy interface left domain or became unidentifiable")
        j=crossing[-1]
        fraction=(row[j]-.6)/(row[j]-row[j+1])
        interface.append(axes[0][j]+fraction*(axes[0][j+1]-axes[0][j]))
    return np.array(interface)


def check_initial(case,cfg,axes,f):
    x=axes[0][None,None,:];y=axes[1][None,:,None];z=axes[2][:,None,None]
    lo=[float(cfg[f"{a}_lower"]) for a in "xyz"[:(3 if case=="cloud" else 2)]]
    length=[float(cfg[f"{a}_upper"])-lo[d] for d,a in enumerate("xyz"[:len(lo)])]
    if case=="acoustic":
        kx,ky=2*np.pi/length[0],2*np.pi/length[1]
        delta=1e-4*np.sin(kx*(x-lo[0])+ky*(y-lo[1]))*np.sinc(1/len(axes[0]))*np.sinc(1/len(axes[1]))
        expected={"rho":1+delta,"p":1/1.4+delta,"vx":delta*kx/np.hypot(kx,ky),"vy":delta*ky/np.hypot(kx,ky),"vz":0}
    elif case=="sod":
        _,rho,u,p=sod_average(axes,cfg,0)
        expected={"rho":rho[None],"p":p[None],"vx":0,"vy":0,"vz":0}
    else:
        shock=x<lo[0]+(.175 if case=="rm" else .075)*length[0]
        if case=="rm":
            interface=lo[0]+.35*length[0]+.02*length[1]*np.cos(2*np.pi*(y-lo[1])/length[1])
            density=2+np.tanh((x-interface)/(.015*length[1]))
        else:
            r2=(x-lo[0]-.4*length[0])**2+(y-lo[1]-.5*length[1])**2+(z-lo[2]-.5*length[2])**2
            density=np.where(r2<(.18*min(length[1:]))**2,5.,1.)
        expected={"rho":np.where(shock,8/3,density),"p":np.where(shock,4.5,1.),
                  "vx":np.where(shock,1.25*np.sqrt(1.4),0.),"vy":0,"vz":0}
    for key,value in expected.items():
        require(np.allclose(f[key],value,rtol=0,atol=2e-13),f"Incorrect initial {case} {key}")


def analyze(case, cfg, axes, frames):
    g=float(cfg["gamma"]);dimension=3 if case=="cloud" else 2
    dv=np.prod([axis[1]-axis[0] for axis in axes[:dimension]])
    records=[]
    phase=None
    if case=="acoustic":
        k=np.array([2*np.pi/(float(cfg[f"{d}_upper"])-float(cfg[f"{d}_lower"])) for d in "xy"])
        x,y=np.meshgrid(axes[0],axes[1]);phase=k[0]*(x-float(cfg["x_lower"]))+k[1]*(y-float(cfg["y_lower"]))
        amplitude=1e-4*np.sinc(k[0]*(axes[0][1]-axes[0][0])/(2*np.pi))*np.sinc(k[1]*(axes[1][1]-axes[1][0])/(2*np.pi))
        omega=np.linalg.norm(k)
    for f in frames:
        r,p=f["rho"],f["p"]
        u,v,w=(f[key] for key in ("vx","vy","vz"))
        record=dict(time=f["time"],step=f["step"],mass=float(r.sum()*dv),
                    energy=float(np.sum(p/(g-1)+.5*r*(u*u+v*v+w*w))*dv),
                    rho_min=float(r.min()),rho_max=float(r.max()),p_min=float(p.min()),p_max=float(p.max()))
        if case=="acoustic":
            ref=amplitude*np.sin(phase-omega*f["time"])
            coef=2j*np.mean((r[0]-1)*np.exp(-1j*phase))/amplitude
            record.update(amplitude_ratio=float(abs(coef)),phase_error=float(np.angle(coef*np.exp(1j*omega*f["time"]))),
                          rho_L1=float(np.mean(abs(r[0]-1-ref))),rho_L1_over_amplitude=float(np.mean(abs(r[0]-1-ref))/1e-4))
        elif case=="rm":
            interface=rm_interface(axes,f)
            ly=float(cfg["y_upper"])-float(cfg["y_lower"])
            phase_y=2*np.pi*(axes[1]-float(cfg["y_lower"]))/ly
            coef=2*np.mean(interface*np.exp(-1j*phase_y))
            vort=np.gradient(v[0],axes[0],axis=1)-np.gradient(u[0],axes[1],axis=0)
            record.update(interface_mean=float(interface.mean()),interface_amplitude=float(abs(coef)),
                          vorticity_rms=float(np.sqrt(np.mean(vort*vort))))
        elif case=="cloud":
            weight=np.maximum(r-3,0)
            require(weight.sum()>0,"High-density core disappeared at this resolution")
            record.update(dense_core_x=float(np.sum(weight*axes[0][None,None,:])/weight.sum()),
                          dense_core_vx=float(np.sum(weight*u)/weight.sum()))
            symmetry=0.
            for key in ("rho","p","vx","vy","vz"):
                for axis,component in ((0,"vz"),(1,"vy")):
                    sign=-1 if key==component else 1
                    symmetry=max(symmetry,float(np.max(abs(f[key]-sign*np.flip(f[key],axis)))/max(1.,np.max(abs(f[key])))))
            record["transverse_symmetry_error"]=symmetry
        records.append(record)
    initial,final=records[0],records[-1]
    summary=dict(case=case,shape=[len(a) for a in axes[:dimension]],steps=final["step"],time=final["time"],snapshots=len(frames))
    if case=="acoustic":
        for key in ("mass","energy"):
            drift=max(abs(row[key]/initial[key]-1) for row in records)
            require(drift<1e-10,f"Periodic {key} is not conserved")
            summary[key+"_relative_drift"]=drift
        require(max(abs(row["phase_error"]) for row in records)<.08,"Acoustic phase propagation failed")
        require(.3<final["amplitude_ratio"]<1.02,"Excessive acoustic damping/growth")
        require(final["rho_L1_over_amplitude"]<.5,"Acoustic error exceeds coarse first-order tolerance")
        summary.update({key:final[key] for key in ("amplitude_ratio","phase_error","rho_L1","rho_L1_over_amplitude")})
    elif case=="sod":
        _,r,u,p=sod_average(axes,cfg,frames[-1]["time"])
        f=frames[-1];normal=(2*f["vx"][0]+f["vy"][0])/np.sqrt(5)
        tangent=(-f["vx"][0]+2*f["vy"][0])/np.sqrt(5)
        for key,num,ref,limit in (("rho",f["rho"][0],r,.03),("normal_velocity",normal,u,.07),("pressure",f["p"][0],p,.04)):
            error=float(np.mean(abs(num-ref)));summary[key+"_L1"]=error
            require(error<limit,f"Oblique Sod {key} L1 error {error:g} exceeds {limit}")
        summary["tangential_velocity_rms"]=float(np.sqrt(np.mean(tangent*tangent)))
        require(summary["tangential_velocity_rms"]<.03,"Excessive spurious tangential velocity")
    elif case=="rm":
        # The shock initially compresses the corrugation. Growth must be measured
        # relative to its POSTSHOCK amplitude, not to the uncompressed initial one.
        lx=float(cfg["x_upper"])-float(cfg["x_lower"])
        ly=float(cfg["y_upper"])-float(cfg["y_lower"])
        impact=.175*lx/(2*np.sqrt(1.4))
        window=[row for row in records if impact<=row["time"]<=impact+.2*ly/(2*np.sqrt(1.4))]
        require(len(window)>=2,"Need snapshots resolving shock-interface passage")
        compressed=min(window,key=lambda row:row["interface_amplitude"])
        late=[row for row in records if row["time"]>=max(impact+.2*ly/(2*np.sqrt(1.4)),.7*final["time"])]
        require(len(late)>=3,"Need several postshock snapshots to measure RM growth")
        slope=float(np.polyfit([row["time"] for row in late],[row["interface_amplitude"] for row in late],1)[0])
        summary.update(interface_amplitude_initial=initial["interface_amplitude"],interface_amplitude_final=final["interface_amplitude"],
                       interface_amplitude_postshock=compressed["interface_amplitude"],postshock_reference_time=compressed["time"],
                       late_amplitude_growth_rate=slope,vorticity_rms_final=final["vorticity_rms"])
        require(final["interface_amplitude"]>1.2*compressed["interface_amplitude"] and slope>0,"No sustained postshock RM interface growth detected")
        require(final["vorticity_rms"]>.05,"No baroclinic vorticity generated")
    elif case=="cloud":
        summary.update(dense_core_displacement=final["dense_core_x"]-initial["dense_core_x"],
                       dense_core_vx=final["dense_core_vx"],peak_density=max(row["rho_max"] for row in records),
                       transverse_symmetry_error=max(row["transverse_symmetry_error"] for row in records))
        require(summary["peak_density"]>6,"No cloud compression")
        require(summary["dense_core_displacement"]>.03 and final["dense_core_vx"]>.1,"No cloud acceleration")
        require(summary["transverse_symmetry_error"]<1e-8,"Cloud transverse symmetry failed")
    return records,summary


def plot(case,cfg,axes,frames,records,output,animate):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.animation import PillowWriter,FFMpegWriter,writers
    x,y,z=axes;iz=(len(z)-1)//2
    fig,axs=plt.subplots(1,3,figsize=(13.2,4.4),layout="constrained")
    extent=(float(cfg["x_lower"]),float(cfg["x_upper"]),float(cfg["y_lower"]),float(cfg["y_upper"]))
    if case=="acoustic":
        fields=[f["rho"][0]-1 for f in frames];vmin,vmax=-1e-4,1e-4
        second=[f["p"][0]-1/1.4 for f in frames];labels=("Density perturbation","Pressure perturbation");cmap="RdBu_r"
    elif case=="rm":
        fields=[f["rho"][0] for f in frames];vmin=min(a.min() for a in fields);vmax=max(a.max() for a in fields)
        second=[np.gradient(f["vy"][0],x,axis=1)-np.gradient(f["vx"][0],y,axis=0) for f in frames]
        labels=("Density / entropy-proxy interface","Vorticity");cmap="inferno"
    else:
        fields=[f["rho"][iz] for f in frames];second=[f["p"][iz] for f in frames]
        vmin=min(a.min() for a in fields);vmax=max(a.max() for a in fields)
        labels=("Density"+(f" at z={z[iz]:.4f}" if case=="cloud" else ""),"Pressure");cmap="inferno"
    im=axs[0].imshow(fields[0],origin="lower",extent=extent,vmin=vmin,vmax=vmax,cmap=cmap,interpolation="nearest")
    lo=min(a.min() for a in second);hi=max(a.max() for a in second)
    if case in ("rm","acoustic"): hi=max(abs(lo),abs(hi),1e-12);lo=-hi
    im2=axs[1].imshow(second[0],origin="lower",extent=extent,vmin=lo,vmax=hi,
                       cmap="RdBu_r" if case in ("rm","acoustic") else "magma",interpolation="nearest")
    fig.colorbar(im,ax=axs[0],shrink=.75);fig.colorbar(im2,ax=axs[1],shrink=.75)
    for ax,label in zip(axs[:2],labels):ax.set(title=label,xlabel="x",ylabel="y")
    times=np.array([f["time"] for f in frames]);overlay=None
    if case=="sod":
        cx=(float(cfg["x_lower"])+float(cfg["x_upper"]))/2
        cy=(float(cfg["y_lower"])+float(cfg["y_upper"]))/2
        position=(2*(x[None,:]-cx)+(y[:,None]-cy))/np.sqrt(5)
        scatter=axs[2].scatter(position.ravel()[::12],fields[0].ravel()[::12],s=2,alpha=.35,label="Numerical cells")
        ss=np.linspace(position.min(),position.max(),500)
        curve,=axs[2].plot(ss,sod_points(ss,0)[0],color="black",lw=1.2,label="Planar exact solution")
        axs[2].set(xlabel="Coordinate along wave normal",ylabel="Density",ylim=(.08,1.05),title="Collapse onto Sod profile")
        axs[2].legend(fontsize=8)
    else:
        key={"acoustic":"amplitude_ratio","rm":"interface_amplitude","cloud":"dense_core_x"}[case]
        values=[r[key] for r in records]
        axs[2].plot(times,values,color="#176b87")
        marker,=axs[2].plot(times[:1],values[:1],"o",color="#c24d2c")
        axs[2].set(xlabel="Time",ylabel=key,title={"acoustic":"Wave amplitude / initial amplitude","rm":"Interface mode amplitude (proxy)","cloud":"High-density core position (proxy)"}[case])
        if case=="rm":overlay,=axs[0].plot(rm_interface(axes,frames[0]),y,color="cyan",lw=.8)
    axs[2].grid(alpha=.2)
    title=fig.suptitle("")
    names={"acoustic":"Periodic diagonal acoustic wave","sod":"Oblique Sod, 26.565 degrees","rm":"Richtmyer-Meshkov, incident Mach 2","cloud":"3D shock-cloud, incident Mach 2 (central slice)"}
    def update(i):
        f=frames[i];im.set_data(fields[i]);im2.set_data(second[i]);title.set_text(f"{names[case]} | t={f['time']:.5f}")
        if case=="sod":
            scatter.set_offsets(np.column_stack((position.ravel()[::12],fields[i].ravel()[::12])))
            curve.set_ydata(sod_points(ss,f["time"])[0])
        else:
            marker.set_data([times[i]],[values[i]])
            if overlay is not None:overlay.set_xdata(rm_interface(axes,f))
    diagnostics=output/"diagnostics";update(len(frames)-1)
    fig.savefig(diagnostics/"final.png",dpi=120)
    if animate:
        directory=output/"animations";directory.mkdir(exist_ok=True)
        # Limit movies, retaining endpoints; all snapshots still enter diagnostics.
        chosen=np.unique(np.linspace(0,len(frames)-1,min(len(frames),60),dtype=int))
        exports=[(PillowWriter(fps=10),directory/f"{case}.gif")]
        if writers.is_available("ffmpeg"):
            exports.append((FFMpegWriter(fps=10,codec="libx264",extra_args=["-pix_fmt","yuv420p"]),directory/f"{case}.mp4"))
        with ExitStack() as stack:
            for writer,path in exports:stack.enter_context(writer.saving(fig,str(path),dpi=100))
            for k,i in enumerate(chosen):
                update(i)
                for writer,_ in exports:writer.grab_frame()
                if k%20==0:print(f"{case}: rendered {k+1}/{len(chosen)}",flush=True)
    plt.close(fig)


def process(case,path,executable,run=False,animate=False):
    cfg=parameters(path);output=path.parent/cfg["output_directory"];log=output/"run.log"
    require(float(cfg["gamma"])==1.4 and cfg.get("constant_gravity_enabled","false").lower()=="false","Requires gamma=1.4, no gravity")
    if run:
        output.mkdir(parents=True,exist_ok=True)
        with log.open("w") as stream:
            result=subprocess.run([str(executable),"-par",str(path)],cwd=path.parent,stdout=stream,stderr=subprocess.STDOUT)
        require(result.returncode==0,f"Simulation failed: {log}")
    paths=snapshot_files(output,"mesh_"+cfg["output_suffix"],log)
    frames=[]
    for source in paths:
        axes,frame=read_frame(source,cfg);frames.append(frame)
    require(frames[0]["time"]==0 and frames[0]["step"]==0,"Missing initial frame")
    require(np.isclose(frames[-1]["time"],float(cfg["t_end"]),rtol=0,atol=1e-12),"Incomplete run")
    check_initial(case,cfg,axes,frames[0])
    records,summary=analyze(case,cfg,axes,frames)
    diagnostics=output/"diagnostics";diagnostics.mkdir(exist_ok=True)
    with (diagnostics/"metrics.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    (diagnostics/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("PASS "+json.dumps(summary),flush=True)
    plot(case,cfg,axes,frames,records,output,animate)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",choices=("acoustic","sod","rm","cloud"))
    parser.add_argument("--parameters",type=Path,required=True)
    parser.add_argument("--executable",type=Path)
    parser.add_argument("--run",action="store_true")
    parser.add_argument("--animate",action="store_true")
    parser.add_argument("--convergence",action="store_true")
    args=parser.parse_args();path=args.parameters.resolve();exe=(args.executable or path.parent/"flare").resolve()
    if args.convergence:
        require(args.case=="acoustic","Convergence target currently supports acoustics only")
        cfg=parameters(path)
        base=process(args.case,path,exe,run=not (path.parent/cfg["output_directory"]/"run.log").exists())
        directory=(path.parent/cfg["output_directory"] / "convergence").resolve();directory.mkdir(parents=True,exist_ok=True)
        cfg["nxa"]=str(2*int(cfg["nxa"]));cfg["nya"]=str(2*int(cfg["nya"]))
        cfg["output_directory"]=str(directory/"fine");cfg["output_step_every"]="40"
        fine_file=directory/"parameters.txt";fine_file.write_text("\n".join(f"{k} = {v}" for k,v in cfg.items())+"\n")
        fine=process(args.case,fine_file,exe,run=True)
        ratio=fine["rho_L1"]/base["rho_L1"]
        require(ratio<.8,f"Acoustic error did not decrease sufficiently on refinement: {ratio}")
        report=dict(coarse=base,fine=fine,error_ratio=ratio,observed_order=float(-np.log2(ratio)))
        (directory/"convergence.json").write_text(json.dumps(report,indent=2)+"\n")
        print("PASS acoustic convergence: "+json.dumps(report),flush=True)
    else:
        process(args.case,path,exe,args.run,args.animate)


if __name__=="__main__":main()
