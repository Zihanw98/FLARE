#!/usr/bin/env python3
"""Render actual FLARE isosurfaces, with lighting, camera orbit and OBJ export.

python3 render3d.py tests/Spherical_Blast_3D/output mesh_blast3d --obj-sequence

Requires numpy, matplotlib, Pillow. MP4 additionally needs ffmpeg/libx264.
Reuses animate3d's validated CSV reader. No Blender installation is required.
Uses marching tetrahedra (six consistently connected tetrahedra per cube),
linear interpolation between cell centers, and optionally between saved times.
Method: https://paulbourke.net/geometry/polygonise/ (tetrahedra section).
This implementation does not copy the reference's code or lookup tables.
The fixed-threshold surface is not a tracked material interface/shock position.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from animate3d import edges, field_values, read_snapshot, snapshot_files


def isosurface(axes, scalar, color, level):
    """Return oriented triangles [face, vertex, (x,y,z,color)] of scalar=level.

    Normals point out of scalar>=level. No surface is added at the grid boundary.
    Cell-center samples are interpreted as vertices of the rendering lattice.
    """
    offsets = ((0,0,0), (1,0,0), (1,1,0), (0,1,0),
               (0,0,1), (1,0,1), (1,1,1), (0,1,1))
    shape = tuple(n - 1 for n in scalar.shape)
    base = np.indices(shape).reshape(3, -1).T
    points, values, colors = [], [], []
    for offset in offsets:
        idx = base + offset
        points.append(np.stack([axes[d][idx[:, d]] for d in range(3)], axis=-1))
        values.append(scalar[tuple(idx.T)])
        colors.append(color[tuple(idx.T)])
    points = np.stack(points, axis=1)
    values, colors = np.stack(values, axis=1), np.stack(colors, axis=1)
    triangles = []
    for tet in ((0,1,2,6), (0,2,3,6), (0,3,7,6),
                (0,7,4,6), (0,4,5,6), (0,5,1,6)):
        p, v, c = points[:, tet], values[:, tet], colors[:, tet]
        cases = np.sum((v >= level) * np.array([1,2,4,8]), axis=1)
        for case in range(1, 15):
            selected = cases == case
            if not selected.any():
                continue
            inside = [i for i in range(4) if case & (1 << i)]
            outside = [i for i in range(4) if not case & (1 << i)]
            pp, vv, cc = p[selected], v[selected], c[selected]
            outward = pp[:, outside].mean(axis=1) - pp[:, inside].mean(axis=1)

            def crossing(a, b):
                weight = (level - vv[:, a]) / (vv[:, b] - vv[:, a])
                xyz = pp[:, a] + weight[:, None] * (pp[:, b] - pp[:, a])
                rgb = cc[:, a] + weight * (cc[:, b] - cc[:, a])
                return np.column_stack((xyz, rgb))

            if len(inside) == 1:
                faces = [np.stack([crossing(inside[0], b) for b in outside], axis=1)]
            elif len(outside) == 1:
                faces = [np.stack([crossing(a, outside[0]) for a in inside], axis=1)]
            else:
                a, b = inside
                c0, d = outside
                ac, ad, bc, bd = crossing(a, c0), crossing(a, d), crossing(b, c0), crossing(b, d)
                faces = [np.stack((ac, ad, bd), axis=1), np.stack((ac, bd, bc), axis=1)]
            for face in faces:
                normal = np.cross(face[:, 1, :3] - face[:, 0, :3],
                                  face[:, 2, :3] - face[:, 0, :3])
                flip = np.sum(normal * outward, axis=1) < 0
                face[flip] = face[flip][:, [0, 2, 1]]
                triangles.append(face[np.linalg.norm(normal, axis=1) > 1e-15])
    return np.concatenate(triangles) if triangles else np.empty((0, 3, 4))


def clip_halfspace(triangles, axis, position, less=True):
    """Clip triangles, interpolating both coordinates and color on cut edges."""
    if len(triangles) == 0:
        return triangles
    distance = (triangles[:, :, axis] - position) * (1 if less else -1)
    inside = distance <= 0
    whole = triangles[inside.all(axis=1)]
    partial = triangles[inside.any(axis=1) & ~inside.all(axis=1)]
    clipped = []
    for triangle in partial:
        polygon = []
        previous = triangle[-1]
        previous_d = (previous[axis] - position) * (1 if less else -1)
        for current in triangle:
            current_d = (current[axis] - position) * (1 if less else -1)
            if (current_d <= 0) != (previous_d <= 0):
                polygon.append(previous + previous_d / (previous_d - current_d) * (current - previous))
            if current_d <= 0:
                polygon.append(current)
            previous, previous_d = current, current_d
        for j in range(1, len(polygon) - 1):
            clipped.append([polygon[0], polygon[j], polygon[j+1]])
    return np.concatenate((whole, np.asarray(clipped).reshape(-1, 3, 4)))


def cutaway(triangles, center):
    # Remove the front quarter, x>center_x and y<center_y; leave cut edges open.
    left = clip_halfspace(triangles, 0, center[0])
    right = clip_halfspace(triangles, 0, center[0], less=False)
    back = clip_halfspace(right, 1, center[1], less=False)
    return np.concatenate((left, back))


def write_obj(path, triangles):
    """Weld coincident vertices for smooth shading/import in Blender."""
    vertices, inverse = np.unique(np.round(triangles[:, :, :3].reshape(-1, 3), 12),
                                  axis=0, return_inverse=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write("# FLARE fixed-threshold isosurface; physical coordinates\n")
        np.savetxt(stream, vertices, fmt="v %.12g %.12g %.12g")
        np.savetxt(stream, inverse.reshape(-1, 3) + 1, fmt="f %d %d %d")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", type=Path)
    parser.add_argument("prefix")
    parser.add_argument("--field", choices=("p", "rho", "speed"), default="p")
    parser.add_argument("--level", type=float, default=1.12)
    parser.add_argument("--frames", type=int, default=72, help="Uniform physical-time samples; fields interpolate linearly")
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--full", action="store_true", help="Keep the full surface instead of removing its front quarter")
    parser.add_argument("--obj-sequence", action="store_true", help="Export full surfaces at original saved times for Blender")
    parser.add_argument("--preview", action="store_true", help="Render only final PNG and OBJ for quick inspection")
    args = parser.parse_args()
    if args.frames < 2 or args.fps < 1 or not np.isfinite(args.level):
        parser.error("frames>=2, fps>=1 and a finite level are required")
    output = args.output_dir or args.directory / "animations" / "surface3d"
    output.mkdir(parents=True, exist_ok=True)
    paths = snapshot_files(args.directory, args.prefix)
    snapshots = [read_snapshot(path) for path in paths]
    axes = snapshots[0].axes
    for snap in snapshots[1:]:
        if any(not np.array_equal(a, b) for a, b in zip(axes, snap.axes)):
            raise ValueError("Grid changed between saved frames")
    times = np.array([snap.time for snap in snapshots])
    bounds = np.array([(edges(axis)[0], edges(axis)[-1]) for axis in axes])
    center = bounds.mean(axis=1)
    scalars = [field_values(snap, args.field) for snap in snapshots]
    # Fix framing over the whole run; never auto-zoom independently per frame.
    occupied = np.logical_or.reduce([scalar >= args.level for scalar in scalars])
    indices = np.nonzero(occupied)
    view_bounds = bounds.copy()
    if len(indices[0]):
        for d in range(3):
            lo, hi = axes[d][indices[d].min()], axes[d][indices[d].max()]
            pad = max((hi - lo) * .22, np.diff(axes[d]).max() * 2)
            view_bounds[d] = lo - pad, hi + pad
    low = min(snap.fields["rho"].min() for snap in snapshots)
    high = max(snap.fields["rho"].max() for snap in snapshots)
    if low == high:
        low, high = low - 1e-6, high + 1e-6
    final_full = isosurface(axes, scalars[-1], snapshots[-1].fields["rho"], args.level)
    if not len(final_full):
        raise ValueError("No final isosurface at this level; choose a threshold inside the data range")
    write_obj(output / "surface_final.obj", final_full)
    write_obj(output / "surface_final_cutaway.obj", cutaway(final_full, center))
    manifest = dict(field=args.field, level=args.level, color="rho", color_range=[low, high],
                    grid=[len(a) for a in axes], source_snapshots=len(snapshots),
                    time_interpolation="linear between saved scalar/density fields", cutaway=not args.full)
    if args.obj_sequence:
        manifest["meshes"] = []
        for i, snap in enumerate(snapshots):
            name = f"meshes/frame_{i:04d}.obj"
            mesh = isosurface(axes, scalars[i], snap.fields["rho"], args.level)
            write_obj(output / name, mesh)
            manifest["meshes"].append(dict(path=name, time=snap.time, step=snap.step))
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.animation import PillowWriter, FFMpegWriter, writers
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from contextlib import ExitStack

    background = "#081321"
    plt.rcParams.update({"text.color": "#e7eef5", "font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(10, 8), facecolor=background)
    ax = fig.add_axes([.03, .09, .88, .77], projection="3d", facecolor=background)
    ax.set_axis_off()
    ax.set_box_aspect(view_bounds[:, 1] - view_bounds[:, 0])
    for d, setter in enumerate((ax.set_xlim, ax.set_ylim, ax.set_zlim)):
        setter(view_bounds[d])
    cmap, norm = plt.get_cmap("YlOrRd"), Normalize(low, high)
    color_axis = fig.add_axes([.90, .30, .018, .36])
    bar = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=color_axis)
    bar.set_label("Density on surface", color="#cad7e6")
    bar.ax.tick_params(colors="#cad7e6")
    bar.outline.set_edgecolor("#40546b")
    fig.text(.06, .94, "SPHERICAL BLAST", fontsize=24, weight="bold")
    fig.text(.06, .90, f"3D {args.field} = {args.level:g} surface  /  {'full' if args.full else 'quarter cutaway'}", fontsize=12, color="#99b6d0")
    clock = fig.text(.06, .045, "", fontsize=12, family="monospace")
    fig.text(.06, .015, f"{len(axes[0])} x {len(axes[1])} x {len(axes[2])} cells  |  fixed threshold  |  interpolated saved fields", fontsize=9, color="#8fa7bd")
    # A faint base grid anchors the fixed physical scale while the camera turns.
    for x in np.linspace(*view_bounds[0], 6):
        ax.plot([x, x], view_bounds[1], [view_bounds[2,0]]*2, color="#20354a", lw=.6)
    for y in np.linspace(*view_bounds[1], 6):
        ax.plot(view_bounds[0], [y, y], [view_bounds[2,0]]*2, color="#20354a", lw=.6)
    artist = None

    def draw(time, phase):
        nonlocal artist
        j = min(max(np.searchsorted(times, time, side="right") - 1, 0), len(times) - 1)
        k = min(j + 1, len(times) - 1)
        weight = 0 if k == j else (time - times[j]) / (times[k] - times[j])
        field = (1 - weight) * scalars[j] + weight * scalars[k]
        rho = (1 - weight) * snapshots[j].fields["rho"] + weight * snapshots[k].fields["rho"]
        triangles = isosurface(axes, field, rho, args.level)
        if not args.full:
            triangles = cutaway(triangles, center)
        if artist is not None:
            artist.remove()
        if len(triangles):
            normal = np.cross(triangles[:,1,:3] - triangles[:,0,:3], triangles[:,2,:3] - triangles[:,0,:3])
            normal /= np.maximum(np.linalg.norm(normal, axis=1)[:,None], 1e-20)
            light = np.array([-.4, -.6, 1.])
            light /= np.linalg.norm(light)
            shade = .35 + .65 * np.abs(normal @ light)
            colors = cmap(norm(triangles[:,:,3].mean(axis=1)))
            colors[:,:3] *= shade[:,None]
            artist = Poly3DCollection(triangles[:,:,:3], facecolors=colors, edgecolors="none", linewidths=0,
                                      antialiased=False, zsort="average")
            ax.add_collection3d(artist)
        else:
            artist = None
        ax.view_init(elev=23 + 7 * np.sin(phase * np.pi), azim=-65 + 40 * phase)
        clock.set_text(f"t = {time:.5f}     |     surface triangles: {len(triangles):,}")

    if not args.preview:
        gif = output / "blast_surface.gif"
        movie_writers = [(PillowWriter(fps=args.fps), gif)]
        if writers.is_available("ffmpeg"):
            movie_writers.append((FFMpegWriter(fps=args.fps, codec="libx264", bitrate=2600,
                                               extra_args=["-pix_fmt", "yuv420p"]), output / "blast_surface.mp4"))
        with ExitStack() as stack:
            for writer, filename in movie_writers:
                stack.enter_context(writer.saving(fig, str(filename), dpi=100))
            for frame, time in enumerate(np.linspace(times[0], times[-1], args.frames)):
                draw(time, frame / (args.frames - 1))
                for writer, _ in movie_writers:
                    writer.grab_frame()
                if frame % 12 == 0:
                    print(f"Rendered {frame + 1}/{args.frames} frames", flush=True)
    draw(times[-1], 1.)
    fig.savefig(output / "blast_surface_last.png", dpi=140, facecolor=background)
    plt.close(fig)
    print(f"Saved renders and Blender-compatible OBJ models in {output}")


if __name__ == "__main__":
    main()
