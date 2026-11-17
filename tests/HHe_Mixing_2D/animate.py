#!/usr/bin/env python3
"""Render composition and fluid colormap movies from this run's CSV manifest."""
import argparse
from contextlib import ExitStack
from pathlib import Path

import numpy as np

from mixing import (COMPOSITION, HYDRO, check_grid, edges, load_frame, require,
                    run_paths, signature)


def render(parameter_file, output=None, stride=1, max_frames=60, fps=10, dpi=90, formats=("gif",)):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.animation import FFMpegWriter, PillowWriter, writers

    if "mp4" in formats:
        require(writers.is_available("ffmpeg"), "MP4 needs ffmpeg; use --format gif")
    cfg, run_output, paths = run_paths(parameter_file)
    output = output or run_output / "animations"
    selected = paths[::stride]
    if selected[-1] != paths[-1]:
        selected.append(paths[-1])
    if len(selected) > max_frames:
        selected = [selected[i] for i in np.linspace(0, len(selected)-1, max_frames, dtype=int)]
    require(len(selected) >= 2, "Need at least two saved states for an animation")
    reference, _ = load_frame(selected[0])
    limits = {key: [np.inf, -np.inf] for key, _, _ in (*COMPOSITION, *HYDRO)}
    stamps = []
    for path in selected:
        stamp = signature(path)
        frame, fields = load_frame(path)
        check_grid(frame, reference)
        require(signature(path) == stamp, f"Snapshot changed: {path}")
        stamps.append(stamp)
        for key, bounds in limits.items():
            bounds[0] = min(bounds[0], float(fields[key].min()))
            bounds[1] = max(bounds[1], float(fields[key].max()))
    for key, bounds in limits.items():
        if key in ("H", "He", "mixing"):
            bounds[:] = [0, 1]
        elif key in ("vx", "vy"):
            span = max(abs(bounds[0]), abs(bounds[1]), 1e-12)
            bounds[:] = [-span, span]
        elif bounds[0] == bounds[1]:
            pad = max(abs(bounds[0])*.01, 1e-12)
            bounds[:] = [bounds[0]-pad, bounds[1]+pad]
    xedges, yedges = edges(reference.x), edges(reference.y)
    extent = (xedges[0], xedges[-1], yedges[0], yedges[-1])
    output.mkdir(parents=True, exist_ok=True)
    for name, descriptors, shape, size in (("hhe_mixing", COMPOSITION, (1, 3), (12.6, 4.8)),
                                            ("hydro_fields", HYDRO, (2, 4), (14.4, 8.2))):
        fig, axes = plt.subplots(*shape, figsize=size, layout="constrained")
        title = fig.suptitle("H/He Kelvin-Helmholtz mixing", fontsize=14)
        artists = []
        for ax, (key, label, cmap) in zip(np.asarray(axes).flat, descriptors):
            artist = ax.imshow(np.zeros((len(reference.y), len(reference.x))), origin="lower",
                               extent=extent, aspect="equal", interpolation="nearest", cmap=cmap,
                               vmin=limits[key][0], vmax=limits[key][1])
            ax.set(title=label, xlabel="x [code]", ylabel="y [code]")
            fig.colorbar(artist, ax=ax, shrink=.78, pad=.025)
            artists.append(artist)
        fig.supxlabel("Fixed color scales; ghosts excluded; one saved state per frame\n"
                      "T is derived for neutral gas; mixing includes finite-volume numerical diffusion",
                      fontsize=9)
        try:
            with ExitStack() as stack:
                exports = []
                for suffix in formats:
                    # Matplotlib's h264 mode enforces even dimensions and yuv420p
                    # for browser/QuickTime compatibility (libx264 alone does not).
                    writer = PillowWriter(fps=fps) if suffix == "gif" else FFMpegWriter(
                        fps=fps, codec="h264", extra_args=["-movflags", "+faststart"])
                    target = output / f"{name}.{suffix}"
                    stack.enter_context(writer.saving(fig, str(target), dpi=dpi))
                    exports.append(writer)
                for index, (path, stamp) in enumerate(zip(selected, stamps)):
                    require(signature(path) == stamp, f"Snapshot changed since scan: {path}")
                    frame, fields = load_frame(path)
                    check_grid(frame, reference)
                    require(signature(path) == stamp, f"Snapshot changed while reading: {path}")
                    title.set_text(f"H/He Kelvin-Helmholtz | t = {frame.time:.4f} [code] | step {frame.step}")
                    for artist, (key, _, _) in zip(artists, descriptors):
                        artist.set_data(fields[key])
                    for writer in exports:
                        writer.grab_frame()
                    if index == 0:
                        fig.savefig(output / f"{name}_first.png", dpi=dpi)
                    if index % 10 == 0 or index == len(selected)-1:
                        print(f"{name}: {index+1}/{len(selected)} frames", flush=True)
                fig.savefig(output / f"{name}_last.png", dpi=dpi)
        finally:
            plt.close(fig)
    print(f"Saved movies and first/last PNGs to {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=Path(__file__).resolve().with_name("parameters.txt"))
    parser.add_argument("--output", type=Path, help="Default: run output/animations")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=60, help="Bound GIF memory; keep first and last")
    parser.add_argument("--fps", type=float, default=10)
    parser.add_argument("--dpi", type=int, default=90)
    parser.add_argument("--format", choices=("gif", "mp4", "both"), default="gif")
    args = parser.parse_args()
    if args.stride < 1 or args.max_frames < 2 or not np.isfinite(args.fps) or not 0 < args.fps <= 50 or args.dpi < 1:
        parser.error("Require stride >= 1, max-frames >= 2, 0 < fps <= 50 and dpi >= 1")
    try:
        render(args.parameters, args.output, args.stride, args.max_frames, args.fps, args.dpi,
               ("gif", "mp4") if args.format == "both" else (args.format,))
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f"Animation failed: {error}\n")


if __name__ == "__main__":
    main()
