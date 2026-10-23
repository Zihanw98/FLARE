"""Render 1-D snapshots as 2-D strips (not a 2-D simulation).

Examples:
  python3 plot_snapshots.py --step 50
  python3 plot_snapshots.py --animate --fps 10
  python3 plot_snapshots.py --step 0 --field rho --ghost
"""

import argparse
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")  # Works without a graphical desktop.
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np


FIELDS = {"rho": (0, "Density", "viridis"),
          "vx": (1, "Velocity", "coolwarm"),
          "p": (2, "Pressure", "magma")}


def read_snapshot(path):
    with path.open() as stream:
        header = stream.readline().strip()
    if not header.startswith("#"):
        raise ValueError(f"{path}: missing metadata header")
    meta = dict(re.findall(r"(\w+)=([^\s]+)", header))
    try:
        meta = {"time": float(meta["time"]),
                **{key: int(meta[key]) for key in ("timestep", "nx", "nxg", "nxa")}}
    except (KeyError, ValueError) as error:
        raise ValueError(f"{path}: invalid metadata") from error
    if (not np.isfinite(meta["time"]) or meta["timestep"] < 0
            or meta["nxa"] <= 0 or meta["nxg"] < 0
            or meta["nx"] != meta["nxa"] + 2 * meta["nxg"]):
        raise ValueError(f"{path}: invalid time or grid layout")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape != (meta["nx"], 3) or not np.all(np.isfinite(data)):
        raise ValueError(f"{path}: expected {meta['nx']} finite rows of rho vx p")
    return meta, data


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--step", type=int, help="Snapshot timestep, e.g. 50")
    mode.add_argument("--animate", action="store_true", help="All snapshots to GIF")
    parser.add_argument("--input-dir", type=Path, default=Path("output"))
    parser.add_argument("--output", type=Path, help="PNG or GIF destination")
    parser.add_argument("--field", choices=["all", *FIELDS], default="all")
    parser.add_argument("--length", type=float, default=100.0,
                        help="Active length Xa (not stored in txt); default: 100")
    parser.add_argument("--ghost", action="store_true", help="Include ghost cells")
    parser.add_argument("--fps", type=float, default=10.0,
                        help="GIF frames/second; snapshots have equal display duration")
    parser.add_argument("--stride", type=int, default=1,
                        help="Animate every Nth snapshot; always keep first/last")
    parser.add_argument("--modelname",type=str,default="test")
    args = parser.parse_args()
    if (not np.isfinite(args.length) or args.length <= 0
            or not np.isfinite(args.fps) or args.fps <= 0 or args.stride < 1):
        parser.error("length and fps must be finite positive values; stride must be >=1")

    try:
        if args.animate:
            paths = [p for p in args.input_dir.glob(f"{args.modelname}_*.txt")
                     if re.fullmatch(rf"{args.modelname}_\d+\.txt", p.name)]
            paths.sort(key=lambda p: int(p.stem.split("_")[1]))
            if not paths:
                raise ValueError(f"No {args.modelname}_*.txt snapshots in {args.input_dir}")
        else:
            paths = [args.input_dir / f"{args.modelname}_{args.step}.txt"]

        snapshots = [read_snapshot(p) for p in paths]
        for path, (meta, _) in zip(paths, snapshots):
            if meta["timestep"] != int(path.stem.split("_")[1]):
                raise ValueError(f"{path}: filename and header timestep disagree")
        first = snapshots[0][0]
        for i, (meta, _) in enumerate(snapshots):
            if any(meta[k] != first[k] for k in ("nx", "nxg", "nxa")):
                raise ValueError("Snapshots have different grid layouts; use one run per folder")
            if i and meta["time"] <= snapshots[i-1][0]["time"]:
                raise ValueError("Snapshot times are not increasing; check for stale files")

        if not args.ghost:
            snapshots = [(m, d[m["nxg"]:m["nxg"] + m["nxa"]])
                         for m, d in snapshots]
        keys = list(FIELDS) if args.field == "all" else [args.field]
        dx = args.length / first["nxa"]
        padding = first["nxg"] * dx if args.ghost else 0.0
        extent = [-padding, args.length + padding, 0, 1]
        fig, axes = plt.subplots(len(keys), 1, figsize=(11, 2.5 * len(keys)),
                                 squeeze=False, layout="constrained")
        images = []
        for ax, key in zip(axes[:, 0], keys):
            column, label, cmap = FIELDS[key]
            # Fixed color scale over ALL snapshots avoids misleading flicker.
            low = min(float(d[:, column].min()) for _, d in snapshots)
            high = max(float(d[:, column].max()) for _, d in snapshots)
            if low == high:
                delta = max(abs(low) * 0.01, 1e-12)
                low, high = low - delta, high + delta
            strip = np.repeat(snapshots[0][1][:, column][None, :], 32, axis=0)
            image = ax.imshow(strip, origin="lower", aspect="auto", extent=extent,
                              interpolation="nearest", cmap=cmap, vmin=low, vmax=high)
            images.append(image)
            ax.set_title(label)
            ax.set_xlabel("x")
            ax.set_yticks([])
            if args.ghost and padding:
                for boundary in (0, args.length):
                    ax.axvline(boundary, color="white", linestyle="--", linewidth=1)
            fig.colorbar(image, ax=ax, label=label, pad=0.02)
        title = fig.suptitle("")

        def update(index):
            meta, data = snapshots[index]
            for image, key in zip(images, keys):
                image.set_data(np.repeat(data[:, FIELDS[key][0]][None, :], 32, axis=0))
            title.set_text(f"1-D shock tube | step={meta['timestep']} | t={meta['time']:.6g}")
            return [*images, title]

        default_name = (f"tube_{args.field}.gif" if args.animate
                        else f"tube_{args.step}_{args.field}.png")
        destination = args.output or Path("figures") / default_name
        suffix = ".gif" if args.animate else ".png"
        if destination.suffix.lower() != suffix:
            raise ValueError(f"Output must have the {suffix} extension")
        destination.parent.mkdir(parents=True, exist_ok=True)
        if args.animate:
            indices = list(range(0, len(snapshots), args.stride))
            if indices[-1] != len(snapshots) - 1:
                indices.append(len(snapshots) - 1)
            animation = FuncAnimation(fig, update, frames=indices,
                                      interval=1000 / args.fps, blit=False)
            animation.save(destination, writer=PillowWriter(fps=args.fps), dpi=100)
        else:
            update(0)
            fig.savefig(destination, dpi=160)
        plt.close(fig)
        print(f"Saved {destination}")
    except (OSError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
