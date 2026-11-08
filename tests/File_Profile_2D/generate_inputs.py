#!/usr/bin/env python3
"""Generate a nonuniform periodic mesh and an advected density wave."""
from math import cos, pi, sin
from pathlib import Path

ROOT = Path(__file__).resolve().parent
N = 32


def cosine_average(left, right, length, shift=0.0):
    k = 2 * pi / length
    return (sin(k * (right - shift)) - sin(k * (left - shift))) / (k * (right - left))


def main():
    x = [i / N + 0.05 * sin(2 * pi * i / N) for i in range(N + 1)]
    y = [0.75 * (i / N + 0.03 * sin(2 * pi * i / N)) for i in range(N + 1)]
    x[0], x[-1], y[0], y[-1] = 0.0, 1.0, 0.0, 0.75
    (ROOT / "mesh.dat").write_text("x y\n" + "".join(f"{a:.17g} {b:.17g}\n" for a, b in zip(x, y)))
    with (ROOT / "profile.dat").open("w") as stream:
        # Deliberately reorder columns to exercise the header-based reader.
        stream.write("gamma vy rho vz pre vx\n")
        for j in range(N):
            for i in range(N):  # x varies fastest.
                rho = 1 + 0.2 * cosine_average(x[i], x[i + 1], 1) * cosine_average(y[j], y[j + 1], 0.75)
                stream.write(f"1.4 -0.1 {rho:.17g} 0 1 0.2\n")


if __name__ == "__main__":
    main()
