#!/usr/bin/env python3
"""Write profile.dat explicitly; edit initial_state() to customize the initial gas."""
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
RHO_H = 1.0
PRESSURE_H = 1.0
PRESSURE_HE = 0.1


def initial_state(x, y, x_mid):
    """Primitive values at cell center (x, y); HI and He_I are mass fractions."""
    hydrogen = x < x_mid
    # H-1 and He-4: equal nucleus number densities require rho_He = 4*rho_H.
    # Unequal pressures drive the flow; initially T_H/T_He = P_H/P_He = 10.
    return {
        "rho": RHO_H if hydrogen else 4.0 * RHO_H,
        "vx": 0.0,
        "vy": 0.0,
        "vz": 0.0,
        "pre": PRESSURE_H if hydrogen else PRESSURE_HE,
        "HI": 1.0 if hydrogen else 0.0,
        "He_I": 0.0 if hydrogen else 1.0,
    }


def main():
    parameters = {}
    for line in (CASE_DIR / "parameters.txt").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            parameters[key.strip()] = value.strip()

    nx, ny = int(parameters["nxa"]), int(parameters["nya"])
    x0, x1 = float(parameters["x_lower"]), float(parameters["x_upper"])
    y0, y1 = float(parameters["y_lower"]), float(parameters["y_upper"])
    x_mid = 0.5 * (x0 + x1)
    gamma = float(parameters["gamma"])
    columns = ("rho", "vx", "vy", "vz", "pre", "gamma", "HI", "He_I")
    output = CASE_DIR / parameters["init_profile_file_name"]
    with output.open("w") as stream:
        stream.write(" ".join(columns) + "\n")
        for j in range(ny):
            y = y0 + (j + 0.5) * (y1 - y0) / ny
            for i in range(nx):  # x varies fastest; no ghost-cell rows.
                x = x0 + (i + 0.5) * (x1 - x0) / nx
                state = initial_state(x, y, x_mid)
                state["gamma"] = gamma
                stream.write(" ".join(f"{state[key]:.17g}" for key in columns) + "\n")
    print(f"Wrote {nx} x {ny} = {nx * ny} active-cell rows to {output}")


if __name__ == "__main__":
    main()
