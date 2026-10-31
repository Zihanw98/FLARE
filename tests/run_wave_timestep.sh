#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT
for dimension in 1 2 3; do
    "${CXX:-g++}" -std=c++17 -O1 -g -Wall -Wextra -pedantic \
        -DHD_NDIM="$dimension" \
        -Isrc -Isrc/core -Isrc/hydro -Isrc/hydro/riemann -Isrc/user -Isrc/utils \
        tests/test_wave_timestep.cpp src/core/mesh.cpp src/utils/hdmath.cpp \
        src/hydro/flux.cpp src/hydro/facefield.cpp src/hydro/godunov.cpp \
        src/hydro/spatial.cpp src/hydro/timestep_CFL.cpp \
        src/hydro/riemann/exact.cpp src/hydro/riemann/factory.cpp \
        src/time/intergator.cpp src/time/euler.cpp \
        -o "$test_dir/check-$dimension"
    "$test_dir/check-$dimension"
done
