#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
exec make -C unit_test test TESTS=wave_timestep "$@"
