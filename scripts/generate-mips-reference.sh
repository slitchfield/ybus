#!/usr/bin/env bash
# Records MIPS's results on the test problems of MATPOWER's MIPS tests, for ybus_ipm's
#   tests: runs scripts/generate_mips_reference.m under Octave, with the MATPOWER that
#   scripts/fetch-matpower.sh fetches, into ybus_ipm/tests/data/mips_reference.txt. The
#   problems are MATPOWER's tests (BSD-licensed); they hold no case data.
#
# Usage: scripts/generate-mips-reference.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATPOWER="$("${ROOT}/scripts/fetch-matpower.sh")"
OUT="${ROOT}/ybus_ipm/tests/data/mips_reference.txt"

{
    echo "# MIPS's results on MATPOWER's MIPS test problems, from scripts/generate_mips_reference.m"
    echo "# $(octave --no-gui --quiet --eval "addpath(genpath('${MATPOWER}')); v = mipsver('all'); printf('MIPS %s, MATPOWER %s, Octave %s', v.Version, mpver, OCTAVE_VERSION)" 2>/dev/null)"
    octave --no-gui --quiet "${ROOT}/scripts/generate_mips_reference.m" "${MATPOWER}"
} > "${OUT}"
echo "wrote ${OUT}"
