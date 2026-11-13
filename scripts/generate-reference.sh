#!/usr/bin/env bash
# Regenerates the committed DC and AC power flow reference results in
#   ybus_powerflow/tests/data/matpower-reference/ by running MATPOWER's rundcpf and runpf under
#   Octave. Needs Octave on PATH; fetches MATPOWER with scripts/fetch-matpower.sh (or uses
#   MATPOWER_DIR).
#
# Usage: scripts/generate-reference.sh [<output-dir> <case>...]
#   With no arguments, regenerates the committed results for case14, case30 and case118. The
#   output is deterministic for a given MATPOWER and Octave, so rerunning it should leave git
#   status clean; any diff is a real change in the results.
#
#   With arguments, writes the named cases' results into <output-dir> instead, for comparing
#   other cases locally, as scripts/benchmark-dc.sh does. Never commit those: only results from
#   cases with public sources are committed (docs/PowerFlowChecklist.md, Phase 1.5).
#
#   A case MATPOWER fails on is reported and skipped, and the script then exits with an error
#   after writing the others.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATPOWER="$("${ROOT}/scripts/fetch-matpower.sh")"
OUT="${ROOT}/ybus_powerflow/tests/data/matpower-reference"
if [[ $# -gt 0 ]]; then
    if [[ $# -lt 2 ]]; then
        echo "usage: $0 [<output-dir> <case>...]" >&2
        exit 2
    fi
    OUT="$1"
    shift
fi
mkdir -p "${OUT}"

# Single-threaded BLAS: under an address-space limit (ulimit -v), as local runs over every case
#   use, multithreaded OpenBLAS hangs Octave at no CPU (rundcopf on case14 under 4 GB, MATPOWER
#   8.1, Octave 8.4.0); single-threaded, the same call takes a fraction of a second
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}" OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

octave --no-gui --quiet "${ROOT}/scripts/generate_reference.m" "${MATPOWER}" "${OUT}" "$@"
