#!/usr/bin/env bash
# Regenerates the committed DC and AC power flow reference results in
#   powerflow/tests/data/matpower-reference/ by running MATPOWER's rundcpf and runpf under Octave.
#   Needs Octave on PATH; fetches MATPOWER with scripts/fetch-matpower.sh (or uses MATPOWER_DIR).
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
OUT="${ROOT}/powerflow/tests/data/matpower-reference"
if [[ $# -gt 0 ]]; then
    if [[ $# -lt 2 ]]; then
        echo "usage: $0 [<output-dir> <case>...]" >&2
        exit 2
    fi
    OUT="$1"
    shift
fi
mkdir -p "${OUT}"

octave --no-gui --quiet "${ROOT}/scripts/generate_reference.m" "${MATPOWER}" "${OUT}" "$@"
