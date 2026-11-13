#!/usr/bin/env bash
# Writes every MATPOWER case as a literal-only case file the ybus_powerflow reader accepts, into
#   testdata/matpower/<matpower-dir-name>-normalized/. Some MATPOWER cases compute their data
#   with code, which the reader rejects; this evaluates them with MATPOWER's own loadcase under
#   Octave. See scripts/normalize_cases.m.
#
# Usage: scripts/normalize-cases.sh
#   Needs Octave on PATH. Uses scripts/fetch-matpower.sh for MATPOWER (so MATPOWER_DIR is
#   honoured). Prints the output directory. Safe to rerun: output that is complete and was made
#   from the same MATPOWER directory is reused.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATPOWER="$("${ROOT}/scripts/fetch-matpower.sh")"
MATPOWER="$(cd "${MATPOWER}" && pwd -P)"
OUT="${ROOT}/testdata/matpower/$(basename "${MATPOWER}")-normalized"
MARKER=".ybus-normalized"

# Complete output from the same MATPOWER directory: nothing to do
if [[ -f "${OUT}/${MARKER}" ]] && [[ "$(head -n 1 "${OUT}/SOURCE")" == "source ${MATPOWER}" ]]; then
    echo "${OUT}"
    exit 0
fi

if ! command -v octave > /dev/null; then
    echo "error: normalizing MATPOWER cases needs Octave (for example: apt install octave)" >&2
    exit 1
fi

# Build in a scratch directory and move it into place only when complete, so an interrupted
#   run never leaves output that looks finished
PART="${OUT}.part"
rm -rf "${PART}"
mkdir -p "${PART}"
octave --no-gui --quiet "${ROOT}/scripts/normalize_cases.m" "${MATPOWER}" "${PART}" >&2
touch "${PART}/${MARKER}"
rm -rf "${OUT}"
mv "${PART}" "${OUT}"
echo "${OUT}"
