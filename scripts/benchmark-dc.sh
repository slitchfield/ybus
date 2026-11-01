#!/usr/bin/env bash
# Surveys the DC power flow over every MATPOWER case: times each step, records peak memory and
#   compares the results with MATPOWER's. The results in docs/benchmarks/ come from this. Each
#   case runs in its own process (powerflow's dc_survey example) under a time limit and an
#   address-space limit, so a case too large for the current solvers is recorded as such
#   instead of stopping the survey.
#
# Usage: scripts/benchmark-dc.sh [<output.csv>]
#   Writes one CSV line per case, smallest case file first, to <output.csv> (default
#   target/dc-survey.csv) and to stdout as it goes. Needs GNU time (/usr/bin/time) and the
#   normalized cases (scripts/normalize-cases.sh, which needs Octave unless they're already
#   made). MATPOWER's own results are generated into testdata/ by scripts/generate-reference.sh,
#   which needs Octave; without it, the comparison columns stay empty. Never commit those
#   results (docs/PowerFlowChecklist.md, Phase 1.5).
#
#   YBUS_BENCH_TIMEOUT  the time limit per case, in seconds (default 3600)
#   YBUS_BENCH_MEMORY   the address-space limit per case, in GiB (default 10)
#   YBUS_BENCH_CASES    the cases to survey, separated by spaces (default: every case)
#   YBUS_BENCH_ORDERING the factorization's ordering, Cholesky or LU: natural, rcm or amd
#                       (default amd, as DcModel::solve uses)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
REFERENCE="${NORMALIZED%-normalized}-reference"
OUT="${1:-${ROOT}/target/dc-survey.csv}"
TIMEOUT="${YBUS_BENCH_TIMEOUT:-3600}"
MEMORY_GIB="${YBUS_BENCH_MEMORY:-10}"
ORDERING="${YBUS_BENCH_ORDERING:-amd}"

if [[ ! -x /usr/bin/time ]]; then
    echo "error: the survey needs GNU time at /usr/bin/time (for example: apt install time)" >&2
    exit 1
fi

if [[ -n "${YBUS_BENCH_CASES:-}" ]]; then
    read -r -a CASES <<< "${YBUS_BENCH_CASES}"
else
    mapfile -t CASES < <(ls -Sr "${NORMALIZED}"/*.m | xargs -n 1 basename -s .m)
fi

# MATPOWER's results, for the cases that don't have them yet. The gen file is written last.
missing=()
for name in "${CASES[@]}"; do
    [[ -f "${REFERENCE}/${name}_dcpf_gen.csv" ]] || missing+=("${name}")
done
if [[ ${#missing[@]} -gt 0 ]]; then
    if command -v octave > /dev/null; then
        "${ROOT}/scripts/generate-reference.sh" "${REFERENCE}" "${missing[@]}" >&2 ||
            echo "warning: MATPOWER failed on some cases, so they aren't compared" >&2
    else
        echo "warning: without Octave, results aren't compared with MATPOWER's" >&2
    fi
fi

cargo build --release --quiet --manifest-path "${ROOT}/Cargo.toml" -p powerflow \
    --example dc_survey
SURVEY="${ROOT}/target/release/examples/dc_survey"

HEADER="$("${SURVEY}" --header)"
# A failed case's line: its name, then an empty field for every other column
EMPTY="${HEADER//[^,]/}"
LOG="$(mktemp -d)"
trap 'rm -rf "${LOG}"' EXIT
mkdir -p "$(dirname "${OUT}")"
echo "${HEADER},status,peak_mib" | tee "${OUT}"
for name in "${CASES[@]}"; do
    # exec replaces the shell, so the limits and GNU time's measurement apply to the survey
    if line="$(/usr/bin/time -f %M -o "${LOG}/rss" timeout "${TIMEOUT}" \
        bash -c 'ulimit -v "$1" && exec "${@:2}"' _ "$((MEMORY_GIB * 1024 * 1024))" \
        "${SURVEY}" --ordering "${ORDERING}" "${NORMALIZED}/${name}.m" "${REFERENCE}" \
        2> "${LOG}/err")"; then
        status=ok
    else
        code=$?
        line="${name}${EMPTY}"
        if [[ ${code} -eq 124 ]]; then
            status="over the ${TIMEOUT} s limit"
        elif grep -q "memory allocation of" "${LOG}/err"; then
            status="over the ${MEMORY_GIB} GiB limit"
        else
            status="error: $(head -n 1 "${LOG}/err" | tr ',' ';')"
        fi
    fi
    # GNU time writes the peak resident set in KiB, after any note about the exit status
    echo "${line},${status},$(($(tail -n 1 "${LOG}/rss") / 1024))" | tee -a "${OUT}"
done
