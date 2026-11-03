#!/usr/bin/env bash
# Compares LU orderings over every MATPOWER case: factors each case's AC Newton Jacobian (at a
#   flat start) and DC B′ under AMD and COLAMD, recording fill and time. The ordering results in
#   docs/benchmarks/ come from this. Each case runs in its own process (powerflow's
#   ordering_survey example) under a time limit and an address-space limit, so a case that fails
#   is recorded as such instead of stopping the survey.
#
# Usage: scripts/benchmark-ordering.sh [<output.csv>]
#   Writes one CSV line per case, smallest case file first, to <output.csv> (default
#   target/ordering-survey.csv) and to stdout as it goes. Needs GNU time (/usr/bin/time) and the
#   normalized cases (scripts/normalize-cases.sh, which needs Octave unless they're already
#   made). It needs no MATPOWER results.
#
#   YBUS_BENCH_TIMEOUT  the time limit per case, in seconds (default 3600)
#   YBUS_BENCH_MEMORY   the address-space limit per case, in GiB (default 10)
#   YBUS_BENCH_CASES    the cases to survey, separated by spaces (default: every case)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
OUT="${1:-${ROOT}/target/ordering-survey.csv}"
TIMEOUT="${YBUS_BENCH_TIMEOUT:-3600}"
MEMORY_GIB="${YBUS_BENCH_MEMORY:-10}"

if [[ ! -x /usr/bin/time ]]; then
    echo "error: the survey needs GNU time at /usr/bin/time (for example: apt install time)" >&2
    exit 1
fi

if [[ -n "${YBUS_BENCH_CASES:-}" ]]; then
    read -r -a CASES <<< "${YBUS_BENCH_CASES}"
else
    mapfile -t CASES < <(ls -Sr "${NORMALIZED}"/*.m | xargs -n 1 basename -s .m)
fi

cargo build --release --quiet --manifest-path "${ROOT}/Cargo.toml" -p powerflow \
    --example ordering_survey
SURVEY="${ROOT}/target/release/examples/ordering_survey"

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
        "${SURVEY}" "${NORMALIZED}/${name}.m" \
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
