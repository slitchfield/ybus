#!/usr/bin/env bash
# Surveys the AC power flow over every MATPOWER case: times each step, records convergence,
#   iterations and peak memory, and compares the results with MATPOWER's runpf. The AC results in
#   docs/benchmarks/ come from this. Each case runs in its own process (powerflow's ac_survey
#   example) under a time limit and an address-space limit, so a case that fails is recorded as
#   such instead of stopping the survey.
#
# Usage: scripts/benchmark-ac.sh [<output.csv>]
#   Writes one CSV line per case, smallest case file first, to <output.csv> (default
#   target/ac-survey.csv) and to stdout as it goes. Needs GNU time (/usr/bin/time) and the
#   normalized cases (scripts/normalize-cases.sh, which needs Octave unless they're already
#   made). MATPOWER's own results are generated into testdata/ by scripts/generate-reference.sh,
#   which needs Octave; without it, the comparison columns stay empty. Never commit those
#   results (docs/PowerFlowChecklist.md, Phase 1.5).
#
#   YBUS_BENCH_TIMEOUT  the time limit per case, in seconds (default 3600)
#   YBUS_BENCH_MEMORY   the address-space limit per case, in GiB (default 10)
#   YBUS_BENCH_CASES    the cases to survey, separated by spaces (default: every case)
#   YBUS_BENCH_METHOD   the solution method: nr (Newton's method, the default), fdxb or fdbx
#   YBUS_BENCH_QLIM     reactive power limits, as runpf's pf.enforce_q_lims: 0 (ignored, the
#                       default), 1 (all at once) or 2 (one at a time); 1 and 2 need nr
#
#   MATPOWER's results are generated for any case that lacks the results being compared (its
#   Newton results, or with limits, its results with them), including a case MATPOWER fails on,
#   which is retried on each run.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
REFERENCE="${NORMALIZED%-normalized}-reference"
OUT="${1:-${ROOT}/target/ac-survey.csv}"
TIMEOUT="${YBUS_BENCH_TIMEOUT:-3600}"
MEMORY_GIB="${YBUS_BENCH_MEMORY:-10}"
METHOD="${YBUS_BENCH_METHOD:-nr}"
QLIM="${YBUS_BENCH_QLIM:-0}"
case "${QLIM}" in
    0) RUN=acpf; FLAGS=() ;;
    1 | 2) RUN="qlim${QLIM}"; FLAGS=(--qlim "${QLIM}") ;;
    *) echo "error: YBUS_BENCH_QLIM must be 0, 1 or 2" >&2; exit 2 ;;
esac

if [[ ! -x /usr/bin/time ]]; then
    echo "error: the survey needs GNU time at /usr/bin/time (for example: apt install time)" >&2
    exit 1
fi

if [[ -n "${YBUS_BENCH_CASES:-}" ]]; then
    read -r -a CASES <<< "${YBUS_BENCH_CASES}"
else
    mapfile -t CASES < <(ls -Sr "${NORMALIZED}"/*.m | xargs -n 1 basename -s .m)
fi

# MATPOWER's results, for the cases that don't have them yet. The summary is written last.
missing=()
for name in "${CASES[@]}"; do
    [[ -f "${REFERENCE}/${name}_${RUN}_summary.csv" ]] || missing+=("${name}")
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
    --example ac_survey
SURVEY="${ROOT}/target/release/examples/ac_survey"

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
        "${SURVEY}" --method "${METHOD}" "${FLAGS[@]}" "${NORMALIZED}/${name}.m" "${REFERENCE}" \
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
