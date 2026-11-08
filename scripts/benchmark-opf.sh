#!/usr/bin/env bash
# Surveys the DC or AC optimal power flow over every MATPOWER case with costs: times it, records
#   its outcome, iterations and peak memory, and compares it with MATPOWER's results, from MP-Core
#   and from the legacy implementation. The OPF results in docs/benchmarks/ come from this. Each
#   case runs in its own process (powerflow's opf_survey example) under a time limit and an
#   address-space limit, so a case that fails is recorded as such instead of stopping the survey.
#
# Usage: scripts/benchmark-opf.sh [<output.csv>]
#   Writes one CSV line per case, smallest case file first, to <output.csv> (default
#   target/opf-survey-<model>-<solver>-<kkt>.csv) and to stdout as it goes. Needs GNU time (/usr/bin/time) and the
#   normalized cases (scripts/normalize-cases.sh, which needs Octave unless they're already made).
#   MATPOWER's results are generated into testdata/ by scripts/generate-reference.sh, which needs
#   Octave; without it, the comparison columns stay empty. Never commit those results
#   (docs/PowerFlowChecklist.md, Phase 1.5).
#
#   YBUS_BENCH_MODEL    dc or ac (default dc)
#   YBUS_BENCH_SOLVER   mips or mehrotra (default mips, MATPOWER's): the OPF's solver, at
#                       MATPOWER's OPF settings
#   YBUS_BENCH_KKT      lu, ldlt or ldlt-uncorrected (default lu): how the solver factors the
#                       KKT matrix; ldlt-uncorrected is LDLᵀ without the predictor-corrector
#                       method's inertia correction
#   YBUS_BENCH_TIMEOUT  the time limit per case, in seconds (default 3600)
#   YBUS_BENCH_MEMORY   the address-space limit per case, in GiB (default 10)
#   YBUS_BENCH_CASES    the cases to survey, separated by spaces (default: every case with costs)
#
#   MATPOWER's results (MP-Core, legacy, and its timings) are generated for any case that lacks
#   them, one case per process under the same limits, with BLAS single-threaded, which Octave
#   needs under an address-space limit.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
REFERENCE="${NORMALIZED%-normalized}-reference"
MODEL="${YBUS_BENCH_MODEL:-dc}"
SOLVER="${YBUS_BENCH_SOLVER:-mips}"
case "${SOLVER}" in
    mips|mehrotra) ;;
    *) echo "error: YBUS_BENCH_SOLVER must be mips or mehrotra" >&2; exit 2 ;;
esac
KKT="${YBUS_BENCH_KKT:-lu}"
case "${KKT}" in
    lu|ldlt|ldlt-uncorrected) ;;
    *) echo "error: YBUS_BENCH_KKT must be lu, ldlt or ldlt-uncorrected" >&2; exit 2 ;;
esac
OUT="${1:-${ROOT}/target/opf-survey-${MODEL}-${SOLVER}-${KKT}.csv}"
TIMEOUT="${YBUS_BENCH_TIMEOUT:-3600}"
MEMORY_GIB="${YBUS_BENCH_MEMORY:-10}"
case "${MODEL}" in
    dc) RUN=dcopf ;;
    ac) RUN=opf ;;
    *) echo "error: YBUS_BENCH_MODEL must be dc or ac" >&2; exit 2 ;;
esac

if [[ ! -x /usr/bin/time ]]; then
    echo "error: the survey needs GNU time at /usr/bin/time (for example: apt install time)" >&2
    exit 1
fi

if [[ -n "${YBUS_BENCH_CASES:-}" ]]; then
    read -r -a CASES <<< "${YBUS_BENCH_CASES}"
else
    # Every case with costs, which the OPF needs
    mapfile -t CASES < <(ls -Sr "${NORMALIZED}"/*.m | xargs grep -l "mpc.gencost" |
        xargs -n 1 basename -s .m)
fi

# MATPOWER's results, for the cases that don't have them all yet (the timing file is written
#   before the summary, and the legacy run's after MP-Core's), one case per process
if command -v octave > /dev/null; then
    for name in "${CASES[@]}"; do
        [[ -f "${REFERENCE}/${name}_${RUN}_legacy_summary.csv" &&
            -f "${REFERENCE}/${name}_${RUN}_time.csv" ]] && continue
        (
            ulimit -v "$((MEMORY_GIB * 1024 * 1024))"
            YBUS_REFERENCE_RUNS="${RUN},${RUN}_legacy" timeout -s KILL "${TIMEOUT}" \
                "${ROOT}/scripts/generate-reference.sh" "${REFERENCE}" "${name}" < /dev/null >&2
        ) || echo "warning: MATPOWER failed on ${name}, so it isn't compared" >&2
    done
else
    echo "warning: without Octave, results aren't compared with MATPOWER's" >&2
fi

cargo build --release --quiet --manifest-path "${ROOT}/Cargo.toml" -p powerflow \
    --example opf_survey
SURVEY="${ROOT}/target/release/examples/opf_survey"

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
        "${SURVEY}" --model "${MODEL}" --solver "${SOLVER}" --kkt "${KKT}" \
        "${NORMALIZED}/${name}.m" "${REFERENCE}" \
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
