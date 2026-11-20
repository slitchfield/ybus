#!/usr/bin/env bash
# Surveys warm starts on OPF sequences: for each case, model, solver, sequence and warm shift,
#   ybus_powerflow's opf_sequence example solves the sequence's load levels with a kept model
#   three ways per step (cold, primal, warm) and compares the primal and warm steps with the cold
#   one by O3's bar. Each run is its own process under a time limit and an address-space limit,
#   so a run that fails is recorded as such instead of stopping the survey.
#
# Usage: scripts/benchmark-opf-sequence.sh [<output.csv>]
#   Writes opf_sequence's lines, with a status column, to <output.csv> (default
#   target/opf-sequence.csv) and to stdout as it goes. Needs the normalized cases
#   (scripts/normalize-cases.sh, which needs Octave unless they're already made).
#
#   YBUS_SEQUENCE_CASES      the cases, separated by spaces (default: case9 case30 case118)
#   YBUS_SEQUENCE_MODELS     dc, ac or both (default: dc ac)
#   YBUS_SEQUENCE_SOLVERS    mips, mehrotra or both (default: mips mehrotra)
#   YBUS_SEQUENCE_SEQUENCES  any of sweep, walk, branch-outages and generator-outages
#                            (default: sweep walk)
#   YBUS_SEQUENCE_SHIFTS     the warm shifts, separated by commas, one warm chain each, in one
#                            run with the cold and primal chains (default: the solvers' own)
#   YBUS_SEQUENCE_SKIP       runs to leave out, as case:model:solver words separated by spaces
#                            (for example case_ACTIVSg25k:ac:mips)
#   YBUS_BENCH_TIMEOUT       the time limit per run, in seconds (default 3600)
#   YBUS_BENCH_MEMORY        the address-space limit per run, in GiB (default 10)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
OUT="${1:-${ROOT}/target/opf-sequence.csv}"
read -r -a CASES <<< "${YBUS_SEQUENCE_CASES:-case9 case30 case118}"
read -r -a MODELS <<< "${YBUS_SEQUENCE_MODELS:-dc ac}"
read -r -a SOLVERS <<< "${YBUS_SEQUENCE_SOLVERS:-mips mehrotra}"
read -r -a SEQUENCES <<< "${YBUS_SEQUENCE_SEQUENCES:-sweep walk}"
SHIFTS="${YBUS_SEQUENCE_SHIFTS:-}"
shift_args=()
[[ -n "${SHIFTS}" ]] && shift_args=(--warm-shifts "${SHIFTS}")
SKIP=" ${YBUS_SEQUENCE_SKIP:-} "
TIMEOUT="${YBUS_BENCH_TIMEOUT:-3600}"
MEMORY_GIB="${YBUS_BENCH_MEMORY:-10}"

cargo build --release --quiet --manifest-path "${ROOT}/Cargo.toml" -p ybus_powerflow \
    --example opf_sequence
SEQUENCE="${ROOT}/target/release/examples/opf_sequence"

HEADER="$("${SEQUENCE}" --header)"
# A failed run's line: its key (case, model, solver, sequence), then an empty field for every
#   other column
EMPTY="$(cut -d, -f5- <<< "${HEADER}")"
EMPTY=",${EMPTY//[^,]/}"
LOG="$(mktemp -d)"
trap 'rm -rf "${LOG}"' EXIT
mkdir -p "$(dirname "${OUT}")"
echo "${HEADER},status" | tee "${OUT}"
for name in "${CASES[@]}"; do
    for model in "${MODELS[@]}"; do
        for solver in "${SOLVERS[@]}"; do
            [[ "${SKIP}" == *" ${name}:${model}:${solver} "* ]] && continue
            for sequence in "${SEQUENCES[@]}"; do
                # exec replaces the shell, so the limits apply to the example itself
                if timeout "${TIMEOUT}" \
                    bash -c 'ulimit -v "$1" && exec "${@:2}"' _ \
                    "$((MEMORY_GIB * 1024 * 1024))" "${SEQUENCE}" --model "${model}" \
                    --solver "${solver}" --sequence "${sequence}" "${shift_args[@]}" \
                    "${NORMALIZED}/${name}.m" > "${LOG}/out" 2> "${LOG}/err"; then
                    sed 's/$/,ok/' "${LOG}/out" | tee -a "${OUT}"
                else
                    code=$?
                    if [[ ${code} -eq 124 ]]; then
                        status="over the ${TIMEOUT} s limit"
                    elif grep -q "memory allocation of" "${LOG}/err"; then
                        status="over the ${MEMORY_GIB} GiB limit"
                    else
                        status="error: $(head -n 1 "${LOG}/err" | tr ',' ';')"
                    fi
                    # The steps finished before the failure, then the failure
                    sed 's/$/,ok/' "${LOG}/out" | tee -a "${OUT}"
                    echo "${name},${model},${solver},${sequence}${EMPTY},${status}" |
                        tee -a "${OUT}"
                fi
            done
        done
    done
done
