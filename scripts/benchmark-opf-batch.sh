#!/usr/bin/env bash
# Surveys batches of OPF solves serially and in parallel (phase O5, target 1): for each case,
#   model, solver and batch, ybus_powerflow's opf_batch example times the batch cold and warm,
#   serially once, then with its par_ form on rayon pools of each thread count, checking every
#   parallel result against the serial one bit for bit. Each run is its own process under a time
#   limit and an address-space limit, so a run that fails is recorded as such instead of stopping
#   the survey.
#
# Usage: scripts/benchmark-opf-batch.sh [<output.csv>]
#   Writes opf_batch's lines, with a status column, to <output.csv> (default
#   target/opf-batch.csv) and to stdout as it goes. Needs the normalized cases
#   (scripts/normalize-cases.sh, which needs Octave unless they're already made).
#
#   YBUS_BATCH_CASES     the cases, separated by spaces (default: case2383wp case3120sp
#                        case9241pegase)
#   YBUS_BATCH_MODELS    dc, ac or both (default: dc ac)
#   YBUS_BATCH_SOLVERS   mips, mehrotra or both (default: mips mehrotra)
#   YBUS_BATCH_BATCHES   outages, hours or both (default: outages hours)
#   YBUS_BATCH_OUTAGES   the branch and the generator outages each (default 14)
#   YBUS_BATCH_THREADS   the thread counts, separated by commas (default 1,2,4,8,14,28)
#   YBUS_BATCH_SKIP      runs to leave out, as case:model:solver words separated by spaces
#   YBUS_BENCH_TIMEOUT   the time limit per run, in seconds (default 7200)
#   YBUS_BENCH_MEMORY    the address-space limit per run, in GiB (default 12). Every thread's
#                        solve holds its own model and factorization, so many threads on a large
#                        case need more; keep the limit below the machine's physical memory, so
#                        a run that needs more fails cleanly instead of exhausting the machine
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
OUT="${1:-${ROOT}/target/opf-batch.csv}"
read -r -a CASES <<< "${YBUS_BATCH_CASES:-case2383wp case3120sp case9241pegase}"
read -r -a MODELS <<< "${YBUS_BATCH_MODELS:-dc ac}"
read -r -a SOLVERS <<< "${YBUS_BATCH_SOLVERS:-mips mehrotra}"
read -r -a BATCHES <<< "${YBUS_BATCH_BATCHES:-outages hours}"
OUTAGES="${YBUS_BATCH_OUTAGES:-14}"
THREADS="${YBUS_BATCH_THREADS:-1,2,4,8,14,28}"
SKIP=" ${YBUS_BATCH_SKIP:-} "
TIMEOUT="${YBUS_BENCH_TIMEOUT:-7200}"
MEMORY_GIB="${YBUS_BENCH_MEMORY:-12}"

cargo build --release --quiet --manifest-path "${ROOT}/Cargo.toml" -p ybus_powerflow \
    --features rayon --example opf_batch
BATCH="${ROOT}/target/release/examples/opf_batch"

HEADER="$("${BATCH}" --header)"
# A failed run's line: its key (case, model, solver, batch), then an empty field for every
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
            for batch in "${BATCHES[@]}"; do
                # exec replaces the shell, so the limits apply to the example itself
                if timeout "${TIMEOUT}" \
                    bash -c 'ulimit -v "$1" && exec "${@:2}"' _ \
                    "$((MEMORY_GIB * 1024 * 1024))" "${BATCH}" --model "${model}" \
                    --solver "${solver}" --batch "${batch}" --outages "${OUTAGES}" \
                    --threads "${THREADS}" "${NORMALIZED}/${name}.m" \
                    > "${LOG}/out" 2> "${LOG}/err"; then
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
                    # The lines finished before the failure, then the failure
                    sed 's/$/,ok/' "${LOG}/out" | tee -a "${OUT}"
                    echo "${name},${model},${solver},${batch}${EMPTY},${status}" |
                        tee -a "${OUT}"
                fi
            done
        done
    done
done
