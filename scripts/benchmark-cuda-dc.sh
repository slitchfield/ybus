#!/usr/bin/env bash
# Surveys the DC power flow's GPU path against the CPU's, for the CUDA spike: batched DC
#   solves, PTDF columns, and single-branch outages solved whole and screened, on each case
#   (ybus_powerflow's cuda_dc_survey example, which documents the columns). Each case runs in its
#   own process under a time limit and a memory limit, so a case that fails is recorded as such
#   instead of stopping the survey.
#
# Usage: scripts/benchmark-cuda-dc.sh [<output.csv>]
#   Appends the survey's CSV lines, as each case finishes, to <output.csv> (default
#   target/cuda-dc-survey.csv, kept out of /tmp so it survives a crash), writing the header
#   first if the file is new. Needs a CUDA device, the normalized cases
#   (scripts/normalize-cases.sh) and systemd-run for the memory limit: CUDA reserves far more
#   address space than it uses, so `ulimit -v` can't limit it, while a cgroup limits what's
#   resident, page-locked memory included.
#
#   YBUS_BENCH_TIMEOUT  the time limit per case, in seconds (default 1800)
#   YBUS_BENCH_MEMORY   the memory limit per case, in GiB (default 10); no swap
#   YBUS_BENCH_CASES    the cases to survey, separated by spaces (default: a spread of sizes
#                       from case118 to case_SyntheticUSA)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
OUT="${1:-${ROOT}/target/cuda-dc-survey.csv}"
TIMEOUT="${YBUS_BENCH_TIMEOUT:-1800}"
MEMORY_GIB="${YBUS_BENCH_MEMORY:-10}"
read -r -a CASES <<< "${YBUS_BENCH_CASES:-case118 case300 case1354pegase case2383wp \
    case_ACTIVSg2000 case9241pegase case13659pegase case_ACTIVSg25k case_ACTIVSg70k \
    case_SyntheticUSA}"

cargo build --release --quiet --manifest-path "${ROOT}/Cargo.toml" -p ybus_powerflow \
    --example cuda_dc_survey --features cuda,rayon
SURVEY="${ROOT}/target/release/examples/cuda_dc_survey"

mkdir -p "$(dirname "${OUT}")"
if [[ ! -s "${OUT}" ]]; then
    "${SURVEY}" --header > "${OUT}"
fi
for name in "${CASES[@]}"; do
    echo "== ${name}" >&2
    status=0
    timeout "${TIMEOUT}" systemd-run --user --scope --quiet \
        -p MemoryMax="${MEMORY_GIB}G" -p MemorySwapMax=0 \
        "${SURVEY}" "${NORMALIZED}/${name}.m" | tee -a "${OUT}" || status=$?
    if [[ "${status}" -ne 0 ]]; then
        echo "${name},failed with status ${status}" | tee -a "${OUT}"
    fi
done
