#!/usr/bin/env bash
# Times dumped KKT matrices on the CPU's pivoted LDLᵀ and on cuDSS (ybus_cuda::GpuLdlt), for O5
# target 4's second gate: a refactorization and two refined solves each, at each thread count.
#
# Usage: scripts/benchmark-cudss-kkt.sh <output.csv> <dir>...
#   Each <dir> holds one case's KKT matrices (kkt-*.mtx) and names its label. The matrices are
#   dumped from interior-point solves (docs/benchmarks/o5-cudss-spike-*.md says how).
#
# Environment:
#   YBUS_KKT_THREADS   the CPU's thread counts, separated by spaces (default: 1 14)
#   YBUS_KKT_TIMEOUT   seconds allowed per case and thread count (default: 900)
#
# GPU runs can bring down WSL: commit first, and give the smallest case first. One case runs at
# a time, under a timeout.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:?usage: benchmark-cudss-kkt.sh <output.csv> <dir>...}"
shift
read -r -a THREADS <<< "${YBUS_KKT_THREADS:-1 14}"
TIMEOUT="${YBUS_KKT_TIMEOUT:-900}"
cargo build --release --quiet --manifest-path "${ROOT}/Cargo.toml" -p ybus_powerflow \
    --features cuda,rayon --example cudss_kkt
BIN="${ROOT}/target/release/examples/cudss_kkt"
"${BIN}" --header > "${OUT}"
for dir in "$@"; do
    label="$(basename "${dir}")"
    for t in "${THREADS[@]}"; do
        RAYON_NUM_THREADS="${t}" timeout "${TIMEOUT}" "${BIN}" "${label}" "${dir}"/kkt-*.mtx \
            | tee -a "${OUT}"
    done
done
