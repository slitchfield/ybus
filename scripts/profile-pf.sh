#!/usr/bin/env bash
# Profiles the power flows with samply, for phase O2's profiles: each workload on each case runs
#   in ybus_powerflow's profile_pf example (the whole pipeline repeated, so each stage weighs as in
#   one run), built with the profiling profile (release with debug line tables), and
#   scripts/profile_summary.py folds the samples into stages and top functions.
#
# Usage: scripts/profile-pf.sh [<output.md>]
#   Writes the summaries to <output.md> (default target/profiles/summary.md); the raw profiles
#   stay in target/profiles/, where `samply load <profile>.json` opens one in the Firefox
#   Profiler. Needs samply (cargo install samply), python3, the normalized cases
#   (scripts/normalize-cases.sh), and /proc/sys/kernel/perf_event_paranoid at 1 or lower.
#
#   YBUS_PROFILE_CASES      the cases, separated by spaces (default: case2383wp case9241pegase
#                           case_ACTIVSg70k case_SyntheticUSA)
#   YBUS_PROFILE_WORKLOADS  the workloads (default: dc nr fdxb fdbx qlim2 n1 n1-parallel)
#   YBUS_PROFILE_SECONDS    how long each run repeats its pipeline, at least once (default 5)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NORMALIZED="$("${ROOT}/scripts/normalize-cases.sh")"
OUT_DIR="${ROOT}/target/profiles"
OUT="${1:-${OUT_DIR}/summary.md}"
SAMPLY="$(command -v samply || echo "${HOME}/.cargo/bin/samply")"
read -r -a CASES <<< "${YBUS_PROFILE_CASES:-case2383wp case9241pegase case_ACTIVSg70k case_SyntheticUSA}"
read -r -a WORKLOADS <<< "${YBUS_PROFILE_WORKLOADS:-dc nr fdxb fdbx qlim2 n1 n1-parallel}"
SECONDS_PER_RUN="${YBUS_PROFILE_SECONDS:-5}"

if [[ ! -x "${SAMPLY}" ]]; then
    echo "error: needs samply (cargo install samply)" >&2
    exit 1
fi
if [[ "$(cat /proc/sys/kernel/perf_event_paranoid)" -gt 1 ]]; then
    echo "error: samply needs /proc/sys/kernel/perf_event_paranoid at 1 or lower" >&2
    exit 1
fi

cargo build --profile profiling --quiet --manifest-path "${ROOT}/Cargo.toml" -p ybus_powerflow \
    --features rayon --example profile_pf
DRIVER="${ROOT}/target/profiling/examples/profile_pf"

mkdir -p "${OUT_DIR}"
profiles=()
for name in "${CASES[@]}"; do
    for workload in "${WORKLOADS[@]}"; do
        base="${OUT_DIR}/${workload}-${name}"
        echo "${workload} ${name}" >&2
        "${SAMPLY}" record --save-only --unstable-presymbolicate -o "${base}.json" -- \
            "${DRIVER}" "${workload}" "${NORMALIZED}/${name}.m" "${SECONDS_PER_RUN}" \
            > /dev/null 2> "${base}.log" || { cat "${base}.log" >&2; continue; }
        profiles+=("${base}.json")
    done
done
python3 "${ROOT}/scripts/profile_summary.py" "${profiles[@]}" > "${OUT}"
echo "wrote ${OUT}" >&2
