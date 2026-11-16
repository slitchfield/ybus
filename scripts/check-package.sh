#!/usr/bin/env bash
# Dry-runs the publish of every publishable crate, as crates.io would see them: each crate is
#   packaged, then built from its package against the packages of the crates before it, in
#   dependency order (docs/PublishingChecklist.md, P4). ybus_ffi is skipped by its
#   publish = false. Nothing is uploaded.
#
# Usage: scripts/check-package.sh
#   Fails if cargo does, and exits 1 if cargo warns about anything but the accepted kinds below.
#
# The accepted warnings (docs/PublishingChecklist.md, P3, item 1):
#   - "ignoring test|example|benchmark ... not included in the published package": the packages
#     hold only src/, the README and the license files, so cargo drops every other target from
#     the published manifest, and says so once for each. There's no way to silence them.
#   - "aborting upload due to dry run": the dry run's own notice, once for each crate.
#   Any other warning, such as "manifest has no description", fails the check. (A license file
#   left out of an include list isn't one: with an SPDX license, cargo doesn't look for it.)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

LOG="$(mktemp)"
trap 'rm -f "${LOG}"' EXIT

# cargo's progress and warnings go to stderr: show them as they come, and keep a copy. With
#   pipefail, cargo's own failure fails the script.
cargo publish --workspace --dry-run --locked 2>&1 | tee "${LOG}"

ACCEPTED='^warning: (ignoring (test|example|benchmark) `[^`]+` as `[^`]+` is not included in the published package|aborting upload due to dry run)$'
UNEXPECTED="$(grep '^warning:' "${LOG}" | grep -Ev "${ACCEPTED}" || true)"
if [[ -n "${UNEXPECTED}" ]]; then
    echo >&2
    echo "error: cargo warned about something other than the accepted warnings:" >&2
    echo "${UNEXPECTED}" >&2
    exit 1
fi
echo "Packaged and verified every publishable crate," \
    "with $(grep -c '^warning: ignoring' "${LOG}") accepted warnings for targets left out"
