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
#   - "crate <name>@<version> already exists on crates.io index": once a version is published,
#     main keeps it until the next release bumps it, so between releases every crate warns.
#     The dry run still packages and builds each one; only a real publish would refuse it,
#     and that runs only after a release's version bump.
#   Any other warning, such as "manifest has no description", fails the check. (A license file
#   left out of an include list isn't one: with an SPDX license, cargo doesn't look for it.)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

LOG="$(mktemp)"
trap 'rm -f "${LOG}"' EXIT

# cargo's progress and warnings go to stderr: show them as they come, and keep a copy. With
#   pipefail, cargo's own failure fails the script. --color never, since CI sets
#   CARGO_TERM_COLOR=always, and colour codes before "warning:" would hide every line from the
#   checks below.
cargo publish --workspace --dry-run --locked --color never 2>&1 | tee "${LOG}"

# The dry run says it aborted each upload. If no such line is found, the output wasn't read as
#   expected, and the warning check below would pass whatever cargo said.
if ! grep -q '^warning: aborting upload due to dry run$' "${LOG}"; then
    echo "error: found no 'aborting upload due to dry run' in cargo's output; can't check its warnings" >&2
    exit 1
fi

ACCEPTED='^warning: (ignoring (test|example|benchmark) `[^`]+` as `[^`]+` is not included in the published package|aborting upload due to dry run|crate ybus_[a-z]+@[0-9]+\.[0-9]+\.[0-9]+ already exists on crates\.io index)$'
UNEXPECTED="$(grep '^warning:' "${LOG}" | grep -Ev "${ACCEPTED}" || true)"
if [[ -n "${UNEXPECTED}" ]]; then
    echo >&2
    echo "error: cargo warned about something other than the accepted warnings:" >&2
    echo "${UNEXPECTED}" >&2
    exit 1
fi
echo "Packaged and verified every publishable crate," \
    "with $(grep -c '^warning: ignoring' "${LOG}") accepted warnings for targets left out"
