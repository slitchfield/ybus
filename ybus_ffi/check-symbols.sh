#!/usr/bin/env bash
# Checks the ybus library's exported functions against the committed record, symbols.txt, and
#   against the header, so that an export added, removed or renamed is a visible change in
#   review: before 1.0 such a change updates symbols.txt in the same PR, and from 1.0 only
#   additions are allowed (README, ABI policy). The record catches names; the committed,
#   verified header catches signatures, and the check that every export is declared in it, and
#   every declaration exported, keeps the two in step. check.sh and CI run this.
#
# Usage: ybus_ffi/check-symbols.sh [--update]
#   --update  rewrites symbols.txt from the built library instead of checking it
set -euo pipefail

crate="$(cd "$(dirname "$0")" && pwd)"
workspace="$(dirname "$crate")"
target="${CARGO_TARGET_DIR:-$workspace/target}"
record="$crate/symbols.txt"

cargo build --quiet -p ybus_ffi
# The functions the shared library defines and exports
exported="$(nm -D --defined-only "$target/debug/libybus.so" | awk '$2 == "T" { print $3 }' |
    LC_ALL=C sort)"

# Every export is the library's own, prefixed ybus_
if strays="$(grep -v '^ybus_' <<< "$exported")"; then
    echo "error: the library exports symbols without the ybus_ prefix:" >&2
    echo "$strays" >&2
    exit 1
fi

if [ "${1:-}" = --update ]; then
    {
        echo "# The ybus library's exported functions, one per line, sorted (LC_ALL=C)."
        echo "# Checked by ybus_ffi/check-symbols.sh; regenerate with --update, in the PR that"
        echo "#   changes the exports."
        echo "$exported"
    } > "$record"
    echo "Updated $record: $(wc -l <<< "$exported") exports"
    exit 0
elif [ -n "${1:-}" ]; then
    echo "usage: $0 [--update]" >&2
    exit 2
fi

# The record, as built
if ! difference="$(diff <(grep -v '^#' "$record") <(echo "$exported"))"; then
    echo "error: the library's exports differ from $record (< recorded, > built):" >&2
    echo "$difference" >&2
    echo "If the change is intended, run ybus_ffi/check-symbols.sh --update." >&2
    exit 1
fi

# The header, as built: every declared function exported, and every export declared
declared="$(grep -v '^[[:space:]]*//' "$crate/include/ybus.h" | grep -oE '\bybus_[a-z0-9_]+\(' |
    tr -d '(' | LC_ALL=C sort -u)"
if ! difference="$(diff <(echo "$declared") <(echo "$exported"))"; then
    echo "error: include/ybus.h and the library disagree (< declared, > exported):" >&2
    echo "$difference" >&2
    exit 1
fi

echo "$(wc -l <<< "$exported") exports match symbols.txt and include/ybus.h"
