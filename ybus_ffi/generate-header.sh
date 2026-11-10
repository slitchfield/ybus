#!/usr/bin/env bash
# Generates ybus_ffi/include/ybus.h from the ybus_ffi crate with cbindgen (the command-line tool,
#   pinned to the version CI installs: cargo install cbindgen --version 0.29.4 --locked). With
#   --verify, it writes nothing and fails if the committed header differs from what it would
#   generate, as check.sh and CI run it.
set -euo pipefail

crate="$(cd "$(dirname "$0")" && pwd)"
cbindgen --config "$crate/cbindgen.toml" --crate ybus_ffi --output "$crate/include/ybus.h" \
    "$@" "$crate"
