#!/usr/bin/env bash
# The C tests: builds the ybus library, then compiles each tests/c/test_*.c against
#   include/ybus.h, links it with the static library as a C program would be linked, and runs
#   it. check.sh runs this once. It needs a C compiler: cc, or the one $CC names.
set -euo pipefail

crate="$(cd "$(dirname "$0")/../.." && pwd)"
workspace="$(dirname "$crate")"
target="${CARGO_TARGET_DIR:-$workspace/target}"
cc="${CC:-cc}"

cargo build --quiet -p ybus_ffi
library="$target/debug/libybus.a"

# The system libraries a Rust static library needs on Linux, as
#   cargo rustc -p ybus_ffi --lib --crate-type staticlib -- --print native-static-libs
#   lists them
system_libraries=(-lgcc_s -lutil -lrt -lpthread -lm -ldl -lc)

# The test case files, such as three_bus.m, which cases.h loads
export YBUS_TEST_DATA="$crate/tests/c"

build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

count=0
for test in "$crate"/tests/c/test_*.c; do
    name="$(basename "$test" .c)"
    "$cc" -std=c11 -Wall -Wextra -Wpedantic -Werror -I "$crate/include" -I "$crate/tests/c" \
        "$test" "$library" "${system_libraries[@]}" -o "$build/$name"
    "$build/$name"
    echo "C test $name passed"
    count=$((count + 1))
done
echo "$count C tests passed"
