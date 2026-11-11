#!/usr/bin/env bash
# The C tests: builds the ybus library, then compiles each tests/c/test_*.c against
#   include/ybus.h, links it with the static library as a C program would be linked, and runs
#   it, then builds the C example (examples/c/solve_case.c) and runs it on case14. check.sh
#   runs this plainly and with --sanitize. It needs a C compiler: cc, or the one $CC
#   names.
#
# Usage: ybus_ffi/tests/c/run.sh [--sanitize]
#   --sanitize  builds and runs every test under AddressSanitizer (with LeakSanitizer) and
#               UndefinedBehaviorSanitizer, failing on the first report. The Rust library isn't
#               instrumented, since stable Rust can't build with -Zsanitizer, but its
#               allocations go through the malloc ASan replaces: a handle leaked or freed twice
#               fails, as does a C buffer overrun or undefined behaviour in C. A write by Rust
#               past a C buffer doesn't show; the exports' length checks guard against that.
set -euo pipefail

mode=plain
flags=()
case "${1:-}" in
    "") ;;
    --sanitize)
        mode=sanitized
        flags=(-fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer -g)
        export ASAN_OPTIONS="detect_leaks=1:abort_on_error=1"
        export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1"
        ;;
    *)
        echo "usage: $0 [--sanitize]" >&2
        exit 2
        ;;
esac

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
# MATPOWER's case files (fetched by scripts/fetch-matpower.sh, not committed), from the
#   MATPOWER_DIR installation if set, as powerflow's tests take them, and the committed
#   references for case14, case30 and case118. A test that needs a case that isn't there skips
#   it, unless YBUS_REQUIRE_CASES=1, as CI sets
if [ -n "${MATPOWER_DIR:-}" ]; then
    matpower="$MATPOWER_DIR"
else
    version="$(sed -n 's/^VERSION="\(.*\)"$/\1/p' "$workspace/scripts/fetch-matpower.sh")"
    matpower="$workspace/testdata/matpower/matpower$version"
fi
export YBUS_MATPOWER_DATA="$matpower/data"
export YBUS_REFERENCE_DATA="$workspace/powerflow/tests/data/matpower-reference"

build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

count=0
for test in "$crate"/tests/c/test_*.c; do
    name="$(basename "$test" .c)"
    "$cc" -std=c11 -Wall -Wextra -Wpedantic -Werror "${flags[@]}" -I "$crate/include" \
        -I "$crate/tests/c" "$test" "$library" "${system_libraries[@]}" -o "$build/$name"
    "$build/$name"
    echo "C test $name passed ($mode)"
    count=$((count + 1))
done
echo "$count C tests passed ($mode)"

# The C example, built as the tests are and run on MATPOWER's case14, when it's fetched (or
#   always, with YBUS_REQUIRE_CASES=1)
"$cc" -std=c11 -Wall -Wextra -Wpedantic -Werror "${flags[@]}" -I "$crate/include" \
    "$crate/examples/c/solve_case.c" "$library" "${system_libraries[@]}" -o "$build/solve_case"
case14="$YBUS_MATPOWER_DATA/case14.m"
if [ -f "$case14" ]; then
    "$build/solve_case" "$case14" > "$build/solve_case.out"
    echo "C example solve_case ran on case14 ($mode)"
elif [ "${YBUS_REQUIRE_CASES:-}" = 1 ]; then
    echo "MATPOWER's case14.m is missing: run scripts/fetch-matpower.sh" >&2
    exit 1
else
    echo "C example solve_case built; skipped running it: MATPOWER's case14.m isn't fetched"
fi
