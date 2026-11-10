# ybus_ffi

The C interface to [`powerflow`](../powerflow/): the `ybus` library, with the header
[`include/ybus.h`](include/). It's being built phase by phase from
[`docs/FfiChecklist.md`](../docs/FfiChecklist.md), and isn't released yet.

## Building and linking

```sh
cargo build -p ybus_ffi --release
```

This builds `target/release/libybus.so` and `target/release/libybus.a`. Compile against the
header and link with `-lybus`. A program linking the static library on Linux also needs the
system libraries Rust's standard library uses, which
`cargo rustc -p ybus_ffi --lib --crate-type staticlib -- --print native-static-libs` lists:

```sh
cc -I ybus_ffi/include program.c target/release/libybus.a \
    -lgcc_s -lutil -lrt -lpthread -lm -ldl -lc
```

## Tests

The Rust tests run with `cargo test -p ybus_ffi`. The C tests in [`tests/c/`](tests/c/) are
small programs, compiled against the header and linked with the static library as a C program
would be, then run, by `tests/c/run.sh`, which `check.sh` runs once.

## Cases

A `ybus_case` is an opaque handle to a case, made by one of three functions and freed by
`ybus_case_free`:

- **`ybus_case_from_path(path, &case)`** reads a MATPOWER file. The path is a NUL-terminated
  UTF-8 string.
- **`ybus_case_from_bytes(bytes, len, &case)`** parses the file's text from memory.
- **`ybus_case_clone(case, &copy)`** copies a case, so the copy can change without changing the
  original.

Each sets `*out` to `NULL` first, and leaves it `NULL` on an error. `ybus_case_free(NULL)` does
nothing. Rows are numbered from 0, as in C.

Every getter returns a status and writes its value through an out-parameter, such as
`ybus_case_bus_count(case, &buses)`. The sizes and scalars are:
- `ybus_case_bus_count`, `ybus_case_gen_count` and `ybus_case_branch_count`
- `ybus_case_base_mva`, and `ybus_case_set_base_mva`
- `ybus_case_has_costs` and `ybus_case_has_bus_names`

### Columns

Each column of the bus, generator and branch tables is read whole, one entry per row in case
order, by a field code. The code is the field's 1-based column number in MATPOWER's matrix for
that table. For example, `YBUS_BUS_PD_MW` is 3, `YBUS_GEN_STATUS` is 8 and `YBUS_BRANCH_X` is
4; the header lists them all.

There's one getter per table and value type:

```c
double pd[3];
ybus_case_get_bus_f64(case, YBUS_BUS_PD_MW, pd, 3);
int64_t status[3];
ybus_case_get_branch_i64(case, YBUS_BRANCH_STATUS, status, 3);
```

- **`int64_t` columns:** labels, areas, zones, bus types (MATPOWER's codes 1 to 4) and statuses
  (1 in service, 0 not). Every other field is a `double` column.
- **A code that's wrong for the function** gives `YBUS_ERROR_INVALID_FIELD`, whether the table
  doesn't have it or it's of the other type.
- **Too short a buffer** gives `YBUS_ERROR_BUFFER_TOO_SMALL` and is left untouched. A longer
  one keeps its tail.

Each getter has a setter, such as `ybus_case_set_bus_f64(case, YBUS_BUS_PD_MW, values, 3)`.
- **Length:** the array must have exactly one entry per row, or it's
  `YBUS_ERROR_LENGTH_MISMATCH`.
- **`int64_t` columns** are checked whole before anything is written. A negative label, area or
  zone, a bus type outside 1 to 4, or a status other than 0 or 1 is
  `YBUS_ERROR_INVALID_VALUE`, and the column is left as it was.
- **`double` columns** aren't checked. Data a solve can't use, such as a branch to a bus that
  doesn't exist, is reported when the case is solved.

`ybus_case_bus_name(case, row, buffer, len, &needed)` reads one bus name, as a NUL-terminated
string. Call it with `NULL, 0` to learn the size. A case without names gives
`YBUS_ERROR_MISSING_DATA`.

## ABI policy

- Every symbol is prefixed `ybus_`, and every macro `YBUS_`.
- `ybus_version()` returns the version of the library linked at run time, encoded as the
  header's `YBUS_VERSION` is (major × 1 000 000 + minor × 1 000 + patch). The header also
  defines `YBUS_VERSION_MAJOR`, `YBUS_VERSION_MINOR` and `YBUS_VERSION_PATCH`.
- Before 1.0, any minor version may break the ABI.
- From 1.0, minor versions only add functions and status codes.

## Errors

- **Status codes.** Every fallible function returns a `ybus_status` (an `int32_t`): `YBUS_OK`,
  which is 0, or the code of what failed. The header lists the codes. A code is never
  renumbered, and a failure without a code of its own, such as one added in a later version,
  is `YBUS_ERROR_OTHER`.
  - 1 to 6 are the boundary's own failures, such as a null argument or a buffer too short.
  - From 100, each `PowerflowError` variant has its own code, in declaration order, named
    after it: `YBUS_ERROR_NO_CASE_FILE` is 100, and `YBUS_ERROR_TOPOLOGY_MISMATCH` is 116.
- **The last error.** Each thread keeps the message of its last failed call, which
  `ybus_last_error_message(buffer, len)` copies out. It returns the bytes the message needs,
  counting the NUL, or 0 if no call on the thread has failed. It writes nothing unless the
  buffer is big enough, so call it with `NULL, 0` first to learn the size. The message is set
  only on failure, so it describes the most recent failure, not the most recent call.
- **Panics.** A panic in the library is a bug. It's caught at the boundary and returned as
  `YBUS_ERROR_PANIC`, rather than aborting the process. Afterwards, only free the handles that
  were passed to that call. Rust's default panic hook prints the panic to stderr first; the
  library doesn't replace the process-wide hook.

## Thread safety

- Any handle can move between threads.
- A handle that isn't being changed (a case, a model, a solution) can be read from many threads
  at once.
- A case being changed, or an options handle being set, needs the caller's own lock.

The Rust types behind the handles are checked to be `Send` and `Sync` by
`powerflow/tests/thread_safety_test.rs`.
