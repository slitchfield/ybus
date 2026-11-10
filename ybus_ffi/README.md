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
- **The last error.** Each thread keeps the message of its last failed call. It's set only on
  failure, so it describes the most recent failure, not the most recent call.
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
