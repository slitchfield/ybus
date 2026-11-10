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
system libraries Rust's standard library uses:

```sh
cc -I ybus_ffi/include program.c target/release/libybus.a -lpthread -ldl -lm
```

## ABI policy

- Every symbol is prefixed `ybus_`, and every macro `YBUS_`.
- `ybus_version()` returns the version of the library linked at run time, encoded as the
  header's `YBUS_VERSION` is (major × 1 000 000 + minor × 1 000 + patch). The header also
  defines `YBUS_VERSION_MAJOR`, `YBUS_VERSION_MINOR` and `YBUS_VERSION_PATCH`.
- Before 1.0, any minor version may break the ABI.
- From 1.0, minor versions only add functions and status codes.

## Thread safety

- Any handle can move between threads.
- A handle that isn't being changed (a case, a model, a solution) can be read from many threads
  at once.
- A case being changed, or an options handle being set, needs the caller's own lock.

The Rust types behind the handles are checked to be `Send` and `Sync` by
`powerflow/tests/thread_safety_test.rs`.
