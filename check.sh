#!/usr/bin/env bash
# This scripts runs various CI-like checks in a convenient way.
set -eux

cargo check --workspace --all-targets
cargo check --workspace --all-features --lib 
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features --  -D warnings -W clippy::all
cargo test --workspace --all-targets
# The other test passes leave out the benchmarks (--all-targets without --benches). Their
#   test runs, one iteration each, take seconds in debug builds and don't depend on the
#   features, and the clippy runs still build them with every feature set.
targets=(--lib --bins --tests --examples)
cargo test --workspace "${targets[@]}" --all-features
# Each optional feature on its own, to catch code that only builds or lints cleanly when
#   another feature is also enabled. Keep this list in sync with [features] in
#   sparse_mat/Cargo.toml. The workspace root has no features of its own, so each one is
#   named with its crate.
for feature in sparse_mat/serde sparse_mat/rayon; do
    cargo clippy --workspace --all-targets --features "$feature" -- -D warnings -W clippy::all
    cargo test --workspace "${targets[@]}" --features "$feature"
done
cargo test --workspace --doc