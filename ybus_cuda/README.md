# ybus_cuda

**Experimental.** The GPU path of
[`ybus_powerflow`](https://crates.io/crates/ybus_powerflow): sparse triangular solves with
many right-hand sides on CUDA, against a Cholesky or LU factorization made on the CPU by
[`ybus_sparse`](https://crates.io/crates/ybus_sparse).

- **`Gpu`:** a CUDA device, its stream and a cuSPARSE handle, or `None` on a machine without
  CUDA.
- **`GpuFactor`:** a factorization uploaded once, which solves blocks of right-hand sides with
  cuSPARSE's SpSM and reports the time spent on transfers.
- **`GpuFactor::pair_solves` and `pair_updates`:** blocks whose right-hand sides are built on
  the device, optionally with a rank-one update of a base solution per column, as a DC study's
  PTDF columns and outages need. They stay on the device as a `GpuBlock`, so only what's asked
  for is downloaded.

Most users want `ybus_powerflow`'s `cuda` feature, which drives this crate, rather than the
crate itself.

## Status

Experimental, at 0.1: the API will change, and each breaking change is a minor version bump.
The GPU's results agree with the CPU's to within rounding, not bit for bit, since the GPU
orders its arithmetic differently.

## Requirements

- **To build:** nothing beyond Rust. CUDA is loaded at run time through
  [`cudarc`](https://crates.io/crates/cudarc), so the crate builds without the CUDA toolkit,
  docs.rs included.
- **To run:** a CUDA 13 driver, cuSPARSE and NVRTC, which compiles the crate's kernels when a
  `Gpu` is made. Without them, `Gpu::new` returns `Ok(None)`.

## Quick start

```toml
[dependencies]
ybus_cuda = "0.1"
ybus_sparse = "0.1"
```

```rust,no_run
use ybus_cuda::{CudaError, Gpu, GpuFactor};
use ybus_sparse::{CsrMatrix, Factorization, SparseMatrix};

fn main() -> Result<(), CudaError> {
    let Some(gpu) = Gpu::new()? else {
        println!("no CUDA on this machine");
        return Ok(());
    };
    let a: CsrMatrix = SparseMatrix::from_triplets(
        3,
        3,
        vec![
            (0, 0, 4.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 4.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 4.0),
        ],
    )
    .unwrap()
    .to_csr();
    let cholesky = a.cholesky().unwrap();

    // Upload the factor once, then solve many right-hand sides in blocks
    let mut factor = GpuFactor::from_cholesky(&gpu, &cholesky)?;
    let rhs: Vec<Vec<f64>> = (0..100).map(|k| vec![1.0, k as f64, 2.0]).collect();
    let on_gpu = factor.solve_many(&rhs)?;
    let on_cpu = cholesky.solve_many(&rhs).unwrap();
    for (gpu, cpu) in on_gpu.iter().zip(&on_cpu) {
        assert!(gpu.iter().zip(cpu).all(|(g, c)| (g - c).abs() < 1e-12));
    }
    println!("{:?}", factor.timings());
    Ok(())
}
```

## Safety

With `ybus_ffi`, this is one of the two crates in the `ybus` workspace with `unsafe` code:
cudarc exposes cuSPARSE only as raw FFI. Each `unsafe` block says why it's sound, and nothing
unsafe crosses the public API.

## Development

- **MSRV:** Rust 1.88, set for the whole workspace in the root `Cargo.toml`.
- **Tests:** in the repository's `ybus_cuda/tests/`. They compare the GPU's solves with the
  CPU's, and skip on a machine without CUDA unless `YBUS_REQUIRE_CUDA=1` is set.
- **API documentation:** `cargo doc -p ybus_cuda --open`.

See the [workspace README](https://github.com/slitchfield/ybus#readme) for `./check.sh`, CI and
the rest of the docs.

## License

Licensed under either of the Apache License, Version 2.0
([`LICENSE-APACHE`](https://github.com/slitchfield/ybus/blob/main/LICENSE-APACHE)) or the MIT
license ([`LICENSE-MIT`](https://github.com/slitchfield/ybus/blob/main/LICENSE-MIT)), at your
option. See
[`AI-DISCLAIMER.md`](https://github.com/slitchfield/ybus/blob/main/AI-DISCLAIMER.md) for how
this licensing relates to the use of AI in developing the crate.
