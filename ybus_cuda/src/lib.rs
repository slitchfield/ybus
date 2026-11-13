//! The GPU path for the power flows, on CUDA through [`cudarc`]: sparse triangular solves with
//! many right-hand sides against a factorization made on the CPU by [`ybus_sparse`].
//!
//! - [`Gpu`]: a CUDA device, its stream and a cuSPARSE handle, or `None` on a machine without
//!   CUDA.
//! - [`GpuFactor`]: a [`ybus_sparse`] Cholesky or LU factorization uploaded once, which solves
//!   blocks of right-hand sides with cuSPARSE's SpSM and reports the time spent on transfers.
//! - [`GpuFactor::pair_solves`] and [`GpuFactor::pair_updates`]: blocks whose right-hand sides
//!   are built on the device (`e_plus − e_minus`), optionally with a rank-one update of a base
//!   solution per column, kept on the device as a [`GpuBlock`] for [`GpuOperator`] products and
//!   per-column maxima, so only what's asked for is downloaded.
//!
//! With `ybus_ffi`, this is one of the two crates in the workspace with `unsafe` code: cudarc
//! exposes cuSPARSE only as raw FFI. Each `unsafe` block says why it's sound, and nothing
//! unsafe crosses the public API.
//! Results are not bit for bit the CPU's, since the GPU orders its arithmetic differently; they
//! agree to within rounding.

use std::fmt;

use cudarc::cusparse::result::CusparseError;
use cudarc::driver::DriverError;

mod factor;
mod kernels;
mod pairs;
mod sparse;

pub use factor::{GpuFactor, Timings};
pub use pairs::{ColumnMax, GpuBlock, GpuOperator, Pair, PairUpdate};

use std::sync::Arc;

use cudarc::driver::{CudaContext, CudaStream};

/// An error from the GPU path.
#[derive(Debug)]
pub enum CudaError {
    /// The CUDA driver failed.
    Driver(DriverError),
    /// cuSPARSE failed.
    Sparse(CusparseError),
    /// A factor has more entries than cuSPARSE's 32-bit indices can address.
    TooLarge {
        /// The number of stored entries.
        entries: usize,
    },
    /// NVRTC couldn't compile this crate's kernels.
    Compile(String),
    /// A right-hand side doesn't have one entry per row of the factorization.
    ShapeMismatch {
        /// The number of rows of the factorization.
        expected: usize,
        /// The length found.
        found: usize,
    },
}

impl fmt::Display for CudaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CudaError::Driver(err) => write!(f, "CUDA driver error: {err:?}"),
            CudaError::Sparse(err) => write!(f, "cuSPARSE error: {:?}", err.0),
            CudaError::TooLarge { entries } => write!(
                f,
                "a factor with {entries} entries is too large for 32-bit indices"
            ),
            CudaError::Compile(log) => write!(f, "NVRTC couldn't compile the kernels: {log}"),
            CudaError::ShapeMismatch { expected, found } => write!(
                f,
                "a right-hand side has {found} entries, but the factorization has {expected} rows"
            ),
        }
    }
}

impl std::error::Error for CudaError {}

impl From<DriverError> for CudaError {
    fn from(err: DriverError) -> Self {
        CudaError::Driver(err)
    }
}

impl From<CusparseError> for CudaError {
    fn from(err: CusparseError) -> Self {
        CudaError::Sparse(err)
    }
}

/// A CUDA device (the first), with the stream all of this crate's work runs on in order, and a
/// cuSPARSE handle bound to that stream.
pub struct Gpu {
    context: Arc<CudaContext>,
    stream: Arc<CudaStream>,
    sparse: sparse::Handle,
    kernels: kernels::Kernels,
}

impl Gpu {
    /// The first CUDA device, or `None` if this machine has no CUDA: the driver, cuSPARSE or
    /// NVRTC library can't be loaded, or the driver reports no device. Other failures are
    /// errors. It compiles this crate's kernels, which takes a fraction of a second.
    pub fn new() -> Result<Option<Gpu>, CudaError> {
        // SAFETY: these only try to load the libraries, which runs their initializers: the
        //   system's CUDA driver, cuSPARSE and NVRTC, which cudarc would load on first use
        //   anyway (and panic if it couldn't, which this check avoids)
        let present = unsafe {
            cudarc::driver::sys::is_culib_present()
                && cudarc::cusparse::sys::is_culib_present()
                && cudarc::nvrtc::sys::is_culib_present()
        };
        if !present {
            return Ok(None);
        }
        let context = match CudaContext::new(0) {
            Ok(context) => context,
            Err(err) if err.0 == cudarc::driver::sys::CUresult::CUDA_ERROR_NO_DEVICE => {
                return Ok(None)
            }
            Err(err) => return Err(err.into()),
        };
        let stream = context.default_stream();
        let sparse = sparse::Handle::new(&context, &stream)?;
        let kernels = kernels::Kernels::compile(&context)?;
        Ok(Some(Gpu {
            context,
            stream,
            sparse,
            kernels,
        }))
    }

    /// The device's name, as the driver reports it.
    pub fn name(&self) -> Result<String, CudaError> {
        Ok(self.context.name()?)
    }

    /// The device memory free and in total, in bytes.
    pub fn memory(&self) -> Result<(usize, usize), CudaError> {
        self.context.bind_to_thread()?;
        Ok(self.context.mem_get_info()?)
    }
}
