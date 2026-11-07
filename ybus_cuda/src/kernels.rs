//! The kernels in `kernels.cu`, compiled once per [`Gpu`](crate::Gpu) by NVRTC, each behind a
//! safe launch that checks its buffers' lengths, so no launch can read or write past one.

use std::ffi::{c_char, CString};
use std::sync::Arc;

use cudarc::driver::{
    CudaContext, CudaFunction, CudaSlice, CudaStream, LaunchConfig, PushKernelArg,
};
use cudarc::nvrtc::{result as nvrtc, sys as nvrtc_sys, Ptx};

use crate::CudaError;

const SOURCE: &str = include_str!("kernels.cu");

/// Threads per block for `column_max`, which its shared arrays are sized for
const MAX_THREADS: u32 = 256;

pub(crate) struct Kernels {
    scatter_pairs: CudaFunction,
    pair_update: CudaFunction,
    zero_entries: CudaFunction,
    broadcast: CudaFunction,
    column_max: CudaFunction,
}

/// An NVRTC program, destroyed on drop.
struct Program(nvrtc_sys::nvrtcProgram);

impl Drop for Program {
    fn drop(&mut self) {
        // SAFETY: the program was created by `create_program` and is destroyed only here
        let _ = unsafe { nvrtc::destroy_program(self.0) };
    }
}

fn compile_error(err: impl std::fmt::Debug) -> CudaError {
    CudaError::Compile(format!("{err:?}"))
}

/// `kernels.cu` compiled for the device's architecture, as a CUBIN rather than PTX: the driver
/// must be at least as new as NVRTC to compile PTX itself (this machine's CUDA 13.2 driver
/// rejects CUDA 13.3's PTX), while a CUBIN loads on any driver of the same major version.
fn cubin(major: i32, minor: i32) -> Result<Vec<u8>, CudaError> {
    let source = CString::new(SOURCE).expect("no NUL in kernels.cu");
    let program = Program(nvrtc::create_program(&source, None).map_err(compile_error)?);
    // Without fused multiply-adds, so the update rounds as the CPU's separate multiply and add
    //   do, and with IEEE division and denormals
    let options = [
        format!("--gpu-architecture=sm_{major}{minor}"),
        "--fmad=false".to_string(),
        "--prec-div=true".to_string(),
        "--ftz=false".to_string(),
    ];
    // SAFETY: `program` is live, and the options are NUL-free strings
    if let Err(err) = unsafe { nvrtc::compile_program(program.0, &options) } {
        // SAFETY: `program` is live; the log is a NUL-terminated string of the returned size
        let log = unsafe { nvrtc::get_program_log(program.0) }
            .map(|log| {
                let bytes: Vec<u8> = log.iter().map(|&c| c as u8).collect();
                String::from_utf8_lossy(&bytes)
                    .trim_end_matches('\0')
                    .to_string()
            })
            .unwrap_or_default();
        return Err(CudaError::Compile(format!("{err:?}: {log}")));
    }
    let mut size = 0usize;
    // SAFETY: `program` is live and compiled for a real architecture, so it has a CUBIN, and
    //   the buffer is the size NVRTC reports for it
    unsafe {
        nvrtc_sys::nvrtcGetCUBINSize(program.0, &mut size)
            .result()
            .map_err(compile_error)?;
        let mut cubin = vec![0u8; size];
        nvrtc_sys::nvrtcGetCUBIN(program.0, cubin.as_mut_ptr() as *mut c_char)
            .result()
            .map_err(compile_error)?;
        Ok(cubin)
    }
}

/// `value` as a kernel's `int`
fn int(value: usize) -> Result<i32, CudaError> {
    i32::try_from(value).map_err(|_| CudaError::TooLarge { entries: value })
}

impl Kernels {
    pub(crate) fn compile(context: &Arc<CudaContext>) -> Result<Self, CudaError> {
        let (major, minor) = context.compute_capability()?;
        let module = context.load_module(Ptx::from_binary(cubin(major, minor)?))?;
        Ok(Kernels {
            scatter_pairs: module.load_function("scatter_pairs")?,
            pair_update: module.load_function("pair_update")?,
            zero_entries: module.load_function("zero_entries")?,
            broadcast: module.load_function("broadcast")?,
            column_max: module.load_function("column_max")?,
        })
    }

    /// Sets column `k` of the `n`×`width` block (zeroed beforehand) to `e_plus[k] − e_minus[k]`,
    /// for rows below `n` or −1 for none.
    pub(crate) fn scatter_pairs(
        &self,
        stream: &Arc<CudaStream>,
        block: &mut CudaSlice<f64>,
        n: usize,
        width: usize,
        plus: &CudaSlice<i32>,
        minus: &CudaSlice<i32>,
    ) -> Result<(), CudaError> {
        assert!(block.len() >= n * width && plus.len() >= width && minus.len() >= width);
        let (n, w) = (int(n)?, int(width)?);
        let mut launch = stream.launch_builder(&self.scatter_pairs);
        launch.arg(block).arg(&n).arg(&w).arg(plus).arg(minus);
        // SAFETY: the kernel's arguments match `scatter_pairs` in kernels.cu; the block holds
        //   `n * width` values and the row lists `width` each (checked above), and the callers
        //   give rows below `n` or −1, so every write is inside the block
        unsafe { launch.launch(LaunchConfig::for_num_elems(w as u32))? };
        Ok(())
    }

    /// Each column's rank-one update: see `pair_update` in kernels.cu. `plus` and `minus` are
    /// rows below `n` or −1.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn pair_update(
        &self,
        stream: &Arc<CudaStream>,
        y: &CudaSlice<f64>,
        out: &mut CudaSlice<f64>,
        n: usize,
        width: usize,
        base: &CudaSlice<f64>,
        plus: &CudaSlice<i32>,
        minus: &CudaSlice<i32>,
        alpha: &CudaSlice<f64>,
        inverse_weight: &CudaSlice<f64>,
        singular: &mut CudaSlice<i32>,
    ) -> Result<(), CudaError> {
        assert!(y.len() >= n * width && out.len() >= n * width && base.len() >= n);
        assert!(plus.len() >= width && minus.len() >= width && alpha.len() >= width);
        assert!(inverse_weight.len() >= width && singular.len() >= width);
        let (n_i, w) = (int(n)?, int(width)?);
        let total =
            u32::try_from(n * width).map_err(|_| CudaError::TooLarge { entries: n * width })?;
        let mut launch = stream.launch_builder(&self.pair_update);
        launch
            .arg(y)
            .arg(out)
            .arg(&n_i)
            .arg(&w)
            .arg(base)
            .arg(plus)
            .arg(minus)
            .arg(alpha)
            .arg(inverse_weight)
            .arg(singular);
        // SAFETY: the arguments match `pair_update` in kernels.cu; every buffer holds what the
        //   kernel indexes (checked above), the rows read are below `n` or −1, and `y` and `out`
        //   are distinct buffers (one borrowed shared, one exclusive)
        unsafe { launch.launch(LaunchConfig::for_num_elems(total))? };
        Ok(())
    }

    /// Sets column `k`'s entry at `rows[k]` (below `n`, or −1 for none) to zero.
    pub(crate) fn zero_entries(
        &self,
        stream: &Arc<CudaStream>,
        block: &mut CudaSlice<f64>,
        n: usize,
        width: usize,
        rows: &CudaSlice<i32>,
    ) -> Result<(), CudaError> {
        assert!(block.len() >= n * width && rows.len() >= width);
        let (n, w) = (int(n)?, int(width)?);
        let mut launch = stream.launch_builder(&self.zero_entries);
        launch.arg(block).arg(&n).arg(&w).arg(rows);
        // SAFETY: the arguments match `zero_entries` in kernels.cu, the buffers hold what it
        //   indexes (checked above), and the rows are below `n` or −1
        unsafe { launch.launch(LaunchConfig::for_num_elems(w as u32))? };
        Ok(())
    }

    /// Sets every column of the `n`×`width` block to `offset`.
    pub(crate) fn broadcast(
        &self,
        stream: &Arc<CudaStream>,
        block: &mut CudaSlice<f64>,
        n: usize,
        width: usize,
        offset: &CudaSlice<f64>,
    ) -> Result<(), CudaError> {
        assert!(block.len() >= n * width && offset.len() >= n);
        let (n_i, w) = (int(n)?, int(width)?);
        let total =
            u32::try_from(n * width).map_err(|_| CudaError::TooLarge { entries: n * width })?;
        let mut launch = stream.launch_builder(&self.broadcast);
        launch.arg(block).arg(&n_i).arg(&w).arg(offset);
        // SAFETY: the arguments match `broadcast` in kernels.cu, and the buffers hold what it
        //   indexes (checked above)
        unsafe { launch.launch(LaunchConfig::for_num_elems(total))? };
        Ok(())
    }

    /// Each column's largest `|value| · weight[row]`, its row (−1 if all are zero) and whether
    /// any value isn't finite, into `best`, `best_row` and `not_finite`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn column_max(
        &self,
        stream: &Arc<CudaStream>,
        block: &CudaSlice<f64>,
        n: usize,
        width: usize,
        weight: &CudaSlice<f64>,
        best: &mut CudaSlice<f64>,
        best_row: &mut CudaSlice<i32>,
        not_finite: &mut CudaSlice<i32>,
    ) -> Result<(), CudaError> {
        assert!(block.len() >= n * width && weight.len() >= n);
        assert!(best.len() >= width && best_row.len() >= width && not_finite.len() >= width);
        let (n_i, w) = (int(n)?, int(width)?);
        let config = LaunchConfig {
            grid_dim: (w as u32, 1, 1),
            block_dim: (MAX_THREADS, 1, 1),
            shared_mem_bytes: 0,
        };
        let mut launch = stream.launch_builder(&self.column_max);
        launch
            .arg(block)
            .arg(&n_i)
            .arg(&w)
            .arg(weight)
            .arg(best)
            .arg(best_row)
            .arg(not_finite);
        // SAFETY: the arguments match `column_max` in kernels.cu, the buffers hold what it
        //   indexes (checked above), and the block size is the 256 its shared arrays hold
        unsafe { launch.launch(config)? };
        Ok(())
    }
}
