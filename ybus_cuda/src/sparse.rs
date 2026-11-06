//! Safe owners for cuSPARSE's objects: the library handle, a triangular matrix on the device,
//! and a block of right-hand sides with the triangular solves analysed for it. cudarc exposes
//! cuSPARSE only as raw FFI, so each object here holds the raw descriptor together with the
//! device memory it points into, and destroys the descriptor on drop.

use std::ffi::c_void;
use std::sync::Arc;

use cudarc::cusparse::result;
use cudarc::cusparse::sys;
use cudarc::driver::{CudaContext, CudaSlice, CudaStream, DevicePtr, DevicePtrMut};

use crate::CudaError;

/// Scalars are passed to cuSPARSE by host pointer (its default pointer mode)
const ONE: f64 = 1.0;

/// The cuSPARSE library handle, with its work ordered on one stream.
pub(crate) struct Handle {
    raw: sys::cusparseHandle_t,
}

impl Handle {
    pub(crate) fn new(context: &Arc<CudaContext>, stream: &CudaStream) -> Result<Self, CudaError> {
        context.bind_to_thread()?;
        let raw = result::create()?;
        let handle = Handle { raw };
        // SAFETY: `raw` is a live handle, and the stream outlives it (both belong to the same
        //   `Gpu`, which drops the handle first); the driver's and the runtime's stream types are
        //   the same object behind different names
        unsafe {
            sys::cusparseSetStream(handle.raw, stream.cu_stream() as sys::cudaStream_t).result()?;
        }
        Ok(handle)
    }
}

impl Drop for Handle {
    fn drop(&mut self) {
        // SAFETY: `raw` was created by `new` and is destroyed only here
        let _ = unsafe { result::destroy(self.raw) };
    }
}

/// Which triangle a matrix holds.
#[derive(Clone, Copy)]
pub(crate) enum Fill {
    Lower,
    Upper,
}

/// A square triangular matrix in CSR form on the device, with its diagonal stored (a unit
/// diagonal stored as ones). It owns the device arrays its descriptor points into, so they
/// outlive the descriptor.
pub(crate) struct Triangle {
    descr: sys::cusparseSpMatDescr_t,
    _pointers: CudaSlice<i32>,
    _indices: CudaSlice<i32>,
    _values: CudaSlice<f64>,
}

/// `values` as cuSPARSE's 32-bit indices, or `TooLarge` if one doesn't fit
fn to_i32(values: &[usize], entries: usize) -> Result<Vec<i32>, CudaError> {
    values
        .iter()
        .map(|&v| i32::try_from(v).map_err(|_| CudaError::TooLarge { entries }))
        .collect()
}

impl Triangle {
    /// Uploads the `n`×`n` CSR matrix with these arrays: row `i`'s entries are at
    /// `pointers[i]..pointers[i + 1]` of `indices` (columns, ascending within each row) and
    /// `values`. The caller ensures the arrays describe a valid triangle of `fill`'s kind.
    pub(crate) fn upload(
        context: &Arc<CudaContext>,
        stream: &Arc<CudaStream>,
        n: usize,
        pointers: &[usize],
        indices: &[usize],
        values: &[f64],
        fill: Fill,
    ) -> Result<Self, CudaError> {
        assert_eq!(
            pointers.len(),
            n + 1,
            "one pointer per row, and one past the end"
        );
        assert_eq!(indices.len(), values.len(), "one index per value");
        let entries = values.len();
        let pointers = stream.clone_htod(&to_i32(pointers, entries)?)?;
        let indices = stream.clone_htod(&to_i32(indices, entries)?)?;
        let values = stream.clone_htod(values)?;
        context.bind_to_thread()?;
        let mut descr = std::ptr::null_mut();
        {
            let (p, _p) = pointers.device_ptr(stream);
            let (i, _i) = indices.device_ptr(stream);
            let (v, _v) = values.device_ptr(stream);
            // SAFETY: the three device arrays hold `n + 1` row offsets and `entries` column
            //   indices and values, as declared, and are owned by the returned `Triangle`, which
            //   destroys the descriptor before dropping them; cuSPARSE only reads them
            unsafe {
                sys::cusparseCreateCsr(
                    &mut descr,
                    n as i64,
                    n as i64,
                    entries as i64,
                    p as *mut c_void,
                    i as *mut c_void,
                    v as *mut c_void,
                    sys::cusparseIndexType_t::CUSPARSE_INDEX_32I,
                    sys::cusparseIndexType_t::CUSPARSE_INDEX_32I,
                    sys::cusparseIndexBase_t::CUSPARSE_INDEX_BASE_ZERO,
                    sys::cudaDataType::CUDA_R_64F,
                )
                .result()?;
            }
        }
        let triangle = Triangle {
            descr,
            _pointers: pointers,
            _indices: indices,
            _values: values,
        };
        let mut fill = match fill {
            Fill::Lower => sys::cusparseFillMode_t::CUSPARSE_FILL_MODE_LOWER,
            Fill::Upper => sys::cusparseFillMode_t::CUSPARSE_FILL_MODE_UPPER,
        };
        let mut diagonal = sys::cusparseDiagType_t::CUSPARSE_DIAG_TYPE_NON_UNIT;
        // SAFETY: `descr` is live, and each attribute is passed as a pointer to a value of the
        //   type cuSPARSE expects for it, with that type's size
        unsafe {
            sys::cusparseSpMatSetAttribute(
                triangle.descr,
                sys::cusparseSpMatAttribute_t::CUSPARSE_SPMAT_FILL_MODE,
                &mut fill as *mut _ as *mut c_void,
                std::mem::size_of_val(&fill),
            )
            .result()?;
            sys::cusparseSpMatSetAttribute(
                triangle.descr,
                sys::cusparseSpMatAttribute_t::CUSPARSE_SPMAT_DIAG_TYPE,
                &mut diagonal as *mut _ as *mut c_void,
                std::mem::size_of_val(&diagonal),
            )
            .result()?;
        }
        Ok(triangle)
    }
}

impl Drop for Triangle {
    fn drop(&mut self) {
        // SAFETY: `descr` was created in `upload` and is destroyed only here, before the
        //   arrays it points into (fields drop after `drop` runs)
        let _ = unsafe { sys::cusparseDestroySpMat(self.descr) };
    }
}

/// A dense column-major `rows`×`columns` matrix on the device, with its descriptor.
struct Dense {
    descr: sys::cusparseDnMatDescr_t,
    data: CudaSlice<f64>,
}

impl Dense {
    fn new(
        context: &Arc<CudaContext>,
        stream: &Arc<CudaStream>,
        rows: usize,
        columns: usize,
    ) -> Result<Self, CudaError> {
        let mut data = stream.alloc_zeros::<f64>(rows * columns)?;
        context.bind_to_thread()?;
        let mut descr = std::ptr::null_mut();
        {
            let (d, _d) = data.device_ptr_mut(stream);
            // SAFETY: `data` holds `rows * columns` values, column-major with leading dimension
            //   `rows`, and is owned by the returned `Dense`, which destroys the descriptor
            //   first
            unsafe {
                sys::cusparseCreateDnMat(
                    &mut descr,
                    rows as i64,
                    columns as i64,
                    rows as i64,
                    d as *mut c_void,
                    sys::cudaDataType::CUDA_R_64F,
                    sys::cusparseOrder_t::CUSPARSE_ORDER_COL,
                )
                .result()?;
            }
        }
        Ok(Dense { descr, data })
    }
}

impl Drop for Dense {
    fn drop(&mut self) {
        // SAFETY: `descr` was created in `new` and is destroyed only here, before `data`
        let _ = unsafe { sys::cusparseDestroyDnMat(self.descr) };
    }
}

/// One triangular solve analysed for a matrix and a pair of dense blocks, with the work buffer
/// the analysis refers to.
struct Spsm {
    descr: sys::cusparseSpSMDescr_t,
    _buffer: CudaSlice<u8>,
}

impl Drop for Spsm {
    fn drop(&mut self) {
        // SAFETY: `descr` was created in `Spsm::analyse` and is destroyed only here, before
        //   the buffer it refers to
        let _ = unsafe { sys::cusparseSpSM_destroyDescr(self.descr) };
    }
}

impl Spsm {
    /// Analyses `matrix x = from` into `to`.
    fn analyse(
        handle: &Handle,
        stream: &Arc<CudaStream>,
        matrix: &Triangle,
        from: &Dense,
        to: &Dense,
    ) -> Result<Self, CudaError> {
        let mut descr = std::ptr::null_mut();
        // SAFETY: `descr` is an out-parameter for a new descriptor
        unsafe { sys::cusparseSpSM_createDescr(&mut descr).result()? };
        let mut size = 0usize;
        // SAFETY: every descriptor is live, and `ONE` is a host `f64`, as the pointer mode and
        //   compute type say
        unsafe {
            sys::cusparseSpSM_bufferSize(
                handle.raw,
                sys::cusparseOperation_t::CUSPARSE_OPERATION_NON_TRANSPOSE,
                sys::cusparseOperation_t::CUSPARSE_OPERATION_NON_TRANSPOSE,
                &ONE as *const f64 as *const c_void,
                matrix.descr,
                from.descr,
                to.descr,
                sys::cudaDataType::CUDA_R_64F,
                sys::cusparseSpSMAlg_t::CUSPARSE_SPSM_ALG_DEFAULT,
                descr,
                &mut size,
            )
            .result()
            .inspect_err(|_| {
                // SAFETY: created above and not yet owned by an `Spsm`
                let _ = sys::cusparseSpSM_destroyDescr(descr);
            })?;
        }
        let buffer = match stream.alloc_zeros::<u8>(size.max(1)) {
            Ok(buffer) => buffer,
            Err(err) => {
                // SAFETY: created above and not yet owned by an `Spsm`
                let _ = unsafe { sys::cusparseSpSM_destroyDescr(descr) };
                return Err(err.into());
            }
        };
        let spsm = Spsm {
            descr,
            _buffer: buffer,
        };
        {
            let (b, _b) = spsm._buffer.device_ptr(stream);
            // SAFETY: as for the buffer size, with a work buffer of the size cuSPARSE asked for,
            //   owned by `spsm` and so alive for every later solve with this descriptor
            unsafe {
                sys::cusparseSpSM_analysis(
                    handle.raw,
                    sys::cusparseOperation_t::CUSPARSE_OPERATION_NON_TRANSPOSE,
                    sys::cusparseOperation_t::CUSPARSE_OPERATION_NON_TRANSPOSE,
                    &ONE as *const f64 as *const c_void,
                    matrix.descr,
                    from.descr,
                    to.descr,
                    sys::cudaDataType::CUDA_R_64F,
                    sys::cusparseSpSMAlg_t::CUSPARSE_SPSM_ALG_DEFAULT,
                    spsm.descr,
                    b as *mut c_void,
                )
                .result()?;
            }
        }
        Ok(spsm)
    }
}

/// A block of right-hand sides on the device, `n` rows by `width` columns, with the two
/// triangular solves analysed for it: the lower one from `first` into `second`, and the upper
/// one from `second` back into `first`, so the solution ends where the right-hand sides
/// started. It must only be used with the two triangles it was analysed for.
pub(crate) struct Block {
    // Declared first, so they're destroyed before the blocks they were analysed with
    lower: Spsm,
    upper: Spsm,
    first: Dense,
    second: Dense,
    width: usize,
}

impl Block {
    pub(crate) fn new(
        context: &Arc<CudaContext>,
        stream: &Arc<CudaStream>,
        handle: &Handle,
        lower: &Triangle,
        upper: &Triangle,
        n: usize,
        width: usize,
    ) -> Result<Self, CudaError> {
        let first = Dense::new(context, stream, n, width)?;
        let second = Dense::new(context, stream, n, width)?;
        let lower = Spsm::analyse(handle, stream, lower, &first, &second)?;
        let upper = Spsm::analyse(handle, stream, upper, &second, &first)?;
        Ok(Block {
            lower,
            upper,
            first,
            second,
            width,
        })
    }

    pub(crate) fn width(&self) -> usize {
        self.width
    }

    /// Copies `host` (column-major, `n * width` values) into the block.
    pub(crate) fn upload(
        &mut self,
        stream: &Arc<CudaStream>,
        host: &[f64],
    ) -> Result<(), CudaError> {
        stream.memcpy_htod(host, &mut self.first.data)?;
        Ok(())
    }

    /// Copies the block (column-major, `n * width` values) into `host`.
    pub(crate) fn download(
        &self,
        stream: &Arc<CudaStream>,
        host: &mut [f64],
    ) -> Result<(), CudaError> {
        stream.memcpy_dtoh(&self.first.data, host)?;
        Ok(())
    }

    /// Solves `lower upper x = b` for every column `b` of the block, in place. `lower` and
    /// `upper` must be the triangles the block was analysed with.
    pub(crate) fn solve(
        &mut self,
        handle: &Handle,
        stream: &Arc<CudaStream>,
        lower: &Triangle,
        upper: &Triangle,
    ) -> Result<(), CudaError> {
        let solves = [
            (
                lower.descr,
                self.lower.descr,
                self.first.descr,
                self.second.descr,
            ),
            (
                upper.descr,
                self.upper.descr,
                self.second.descr,
                self.first.descr,
            ),
        ];
        // Taken for their ordering on the stream: the solves read and write both blocks
        let (_f, _first) = self.first.data.device_ptr_mut(stream);
        let (_s, _second) = self.second.data.device_ptr_mut(stream);
        for (matrix, spsm, from, to) in solves {
            // SAFETY: the descriptors are live and were analysed together (the caller passes
            //   the triangles the block was built with), the analysis's work buffer is owned by
            //   `spsm`, and `ONE` is a host `f64`
            unsafe {
                sys::cusparseSpSM_solve(
                    handle.raw,
                    sys::cusparseOperation_t::CUSPARSE_OPERATION_NON_TRANSPOSE,
                    sys::cusparseOperation_t::CUSPARSE_OPERATION_NON_TRANSPOSE,
                    &ONE as *const f64 as *const c_void,
                    matrix,
                    from,
                    to,
                    sys::cudaDataType::CUDA_R_64F,
                    sys::cusparseSpSMAlg_t::CUSPARSE_SPSM_ALG_DEFAULT,
                    spsm,
                )
                .result()?;
            }
        }
        Ok(())
    }
}
