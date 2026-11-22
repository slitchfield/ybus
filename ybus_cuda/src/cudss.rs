//! A symmetric indefinite `L D Lᵀ` on the device through NVIDIA's cuDSS, for KKT systems that
//! are refactored with new values on a fixed pattern.
//!
//! cuDSS is loaded at run time (`libcudss.so.1`), as cudarc loads the driver and cuSPARSE, so
//! the crate builds without it and a machine without it finds out when [`GpuLdlt::new`]
//! returns `Ok(None)`. The bindings are hand-written against cuDSS 1.0's header
//! (`cudss.h`, `cudss_data_types.h`), for the few functions used here.

use std::ffi::{c_int, c_void};
use std::sync::OnceLock;
use std::time::{Duration, Instant};

use cudarc::cusparse::sys::cudaDataType;
use cudarc::driver::{CudaSlice, DevicePtr, DevicePtrMut, PinnedHostSlice};
use ybus_sparse::compressed_matrix::Inertia;
use ybus_sparse::CsrMatrix;

use crate::{CudaError, Gpu};

type Status = c_int;
type Opaque = *mut c_void;

/// `cudssConfigParam_t`'s deterministic mode, by its position in the enum
const CONFIG_DETERMINISTIC_MODE: c_int = 24;
/// `cudssDataParam_t`'s error indicator and inertia
const DATA_INFO: c_int = 0;
const DATA_NPIVOTS: c_int = 2;
const DATA_INERTIA: c_int = 3;
/// `cudssPhase_t`
const PHASE_ANALYSIS: c_int = 0b11;
const PHASE_FACTORIZATION: c_int = 1 << 2;
const PHASE_REFACTORIZATION: c_int = 1 << 3;
const PHASE_SOLVE: c_int = 0b11_1111_0000;
/// `cudssMatrixType_t`, `cudssMatrixViewType_t`, `cudssIndexBase_t`, `cudssLayout_t`
const MTYPE_SYMMETRIC: c_int = 1;
const MVIEW_LOWER: c_int = 1;
const BASE_ZERO: c_int = 0;
const LAYOUT_COL_MAJOR: c_int = 0;

/// The cuDSS functions used, from the loaded library, which they keep alive
struct Api {
    create: unsafe extern "C" fn(*mut Opaque) -> Status,
    destroy: unsafe extern "C" fn(Opaque) -> Status,
    set_stream: unsafe extern "C" fn(Opaque, *mut c_void) -> Status,
    config_create: unsafe extern "C" fn(*mut Opaque) -> Status,
    config_destroy: unsafe extern "C" fn(Opaque) -> Status,
    config_set: unsafe extern "C" fn(Opaque, c_int, *const c_void, usize) -> Status,
    data_create: unsafe extern "C" fn(Opaque, *mut Opaque) -> Status,
    data_destroy: unsafe extern "C" fn(Opaque, Opaque) -> Status,
    data_get: unsafe extern "C" fn(Opaque, Opaque, c_int, *mut c_void, usize, *mut usize) -> Status,
    execute: unsafe extern "C" fn(Opaque, c_int, Opaque, Opaque, Opaque, Opaque, Opaque) -> Status,
    #[allow(clippy::type_complexity)]
    matrix_create_csr: unsafe extern "C" fn(
        *mut Opaque,
        i64,
        i64,
        i64,
        *const c_void,
        *const c_void,
        *const c_void,
        *const c_void,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
        c_int,
    ) -> Status,
    matrix_create_dn:
        unsafe extern "C" fn(*mut Opaque, i64, i64, i64, *const c_void, c_int, c_int) -> Status,
    matrix_destroy: unsafe extern "C" fn(Opaque) -> Status,
    _library: libloading::Library,
}

impl Api {
    /// cuDSS, loaded once per process, or `None` if it can't be
    fn get() -> Option<&'static Api> {
        static API: OnceLock<Option<Api>> = OnceLock::new();
        API.get_or_init(|| Api::load().ok()).as_ref()
    }

    fn load() -> Result<Api, libloading::Error> {
        // SAFETY: loading runs the library's initializers: NVIDIA's cuDSS, which links the CUDA
        //   runtime and cuBLAS this crate's users already run; each symbol is declared with the
        //   signature cuDSS 1.0's header gives it, with enums as C `int` and handles as pointers
        unsafe {
            let library = libloading::Library::new("libcudss.so.1")?;
            Ok(Api {
                create: *library.get(b"cudssCreate\0")?,
                destroy: *library.get(b"cudssDestroy\0")?,
                set_stream: *library.get(b"cudssSetStream\0")?,
                config_create: *library.get(b"cudssConfigCreate\0")?,
                config_destroy: *library.get(b"cudssConfigDestroy\0")?,
                config_set: *library.get(b"cudssConfigSet\0")?,
                data_create: *library.get(b"cudssDataCreate\0")?,
                data_destroy: *library.get(b"cudssDataDestroy\0")?,
                data_get: *library.get(b"cudssDataGet\0")?,
                execute: *library.get(b"cudssExecute\0")?,
                matrix_create_csr: *library.get(b"cudssMatrixCreateCsr\0")?,
                matrix_create_dn: *library.get(b"cudssMatrixCreateDn\0")?,
                matrix_destroy: *library.get(b"cudssMatrixDestroy\0")?,
                _library: library,
            })
        }
    }
}

/// A cuDSS call's status as a result
fn check(call: &'static str, status: Status) -> Result<(), CudaError> {
    if status == 0 {
        Ok(())
    } else {
        Err(CudaError::Cudss { call, status })
    }
}

/// Time spent by a [`GpuLdlt`], split into its stages, each up to the stream's completion.
/// Gathering a matrix's lower triangle on the host is in none of them.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LdltTimings {
    /// Copying values and right-hand sides to the device.
    pub upload: Duration,
    /// The analysis: cuDSS's reordering and symbolic factorization, once.
    pub analysis: Duration,
    /// The numeric factorizations and refactorizations.
    pub factorization: Duration,
    /// The solves.
    pub solve: Duration,
    /// Copying solutions back.
    pub download: Duration,
    /// Factorizations, the first included.
    pub factorizations: usize,
    /// Solves.
    pub solves: usize,
}

/// A sparse symmetric indefinite matrix factored on the device by cuDSS as `L D Lᵀ`, with its
/// inertia: analyzed once for a pattern, then refactored for each new set of values on it, as an
/// interior-point method's KKT systems are.
///
/// It runs on [`Gpu`]'s stream, in cuDSS's deterministic mode, so the same matrix gives the same
/// factors and solutions every time. They're not bit for bit `ybus_sparse`'s [`Ldlt`]'s: cuDSS
/// orders and pivots its arithmetic its own way.
///
/// **Pivoting:** cuDSS's default, which pivots within its supernodes and perturbs the pivots
/// that are too small, rather than delaying them to a parent as [`Ldlt`] does. Measured on
/// saddle-point matrices (`tests/gpu_ldlt_test.rs`):
/// - **Solves need refinement:** an unrefined solve's relative residual is about 1e-4, and
///   three corrections against the matrix bring it to about 1e-13.
/// - **The inertia matched [`Ldlt`]'s** on the tests' matrices. Turning the perturbation off made
///   it wrong.
/// - **A singular matrix factors:** its zero pivot is perturbed and no error is reported, so the
///   inertia is a nearby nonsingular matrix's. Don't rely on [`GpuLdlt::factor`] to detect
///   singularity.
///
/// **Not for interior-point KKT matrices:** an OPF's KKT matrix `[W Jᵀ; J 0]` has a zero
/// constraint block, and cuDSS's ordering eliminates most constraint rows before the variables
/// they couple to, so their pivots are exactly zero. A 2×2 pivot with a coupled variable would
/// take them, but cuDSS pivots only within a supernode and can't delay a pivot to its parent, so
/// it perturbs them to its epsilon (1e-13): a third of all pivots on the OPF's matrices. Dividing
/// by those makes the factors grow (to 1e21 and more in `ybus_sparse`'s static LDLᵀ under the
/// same rule), and the rounding in later updates swamps the matrix: the inertia is wrong, and
/// refinement doesn't reduce the solves' residuals
/// (`docs/benchmarks/o5-cudss-spike-2026-11-22.md`). [`GpuLdlt::perturbed_pivots`] reports it.
/// Regularizing the matrix enough to be quasi-definite removes the perturbations.
///
/// [`Ldlt`]: ybus_sparse::compressed_matrix::Ldlt
pub struct GpuLdlt<'g> {
    gpu: &'g Gpu,
    api: &'static Api,
    handle: Opaque,
    config: Opaque,
    data: Opaque,
    matrix: Opaque,
    solution: Opaque,
    rhs: Opaque,
    n: usize,
    /// The full pattern, which every factored matrix must have, and the positions of its lower
    ///   triangle's values among the full matrix's, in the lower triangle's CSR order
    ptr: Vec<usize>,
    idx: Vec<usize>,
    lower: Vec<usize>,
    /// The device's copies of the lower triangle, and of one right-hand side and its solution;
    ///   the descriptors point into them, so they outlive the descriptors
    _pointers: CudaSlice<i32>,
    _indices: CudaSlice<i32>,
    values: CudaSlice<f64>,
    x: CudaSlice<f64>,
    b: CudaSlice<f64>,
    /// Page-locked staging for the values and for one vector
    staged_values: PinnedHostSlice<f64>,
    staged_vector: PinnedHostSlice<f64>,
    factored: bool,
    timings: LdltTimings,
}

impl<'g> GpuLdlt<'g> {
    /// Analyzes `a`, a symmetric matrix stored in full (both triangles, as
    /// [`CsrMatrix::ldlt`] takes it), for factoring on `gpu`, or `None` if cuDSS isn't
    /// installed. Only `a`'s lower triangle is uploaded, values included, since cuDSS's matching
    /// uses them; [`GpuLdlt::factor`] then factors it, or another matrix on the pattern.
    ///
    /// Returns `CudaError::ShapeMismatch` for a non-square matrix, `CudaError::TooLarge` past
    /// 32-bit indices, and cuDSS's or the driver's errors.
    pub fn new(gpu: &'g Gpu, a: &CsrMatrix<f64>) -> Result<Option<GpuLdlt<'g>>, CudaError> {
        let Some(api) = Api::get() else {
            return Ok(None);
        };
        let (rows, columns) = a.shape();
        if rows != columns {
            return Err(CudaError::ShapeMismatch {
                expected: rows,
                found: columns,
            });
        }
        let n = rows;
        if n > i32::MAX as usize {
            return Err(CudaError::TooLarge { entries: n });
        }
        let (ptr, idx) = (a.pointers().to_vec(), a.indices().to_vec());
        let mut lower_ptr = Vec::with_capacity(n + 1);
        let mut lower_idx = Vec::new();
        let mut lower = Vec::new();
        lower_ptr.push(0i32);
        for i in 0..n {
            for (p, &j) in idx.iter().enumerate().take(ptr[i + 1]).skip(ptr[i]) {
                if j <= i {
                    lower_idx.push(j as i32);
                    lower.push(p);
                }
            }
            lower_ptr.push(i32::try_from(lower.len()).map_err(|_| CudaError::TooLarge {
                entries: lower.len(),
            })?);
        }
        let stream = &gpu.stream;
        gpu.context.bind_to_thread()?;
        let start = Instant::now();
        let pointers = stream.clone_htod(&lower_ptr)?;
        let indices = stream.clone_htod(&lower_idx)?;
        let values = stream.alloc_zeros::<f64>(lower.len().max(1))?;
        let x = stream.alloc_zeros::<f64>(n.max(1))?;
        let b = stream.alloc_zeros::<f64>(n.max(1))?;
        // SAFETY: page-locked host buffers, written in full before each copy that reads them
        let staged_values = unsafe { gpu.context.alloc_pinned::<f64>(lower.len().max(1))? };
        let staged_vector = unsafe { gpu.context.alloc_pinned::<f64>(n.max(1))? };
        stream.synchronize()?;
        let upload = start.elapsed();
        let mut ldlt = GpuLdlt {
            gpu,
            api,
            handle: std::ptr::null_mut(),
            config: std::ptr::null_mut(),
            data: std::ptr::null_mut(),
            matrix: std::ptr::null_mut(),
            solution: std::ptr::null_mut(),
            rhs: std::ptr::null_mut(),
            n,
            ptr,
            idx,
            lower,
            _pointers: pointers,
            _indices: indices,
            values,
            x,
            b,
            staged_values,
            staged_vector,
            factored: false,
            timings: LdltTimings {
                upload,
                ..LdltTimings::default()
            },
        };
        // Each object is stored as soon as it exists, so `Drop` destroys whatever was made if a
        //   later step fails
        // SAFETY: each call gets out-parameters for new objects, or objects created just above
        //   and still live; the stream is `gpu`'s, which outlives this `GpuLdlt` (it borrows
        //   `gpu`), and the driver's and the runtime's stream types are the same object behind
        //   different names
        unsafe {
            check("cudssCreate", (api.create)(&mut ldlt.handle))?;
            check(
                "cudssSetStream",
                (api.set_stream)(ldlt.handle, stream.cu_stream() as *mut c_void),
            )?;
            check("cudssConfigCreate", (api.config_create)(&mut ldlt.config))?;
            let deterministic: c_int = 1;
            check(
                "cudssConfigSet",
                (api.config_set)(
                    ldlt.config,
                    CONFIG_DETERMINISTIC_MODE,
                    (&deterministic as *const c_int).cast(),
                    std::mem::size_of::<c_int>(),
                ),
            )?;
            check(
                "cudssDataCreate",
                (api.data_create)(ldlt.handle, &mut ldlt.data),
            )?;
        }
        ldlt.describe(lower_ptr.len())?;
        // `a`'s values before the analysis: cuDSS 1.0's default matching (`AUTO`) computes its
        //   permutation and scaling from them, and from zeros it does nothing useful
        let start = Instant::now();
        ldlt.stage_values(a)?;
        ldlt.timings.upload += start.elapsed();
        let start = Instant::now();
        ldlt.execute(PHASE_ANALYSIS, "cudssExecute (analysis)")?;
        stream.synchronize()?;
        ldlt.timings.analysis = start.elapsed();
        Ok(Some(ldlt))
    }

    /// The descriptors of the matrix and of one right-hand side and solution, over the device
    ///   arrays this `GpuLdlt` owns
    fn describe(&mut self, pointers: usize) -> Result<(), CudaError> {
        let stream = &self.gpu.stream;
        let (p, _p) = self._pointers.device_ptr(stream);
        let (i, _i) = self._indices.device_ptr(stream);
        let (v, _v) = self.values.device_ptr(stream);
        let (x, _x) = self.x.device_ptr_mut(stream);
        let (b, _b) = self.b.device_ptr(stream);
        let n = self.n as i64;
        debug_assert_eq!(pointers, self.n + 1);
        // SAFETY: the device arrays hold `n + 1` row offsets, and `lower.len()` column indices
        //   and values, of the lower triangle in zero-based CSR, and `n` values for each vector,
        //   as declared; this `GpuLdlt` owns them and destroys the descriptors first
        unsafe {
            check(
                "cudssMatrixCreateCsr",
                (self.api.matrix_create_csr)(
                    &mut self.matrix,
                    n,
                    n,
                    self.lower.len() as i64,
                    p as *const c_void,
                    std::ptr::null(),
                    i as *const c_void,
                    v as *const c_void,
                    cudaDataType::CUDA_R_32I as c_int,
                    cudaDataType::CUDA_R_32I as c_int,
                    cudaDataType::CUDA_R_64F as c_int,
                    MTYPE_SYMMETRIC,
                    MVIEW_LOWER,
                    BASE_ZERO,
                ),
            )?;
            check(
                "cudssMatrixCreateDn",
                (self.api.matrix_create_dn)(
                    &mut self.solution,
                    n,
                    1,
                    n,
                    x as *const c_void,
                    cudaDataType::CUDA_R_64F as c_int,
                    LAYOUT_COL_MAJOR,
                ),
            )?;
            check(
                "cudssMatrixCreateDn",
                (self.api.matrix_create_dn)(
                    &mut self.rhs,
                    n,
                    1,
                    n,
                    b as *const c_void,
                    cudaDataType::CUDA_R_64F as c_int,
                    LAYOUT_COL_MAJOR,
                ),
            )?;
        }
        Ok(())
    }

    /// One `cudssExecute` phase, and the factorization's error indicator after it
    fn execute(&mut self, phase: c_int, call: &'static str) -> Result<(), CudaError> {
        // SAFETY: every object was created in `new` and is live
        unsafe {
            check(
                call,
                (self.api.execute)(
                    self.handle,
                    phase,
                    self.config,
                    self.data,
                    self.matrix,
                    self.solution,
                    self.rhs,
                ),
            )
        }
    }

    /// An integer the data object reports, read once the stream is done
    fn data_ints<const N: usize>(&self, param: c_int) -> Result<[c_int; N], CudaError> {
        self.gpu.stream.synchronize()?;
        let mut value = [0 as c_int; N];
        let mut written = 0usize;
        // SAFETY: `value` is a host buffer of the size passed, which cuDSS writes at most
        unsafe {
            check(
                "cudssDataGet",
                (self.api.data_get)(
                    self.handle,
                    self.data,
                    param,
                    value.as_mut_ptr().cast(),
                    std::mem::size_of_val(&value),
                    &mut written,
                ),
            )?;
        }
        Ok(value)
    }

    /// `a`'s lower triangle's values, gathered into the page-locked staging and copied to the
    ///   device's array, which the matrix's descriptor points into
    fn stage_values(&mut self, a: &CsrMatrix<f64>) -> Result<(), CudaError> {
        {
            let staged = self.staged_values.as_mut_slice()?;
            let data = a.values();
            for (slot, &p) in staged.iter_mut().zip(&self.lower) {
                *slot = data[p];
            }
        }
        self.gpu
            .stream
            .memcpy_htod(&self.staged_values, &mut self.values)?;
        self.gpu.stream.synchronize()?;
        Ok(())
    }

    /// The order of the matrix.
    pub fn dim(&self) -> usize {
        self.n
    }

    /// Factors `a`, which must have the analyzed pattern: the first time a factorization,
    /// then refactorizations, which reuse the analysis, including the matching and scaling
    /// cuDSS chose from the values [`GpuLdlt::new`] was given.
    ///
    /// Returns `CudaError::PatternMismatch` for another pattern, `CudaError::Singular` if cuDSS
    /// reports a failure in its error indicator, and cuDSS's or the driver's errors. A singular
    /// matrix usually factors anyway, perturbed (see [`GpuLdlt`]).
    pub fn factor(&mut self, a: &CsrMatrix<f64>) -> Result<(), CudaError> {
        if a.shape() != (self.n, self.n)
            || a.pointers() != self.ptr.as_slice()
            || a.indices() != self.idx.as_slice()
        {
            return Err(CudaError::PatternMismatch);
        }
        self.gpu.context.bind_to_thread()?;
        let start = Instant::now();
        self.stage_values(a)?;
        self.timings.upload += start.elapsed();
        let start = Instant::now();
        if self.factored {
            self.execute(PHASE_REFACTORIZATION, "cudssExecute (refactorization)")?;
        } else {
            self.execute(PHASE_FACTORIZATION, "cudssExecute (factorization)")?;
        }
        let [info] = self.data_ints::<1>(DATA_INFO)?;
        self.timings.factorization += start.elapsed();
        self.timings.factorizations += 1;
        if info != 0 {
            self.factored = false;
            return Err(CudaError::Singular { info });
        }
        self.factored = true;
        Ok(())
    }

    /// Solves `A x = b` with the last factorization.
    ///
    /// Returns `CudaError::ShapeMismatch` unless `b` has one entry per row, and
    /// `CudaError::NotFactored` before a successful [`GpuLdlt::factor`].
    pub fn solve(&mut self, b: &[f64]) -> Result<Vec<f64>, CudaError> {
        if b.len() != self.n {
            return Err(CudaError::ShapeMismatch {
                expected: self.n,
                found: b.len(),
            });
        }
        if !self.factored {
            return Err(CudaError::NotFactored);
        }
        let stream = &self.gpu.stream;
        self.gpu.context.bind_to_thread()?;
        let start = Instant::now();
        self.staged_vector.as_mut_slice()?[..self.n].copy_from_slice(b);
        stream.memcpy_htod(&self.staged_vector, &mut self.b)?;
        stream.synchronize()?;
        self.timings.upload += start.elapsed();
        let start = Instant::now();
        self.execute(PHASE_SOLVE, "cudssExecute (solve)")?;
        stream.synchronize()?;
        self.timings.solve += start.elapsed();
        self.timings.solves += 1;
        let start = Instant::now();
        stream.memcpy_dtoh(&self.x, &mut self.staged_vector)?;
        stream.synchronize()?;
        let x = self.staged_vector.as_slice()?[..self.n].to_vec();
        self.timings.download += start.elapsed();
        Ok(x)
    }

    /// The last factored matrix's inertia, as cuDSS reports it from `D`.
    ///
    /// Returns `CudaError::NotFactored` before a successful [`GpuLdlt::factor`].
    pub fn inertia(&self) -> Result<Inertia, CudaError> {
        if !self.factored {
            return Err(CudaError::NotFactored);
        }
        let [positive, negative] = self.data_ints::<2>(DATA_INERTIA)?;
        Ok(Inertia {
            positive: positive.max(0) as usize,
            negative: negative.max(0) as usize,
        })
    }

    /// How many 1×1 pivots the last factorization perturbed (`CUDSS_DATA_NPIVOTS`): pivots too
    /// small to take, replaced by cuDSS's signed epsilon. 0 means the factorization is exact up
    /// to rounding; many mean the factors, the solves and the inertia can't be trusted (see
    /// [`GpuLdlt`]).
    ///
    /// Returns `CudaError::NotFactored` before a successful [`GpuLdlt::factor`].
    pub fn perturbed_pivots(&self) -> Result<usize, CudaError> {
        if !self.factored {
            return Err(CudaError::NotFactored);
        }
        let [count] = self.data_ints::<1>(DATA_NPIVOTS)?;
        Ok(count.max(0) as usize)
    }

    /// The time spent so far, by stage.
    pub fn timings(&self) -> LdltTimings {
        self.timings
    }

    /// Restarts the timings.
    pub fn reset_timings(&mut self) {
        self.timings = LdltTimings::default();
    }
}

impl Drop for GpuLdlt<'_> {
    fn drop(&mut self) {
        let api = self.api;
        let _ = self.gpu.context.bind_to_thread();
        let _ = self.gpu.stream.synchronize();
        // SAFETY: each object was created in `new` (or is null, never made) and is destroyed
        //   only here, the descriptors before the data, the config and the handle, and all
        //   before the device arrays they point into, which drop after this
        unsafe {
            for matrix in [self.matrix, self.solution, self.rhs] {
                if !matrix.is_null() {
                    let _ = (api.matrix_destroy)(matrix);
                }
            }
            if !self.data.is_null() {
                let _ = (api.data_destroy)(self.handle, self.data);
            }
            if !self.config.is_null() {
                let _ = (api.config_destroy)(self.config);
            }
            if !self.handle.is_null() {
                let _ = (api.destroy)(self.handle);
            }
        }
    }
}
