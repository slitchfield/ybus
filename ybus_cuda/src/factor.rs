//! A CPU factorization's triangles on the GPU, solving blocks of right-hand sides.

use std::time::{Duration, Instant};

use cudarc::driver::PinnedHostSlice;
use ybus_sparse::compressed_matrix::{Cholesky, Lu};

use crate::sparse::{Block, Csr, Fill};
use crate::{CudaError, Gpu};

/// The widest block [`GpuFactor`] chooses by itself
const MAX_DEFAULT_WIDTH: usize = 1024;

/// The widest block a single-row factor's device paths use: cuSPARSE 12.8's SpSM gives wrong
/// results for a 1×1 matrix with more right-hand sides than this
const SINGLE_ROW_WIDTH: usize = 64;

/// Time spent by [`GpuFactor::solve_many`], split into its stages. The host's packing of the
/// right-hand sides (with the factorization's permutation) and unpacking of the solutions is
/// in none of them.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Timings {
    /// Copying blocks of right-hand sides to the device.
    pub upload: Duration,
    /// The triangular solves.
    pub solve: Duration,
    /// Copying blocks of solutions back.
    pub download: Duration,
    /// Blocks solved.
    pub blocks: usize,
}

/// A [`ybus_sparse`] factorization uploaded to the GPU, for solving many right-hand sides
/// against it: they go in blocks of [`GpuFactor::block_width`] columns, each solved by
/// cuSPARSE's SpSM (forward with the lower triangle, back with the upper). Blocks go to and
/// from the device through page-locked host memory.
pub struct GpuFactor<'g> {
    pub(crate) gpu: &'g Gpu,
    pub(crate) n: usize,
    /// Row `i` of a packed block is entry `gather[i]` of a right-hand side
    gather: Vec<usize>,
    /// Entry `scatter[i]` of a solution is row `i` of a solved block
    pub(crate) scatter: Vec<usize>,
    /// Entry `u` of a right-hand side is row `rhs_row[u]` of a packed block
    pub(crate) rhs_row: Vec<usize>,
    /// Entry `u` of a solution is row `solution_row[u]` of a solved block
    pub(crate) solution_row: Vec<usize>,
    /// For a single row, the two triangles' diagonal entries: `solve_many` solves such a factor
    /// on the host (see `SINGLE_ROW_WIDTH`)
    single: Option<(f64, f64)>,
    // Declared before the triangles, so it's destroyed before them
    pub(crate) block: Option<Block>,
    pub(crate) lower: Csr,
    pub(crate) upper: Csr,
    width: usize,
    /// Page-locked host memory for every copy to and from the device, grown as needed
    pub(crate) staging: Staging,
    timings: Timings,
}

/// Page-locked host memory, grown as needed: copies to and from it run at the bus's full rate,
/// where pageable memory has the driver copy through a buffer of its own.
pub(crate) struct Staging {
    slice: Option<PinnedHostSlice<f64>>,
}

impl Staging {
    /// At least `len` values of page-locked memory, as a slice of exactly `len`.
    pub(crate) fn get(&mut self, gpu: &Gpu, len: usize) -> Result<&mut [f64], CudaError> {
        if self.slice.as_ref().is_none_or(|slice| slice.len() < len) {
            self.slice = None;
            // SAFETY: plain page-locked memory (no flags: not write-combined, which would make
            //   reading it from the host slow), zeroed before any slice of it is made, so it's
            //   never read uninitialized
            let slice = unsafe {
                let mut slice = gpu.context.alloc_pinned_with_flags::<f64>(len.max(1), 0)?;
                std::ptr::write_bytes(slice.as_mut_ptr()?, 0, slice.len());
                slice
            };
            self.slice = Some(slice);
        }
        let slice = self
            .slice
            .as_mut()
            .expect("allocated above")
            .as_mut_slice()?;
        Ok(&mut slice[..len])
    }
}

/// The inverse of a permutation: `inverse[perm[i]] = i`.
fn inverse(perm: &[usize]) -> Vec<usize> {
    let mut inverse = vec![0; perm.len()];
    for (i, &p) in perm.iter().enumerate() {
        inverse[p] = i;
    }
    inverse
}

impl<'g> GpuFactor<'g> {
    /// Uploads a Cholesky factorization `P A Pᵀ = L Lᵀ`: `L` and `Lᵀ` as triangles, `P` as the
    /// packing and unpacking order.
    pub fn from_cholesky(gpu: &'g Gpu, cholesky: &Cholesky<f64>) -> Result<Self, CudaError> {
        let l = cholesky.l();
        // `L` in CSC form is `Lᵀ` in CSR form
        let l_transpose = l.to_csc();
        let perm = cholesky.permutation().to_vec();
        Self::new(
            gpu,
            (l.pointers(), l.indices(), l.values()),
            (
                l_transpose.pointers(),
                l_transpose.indices(),
                l_transpose.values(),
            ),
            perm.clone(),
            perm,
        )
    }

    /// Uploads an LU factorization `P A Q = L U`: `L` (with its unit diagonal) and `U` as
    /// triangles, `P` as the packing order and `Q` as the unpacking order.
    pub fn from_lu(gpu: &'g Gpu, lu: &Lu<f64>) -> Result<Self, CudaError> {
        let l = lu.l().to_csr();
        let u = lu.u().to_csr();
        Self::new(
            gpu,
            (l.pointers(), l.indices(), l.values()),
            (u.pointers(), u.indices(), u.values()),
            lu.permutation().to_vec(),
            lu.column_permutation().to_vec(),
        )
    }

    fn new(
        gpu: &'g Gpu,
        lower: (&[usize], &[usize], &[f64]),
        upper: (&[usize], &[usize], &[f64]),
        gather: Vec<usize>,
        scatter: Vec<usize>,
    ) -> Result<Self, CudaError> {
        let n = gather.len();
        let single = (n == 1).then(|| (lower.2[0], upper.2[0]));
        let upload = |(pointers, indices, values), fill| {
            Csr::upload(
                &gpu.context,
                &gpu.stream,
                n,
                n,
                pointers,
                indices,
                values,
                Some(fill),
            )
        };
        let lower = upload(lower, Fill::Lower)?;
        let upper = upload(upper, Fill::Upper)?;
        // Two blocks of `n` doubles per column, in at most a quarter of the free memory
        let (free, _) = gpu.memory()?;
        let width = (free / 4 / (2 * 8 * n.max(1))).clamp(1, MAX_DEFAULT_WIDTH);
        let mut factor = GpuFactor {
            gpu,
            n,
            rhs_row: inverse(&gather),
            solution_row: inverse(&scatter),
            gather,
            scatter,
            single,
            block: None,
            lower,
            upper,
            width,
            staging: Staging { slice: None },
            timings: Timings::default(),
        };
        factor.set_block_width(width);
        Ok(factor)
    }

    /// The number of rows of the factorization.
    pub fn dim(&self) -> usize {
        self.n
    }

    /// The number of right-hand sides solved together.
    pub fn block_width(&self) -> usize {
        self.width
    }

    /// Sets the number of right-hand sides solved together (at least 1, and at most 64 for a
    /// single-row factor). A block takes two `dim() × width` arrays of doubles on the device.
    pub fn set_block_width(&mut self, width: usize) {
        let most = if self.n == 1 {
            SINGLE_ROW_WIDTH
        } else {
            usize::MAX
        };
        self.width = width.clamp(1, most);
    }

    /// The time spent in [`GpuFactor::solve_many`] since the factor was made or
    /// [`GpuFactor::reset_timings`] was last called.
    pub fn timings(&self) -> Timings {
        self.timings
    }

    /// Sets the timings back to zero.
    pub fn reset_timings(&mut self) {
        self.timings = Timings::default();
    }

    /// Makes sure the block (and its solves' analysis) is for the current width.
    pub(crate) fn ensure_block(&mut self) -> Result<(), CudaError> {
        if self.block.as_ref().map(Block::width) != Some(self.width) {
            // Drop the old block's analysis before making the new one's
            self.block = None;
            self.block = Some(Block::new(
                &self.gpu.context,
                &self.gpu.stream,
                &self.gpu.sparse,
                &self.lower,
                &self.upper,
                self.n,
                self.width,
            )?);
        }
        Ok(())
    }

    /// Solves `A x = b` for each right-hand side `b` in `rhs`, returning the solutions in the
    /// same order, as `ybus_sparse`'s `solve_many` does, to within rounding. The last block is
    /// padded with zero columns. Returns `CudaError::ShapeMismatch` for the first right-hand
    /// side whose length isn't [`GpuFactor::dim`].
    pub fn solve_many(&mut self, rhs: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, CudaError> {
        let n = self.n;
        if let Some(b) = rhs.iter().find(|b| b.len() != n) {
            return Err(CudaError::ShapeMismatch {
                expected: n,
                found: b.len(),
            });
        }
        if n == 0 || rhs.is_empty() {
            return Ok(vec![Vec::new(); rhs.len()]);
        }
        if let Some((lower, upper)) = self.single {
            // Forward, then back substitution, dividing in the CPU's order
            return Ok(rhs.iter().map(|b| vec![b[0] / lower / upper]).collect());
        }
        let width = self.width;
        self.ensure_block()?;
        let block = self.block.as_mut().expect("made above");
        let stream = &self.gpu.stream;
        let host = self.staging.get(self.gpu, n * width)?;
        let mut solutions = Vec::with_capacity(rhs.len());
        for chunk in rhs.chunks(width) {
            for (column, b) in host.chunks_exact_mut(n).zip(chunk) {
                for (packed, &from) in column.iter_mut().zip(&self.gather) {
                    *packed = b[from];
                }
            }
            for column in host.chunks_exact_mut(n).skip(chunk.len()) {
                column.fill(0.0);
            }

            let start = Instant::now();
            block.upload(stream, host)?;
            stream.synchronize()?;
            let uploaded = Instant::now();
            block.solve(&self.gpu.sparse, stream, &self.lower, &self.upper)?;
            stream.synchronize()?;
            let solved = Instant::now();
            block.download(stream, host)?;
            stream.synchronize()?;
            let downloaded = Instant::now();
            self.timings.upload += uploaded - start;
            self.timings.solve += solved - uploaded;
            self.timings.download += downloaded - solved;
            self.timings.blocks += 1;

            for column in host.chunks_exact(n).take(chunk.len()) {
                let mut x = vec![0.0; n];
                for (&to, &value) in self.scatter.iter().zip(column) {
                    x[to] = value;
                }
                solutions.push(x);
            }
        }
        Ok(solutions)
    }
}
