//! A CPU factorization's triangles on the GPU, solving blocks of right-hand sides.

use std::time::{Duration, Instant};

use sparse_mat::compressed_matrix::{Cholesky, Lu};

use crate::sparse::{Block, Fill, Triangle};
use crate::{CudaError, Gpu};

/// The widest block [`GpuFactor`] chooses by itself
const MAX_DEFAULT_WIDTH: usize = 1024;

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

/// A [`sparse_mat`] factorization uploaded to the GPU, for solving many right-hand sides
/// against it: they go in blocks of [`GpuFactor::block_width`] columns, each solved by
/// cuSPARSE's SpSM (forward with the lower triangle, back with the upper), with the
/// factorization's permutations applied on the host as the blocks are packed and unpacked.
pub struct GpuFactor<'g> {
    gpu: &'g Gpu,
    n: usize,
    /// Row `i` of a packed block is entry `gather[i]` of a right-hand side
    gather: Vec<usize>,
    /// Entry `scatter[i]` of a solution is row `i` of a solved block
    scatter: Vec<usize>,
    /// For a single row, the two triangles' diagonal entries: such a factor is solved on the
    /// host, since cuSPARSE 12.8's SpSM gives wrong results for a 1×1 matrix with more than
    /// 64 right-hand sides
    single: Option<(f64, f64)>,
    // Declared before the triangles, so it's destroyed before them
    block: Option<Block>,
    lower: Triangle,
    upper: Triangle,
    width: usize,
    host: Vec<f64>,
    timings: Timings,
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
            Triangle::upload(
                &gpu.context,
                &gpu.stream,
                n,
                pointers,
                indices,
                values,
                fill,
            )
        };
        let lower = upload(lower, Fill::Lower)?;
        let upper = upload(upper, Fill::Upper)?;
        // Two blocks of `n` doubles per column, in at most a quarter of the free memory
        let (free, _) = gpu.memory()?;
        let width = (free / 4 / (2 * 8 * n.max(1))).clamp(1, MAX_DEFAULT_WIDTH);
        Ok(GpuFactor {
            gpu,
            n,
            gather,
            scatter,
            single,
            block: None,
            lower,
            upper,
            width,
            host: Vec::new(),
            timings: Timings::default(),
        })
    }

    /// The number of rows of the factorization.
    pub fn dim(&self) -> usize {
        self.n
    }

    /// The number of right-hand sides solved together.
    pub fn block_width(&self) -> usize {
        self.width
    }

    /// Sets the number of right-hand sides solved together (at least 1). A block takes two
    /// `dim() × width` arrays of doubles on the device.
    pub fn set_block_width(&mut self, width: usize) {
        self.width = width.max(1);
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

    /// Solves `A x = b` for each right-hand side `b` in `rhs`, returning the solutions in the
    /// same order, as `sparse_mat`'s `solve_many` does, to within rounding. The last block is
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
        if self.block.as_ref().map(Block::width) != Some(width) {
            // Drop the old block's analysis before making the new one's
            self.block = None;
            self.block = Some(Block::new(
                &self.gpu.context,
                &self.gpu.stream,
                &self.gpu.sparse,
                &self.lower,
                &self.upper,
                n,
                width,
            )?);
        }
        let block = self.block.as_mut().expect("made above");
        let stream = &self.gpu.stream;
        self.host.resize(n * width, 0.0);
        let mut solutions = Vec::with_capacity(rhs.len());
        for chunk in rhs.chunks(width) {
            for (column, b) in self.host.chunks_exact_mut(n).zip(chunk) {
                for (packed, &from) in column.iter_mut().zip(&self.gather) {
                    *packed = b[from];
                }
            }
            for column in self.host.chunks_exact_mut(n).skip(chunk.len()) {
                column.fill(0.0);
            }

            let start = Instant::now();
            block.upload(stream, &self.host)?;
            stream.synchronize()?;
            let uploaded = Instant::now();
            block.solve(&self.gpu.sparse, stream, &self.lower, &self.upper)?;
            stream.synchronize()?;
            let solved = Instant::now();
            block.download(stream, &mut self.host)?;
            stream.synchronize()?;
            let downloaded = Instant::now();
            self.timings.upload += uploaded - start;
            self.timings.solve += solved - uploaded;
            self.timings.download += downloaded - solved;
            self.timings.blocks += 1;

            for column in self.host.chunks_exact(n).take(chunk.len()) {
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
