//! Solvers on `SparseMatrix`, which forward to `CompressedMatrix` until Phase 2 removes them
//! from `SparseMatrix`.

use num_traits::Float;

use super::SparseMatrix;
use crate::compressed_matrix::Lu;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar + Float> SparseMatrix<T> {
    /// Forward substitution; forwards to `CsrMatrix::solve_lower_triangular`.
    pub fn solve_lower_triangular(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csr().solve_lower_triangular(b)
    }

    /// Back substitution; forwards to `CsrMatrix::solve_upper_triangular`.
    pub fn solve_upper_triangular(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csr().solve_upper_triangular(b)
    }

    /// LU factorization; forwards to `CsrMatrix::lu`.
    pub fn lu(&self) -> Result<Lu<T>, SparseError> {
        self.to_csr().lu()
    }

    /// Solves `A x = b` by LU; forwards to `CsrMatrix::solve`.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csr().solve(b)
    }

    /// The Cholesky factor `L`; forwards to `CsrMatrix::cholesky`.
    pub fn cholesky(&self) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr().cholesky().map(|l| l.to_sparse())
    }

    /// Solves `A x = b` by Cholesky; forwards to `CsrMatrix::solve_spd`.
    pub fn solve_spd(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csr().solve_spd(b)
    }

    /// Conjugate gradient; forwards to `CsrMatrix::conjugate_gradient`.
    pub fn conjugate_gradient(
        &self,
        b: &[T],
        tolerance: T,
        max_iterations: usize,
    ) -> Result<Vec<T>, SparseError> {
        self.to_csr()
            .conjugate_gradient(b, tolerance, max_iterations)
    }
}
