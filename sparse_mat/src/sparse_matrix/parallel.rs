//! Parallel products on `SparseMatrix` (feature `rayon`), which forward to `CompressedMatrix`
//! until Phase 2 removes them from `SparseMatrix`.

use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar + Send + Sync> SparseMatrix<T> {
    /// Parallel version of `mul_vec`; forwards to `CsrMatrix::par_mul_vec`.
    pub fn par_mul_vec(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csr().par_mul_vec(x)
    }

    /// Parallel version of `vec_mul`; forwards to `CscMatrix::par_vec_mul`.
    pub fn par_vec_mul(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csc().par_vec_mul(x)
    }

    /// Parallel version of `matmul`; forwards to `CsrMatrix::par_matmul`, dropping entries
    /// that cancel to zero.
    pub fn par_matmul(&self, other: &SparseMatrix<T>) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .par_matmul(&other.to_csr())
            .map(|product| product.to_sparse())
    }
}
