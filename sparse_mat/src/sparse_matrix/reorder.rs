//! Reordering on `SparseMatrix`, which forwards to `CompressedMatrix` until Phase 2 removes it
//! from `SparseMatrix`.

use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar> SparseMatrix<T> {
    /// Reorders rows and columns; forwards to `CompressedMatrix::permute`.
    pub fn permute(
        &self,
        row_perm: &[usize],
        col_perm: &[usize],
    ) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .permute(row_perm, col_perm)
            .map(|permuted| permuted.to_sparse())
    }

    /// `P A Pᵀ`; forwards to `CompressedMatrix::permute_symmetric`.
    pub fn permute_symmetric(&self, perm: &[usize]) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .permute_symmetric(perm)
            .map(|permuted| permuted.to_sparse())
    }
}
