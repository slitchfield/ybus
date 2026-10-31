//! Block construction on `SparseMatrix`, which forwards to `CompressedMatrix` until Phase 2
//! removes it from `SparseMatrix`.

use std::borrow::Borrow;

use super::SparseMatrix;
use crate::compressed_matrix::CsrMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar> SparseMatrix<T> {
    /// The Kronecker product; forwards to `CompressedMatrix::kron`.
    pub fn kron(&self, other: &SparseMatrix<T>) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .kron(&other.to_csr())
            .map(|product| product.to_sparse())
    }

    /// Joins matrices side by side; forwards to `CompressedMatrix::hstack`.
    pub fn hstack<I>(blocks: I) -> Result<SparseMatrix<T>, SparseError>
    where
        I: IntoIterator,
        I::Item: Borrow<SparseMatrix<T>>,
    {
        let blocks: Vec<CsrMatrix<T>> = blocks.into_iter().map(|b| b.borrow().to_csr()).collect();
        CsrMatrix::hstack(&blocks).map(|stacked| stacked.to_sparse())
    }

    /// Joins matrices top to bottom; forwards to `CompressedMatrix::vstack`.
    pub fn vstack<I>(blocks: I) -> Result<SparseMatrix<T>, SparseError>
    where
        I: IntoIterator,
        I::Item: Borrow<SparseMatrix<T>>,
    {
        let blocks: Vec<CsrMatrix<T>> = blocks.into_iter().map(|b| b.borrow().to_csr()).collect();
        CsrMatrix::vstack(&blocks).map(|stacked| stacked.to_sparse())
    }
}
