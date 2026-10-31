//! Properties on `SparseMatrix`, which forward to `CompressedMatrix` until Phase 2 removes
//! them from `SparseMatrix`, except `is_square`, which stays.

use num_traits::Float;

use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar> SparseMatrix<T> {
    /// Whether the matrix has as many rows as columns.
    pub fn is_square(&self) -> bool {
        self.shape.0 == self.shape.1
    }

    /// The fraction of positions that hold an entry; forwards to `CompressedMatrix::density`.
    pub fn density(&self) -> f64 {
        self.to_csr().density()
    }

    /// The sum of all entries in row-major order; forwards to `CsrMatrix::sum`.
    pub fn sum(&self) -> T {
        self.to_csr().sum()
    }

    /// The sum of the main diagonal; forwards to `CompressedMatrix::trace`.
    pub fn trace(&self) -> Result<T, SparseError> {
        self.to_csr().trace()
    }

    /// Whether the matrix equals its transpose; forwards to `CompressedMatrix::is_symmetric`.
    pub fn is_symmetric(&self) -> bool {
        self.to_csr().is_symmetric()
    }
}

impl<T: Scalar + Float> SparseMatrix<T> {
    /// The Frobenius norm; forwards to `CompressedMatrix::norm_frobenius`.
    pub fn norm_frobenius(&self) -> T {
        self.to_csr().norm_frobenius()
    }

    /// The 1-norm; forwards to `CompressedMatrix::norm_1`.
    pub fn norm_1(&self) -> T {
        self.to_csr().norm_1()
    }

    /// The infinity norm; forwards to `CompressedMatrix::norm_inf`.
    pub fn norm_inf(&self) -> T {
        self.to_csr().norm_inf()
    }
}
