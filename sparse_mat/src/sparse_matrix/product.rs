//! Matrix products on `SparseMatrix`, which forward to `CompressedMatrix` until Phase 2 removes
//! them from `SparseMatrix`.

use std::ops::Mul;

use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar> SparseMatrix<T> {
    /// The matrix-vector product `A x`; forwards to `CsrMatrix::mul_vec`.
    pub fn mul_vec(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csr().mul_vec(x)
    }

    /// The vector-matrix product `x A`; forwards to `CscMatrix::vec_mul`.
    pub fn vec_mul(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        self.to_csc().vec_mul(x)
    }

    /// The product with a dense matrix; forwards to `CsrMatrix::mul_dense`.
    pub fn mul_dense(&self, b: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        self.to_csr().mul_dense(b)
    }

    /// The sparse matrix product; forwards to `CsrMatrix::matmul`. Entries that cancel to
    /// exactly zero are not stored, since `SparseMatrix` never stores a zero.
    pub fn matmul(&self, other: &SparseMatrix<T>) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .matmul(&other.to_csr())
            .map(|product| product.to_sparse())
    }
}

/// Matrix product; see `SparseMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar> Mul<&SparseMatrix<T>> for &SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn mul(self, other: &SparseMatrix<T>) -> SparseMatrix<T> {
        self.matmul(other).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Matrix product; see `SparseMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar> Mul<SparseMatrix<T>> for SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn mul(self, other: SparseMatrix<T>) -> SparseMatrix<T> {
        &self * &other
    }
}

/// Matrix product; see `SparseMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar> Mul<&SparseMatrix<T>> for SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn mul(self, other: &SparseMatrix<T>) -> SparseMatrix<T> {
        &self * other
    }
}

/// Matrix product; see `SparseMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar> Mul<SparseMatrix<T>> for &SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn mul(self, other: SparseMatrix<T>) -> SparseMatrix<T> {
        self * &other
    }
}

/// Matrix-vector product `A x`; see `SparseMatrix::mul_vec`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's column count.
impl<T: Scalar> Mul<&[T]> for &SparseMatrix<T> {
    type Output = Vec<T>;

    fn mul(self, x: &[T]) -> Vec<T> {
        self.mul_vec(x).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Matrix-vector product `A x`; see `SparseMatrix::mul_vec`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's column count.
impl<T: Scalar> Mul<&Vec<T>> for &SparseMatrix<T> {
    type Output = Vec<T>;

    fn mul(self, x: &Vec<T>) -> Vec<T> {
        self.mul_vec(x).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Vector-matrix product `x A`; see `SparseMatrix::vec_mul`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's row count.
impl<T: Scalar> Mul<&SparseMatrix<T>> for &[T] {
    type Output = Vec<T>;

    fn mul(self, matrix: &SparseMatrix<T>) -> Vec<T> {
        matrix.vec_mul(self).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Vector-matrix product `x A`; see `SparseMatrix::vec_mul`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's row count.
impl<T: Scalar> Mul<&SparseMatrix<T>> for &Vec<T> {
    type Output = Vec<T>;

    fn mul(self, matrix: &SparseMatrix<T>) -> Vec<T> {
        matrix.vec_mul(self).unwrap_or_else(|err| panic!("{err}"))
    }
}
