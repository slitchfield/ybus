//! Standard trait implementations for `SparseMatrix`.

use std::fmt;

use super::access::Entries;
use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

/// Shows the shape and the stored entries in row-major order.
impl<T: Scalar + fmt::Debug> fmt::Debug for SparseMatrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let entries: Vec<(usize, usize, T)> = self.iter().collect();
        f.debug_struct("SparseMatrix")
            .field("shape", &self.shape)
            .field("entries", &entries)
            .finish()
    }
}

/// Two matrices are equal when they have the same shape and the same stored entries, whatever
/// order they were inserted in.
///
/// Like floats themselves this is only a partial equivalence: a matrix holding `NaN` is not
/// equal to itself.
impl<T: Scalar> PartialEq for SparseMatrix<T> {
    fn eq(&self, other: &Self) -> bool {
        self.shape == other.shape && self.values == other.values
    }
}

/// Iterates over the stored entries by reference; equivalent to `SparseMatrix::iter`.
impl<T: Scalar> IntoIterator for &SparseMatrix<T> {
    type Item = (usize, usize, T);
    type IntoIter = Entries<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Consumes the matrix, yielding its stored entries as `(row, col, value)` in row-major
/// order.
impl<T: Scalar> IntoIterator for SparseMatrix<T> {
    type Item = (usize, usize, T);
    type IntoIter = Entries<T>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Builds an `R x C` matrix from a fixed-size array of rows, skipping zeros.
impl<T: Scalar, const R: usize, const C: usize> From<[[T; C]; R]> for SparseMatrix<T> {
    fn from(rows: [[T; C]; R]) -> Self {
        let mut values = std::collections::HashMap::new();
        for (row_idx, row) in rows.iter().enumerate() {
            for (col_idx, &val) in row.iter().enumerate() {
                if val != T::zero() {
                    values.insert((row_idx, col_idx), val);
                }
            }
        }
        SparseMatrix::from_values((R, C), values)
    }
}

/// Builds a matrix from dense rows; see `SparseMatrix::from_dense`. Fails with
/// `SparseError::RaggedRows` if the rows differ in length.
impl<T: Scalar> TryFrom<Vec<Vec<T>>> for SparseMatrix<T> {
    type Error = SparseError;

    fn try_from(rows: Vec<Vec<T>>) -> Result<Self, Self::Error> {
        SparseMatrix::from_dense(rows)
    }
}

/// Expands the matrix to dense rows; see `SparseMatrix::to_dense`.
impl<T: Scalar> From<&SparseMatrix<T>> for Vec<Vec<T>> {
    fn from(matrix: &SparseMatrix<T>) -> Self {
        matrix.to_dense()
    }
}

/// Expands the matrix to dense rows; see `SparseMatrix::to_dense`.
impl<T: Scalar> From<SparseMatrix<T>> for Vec<Vec<T>> {
    fn from(matrix: SparseMatrix<T>) -> Self {
        matrix.to_dense()
    }
}
