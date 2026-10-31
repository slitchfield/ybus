//! Read access to a matrix's stored entries, and the sparse vector type that compressed
//! matrices return for rows, columns and diagonals.

use std::ops::Index;

use super::SparseMatrix;
use crate::Scalar;

/// One row or column of a matrix, holding only its stored entries.
///
/// Returned by `CsrMatrix::row`, `CscMatrix::col` and `CompressedMatrix::diagonal`, which
/// include any stored zeros.
#[derive(Clone, Debug, PartialEq)]
pub struct SparseVector<T = f64> {
    pub(crate) len: usize,
    pub(crate) indices: Vec<usize>,
    pub(crate) values: Vec<T>,
}

impl<T: Scalar> SparseVector<T> {
    /// The length of the vector, counting implicit zeros.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the vector has length zero.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The number of stored (non-zero) entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// The positions of the stored entries, in increasing order.
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    /// The stored values, parallel to `indices`.
    pub fn values(&self) -> &[T] {
        &self.values
    }

    /// The element at `index`, or `None` if it is an implicit zero or out of range.
    pub fn get(&self, index: usize) -> Option<T> {
        self.indices
            .binary_search(&index)
            .ok()
            .map(|pos| self.values[pos])
    }

    /// Iterates over the stored entries as `(index, value)`.
    pub fn iter(&self) -> impl Iterator<Item = (usize, T)> + '_ {
        self.indices
            .iter()
            .copied()
            .zip(self.values.iter().copied())
    }

    /// Expands to a dense vector of length `len`.
    pub fn to_dense(&self) -> Vec<T> {
        let mut dense = vec![T::zero(); self.len];
        for (index, value) in self.iter() {
            dense[index] = value;
        }
        dense
    }
}

/// Iterator over a matrix's stored entries as `(row, col, value)`, created by
/// `SparseMatrix::iter`.
///
/// Entries come in row-major order: by row, then by column within each row.
pub struct Entries<T: Scalar = f64> {
    entries: std::vec::IntoIter<(usize, usize, T)>,
}

impl<T: Scalar> Iterator for Entries<T> {
    type Item = (usize, usize, T);

    fn next(&mut self) -> Option<Self::Item> {
        self.entries.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.entries.size_hint()
    }
}

impl<T: Scalar> ExactSizeIterator for Entries<T> {}

impl<T: Scalar> SparseMatrix<T> {
    /// Iterates over the stored entries as `(row, col, value)` triplets in row-major order.
    ///
    /// Collects and sorts the entries when called, so it costs O(nnz log nnz), never
    /// O(rows * cols). For repeated traversal, convert with `to_csr` and iterate that.
    pub fn iter(&self) -> Entries<T> {
        let mut entries: Vec<(usize, usize, T)> = self
            .values
            .iter()
            .map(|(&(row, col), &value)| (row, col, value))
            .collect();
        entries.sort_unstable_by_key(|&(row, col, _)| (row, col));
        Entries {
            entries: entries.into_iter(),
        }
    }
}

/// Reads the element at `(row, col)`, yielding zero for entries that are not stored.
///
/// # Panics
///
/// Panics if the position is outside the matrix. Use `peek_at` for a checked read.
impl<T: Scalar> Index<(usize, usize)> for SparseMatrix<T> {
    type Output = T;

    fn index(&self, (row, col): (usize, usize)) -> &T {
        assert!(
            row < self.shape.0 && col < self.shape.1,
            "index ({row}, {col}) out of bounds for shape {:?}",
            self.shape
        );
        self.values.get(&(row, col)).unwrap_or(&self.zero)
    }
}
