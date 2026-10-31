//! Read access to a matrix's stored entries, rows and columns.

use std::borrow::Cow;
use std::ops::{Index, Range};

use super::{Compressed, SparseMatrix};
use crate::sparse_error::SparseError;
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

/// Iterator over a matrix's columns, created by `SparseMatrix::col_iter`.
///
/// Each item is one dense column of length `shape.0`.
pub struct ColIterator<'a, T: Scalar = f64> {
    shape: (usize, usize),
    /// Borrowed from the matrix's CSC cache when fresh, otherwise built for this iterator.
    csc: Cow<'a, Compressed<T>>,
    col_iter_idx: usize,
}

impl<T: Scalar> Iterator for ColIterator<'_, T> {
    /// One dense column.
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.col_iter_idx >= self.shape.1 {
            return None;
        }
        let start = self.csc.ptr[self.col_iter_idx];
        let end = self.csc.ptr[self.col_iter_idx + 1];

        let mut column = vec![T::zero(); self.shape.0];
        for (row, val) in std::iter::zip(&self.csc.idx[start..end], &self.csc.data[start..end]) {
            column[*row] = *val;
        }

        self.col_iter_idx += 1;
        Some(column)
    }
}

/// Iterator over a matrix's stored entries as `(row, col, value)`, created by
/// `SparseMatrix::iter`.
///
/// Entries come in row-major order: by row, then by column within each row.
pub struct Entries<'a, T: Scalar = f64> {
    /// Borrowed from the matrix's cache when fresh, otherwise built for this iterator.
    csr: Cow<'a, Compressed<T>>,
    row: usize,
    pos: usize,
}

impl<T: Scalar> Iterator for Entries<'_, T> {
    type Item = (usize, usize, T);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.csr.idx.len() {
            return None;
        }
        // Skip past rows with no entries left; ptr's last element equals the entry count,
        //   so this stops before running off the end.
        while self.csr.ptr[self.row + 1] <= self.pos {
            self.row += 1;
        }
        let item = (self.row, self.csr.idx[self.pos], self.csr.data[self.pos]);
        self.pos += 1;
        Some(item)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.csr.idx.len() - self.pos;
        (remaining, Some(remaining))
    }
}

impl<T: Scalar> ExactSizeIterator for Entries<'_, T> {}

impl<T: Scalar> SparseMatrix<T> {
    /// Iterates over the stored entries as `(row, col, value)` triplets in row-major order.
    ///
    /// Unlike `row_iter`, this only visits stored entries, so it costs time proportional to
    /// `num_nonzero` rather than `rows * cols`. Call `explicitly_compress` first when
    /// iterating repeatedly; otherwise each iterator rebuilds the CSR arrays itself.
    pub fn iter(&self) -> Entries<'_, T> {
        Entries {
            csr: self.csr_cow(),
            row: 0,
            pos: 0,
        }
    }
}

impl<T: Scalar> SparseMatrix<T> {
    /// Row `i` as a sparse vector; forwards to `CsrMatrix::row`.
    pub fn row(&self, i: usize) -> Result<SparseVector<T>, SparseError> {
        self.to_csr().row(i)
    }

    /// Column `j` as a sparse vector; forwards to `CscMatrix::col`.
    pub fn col(&self, j: usize) -> Result<SparseVector<T>, SparseError> {
        self.to_csc().col(j)
    }

    /// Iterates over the columns, yielding each as a dense `Vec<T>`.
    pub fn col_iter(&self) -> ColIterator<'_, T> {
        ColIterator {
            shape: self.shape,
            csc: self.csc_cow(),
            col_iter_idx: 0,
        }
    }

    /// The main diagonal; forwards to `CompressedMatrix::diagonal`.
    pub fn diagonal(&self) -> SparseVector<T> {
        self.to_csr().diagonal()
    }

    /// A contiguous block; forwards to `CompressedMatrix::submatrix`.
    pub fn submatrix(
        &self,
        rows: Range<usize>,
        cols: Range<usize>,
    ) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .submatrix(rows, cols)
            .map(|block| block.to_sparse())
    }

    /// Selected rows and columns; forwards to `CompressedMatrix::select`.
    pub fn select(&self, rows: &[usize], cols: &[usize]) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .select(rows, cols)
            .map(|selected| selected.to_sparse())
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
