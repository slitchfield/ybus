//! Read access to rows, columns and the diagonal, and extracting or reordering parts of a
//! matrix.
//!
//! Each access lives on the orientation where it's a slice: rows on `CsrMatrix`, columns on
//! `CscMatrix`. Selections and permutations work in either orientation and keep it. Stored
//! zeros are carried along like any other entry.

use std::ops::Range;

use super::{CompressedMatrix, CscMatrix, CsrMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::sparse_matrix::SparseVector;
use crate::Scalar;

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// Major index `major` as a sparse vector of length `len`, stored zeros included.
    fn segment(&self, major: usize, len: usize) -> SparseVector<T> {
        let range = self.ptr[major]..self.ptr[major + 1];
        SparseVector {
            len,
            indices: self.idx[range.clone()].to_vec(),
            values: self.data[range].to_vec(),
        }
    }

    /// Each major index in turn as dense values of length `len`.
    fn dense_segments(&self, len: usize) -> impl Iterator<Item = Vec<T>> + '_ {
        (0..self.major_dim()).map(move |major| {
            let mut dense = vec![T::zero(); len];
            for p in self.ptr[major]..self.ptr[major + 1] {
                dense[self.idx[p]] = self.data[p];
            }
            dense
        })
    }

    /// The main diagonal as a sparse vector of length `min(rows, cols)`, holding the diagonal
    /// positions in the pattern (stored zeros included). Each is found by binary search.
    pub fn diagonal(&self) -> SparseVector<T> {
        let len = self.shape.0.min(self.shape.1);
        let (indices, values) = (0..len)
            .filter_map(|i| self.stored_value(i, i).map(|value| (i, value)))
            .unzip();
        SparseVector {
            len,
            indices,
            values,
        }
    }

    /// Extracts the block `rows x cols` as a new matrix of shape `(rows.len(), cols.len())`
    /// in the same orientation, with indices shifted so the block starts at `(0, 0)`.
    ///
    /// Returns `SparseError::InvalidRange` if either range is reversed or ends past the
    /// matrix.
    pub fn submatrix(
        &self,
        rows: Range<usize>,
        cols: Range<usize>,
    ) -> Result<CompressedMatrix<T, O>, SparseError> {
        for (range, len) in [(&rows, self.shape.0), (&cols, self.shape.1)] {
            if range.start > range.end || range.end > len {
                return Err(SparseError::InvalidRange {
                    start: range.start,
                    end: range.end,
                    len,
                });
            }
        }
        let (majors, minors) = if O::ROW_MAJOR {
            (rows.clone(), cols.clone())
        } else {
            (cols.clone(), rows.clone())
        };
        let mut ptr = Vec::with_capacity(majors.len() + 1);
        ptr.push(0);
        let (mut idx, mut data) = (vec![], vec![]);
        for major in majors {
            let segment = self.ptr[major]..self.ptr[major + 1];
            let indices = &self.idx[segment.clone()];
            // Indices are sorted, so the block's part of the segment is contiguous
            let start = segment.start + indices.partition_point(|&minor| minor < minors.start);
            let end = segment.start + indices.partition_point(|&minor| minor < minors.end);
            idx.extend(
                self.idx[start..end]
                    .iter()
                    .map(|&minor| minor - minors.start),
            );
            data.extend_from_slice(&self.data[start..end]);
            ptr.push(idx.len());
        }
        Ok(Self::from_parts_unchecked(
            (rows.len(), cols.len()),
            ptr,
            idx,
            data,
        ))
    }

    /// Extracts the rows and columns named by two index lists, as a new matrix of shape
    /// `(rows.len(), cols.len())` in the same orientation: entry `(i, j)` of the result is
    /// entry `(rows[i], cols[j])` of `self`.
    ///
    /// Unlike `submatrix`, the indices need not be contiguous or sorted, so this can pull out
    /// any subset (such as the PV and PQ buses of a power system) in any order. An index may
    /// appear more than once, which repeats that row or column, as multiplying by a selection
    /// matrix would. Costs O(nnz) plus sorting each result row (CSR) or column (CSC). Returns
    /// `SparseError::OutOfRange` for the first index that is not a row or column of the
    /// matrix, checking `rows` first.
    pub fn select(
        &self,
        rows: &[usize],
        cols: &[usize],
    ) -> Result<CompressedMatrix<T, O>, SparseError> {
        for (indices, len) in [(rows, self.shape.0), (cols, self.shape.1)] {
            if let Some(&index) = indices.iter().find(|&&index| index >= len) {
                return Err(SparseError::OutOfRange { index, len });
            }
        }
        let (majors, minors) = if O::ROW_MAJOR {
            (rows, cols)
        } else {
            (cols, rows)
        };
        // For each original minor index, the result positions it appears at
        let mut positions: Vec<Vec<usize>> = vec![vec![]; self.minor_dim()];
        for (position, &minor) in minors.iter().enumerate() {
            positions[minor].push(position);
        }
        let mut ptr = Vec::with_capacity(majors.len() + 1);
        ptr.push(0);
        let (mut idx, mut data) = (vec![], vec![]);
        let mut segment: Vec<(usize, T)> = vec![];
        for &major in majors {
            segment.clear();
            for p in self.ptr[major]..self.ptr[major + 1] {
                for &position in &positions[self.idx[p]] {
                    segment.push((position, self.data[p]));
                }
            }
            segment.sort_unstable_by_key(|&(position, _)| position);
            for &(position, value) in &segment {
                idx.push(position);
                data.push(value);
            }
            ptr.push(idx.len());
        }
        Ok(Self::from_parts_unchecked(
            (rows.len(), cols.len()),
            ptr,
            idx,
            data,
        ))
    }

    /// Reorders rows and columns: entry `(i, j)` of the result is entry
    /// `(row_perm[i], col_perm[j])` of `self`. As matrices, this is `P A Qᵀ` where `P` and `Q`
    /// are the permutation matrices whose row `i` selects `row_perm[i]` and `col_perm[i]`.
    ///
    /// The same result as `select(row_perm, col_perm)`, but both lists are checked to be true
    /// permutations, so no row or column is lost or repeated. Returns
    /// `SparseError::InvalidPermutation` if either list has the wrong length, an index out of
    /// range, or a repeated index.
    pub fn permute(
        &self,
        row_perm: &[usize],
        col_perm: &[usize],
    ) -> Result<CompressedMatrix<T, O>, SparseError> {
        check_permutation(row_perm, self.shape.0)?;
        check_permutation(col_perm, self.shape.1)?;
        self.select(row_perm, col_perm)
    }

    /// Applies the same permutation to rows and columns, `P A Pᵀ`, as fill-reducing orderings
    /// do. This keeps the diagonal on the diagonal, so it preserves symmetry, the trace and
    /// the set of diagonal values.
    ///
    /// Returns `SparseError::NotSquare` for a non-square matrix, and
    /// `SparseError::InvalidPermutation` as for `permute`.
    pub fn permute_symmetric(&self, perm: &[usize]) -> Result<CompressedMatrix<T, O>, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        self.permute(perm, perm)
    }
}

impl<T: Scalar> CsrMatrix<T> {
    /// Row `i` as a sparse vector of length `cols`, copied from its slice of the arrays,
    /// stored zeros included.
    ///
    /// Returns `SparseError::OutOfRange` if `i` is not a row of the matrix.
    pub fn row(&self, i: usize) -> Result<SparseVector<T>, SparseError> {
        if i >= self.shape.0 {
            return Err(SparseError::OutOfRange {
                index: i,
                len: self.shape.0,
            });
        }
        Ok(self.segment(i, self.shape.1))
    }

    /// Iterates over the rows, yielding each as a dense `Vec<T>` of length `cols`.
    pub fn row_iter(&self) -> impl Iterator<Item = Vec<T>> + '_ {
        self.dense_segments(self.shape.1)
    }
}

impl<T: Scalar> CscMatrix<T> {
    /// Column `j` as a sparse vector of length `rows`, copied from its slice of the arrays,
    /// stored zeros included.
    ///
    /// Returns `SparseError::OutOfRange` if `j` is not a column of the matrix.
    pub fn col(&self, j: usize) -> Result<SparseVector<T>, SparseError> {
        if j >= self.shape.1 {
            return Err(SparseError::OutOfRange {
                index: j,
                len: self.shape.1,
            });
        }
        Ok(self.segment(j, self.shape.0))
    }

    /// Iterates over the columns, yielding each as a dense `Vec<T>` of length `rows`.
    pub fn col_iter(&self) -> impl Iterator<Item = Vec<T>> + '_ {
        self.dense_segments(self.shape.0)
    }
}

/// Checks that `perm` is a permutation of `0..len`.
pub(crate) fn check_permutation(perm: &[usize], len: usize) -> Result<(), SparseError> {
    let invalid = |reason| Err(SparseError::InvalidPermutation { reason });
    if perm.len() != len {
        return invalid("length must equal the dimension being permuted");
    }
    let mut seen = vec![false; len];
    for &index in perm {
        if index >= len {
            return invalid("index out of range");
        }
        if seen[index] {
            return invalid("index repeated");
        }
        seen[index] = true;
    }
    Ok(())
}
