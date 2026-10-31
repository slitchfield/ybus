//! Matrix products: matrix-vector, vector-matrix, matrix-dense and matrix-matrix (SpGEMM).
//!
//! Every product walks the compressed arrays, whose entries come in sorted order, so results
//! are deterministic even for floating point sums. Each result entry adds its terms in
//! increasing order of the inner index in either orientation, so a `CsrMatrix` and a
//! `CscMatrix` holding the same matrix give bit-for-bit equal results.

use std::collections::HashMap;
use std::ops::Mul;

use super::{CompressedMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// `y[i] = Σ data[p] · x[idx[p]]` over each major index `i`'s entries: `A x` for CSR arrays,
/// `x A` for CSC arrays. The caller has checked that `x` covers the minor indices.
pub(crate) fn dot_kernel<T: Scalar>(ptr: &[usize], idx: &[usize], data: &[T], x: &[T]) -> Vec<T> {
    ptr.windows(2)
        .map(|bounds| {
            let range = bounds[0]..bounds[1];
            idx[range.clone()]
                .iter()
                .zip(&data[range])
                .fold(T::zero(), |sum, (&minor, &value)| sum + value * x[minor])
        })
        .collect()
}

/// `y[idx[p]] += data[p] · x[i]` over each major index `i`'s entries, into a result of length
/// `len`: `x A` for CSR arrays, `A x` for CSC arrays. The caller has checked that `x` covers the
/// major indices.
fn scatter_kernel<T: Scalar>(
    ptr: &[usize],
    idx: &[usize],
    data: &[T],
    x: &[T],
    len: usize,
) -> Vec<T> {
    let mut y = vec![T::zero(); len];
    for (major, bounds) in ptr.windows(2).enumerate() {
        for p in bounds[0]..bounds[1] {
            y[idx[p]] = y[idx[p]] + data[p] * x[major];
        }
    }
    y
}

/// One row of a Gustavson product of row-compressed `a` and `b`: the sum of `b`'s rows selected
/// and scaled by the entries `a_idx`/`a_data` of one row of `a`, as `(column, value)` pairs in
/// column order. Every reachable column is kept, even where the terms cancel to zero.
/// `accumulator` is scratch space, reused across rows.
pub(crate) fn gustavson_row<T: Scalar>(
    a_idx: &[usize],
    a_data: &[T],
    b: (&[usize], &[usize], &[T]),
    accumulator: &mut HashMap<usize, T>,
) -> Vec<(usize, T)> {
    let (b_ptr, b_idx, b_data) = b;
    accumulator.clear();
    for (&k, &a_value) in a_idx.iter().zip(a_data) {
        for p in b_ptr[k]..b_ptr[k + 1] {
            let entry = accumulator.entry(b_idx[p]).or_insert(T::zero());
            *entry = *entry + a_value * b_data[p];
        }
    }
    let mut row: Vec<(usize, T)> = accumulator.drain().collect();
    row.sort_unstable_by_key(|&(col, _)| col);
    row
}

/// Assembles rows of `(minor, value)` pairs, one per major index, into compressed arrays.
pub(crate) fn assemble_rows<T>(rows: Vec<Vec<(usize, T)>>) -> (Vec<usize>, Vec<usize>, Vec<T>) {
    let mut ptr = Vec::with_capacity(rows.len() + 1);
    ptr.push(0);
    let mut idx = Vec::new();
    let mut data = Vec::new();
    for row in rows {
        for (minor, value) in row {
            idx.push(minor);
            data.push(value);
        }
        ptr.push(idx.len());
    }
    (ptr, idx, data)
}

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// The matrix-vector product `A x` for a column vector `x`, a vector of length `rows`.
    /// CSR computes each entry as a dot product with a row; CSC adds each column in, scaled.
    ///
    /// Returns `SparseError::ShapeMismatch` (with the vector treated as `len x 1`) unless
    /// `x.len()` equals the column count.
    pub fn mul_vec(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        if x.len() != self.shape.1 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: (x.len(), 1),
            });
        }
        Ok(if O::ROW_MAJOR {
            dot_kernel(&self.ptr, &self.idx, &self.data, x)
        } else {
            scatter_kernel(&self.ptr, &self.idx, &self.data, x, self.shape.0)
        })
    }

    /// The vector-matrix product `x A` for a row vector `x`, a vector of length `cols`: `Aᵀ x`
    /// without forming the transpose. CSC computes each entry as a dot product with a column;
    /// CSR adds each row in, scaled.
    ///
    /// Returns `SparseError::ShapeMismatch` (with the vector treated as `1 x len`) unless
    /// `x.len()` equals the row count.
    pub fn vec_mul(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        if x.len() != self.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: (1, x.len()),
                right: self.shape,
            });
        }
        Ok(if O::ROW_MAJOR {
            scatter_kernel(&self.ptr, &self.idx, &self.data, x, self.shape.1)
        } else {
            dot_kernel(&self.ptr, &self.idx, &self.data, x)
        })
    }

    /// The product `A B` with a dense matrix `B` given as rows, as dense rows.
    ///
    /// Returns `SparseError::RaggedRows` if `B`'s rows differ in length, and
    /// `SparseError::ShapeMismatch` unless `B` has as many rows as `A` has columns.
    pub fn mul_dense(&self, b: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        let b_cols = b.first().map_or(0, Vec::len);
        if let Some((row, bad)) = b.iter().enumerate().find(|(_, row)| row.len() != b_cols) {
            return Err(SparseError::RaggedRows {
                row,
                expected: b_cols,
                found: bad.len(),
            });
        }
        if b.len() != self.shape.1 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: (b.len(), b_cols),
            });
        }

        let mut out = vec![vec![T::zero(); b_cols]; self.shape.0];
        for major in 0..self.major_dim() {
            for p in self.ptr[major]..self.ptr[major + 1] {
                let (row, k) = Self::row_col(major, self.idx[p]);
                let value = self.data[p];
                for (out, &b_kj) in out[row].iter_mut().zip(&b[k]) {
                    *out = *out + value * b_kj;
                }
            }
        }
        Ok(out)
    }

    /// The sparse matrix product `A B` (SpGEMM), in the operands' orientation, using
    /// Gustavson's algorithm: for CSR, row `i` of the result is the sum of the rows of `B`
    /// selected and scaled by row `i` of `A`; CSC works by columns, through `(AB)ᵀ = BᵀAᵀ`. Work
    /// is proportional to the number of multiplications between stored entries.
    ///
    /// The result stores every position the operands' patterns reach, even where the terms
    /// cancel to exactly zero, so its pattern doesn't depend on the values. Returns
    /// `SparseError::ShapeMismatch` unless `A`'s column count equals `B`'s row count; `*`
    /// panics instead.
    pub fn matmul(
        &self,
        other: &CompressedMatrix<T, O>,
    ) -> Result<CompressedMatrix<T, O>, SparseError> {
        if self.shape.1 != other.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: other.shape,
            });
        }
        // Row-compressed arrays of the left and right factors: A and B for CSR, and Bᵀ and Aᵀ
        //   (whose CSR arrays are B's and A's CSC arrays) for CSC
        let (left, right) = if O::ROW_MAJOR {
            (self, other)
        } else {
            (other, self)
        };
        let mut accumulator = HashMap::new();
        let rows: Vec<Vec<(usize, T)>> = left
            .ptr
            .windows(2)
            .map(|bounds| {
                let range = bounds[0]..bounds[1];
                gustavson_row(
                    &left.idx[range.clone()],
                    &left.data[range],
                    (&right.ptr, &right.idx, &right.data),
                    &mut accumulator,
                )
            })
            .collect();
        let (ptr, idx, data) = assemble_rows(rows);
        Ok(CompressedMatrix::from_parts_unchecked(
            (self.shape.0, other.shape.1),
            ptr,
            idx,
            data,
        ))
    }
}

/// Matrix product; see `CompressedMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar, O: Orientation> Mul<&CompressedMatrix<T, O>> for &CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn mul(self, other: &CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
        self.matmul(other).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Matrix product; see `CompressedMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar, O: Orientation> Mul<CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn mul(self, other: CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
        &self * &other
    }
}

/// Matrix product; see `CompressedMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar, O: Orientation> Mul<&CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn mul(self, other: &CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
        &self * other
    }
}

/// Matrix product; see `CompressedMatrix::matmul`.
///
/// # Panics
///
/// Panics unless the left matrix has as many columns as the right has rows.
impl<T: Scalar, O: Orientation> Mul<CompressedMatrix<T, O>> for &CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn mul(self, other: CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
        self * &other
    }
}

/// Matrix-vector product `A x`; see `CompressedMatrix::mul_vec`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's column count.
impl<T: Scalar, O: Orientation> Mul<&[T]> for &CompressedMatrix<T, O> {
    type Output = Vec<T>;

    fn mul(self, x: &[T]) -> Vec<T> {
        self.mul_vec(x).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Matrix-vector product `A x`; see `CompressedMatrix::mul_vec`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's column count.
impl<T: Scalar, O: Orientation> Mul<&Vec<T>> for &CompressedMatrix<T, O> {
    type Output = Vec<T>;

    fn mul(self, x: &Vec<T>) -> Vec<T> {
        self.mul_vec(x).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Vector-matrix product `x A`; see `CompressedMatrix::vec_mul`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's row count.
impl<T: Scalar, O: Orientation> Mul<&CompressedMatrix<T, O>> for &[T] {
    type Output = Vec<T>;

    fn mul(self, matrix: &CompressedMatrix<T, O>) -> Vec<T> {
        matrix.vec_mul(self).unwrap_or_else(|err| panic!("{err}"))
    }
}

/// Vector-matrix product `x A`; see `CompressedMatrix::vec_mul`.
///
/// # Panics
///
/// Panics unless the vector's length equals the matrix's row count.
impl<T: Scalar, O: Orientation> Mul<&CompressedMatrix<T, O>> for &Vec<T> {
    type Output = Vec<T>;

    fn mul(self, matrix: &CompressedMatrix<T, O>) -> Vec<T> {
        matrix.vec_mul(self).unwrap_or_else(|err| panic!("{err}"))
    }
}
