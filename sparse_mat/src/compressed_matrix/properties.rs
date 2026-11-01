//! Properties and reductions of a matrix: shape predicates, sums and norms.

use num_traits::Float;

use super::{CompressedMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// The fraction of positions in the pattern: `num_stored / (rows * cols)`, stored zeros
    /// included. An empty matrix (zero rows or columns) has density `0.0`.
    pub fn density(&self) -> f64 {
        let area = self.shape.0 as f64 * self.shape.1 as f64;
        if area == 0.0 {
            0.0
        } else {
            self.num_stored() as f64 / area
        }
    }

    /// The sum of all entries, added in storage order (row-major for CSR, column-major for
    /// CSC), so the result is deterministic for floating point types.
    pub fn sum(&self) -> T {
        self.data
            .iter()
            .fold(T::zero(), |total, &value| total + value)
    }

    /// The sum of the main diagonal, added in index order. Returns `SparseError::NotSquare` for
    /// a non-square matrix.
    pub fn trace(&self) -> Result<T, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        Ok((0..self.shape.0)
            .filter_map(|i| self.stored_value(i, i))
            .fold(T::zero(), |total, value| total + value))
    }

    /// The largest distance `|row − col|` of a stored entry from the diagonal, stored zeros
    /// included: 0 for a diagonal (or empty) matrix, 1 for a tridiagonal one. Orderings such as
    /// reverse Cuthill–McKee aim to make it small.
    pub fn bandwidth(&self) -> usize {
        self.iter()
            .map(|(row, col, _)| row.abs_diff(col))
            .max()
            .unwrap_or(0)
    }

    /// The size of the lower envelope: the sum over rows `i` of `i − f(i)`, where `f(i)` is the
    /// smallest column of a stored entry in row `i` at or left of the diagonal (`f(i) = i` if
    /// there is none). Stored zeros count. For a symmetric matrix, a Cholesky factor's fill stays
    /// inside this envelope, so it bounds the factor's entries below the diagonal.
    pub fn profile(&self) -> usize {
        let mut first = (0..self.shape.0).collect::<Vec<usize>>();
        for (row, col, _) in self.iter() {
            if col < first[row] {
                first[row] = col;
            }
        }
        first.iter().enumerate().map(|(row, &col)| row - col).sum()
    }

    /// Whether the matrix equals its transpose, mathematically: stored zeros don't matter, and
    /// floating point entries must match exactly. Non-square matrices are never symmetric.
    /// Costs O(nnz), comparing with a reoriented copy.
    pub fn is_symmetric(&self) -> bool {
        self.is_square() && *self == self.clone().transpose()
    }

    /// The sums of absolute values along each minor index, accumulated in storage order.
    fn minor_sums_abs(&self) -> Vec<T>
    where
        T: Float,
    {
        let mut sums = vec![T::zero(); self.minor_dim()];
        for (&minor, value) in self.idx.iter().zip(&self.data) {
            sums[minor] = sums[minor] + value.abs();
        }
        sums
    }

    /// The sums of absolute values along each major index.
    fn major_sums_abs(&self) -> Vec<T>
    where
        T: Float,
    {
        self.ptr
            .windows(2)
            .map(|bounds| {
                self.data[bounds[0]..bounds[1]]
                    .iter()
                    .fold(T::zero(), |total, value| total + value.abs())
            })
            .collect()
    }
}

impl<T: Scalar + Float, O: Orientation> CompressedMatrix<T, O> {
    /// The Frobenius norm: the square root of the sum of squared magnitudes.
    pub fn norm_frobenius(&self) -> T {
        self.data
            .iter()
            .fold(T::zero(), |total, &value| total + value * value)
            .sqrt()
    }

    /// The 1-norm: the largest column sum of absolute values. Each column sum adds its entries
    /// in row order in either orientation.
    pub fn norm_1(&self) -> T {
        let sums = if O::ROW_MAJOR {
            self.minor_sums_abs()
        } else {
            self.major_sums_abs()
        };
        sums.into_iter()
            .fold(T::zero(), |largest, sum| largest.max(sum))
    }

    /// The infinity norm: the largest row sum of absolute values. Each row sum adds its entries
    /// in column order in either orientation.
    pub fn norm_inf(&self) -> T {
        let sums = if O::ROW_MAJOR {
            self.major_sums_abs()
        } else {
            self.minor_sums_abs()
        };
        sums.into_iter()
            .fold(T::zero(), |largest, sum| largest.max(sum))
    }
}
