//! Parallel products using `rayon` (feature `rayon`).
//!
//! Each result entry is computed by exactly one thread with the same accumulation order as
//! the sequential version, so `par_*` results are bit-for-bit equal to their sequential
//! counterparts, not merely close. Each product lives on the orientation where it splits into
//! independent dot products: `par_mul_vec` on CSR rows, `par_vec_mul` on CSC columns.

use std::collections::HashMap;

use rayon::prelude::*;

use super::product::{assemble_rows, gustavson_row};
use super::{CompressedMatrix, CscMatrix, CsrMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// Dot products of `x` with each major index's entries, one per thread.
fn par_dot<T: Scalar + Send + Sync>(ptr: &[usize], idx: &[usize], data: &[T], x: &[T]) -> Vec<T> {
    ptr.par_windows(2)
        .map(|bounds| {
            let range = bounds[0]..bounds[1];
            idx[range.clone()]
                .iter()
                .zip(&data[range])
                .fold(T::zero(), |sum, (&minor, &value)| sum + value * x[minor])
        })
        .collect()
}

impl<T: Scalar + Send + Sync> CsrMatrix<T> {
    /// Parallel version of `mul_vec`: entries of the result, one per row, are computed on
    /// separate threads.
    pub fn par_mul_vec(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        if x.len() != self.shape.1 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: (x.len(), 1),
            });
        }
        Ok(par_dot(&self.ptr, &self.idx, &self.data, x))
    }
}

impl<T: Scalar + Send + Sync> CscMatrix<T> {
    /// Parallel version of `vec_mul`: entries of the result, one per column, are computed on
    /// separate threads.
    pub fn par_vec_mul(&self, x: &[T]) -> Result<Vec<T>, SparseError> {
        if x.len() != self.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: (1, x.len()),
                right: self.shape,
            });
        }
        Ok(par_dot(&self.ptr, &self.idx, &self.data, x))
    }
}

impl<T: Scalar + Send + Sync, O: Orientation> CompressedMatrix<T, O> {
    /// Parallel version of `matmul`: each major index of the product (a row for CSR, a column
    /// for CSC) is built on one thread, with the same stored zeros.
    pub fn par_matmul(
        &self,
        other: &CompressedMatrix<T, O>,
    ) -> Result<CompressedMatrix<T, O>, SparseError> {
        if self.shape.1 != other.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: other.shape,
            });
        }
        let (left, right) = if O::ROW_MAJOR {
            (self, other)
        } else {
            (other, self)
        };
        let rows: Vec<Vec<(usize, T)>> = left
            .ptr
            .par_windows(2)
            .map_init(HashMap::new, |accumulator, bounds| {
                let range = bounds[0]..bounds[1];
                gustavson_row(
                    &left.idx[range.clone()],
                    &left.data[range],
                    (&right.ptr, &right.idx, &right.data),
                    accumulator,
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
