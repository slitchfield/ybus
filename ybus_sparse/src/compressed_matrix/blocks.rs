//! Building larger matrices from smaller ones: Kronecker product and stacking, for operands of
//! one orientation.

use std::borrow::Borrow;

use super::{CompressedMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// The Kronecker product `A ⊗ B`: an `(m·p) x (n·q)` matrix, for `A` of shape `m x n` and
    /// `B` of shape `p x q`, made of a copy of `B` scaled by each entry of `A`. The pattern is
    /// every pairing of stored entries, zero products included.
    ///
    /// Costs `A.num_stored() * B.num_stored()`. Returns `SparseError::ShapeOverflow` if a
    /// dimension of the result does not fit in a `usize`.
    pub fn kron(&self, other: &CompressedMatrix<T, O>) -> Result<Self, SparseError> {
        let rows = self
            .shape
            .0
            .checked_mul(other.shape.0)
            .ok_or(SparseError::ShapeOverflow)?;
        let cols = self
            .shape
            .1
            .checked_mul(other.shape.1)
            .ok_or(SparseError::ShapeOverflow)?;
        let other_minor = other.minor_dim();
        let mut ptr = Vec::with_capacity(self.major_dim() * other.major_dim() + 1);
        ptr.push(0);
        let capacity = self.idx.len() * other.idx.len();
        let (mut idx, mut data) = (Vec::with_capacity(capacity), Vec::with_capacity(capacity));
        // Result major index a * other_major + b, minor index a' * other_minor + b'; visiting
        //   both operands' entries in order keeps each result segment sorted
        for a_major in 0..self.major_dim() {
            for b_major in 0..other.major_dim() {
                for p in self.ptr[a_major]..self.ptr[a_major + 1] {
                    for q in other.ptr[b_major]..other.ptr[b_major + 1] {
                        idx.push(self.idx[p] * other_minor + other.idx[q]);
                        data.push(self.data[p] * other.data[q]);
                    }
                }
                ptr.push(idx.len());
            }
        }
        Ok(Self::from_parts_unchecked((rows, cols), ptr, idx, data))
    }

    /// Joins matrices side by side into one matrix with the same number of rows.
    ///
    /// Accepts anything that iterates over matrices or references to them. No blocks gives a
    /// `0 x 0` matrix. Returns `SparseError::ShapeMismatch` (naming the first block and the
    /// first one that differs) unless every block has the same number of rows, and
    /// `SparseError::ShapeOverflow` if the columns do not fit in a `usize`. For CSC this
    /// concatenates the arrays; for CSR it interleaves each row's segments.
    pub fn hstack<I>(blocks: I) -> Result<Self, SparseError>
    where
        I: IntoIterator,
        I::Item: Borrow<CompressedMatrix<T, O>>,
    {
        Self::stack(blocks, true)
    }

    /// Joins matrices top to bottom into one matrix with the same number of columns.
    ///
    /// The vertical counterpart of `hstack`, with the same errors for column counts. For CSR
    /// this concatenates the arrays; for CSC it interleaves each column's segments.
    pub fn vstack<I>(blocks: I) -> Result<Self, SparseError>
    where
        I: IntoIterator,
        I::Item: Borrow<CompressedMatrix<T, O>>,
    {
        Self::stack(blocks, false)
    }

    /// Shared implementation of `hstack` (`horizontal`) and `vstack`.
    fn stack<I>(blocks: I, horizontal: bool) -> Result<Self, SparseError>
    where
        I: IntoIterator,
        I::Item: Borrow<CompressedMatrix<T, O>>,
    {
        let blocks: Vec<I::Item> = blocks.into_iter().collect();
        let blocks: Vec<&CompressedMatrix<T, O>> = blocks.iter().map(Borrow::borrow).collect();
        let Some(first) = blocks.first() else {
            return Ok(Self::default());
        };
        // Lengths along the axis being stacked, which must fit in a usize, and across it
        let along = |block: &CompressedMatrix<T, O>| {
            if horizontal {
                block.shape.1
            } else {
                block.shape.0
            }
        };
        let across = |block: &CompressedMatrix<T, O>| {
            if horizontal {
                block.shape.0
            } else {
                block.shape.1
            }
        };
        let mut stacked_len = 0usize;
        for block in &blocks {
            if across(block) != across(first) {
                return Err(SparseError::ShapeMismatch {
                    left: first.shape,
                    right: block.shape,
                });
            }
            stacked_len = stacked_len
                .checked_add(along(block))
                .ok_or(SparseError::ShapeOverflow)?;
        }
        let shape = if horizontal {
            (first.shape.0, stacked_len)
        } else {
            (stacked_len, first.shape.1)
        };

        let mut ptr = vec![0];
        let (mut idx, mut data) = (vec![], vec![]);
        // Stacking along the compressed axis appends each block's segments in turn; across it,
        //   each segment is the blocks' segments side by side, minor indices offset
        if horizontal != O::ROW_MAJOR {
            for block in &blocks {
                for major in 0..block.major_dim() {
                    let range = block.ptr[major]..block.ptr[major + 1];
                    idx.extend_from_slice(&block.idx[range.clone()]);
                    data.extend_from_slice(&block.data[range]);
                    ptr.push(idx.len());
                }
            }
        } else {
            for major in 0..first.major_dim() {
                let mut offset = 0;
                for block in &blocks {
                    let range = block.ptr[major]..block.ptr[major + 1];
                    idx.extend(block.idx[range.clone()].iter().map(|&minor| minor + offset));
                    data.extend_from_slice(&block.data[range]);
                    offset += block.minor_dim();
                }
                ptr.push(idx.len());
            }
        }
        Ok(Self::from_parts_unchecked(shape, ptr, idx, data))
    }
}
