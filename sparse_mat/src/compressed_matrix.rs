//! Compressed sparse matrices with a fixed pattern, for numerical work.
//!
//! A [`CompressedMatrix`] stores its entries compressed along one axis: by rows in a
//! [`CsrMatrix`] (compressed sparse row) or by columns in a [`CscMatrix`] (compressed sparse
//! column). The orientation is part of the type, so an algorithm that needs one orientation
//! can't be handed the other, and converting between them is always an explicit call. Build a
//! matrix with [`SparseMatrix`], which supports arbitrary edits, and convert it with
//! [`SparseMatrix::to_csr`] or [`SparseMatrix::to_csc`], or build one directly from triplets or
//! raw arrays.
//!
//! The pattern (which positions are stored) is fixed when the matrix is built. Values can
//! change in place through [`CompressedMatrix::values_mut`], but entries can't be added or
//! removed, and a stored entry keeps its place even when its value is zero. That is what lets a
//! matrix be refilled with new values, or refactored, without its structure changing.
//! [`CompressedMatrix::num_stored`] counts the pattern and [`CompressedMatrix::num_nonzero`]
//! the entries whose value isn't zero.
//!
//! Equality is mathematical: two matrices are equal when they have the same shape and the same
//! nonzero entries, whatever their orientation and stored zeros.

use std::fmt;
use std::marker::PhantomData;
use std::sync::OnceLock;

use crate::sparse_error::SparseError;
use crate::sparse_matrix::{write_dense_grid, SparseMatrix};
use crate::Scalar;

mod access;
mod assembly;
mod blocks;
mod cholesky;
mod lu;
mod ops;
mod ordering;
#[cfg(feature = "rayon")]
mod parallel;
pub(crate) mod product;
mod properties;
#[cfg(feature = "serde")]
mod serde_impl;
mod solvers;

pub use assembly::AssemblyMap;
pub use cholesky::{Cholesky, SymbolicCholesky};
pub use lu::{Lu, LuOptions, Refactor, SymbolicLu};
pub use ordering::Ordering;

mod sealed {
    pub trait Sealed {}
}

/// The orientation of a [`CompressedMatrix`]: which axis is compressed. Implemented only by
/// [`Csr`] and [`Csc`], which hold no data, so a matrix is `Send` and `Sync` when its values are.
pub trait Orientation:
    sealed::Sealed + Copy + fmt::Debug + Default + Send + Sync + 'static
{
    /// The other orientation, which transposing produces.
    type Transposed: Orientation<Transposed = Self>;
    /// Whether the rows are compressed (CSR) rather than the columns (CSC).
    const ROW_MAJOR: bool;
    /// `"CSR"` or `"CSC"`.
    const NAME: &'static str;
}

/// Compressed sparse row orientation: each row's entries are stored together, sorted by column.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Csr;

/// Compressed sparse column orientation: each column's entries are stored together, sorted by
/// row.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Csc;

impl sealed::Sealed for Csr {}
impl sealed::Sealed for Csc {}

impl Orientation for Csr {
    type Transposed = Csc;
    const ROW_MAJOR: bool = true;
    const NAME: &'static str = "CSR";
}

impl Orientation for Csc {
    type Transposed = Csr;
    const ROW_MAJOR: bool = false;
    const NAME: &'static str = "CSC";
}

/// A sparse matrix of `T` values compressed along the axis `O` chooses, with a fixed pattern.
/// Use the aliases [`CsrMatrix`] and [`CscMatrix`]; see the [module docs](self) for the model.
///
/// The compressed axis is the major axis: rows for CSR, columns for CSC. Entry `p` of major
/// index `i` lies in `pointers()[i]..pointers()[i + 1]`, at minor index `indices()[p]` with value
/// `values()[p]`, and the minor indices of each major index are strictly increasing.
#[derive(Clone)]
pub struct CompressedMatrix<T = f64, O: Orientation = Csr> {
    shape: (usize, usize),
    ptr: Vec<usize>,
    idx: Vec<usize>,
    data: Vec<T>,
    orientation: PhantomData<O>,
    /// The slot of each diagonal entry in `data`, computed on first use by
    /// `diagonal_positions`. The pattern is fixed, so it never goes stale.
    diagonal: OnceLock<Vec<Option<usize>>>,
}

/// A matrix in compressed sparse row form.
pub type CsrMatrix<T = f64> = CompressedMatrix<T, Csr>;

/// A matrix in compressed sparse column form.
pub type CscMatrix<T = f64> = CompressedMatrix<T, Csc>;

/// Allocates the `dim + 1` zeroed pointers of a compressed form.
///
/// Compressed storage costs O(dim) memory however sparse the matrix is, so a dimension near
/// `usize::MAX` cannot have one. Panics with a message naming the dimension, instead of
/// overflowing `dim + 1` (which wraps silently in release builds) or failing inside the
/// allocator.
pub(crate) fn zeroed_pointers(dim: usize) -> Vec<usize> {
    let len = dim.checked_add(1).unwrap_or_else(|| {
        panic!("compressed storage for dimension {dim} needs dim + 1 pointers, overflowing usize")
    });
    let mut ptr = Vec::new();
    ptr.try_reserve_exact(len).unwrap_or_else(|err| {
        panic!("compressed storage for dimension {dim} needs {len} pointers: {err}")
    });
    ptr.resize(len, 0);
    ptr
}

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// Wraps arrays that already satisfy the type's invariants: `ptr` has one entry per major
    /// index plus one, from 0 to `idx.len()`, and each segment of `idx` is strictly increasing
    /// and below the minor dimension. The caller guarantees this; it's checked in debug builds.
    pub(crate) fn from_parts_unchecked(
        shape: (usize, usize),
        ptr: Vec<usize>,
        idx: Vec<usize>,
        data: Vec<T>,
    ) -> Self {
        let local = Self {
            shape,
            ptr,
            idx,
            data,
            orientation: PhantomData,
            diagonal: OnceLock::new(),
        };
        debug_assert_eq!(
            validate(
                local.major_dim(),
                local.minor_dim(),
                O::ROW_MAJOR,
                &local.ptr,
                &local.idx,
                local.data.len()
            ),
            Ok(())
        );
        local
    }

    /// Builds a `rows x cols` matrix from raw compressed arrays: `pointers` has one entry per
    /// major index (row for CSR, column for CSC) plus one, non-decreasing from 0 to the entry
    /// count, and major index `i`'s entries are `indices[pointers[i]..pointers[i + 1]]`, strictly
    /// increasing and in range, with their `values`.
    ///
    /// Zero values are kept as stored entries. Returns `SparseError::InvalidCompressed` naming
    /// the problem if the arrays are inconsistent with each other or with the shape.
    pub fn try_from_raw_parts(
        rows: usize,
        cols: usize,
        pointers: Vec<usize>,
        indices: Vec<usize>,
        values: Vec<T>,
    ) -> Result<Self, SparseError> {
        let (major, minor) = if O::ROW_MAJOR {
            (rows, cols)
        } else {
            (cols, rows)
        };
        validate(
            major,
            minor,
            O::ROW_MAJOR,
            &pointers,
            &indices,
            values.len(),
        )
        .map_err(|reason| SparseError::InvalidCompressed { reason })?;
        Ok(Self {
            shape: (rows, cols),
            ptr: pointers,
            idx: indices,
            data: values,
            orientation: PhantomData,
            diagonal: OnceLock::new(),
        })
    }

    /// Builds a `rows x cols` matrix from `(row, col, value)` triplets using **COO semantics**:
    /// triplets that share a position are summed, in the order given.
    ///
    /// Every position that appears is stored, even when its values are zero or sum to exactly
    /// zero, so the pattern depends only on the positions. (`SparseMatrix::from_triplets`
    /// drops such zeros.) Returns `SparseError::OutOfBounds` if any triplet lies outside
    /// `rows x cols`.
    pub fn from_triplets(
        rows: usize,
        cols: usize,
        // COO semantics: duplicate (row, col) entries are summed
        triplets: impl IntoIterator<Item = (usize, usize, T)>,
    ) -> Result<Self, SparseError> {
        let triplets: Vec<(usize, usize, T)> = triplets.into_iter().collect();
        Self::assemble(rows, cols, &triplets, None)
    }

    /// `from_triplets` on collected triplets. With `slots`, also records each triplet's slot in
    /// `data` and whether it was the first of its position, which a refill writes rather than
    /// adds to (`assembly.rs`).
    pub(crate) fn assemble(
        rows: usize,
        cols: usize,
        triplets: &[(usize, usize, T)],
        mut slots: Option<&mut [(usize, bool)]>,
    ) -> Result<Self, SparseError> {
        if let Some(&(row, col, _)) = triplets
            .iter()
            .find(|&&(row, col, _)| row >= rows || col >= cols)
        {
            return Err(SparseError::OutOfBounds {
                row,
                col,
                shape: (rows, cols),
            });
        }
        let (major_dim, minor_dim) = if O::ROW_MAJOR {
            (rows, cols)
        } else {
            (cols, rows)
        };
        let major_minor =
            |row: usize, col: usize| if O::ROW_MAJOR { (row, col) } else { (col, row) };

        // Group the triplets by major index with a counting sort, keeping their input order and
        //   where each came from
        let mut start = zeroed_pointers(major_dim);
        for &(row, col, _) in triplets {
            start[major_minor(row, col).0 + 1] += 1;
        }
        for i in 0..major_dim {
            start[i + 1] += start[i];
        }
        let mut cursor = start.clone();
        let mut grouped = vec![(0, 0, T::zero()); triplets.len()];
        for (k, &(row, col, value)) in triplets.iter().enumerate() {
            let (major, minor) = major_minor(row, col);
            grouped[cursor[major]] = (minor, k, value);
            cursor[major] += 1;
        }

        // Sort each group by minor index (stably, so duplicates are summed in input order)
        //   and merge the duplicates
        let mut ptr = zeroed_pointers(major_dim);
        let mut idx = Vec::with_capacity(triplets.len());
        let mut data: Vec<T> = Vec::with_capacity(triplets.len());
        for major in 0..major_dim {
            let group = &mut grouped[start[major]..start[major + 1]];
            group.sort_by_key(|entry| entry.0);
            for &(minor, k, value) in group.iter() {
                // A duplicate of the previous entry in this group: add it in
                let first = !(idx.len() > ptr[major] && idx.last() == Some(&minor));
                if first {
                    idx.push(minor);
                    data.push(value);
                } else {
                    let sum = data.last_mut().expect("idx and data have the same length");
                    *sum = *sum + value;
                }
                if let Some(slots) = slots.as_deref_mut() {
                    slots[k] = (data.len() - 1, first);
                }
            }
            ptr[major + 1] = idx.len();
        }
        debug_assert!(idx.iter().all(|&minor| minor < minor_dim));
        Ok(Self::from_parts_unchecked((rows, cols), ptr, idx, data))
    }

    /// Builds a matrix from dense rows, like `SparseMatrix::from_dense`: zeros are not stored,
    /// and the shape is inferred from the rows. Returns `SparseError::RaggedRows` if the rows'
    /// lengths differ.
    pub fn from_dense<R: AsRef<[T]>>(rows: &[R]) -> Result<Self, SparseError> {
        SparseMatrix::from_dense(rows).map(|matrix| Self::from_sparse(&matrix))
    }

    /// The `n x n` identity matrix.
    pub fn identity(n: usize) -> Self {
        Self::from_sparse(&SparseMatrix::identity(n))
    }

    /// A square matrix with `diagonal` on its main diagonal. Zeros are not stored.
    pub fn from_diagonal(diagonal: impl IntoIterator<Item = T>) -> Self {
        Self::from_sparse(&SparseMatrix::from_diagonal(diagonal))
    }

    /// The same matrix in compressed form, storing exactly `matrix`'s entries (none of them
    /// zero). Equivalent to `matrix.to_csr()` or `matrix.to_csc()`.
    pub fn from_sparse(matrix: &SparseMatrix<T>) -> Self {
        matrix.to_compressed()
    }

    /// The same matrix as a `SparseMatrix`, which can be edited freely. Stored zeros are
    /// dropped, since `SparseMatrix` never stores a zero.
    pub fn to_sparse(&self) -> SparseMatrix<T> {
        let values = self
            .iter()
            .filter(|&(_, _, value)| value != T::zero())
            .map(|(row, col, value)| ((row, col), value))
            .collect();
        SparseMatrix::from_values(self.shape, values)
    }

    /// The transpose, which is the same arrays in the other orientation, so it costs O(1): the
    /// transpose of a `CsrMatrix` is a `CscMatrix` and vice versa. Use `to_csr`/`to_csc` on the
    /// result to get it back in the original orientation, at O(nnz).
    pub fn transpose(self) -> CompressedMatrix<T, O::Transposed> {
        CompressedMatrix {
            shape: (self.shape.1, self.shape.0),
            ptr: self.ptr,
            idx: self.idx,
            data: self.data,
            orientation: PhantomData,
            // Entry (i, i) of the transpose is the same entry, in the same slot
            diagonal: self.diagonal,
        }
    }

    /// The same matrix in compressed sparse row form, with the same pattern, stored zeros
    /// included. A copy if it's already CSR; otherwise a reorientation in O(nnz + rows + cols).
    pub fn to_csr(&self) -> CsrMatrix<T> {
        self.reoriented()
    }

    /// The same matrix in compressed sparse column form, with the same pattern, stored zeros
    /// included. A copy if it's already CSC; otherwise a reorientation in O(nnz + rows + cols).
    pub fn to_csc(&self) -> CscMatrix<T> {
        self.reoriented()
    }

    /// The same matrix compressed along `P`'s axis.
    pub(crate) fn reoriented<P: Orientation>(&self) -> CompressedMatrix<T, P> {
        if O::ROW_MAJOR == P::ROW_MAJOR {
            return CompressedMatrix {
                shape: self.shape,
                ptr: self.ptr.clone(),
                idx: self.idx.clone(),
                data: self.data.clone(),
                orientation: PhantomData,
                diagonal: self.diagonal.clone(),
            };
        }
        // A counting sort by minor index. Visiting the major indices in order leaves each new
        //   segment sorted, so no further sort is needed.
        let minor_dim = self.minor_dim();
        let mut ptr = zeroed_pointers(minor_dim);
        for &minor in &self.idx {
            ptr[minor + 1] += 1;
        }
        for i in 0..minor_dim {
            ptr[i + 1] += ptr[i];
        }
        let mut cursor = ptr.clone();
        let mut idx = vec![0; self.idx.len()];
        let mut data = vec![T::zero(); self.data.len()];
        for major in 0..self.major_dim() {
            for p in self.ptr[major]..self.ptr[major + 1] {
                let minor = self.idx[p];
                idx[cursor[minor]] = major;
                data[cursor[minor]] = self.data[p];
                cursor[minor] += 1;
            }
        }
        CompressedMatrix {
            shape: self.shape,
            ptr,
            idx,
            data,
            orientation: PhantomData,
            diagonal: OnceLock::new(),
        }
    }

    /// The shape `(rows, cols)`.
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    /// Whether the matrix has as many rows as columns.
    pub fn is_square(&self) -> bool {
        self.shape.0 == self.shape.1
    }

    /// The number of stored entries: the size of the pattern, stored zeros included.
    pub fn num_stored(&self) -> usize {
        self.idx.len()
    }

    /// The number of stored entries whose value isn't zero.
    pub fn num_nonzero(&self) -> usize {
        self.data
            .iter()
            .filter(|&&value| value != T::zero())
            .count()
    }

    /// The pointers: major index `i`'s entries are at `pointers()[i]..pointers()[i + 1]`.
    pub fn pointers(&self) -> &[usize] {
        &self.ptr
    }

    /// The minor index of each stored entry (a column for CSR, a row for CSC), strictly
    /// increasing within each major index.
    pub fn indices(&self) -> &[usize] {
        &self.idx
    }

    /// The value of each stored entry, parallel to `indices()`.
    pub fn values(&self) -> &[T] {
        &self.data
    }

    /// The values, to change in place. The pattern stays as it is: setting a value to zero
    /// keeps the entry stored.
    pub fn values_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// The arrays, as `(shape, pointers, indices, values)`, which `try_from_raw_parts` accepts
    /// back.
    pub fn into_raw_parts(self) -> ((usize, usize), Vec<usize>, Vec<usize>, Vec<T>) {
        (self.shape, self.ptr, self.idx, self.data)
    }

    /// Iterates over the stored `(row, col, value)` entries in storage order: row-major for
    /// CSR, column-major for CSC. Stored zeros are included.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, T)> + '_ {
        (0..self.major_dim()).flat_map(move |major| {
            (self.ptr[major]..self.ptr[major + 1]).map(move |p| {
                let (row, col) = Self::row_col(major, self.idx[p]);
                (row, col, self.data[p])
            })
        })
    }

    /// The matrix as dense rows. Allocates `rows * cols` elements, so this is only suitable for
    /// shapes that fit in memory.
    pub fn to_dense(&self) -> Vec<Vec<T>> {
        let mut dense = vec![vec![T::zero(); self.shape.1]; self.shape.0];
        for (row, col, value) in self.iter() {
            dense[row][col] = value;
        }
        dense
    }

    /// The length of the compressed axis: rows for CSR, columns for CSC.
    pub(crate) fn major_dim(&self) -> usize {
        if O::ROW_MAJOR {
            self.shape.0
        } else {
            self.shape.1
        }
    }

    /// The length of the other axis.
    pub(crate) fn minor_dim(&self) -> usize {
        if O::ROW_MAJOR {
            self.shape.1
        } else {
            self.shape.0
        }
    }

    /// The `(row, col)` of the entry at `(major, minor)`.
    pub(crate) fn row_col(major: usize, minor: usize) -> (usize, usize) {
        if O::ROW_MAJOR {
            (major, minor)
        } else {
            (minor, major)
        }
    }

    /// The stored entries whose value isn't zero, in storage order.
    fn nonzeros(&self) -> impl Iterator<Item = (usize, usize, T)> + '_ {
        self.iter().filter(|&(_, _, value)| value != T::zero())
    }

    /// The slot in `data` of the entry at `(row, col)`, or `None` if the position isn't in the
    /// pattern, found by binary search within its major index. The caller has checked the
    /// bounds.
    pub(crate) fn slot(&self, row: usize, col: usize) -> Option<usize> {
        let (major, minor) = if O::ROW_MAJOR { (row, col) } else { (col, row) };
        let segment = self.ptr[major]..self.ptr[major + 1];
        self.idx[segment.clone()]
            .binary_search(&minor)
            .ok()
            .map(|offset| segment.start + offset)
    }

    /// The value stored at `(row, col)`, or `None` if the position isn't in the pattern. The
    /// caller has checked the bounds.
    pub(crate) fn stored_value(&self, row: usize, col: usize) -> Option<T> {
        self.slot(row, col).map(|p| self.data[p])
    }
}

/// Checks raw compressed arrays against a `major x minor` shape (in compressed terms), naming
/// the first problem in the words of `try_from_raw_parts`.
fn validate(
    major: usize,
    minor: usize,
    row_major: bool,
    ptr: &[usize],
    idx: &[usize],
    data_len: usize,
) -> Result<(), &'static str> {
    // Checked: major + 1 overflows for major == usize::MAX, and no pointer array can be that long
    if major.checked_add(1) != Some(ptr.len()) {
        return Err(if row_major {
            "pointers must have rows + 1 entries"
        } else {
            "pointers must have cols + 1 entries"
        });
    }
    if idx.len() != data_len {
        return Err("indices and values must have the same length");
    }
    if ptr[0] != 0 || ptr[major] != idx.len() {
        return Err("pointers must start at 0 and end at the entry count");
    }
    for bounds in ptr.windows(2) {
        if bounds[0] > bounds[1] {
            return Err("pointers must be non-decreasing");
        }
        let mut previous: Option<usize> = None;
        for &index in &idx[bounds[0]..bounds[1]] {
            if index >= minor {
                return Err(if row_major {
                    "column index out of range"
                } else {
                    "row index out of range"
                });
            }
            if previous.is_some_and(|previous| previous >= index) {
                return Err(if row_major {
                    "column indices must be strictly increasing within a row"
                } else {
                    "row indices must be strictly increasing within a column"
                });
            }
            previous = Some(index);
        }
    }
    Ok(())
}

/// An empty `0 x 0` matrix.
impl<T: Scalar, O: Orientation> Default for CompressedMatrix<T, O> {
    fn default() -> Self {
        Self::from_parts_unchecked((0, 0), vec![0], vec![], vec![])
    }
}

/// Shows the orientation, the shape and the stored entries in storage order, stored zeros
/// included.
impl<T: Scalar + fmt::Debug, O: Orientation> fmt::Debug for CompressedMatrix<T, O> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let entries: Vec<(usize, usize, T)> = self.iter().collect();
        f.debug_struct(if O::ROW_MAJOR {
            "CsrMatrix"
        } else {
            "CscMatrix"
        })
        .field("shape", &self.shape)
        .field("entries", &entries)
        .finish()
    }
}

/// Formats the matrix as a bracketed grid of dense rows, one per line, as `SparseMatrix` does.
impl<T: Scalar + fmt::Display, O: Orientation> fmt::Display for CompressedMatrix<T, O> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write_dense_grid(f, self.shape.1, self.to_dense())
    }
}

/// Mathematical equality: the same shape and the same nonzero entries. Orientation and stored
/// zeros don't matter, so a `CsrMatrix` can equal a `CscMatrix`; comparing two orientations
/// reorients one of them, at O(nnz).
///
/// Like floats themselves this is only a partial equivalence: a matrix holding `NaN` is not
/// equal to itself.
impl<T: Scalar, O: Orientation, P: Orientation> PartialEq<CompressedMatrix<T, P>>
    for CompressedMatrix<T, O>
{
    fn eq(&self, other: &CompressedMatrix<T, P>) -> bool {
        if self.shape != other.shape {
            return false;
        }
        if O::ROW_MAJOR == P::ROW_MAJOR {
            self.nonzeros().eq(other.nonzeros())
        } else {
            self.nonzeros().eq(other.reoriented::<O>().nonzeros())
        }
    }
}

/// Mathematical equality with a `SparseMatrix`: the same shape and the same nonzero entries.
impl<T: Scalar, O: Orientation> PartialEq<SparseMatrix<T>> for CompressedMatrix<T, O> {
    fn eq(&self, other: &SparseMatrix<T>) -> bool {
        self.shape == other.shape
            && self.num_nonzero() == other.num_nonzero()
            && self
                .nonzeros()
                .all(|(row, col, value)| other.peek_at(row, col) == Ok(Some(value)))
    }
}

/// Mathematical equality with a `CompressedMatrix`: the same shape and the same nonzero
/// entries.
impl<T: Scalar, O: Orientation> PartialEq<CompressedMatrix<T, O>> for SparseMatrix<T> {
    fn eq(&self, other: &CompressedMatrix<T, O>) -> bool {
        other == self
    }
}

/// Converts with `CompressedMatrix::from_sparse`.
impl<T: Scalar, O: Orientation> From<&SparseMatrix<T>> for CompressedMatrix<T, O> {
    fn from(matrix: &SparseMatrix<T>) -> Self {
        Self::from_sparse(matrix)
    }
}
