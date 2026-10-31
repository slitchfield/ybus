use std::borrow::Cow;
use std::collections::HashMap;

use crate::compressed_matrix::{
    zeroed_pointers, CompressedMatrix, CscMatrix, CsrMatrix, Orientation,
};
use crate::sparse_error::SparseError;
use crate::Scalar;

mod access;
mod blocks;
mod complex;
mod matrix_market;
mod ops;
#[cfg(feature = "rayon")]
mod parallel;
mod product;
mod properties;
mod reorder;
#[cfg(feature = "serde")]
mod serde_impl;
mod solvers;
mod traits;

pub use solvers::Lu;

pub use access::{ColIterator, Entries, SparseVector};

/// A sparse matrix of `T` values (`f64` by default; any [`Scalar`] type works).
///
/// Entries live in a dictionary of keys (a `HashMap` from `(row, col)` to value), which
/// makes single-element reads and writes cheap. A zero is never stored: inserting a zero
/// removes any existing entry, so `num_nonzero` is the number of stored entries.
///
/// Compressed sparse row (CSR) and compressed sparse column (CSC) copies are built on demand
/// by `explicitly_compress` and `explicitly_compress_csc` and exposed through `compressed`
/// and `compressed_csc`. Any mutation drops both copies. Operations that need compressed
/// storage (row and column access, products) read a fresh cached copy when there is one and
/// otherwise build a temporary one, so they are always correct but cheaper after compressing.
///
/// A compressed copy needs one pointer per row (CSR) or column (CSC) plus one, however few
/// entries are stored. Shapes too large for that (for example a dimension of `usize::MAX`)
/// still work for element access and construction, but operations that need the compressed
/// form panic with a message saying so.
#[derive(Clone)]
pub struct SparseMatrix<T = f64> {
    /// The matrix dimensions as `(rows, cols)`.
    pub shape: (usize, usize),
    /// The stored (non-zero) entries, keyed by `(row, col)`.
    values: HashMap<(usize, usize), T>,

    /// Cached CSR copy of `values`; `None` once any mutation makes it stale.
    csr: Option<Compressed<T>>,
    /// Cached CSC copy of `values`; `None` once any mutation makes it stale.
    csc: Option<Compressed<T>>,
    /// A zero to hand out by reference for entries that are not stored (see `Index`).
    zero: T,
}

/// Borrowed view of the matrix's cached CSR arrays, as returned by
/// `SparseMatrix::compressed`.
pub struct CsrView<'a, T = f64> {
    /// Row pointers: row `i`'s entries are `cols[rows[i]..rows[i + 1]]`.
    pub rows: &'a [usize],
    /// Column index of each stored entry, sorted within each row.
    pub cols: &'a [usize],
    /// Value of each stored entry, parallel to `cols`.
    pub data: &'a [T],
}

/// Borrowed view of the matrix's cached CSC arrays, as returned by
/// `SparseMatrix::compressed_csc`.
pub struct CscView<'a, T = f64> {
    /// Column pointers: column `j`'s entries are `rows[cols[j]..cols[j + 1]]`.
    pub cols: &'a [usize],
    /// Row index of each stored entry, sorted within each column.
    pub rows: &'a [usize],
    /// Value of each stored entry, parallel to `rows`.
    pub data: &'a [T],
}

/// Compressed sparse storage along one axis.
///
/// For CSR, `ptr` holds row pointers and `idx` column indices; for CSC the roles swap. `idx`
/// is sorted within each segment of `ptr`, and `data` is parallel to `idx`.
#[derive(Clone)]
pub(crate) struct Compressed<T> {
    pub(crate) ptr: Vec<usize>,
    pub(crate) idx: Vec<usize>,
    pub(crate) data: Vec<T>,
}

/// Iterator over a matrix's rows, created by `SparseMatrix::row_iter`.
///
/// Each item is one dense row of length `shape.1`.
pub struct RowIterator<'a, T: Scalar = f64> {
    shape: (usize, usize),
    /// Borrowed from the matrix's cache when fresh, otherwise built for this iterator.
    csr: Cow<'a, Compressed<T>>,
    row_iter_idx: usize,
}

impl<T: Scalar> Iterator for RowIterator<'_, T> {
    /// One dense row.
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.row_iter_idx < self.shape.0 {
            let start = self.csr.ptr[self.row_iter_idx];
            let end = self.csr.ptr[self.row_iter_idx + 1];

            let colslice = &self.csr.idx[start..end];
            let dataslice = &self.csr.data[start..end];

            let mut retvec: Vec<T> = vec![T::zero(); self.shape.1];

            for (col, val) in std::iter::zip(colslice, dataslice) {
                retvec[*col] = *val;
            }

            self.row_iter_idx += 1;

            Some(retvec)
        } else {
            None
        }
    }
}

impl<T: Scalar> Default for SparseMatrix<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Scalar> SparseMatrix<T> {
    /// Builds compressed storage from `values` without modifying the caches: CSR
    /// (`by_col = false`) or CSC (`by_col = true`).
    fn build_compressed(&self, by_col: bool) -> Compressed<T> {
        let major_dim = if by_col { self.shape.1 } else { self.shape.0 };
        let major_of = |&(row, col): &(usize, usize)| if by_col { col } else { row };

        // Count entries per major index, then prefix-sum into segment pointers
        let mut ptr = zeroed_pointers(major_dim);
        for key in self.values.keys() {
            ptr[major_of(key) + 1] += 1;
        }
        let mut running = 0;
        for p in ptr.iter_mut() {
            running += *p;
            *p = running;
        }

        // Scatter each entry into its segment, using ptr[major] as the segment's write
        //   cursor. Afterwards each ptr[major] has advanced to the next segment's start, so
        //   shifting the array right by one restores the pointers without a second copy.
        let mut entries = vec![(0, T::zero()); self.values.len()];
        for (key, val) in self.values.iter() {
            let major = major_of(key);
            let minor = if by_col { key.0 } else { key.1 };
            entries[ptr[major]] = (minor, *val);
            ptr[major] += 1;
        }
        ptr.copy_within(..major_dim, 1);
        ptr[0] = 0;

        for segment in ptr.windows(2) {
            entries[segment[0]..segment[1]].sort_unstable_by_key(|entry| entry.0);
        }

        let (idx, data) = entries.into_iter().unzip();
        Compressed { ptr, idx, data }
    }

    /// The matrix in compressed sparse row form, storing exactly its entries (none of them
    /// zero), with a fixed pattern for numerical work. Costs O(nnz) plus sorting each row,
    /// unless a fresh CSR cache exists, which is copied.
    pub fn to_csr(&self) -> CsrMatrix<T> {
        self.to_compressed()
    }

    /// The matrix in compressed sparse column form, storing exactly its entries (none of them
    /// zero), with a fixed pattern for numerical work. Costs O(nnz) plus sorting each column,
    /// unless a fresh CSC cache exists, which is copied.
    pub fn to_csc(&self) -> CscMatrix<T> {
        self.to_compressed()
    }

    /// The matrix compressed along `O`'s axis; `to_csr` and `to_csc` name the two cases.
    pub fn to_compressed<O: Orientation>(&self) -> CompressedMatrix<T, O> {
        let compressed = if O::ROW_MAJOR {
            self.csr_cow()
        } else {
            self.csc_cow()
        }
        .into_owned();
        CompressedMatrix::from_parts_unchecked(
            self.shape,
            compressed.ptr,
            compressed.idx,
            compressed.data,
        )
    }

    /// Drops the cached compressed copies; call after any change to `values` or `shape`.
    fn invalidate(&mut self) {
        self.csr = None;
        self.csc = None;
    }

    /// The CSC form: borrowed from the cache when fresh, otherwise built on the spot.
    pub(crate) fn csc_cow(&self) -> Cow<'_, Compressed<T>> {
        match &self.csc {
            Some(csc) => Cow::Borrowed(csc),
            None => Cow::Owned(self.build_compressed(true)),
        }
    }

    /// The CSR form: borrowed from the cache when fresh, otherwise built on the spot.
    pub(crate) fn csr_cow(&self) -> Cow<'_, Compressed<T>> {
        match &self.csr {
            Some(csr) => Cow::Borrowed(csr),
            None => Cow::Owned(self.build_compressed(false)),
        }
    }

    /// Rebuilds the cached CSR arrays so that `compressed` returns `Some` again.
    ///
    /// Worthwhile before iterating repeatedly, since `row_iter` otherwise rebuilds the
    /// arrays for every iterator while the cache is stale.
    pub fn explicitly_compress(&mut self) {
        self.csr = Some(self.build_compressed(false));
    }

    /// Rebuilds the cached CSC arrays so that `compressed_csc` returns `Some` again.
    ///
    /// Worthwhile before repeated column access or transposed products.
    pub fn explicitly_compress_csc(&mut self) {
        self.csc = Some(self.build_compressed(true));
    }

    /// The cached CSR arrays, or `None` if they are stale (or were never built).
    /// Call `explicitly_compress` to refresh them.
    pub fn compressed(&self) -> Option<CsrView<'_, T>> {
        self.csr.as_ref().map(|csr| CsrView {
            rows: &csr.ptr,
            cols: &csr.idx,
            data: &csr.data,
        })
    }

    /// The cached CSC arrays, or `None` if they are stale (or were never built).
    /// Call `explicitly_compress_csc` to refresh them.
    pub fn compressed_csc(&self) -> Option<CscView<'_, T>> {
        self.csc.as_ref().map(|csc| CscView {
            cols: &csc.ptr,
            rows: &csc.idx,
            data: &csc.data,
        })
    }

    /// Iterates over the rows of the matrix, yielding each as a dense `Vec<T>`.
    ///
    /// Always reflects the current contents. Call `explicitly_compress` first when
    /// iterating repeatedly; otherwise each iterator rebuilds the CSR arrays itself.
    #[allow(dead_code)]
    pub fn row_iter(&self) -> RowIterator<'_, T> {
        RowIterator {
            shape: self.shape,
            csr: self.csr_cow(),
            row_iter_idx: 0,
        }
    }

    /// Creates an empty `0 x 0` matrix.
    #[allow(dead_code)]
    pub fn new() -> SparseMatrix<T> {
        Self {
            shape: (0, 0),
            values: HashMap::new(),
            csr: None,
            csc: None,
            zero: T::zero(),
        }
    }

    /// Creates an `n x m` matrix with no stored entries.
    #[allow(dead_code)]
    pub fn empty_with_shape(n: usize, m: usize) -> SparseMatrix<T> {
        // No up-front reservation: n * m says nothing about how many entries a sparse
        //   matrix will hold, and can overflow or be enormous. The map grows as needed.
        Self::from_values((n, m), HashMap::new())
    }

    /// Wraps an existing entry map in a matrix of the given shape, with a stale cache.
    /// The caller guarantees every key is in bounds and no value is zero.
    pub(crate) fn from_values(
        shape: (usize, usize),
        values: HashMap<(usize, usize), T>,
    ) -> SparseMatrix<T> {
        Self {
            shape,
            values,
            csr: None,
            csc: None,
            zero: T::zero(),
        }
    }

    /// Builds a matrix from `(row, col, value)` triplets using **COO (coordinate list)
    /// semantics**: triplets that share a `(row, col)` are *summed*, not overwritten, and
    /// entries that sum to exactly zero are not stored. Compare `insert_triplets`, which
    /// overwrites duplicates (last write wins).
    ///
    /// Returns `SparseError::OutOfBounds` if any triplet lies outside `rows x cols`.
    #[allow(dead_code)]
    pub fn from_triplets(
        rows: usize,
        cols: usize,
        // COO semantics: duplicate (row, col) entries are summed
        triplets: impl IntoIterator<Item = (usize, usize, T)>,
    ) -> Result<SparseMatrix<T>, SparseError> {
        let triplets = triplets.into_iter();
        // Sized by the number of triplets, not rows * cols, which can be enormous
        let mut values = HashMap::with_capacity(triplets.size_hint().0);

        for (row, col, val) in triplets {
            if row >= rows || col >= cols {
                return Err(SparseError::OutOfBounds {
                    row,
                    col,
                    shape: (rows, cols),
                });
            }
            let entry = values.entry((row, col)).or_insert(T::zero());
            *entry = *entry + val;
        }

        let mut local = Self::from_values((rows, cols), values);
        // Zero inputs and exact cancellations leave stored zeros; drop them
        local.prune_zeroes();
        Ok(local)
    }

    /// Builds a matrix from dense rows, e.g. `Vec<Vec<T>>` or `[[T; 3]; 2]`. Zeros are
    /// not stored. The shape is inferred: the number of rows, by the length of the first row.
    /// Empty input gives a 0 x 0 matrix.
    ///
    /// Returns `SparseError::RaggedRows` if any row's length differs from the first row's.
    /// Accepts `row_iter()` output, so `from_dense(m.row_iter())` rebuilds `m` (except a
    /// `0 x m` matrix, which has no rows to carry its column count and rebuilds as `0 x 0`).
    #[allow(dead_code)]
    pub fn from_dense<R: AsRef<[T]>>(
        rows: impl IntoIterator<Item = R>,
    ) -> Result<SparseMatrix<T>, SparseError> {
        let mut values = HashMap::new();
        let mut ncols: Option<usize> = None;
        let mut nrows: usize = 0;

        for (row_idx, row) in rows.into_iter().enumerate() {
            let row = row.as_ref();
            let expected = *ncols.get_or_insert(row.len());
            if row.len() != expected {
                return Err(SparseError::RaggedRows {
                    row: row_idx,
                    expected,
                    found: row.len(),
                });
            }
            for (col_idx, &val) in row.iter().enumerate() {
                // Zeros (including -0.0 for floats) are never stored
                if val != T::zero() {
                    values.insert((row_idx, col_idx), val);
                }
            }
            nrows += 1;
        }

        Ok(Self::from_values((nrows, ncols.unwrap_or(0)), values))
    }

    /// Builds a `rows x cols` matrix from raw CSR arrays: `row_ptr` has `rows + 1`
    /// non-decreasing entries starting at 0, and row `i`'s entries are
    /// `col_idx[row_ptr[i]..row_ptr[i + 1]]` with matching `data`, with strictly increasing
    /// column indices in range.
    ///
    /// Zero values are dropped. When there are none, the arrays become the matrix's CSR cache
    /// directly, so `compressed` is `Some` straight away. Returns
    /// `SparseError::InvalidCompressed` if the arrays are inconsistent, as
    /// `CsrMatrix::try_from_raw_parts` does.
    #[allow(dead_code)]
    pub fn from_csr(
        rows: usize,
        cols: usize,
        row_ptr: Vec<usize>,
        col_idx: Vec<usize>,
        data: Vec<T>,
    ) -> Result<SparseMatrix<T>, SparseError> {
        let csr = CsrMatrix::try_from_raw_parts(rows, cols, row_ptr, col_idx, data)?;
        let mut local = csr.to_sparse();
        if local.values.len() == csr.num_stored() {
            // No zeros were dropped, so the arrays serve as the CSR cache as they are
            let (_, ptr, idx, data) = csr.into_raw_parts();
            local.csr = Some(Compressed { ptr, idx, data });
        }
        Ok(local)
    }

    /// Builds a square matrix with `diagonal` on its main diagonal. The size is the number of
    /// items, and zeros are not stored.
    ///
    /// Named `from_diagonal` because `diagonal` reads a matrix's diagonal back out.
    #[allow(dead_code)]
    pub fn from_diagonal(diagonal: impl IntoIterator<Item = T>) -> SparseMatrix<T> {
        let mut values = HashMap::new();
        let mut size = 0;
        for (i, val) in diagonal.into_iter().enumerate() {
            if val != T::zero() {
                values.insert((i, i), val);
            }
            size = i + 1;
        }
        Self::from_values((size, size), values)
    }

    /// Returns the matrix as dense rows. Allocates `rows * cols` elements, so this is only
    /// suitable for shapes that fit in memory. The inverse of `from_dense`, except that a
    /// `0 x m` matrix has no rows to carry its column count, so it round-trips as `0 x 0`.
    #[allow(dead_code)]
    pub fn to_dense(&self) -> Vec<Vec<T>> {
        let mut dense = vec![vec![T::zero(); self.shape.1]; self.shape.0];
        for ((row, col), val) in self.values.iter() {
            dense[*row][*col] = *val;
        }
        dense
    }

    /// Creates the `n x n` identity matrix.
    #[allow(dead_code)]
    pub fn identity(n: usize) -> SparseMatrix<T> {
        let mut local = Self::empty_with_shape(n, n);
        local.values.reserve(n);
        for diag_idx in 0..n {
            // Unwrap safe by construction of diag_idx
            local.insert(diag_idx, diag_idx, T::one()).unwrap();
        }
        local
    }

    /// Returns the transpose as a new matrix, leaving `self` unchanged. See
    /// `transpose_inplace` to transpose without allocating a second matrix.
    #[allow(dead_code)]
    pub fn create_transpose(&self) -> SparseMatrix<T> {
        let mut local = Self::empty_with_shape(self.shape.1, self.shape.0);
        local.values.reserve(self.values.len());
        for ((row, col), val) in self.values.iter() {
            // Unwrap safe by construction of col, row
            local.insert(*col, *row, *val).unwrap(); // Deref okay due to elementary r, c, v types
        }
        local
    }

    /// Sets the entry at `(row, col)` to `value`, replacing any existing value.
    ///
    /// Inserting a zero removes the entry instead, since zeros are never stored.
    /// Returns `SparseError::OutOfBounds` if the position is outside the matrix.
    #[allow(dead_code)]
    pub fn insert(&mut self, row: usize, col: usize, value: T) -> Result<(), SparseError> {
        if row >= self.shape.0 || col >= self.shape.1 {
            return Err(SparseError::OutOfBounds {
                row,
                col,
                shape: self.shape,
            });
        }

        // If new value is exactly zero, clear any existing entry
        if value == T::zero() {
            // UNWRAP: already validated row, col above
            self.clear_at(row, col).unwrap();
        } else {
            self.values.insert((row, col), value);
        }
        self.invalidate();
        Ok(())
    }

    /// Inserts several `(row, col, value)` triplets, as if by repeated `insert`.
    ///
    /// Duplicate `(row, col)` triplets overwrite (last write wins). See `from_triplets`
    /// for COO semantics, where duplicates are summed.
    ///
    /// Every triplet is bounds-checked before any is applied, so on
    /// `SparseError::OutOfBounds` the matrix is left unchanged.
    #[allow(dead_code)]
    pub fn insert_triplets(&mut self, triplets: Vec<(usize, usize, T)>) -> Result<(), SparseError> {
        // Iterate through triplets first to validate no oob errors before changing
        //   underlying matrix state
        for (row, col, _) in triplets.iter() {
            if *row >= self.shape.0 || *col >= self.shape.1 {
                return Err(SparseError::OutOfBounds {
                    row: *row,
                    col: *col,
                    shape: self.shape,
                });
            }
        }

        for (row, col, val) in triplets.iter() {
            // If new value is exactly zero, clear any existing entry
            if *val == T::zero() {
                // UNWRAP: Already validated row, col above
                self.clear_at(*row, *col).unwrap();
            } else {
                self.values.insert((*row, *col), *val);
            }
        }
        self.invalidate();
        Ok(())
    }

    /// Removes the entry at `(row, col)` and returns it, or `None` if nothing was stored
    /// there.
    ///
    /// Returns `SparseError::OutOfBounds` if the position is outside the matrix.
    #[allow(dead_code)]
    pub fn clear_at(&mut self, row: usize, col: usize) -> Result<Option<T>, SparseError> {
        if row >= self.shape.0 || col >= self.shape.1 {
            return Err(SparseError::OutOfBounds {
                row,
                col,
                shape: self.shape,
            });
        }

        self.invalidate();
        Ok(self.values.remove(&(row, col)))
    }

    /// Returns the entry at `(row, col)`, or `None` if nothing is stored there (an implicit
    /// zero).
    ///
    /// Returns `SparseError::OutOfBounds` if the position is outside the matrix.
    #[allow(dead_code)]
    pub fn peek_at(&self, row: usize, col: usize) -> Result<Option<T>, SparseError> {
        if row >= self.shape.0 || col >= self.shape.1 {
            return Err(SparseError::OutOfBounds {
                row,
                col,
                shape: self.shape,
            });
        }

        Ok(self.values.get(&(row, col)).copied())
    }

    /// Returns the number of stored entries.
    #[allow(dead_code)]
    pub fn num_nonzero(&self) -> usize {
        self.values.len()
    }

    /// Transposes the matrix in place, swapping `shape` and every entry's row and column.
    /// See `create_transpose` for a version that returns a new matrix.
    #[allow(dead_code)]
    pub fn transpose_inplace(&mut self) {
        // Naive impl, could do better
        self.shape = (self.shape.1, self.shape.0);

        let triplets: Vec<((usize, usize), T)> = self.values.drain().collect();

        for ((row, col), val) in triplets {
            self.values.insert((col, row), val);
        }
        self.invalidate();
    }

    /// Removes any stored entries equal to zero.
    ///
    /// Every public method already avoids storing zeros, so this is a safety net that
    /// leaves the cache untouched when there is nothing to remove.
    #[allow(dead_code)]
    pub fn prune_zeroes(&mut self) {
        let init_len = self.values.len();
        self.values.retain(|_, v| *v != T::zero());
        let final_len = self.values.len();
        if init_len != final_len {
            self.invalidate();
        }
    }
}

use std::fmt;
/// Formats the matrix as a bracketed grid of dense rows, one per line.
impl<T: Scalar + fmt::Display> fmt::Display for SparseMatrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write_dense_grid(f, self.shape.1, self.row_iter())
    }
}

/// Writes dense rows of `cols` values as the bracketed grid both matrix types display as.
pub(crate) fn write_dense_grid<T: fmt::Display>(
    f: &mut fmt::Formatter,
    cols: usize,
    rows: impl IntoIterator<Item = Vec<T>>,
) -> fmt::Result {
    let inner_line_width = 8 * cols; // 6 chars per col + comma + space + leading space
    write!(f, "\t/")?;
    for _ in 0..inner_line_width {
        write!(f, " ")?;
    }
    writeln!(f, "\\")?;
    // TODO: Account for variable number of digits in cols
    for row in rows {
        write!(f, "\t| ")?;
        for (idx, elem) in row.iter().enumerate() {
            if idx != 0 {
                write!(f, ", ")?;
            }
            write!(f, "{:>6.2}", elem)?; // TODO: dynamic precision based on longest values
        }
        writeln!(f, " |")?;
    }
    write!(f, "\t\\")?;
    for _ in 0..inner_line_width {
        write!(f, " ")?;
    }
    write!(f, "/")?;
    writeln!(f)
}
