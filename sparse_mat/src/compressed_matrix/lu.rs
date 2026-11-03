//! Sparse LU factorization: Gilbert and Peierls's left-looking algorithm with threshold partial
//! pivoting, for real and complex matrices.
//!
//! [`CscMatrix::analyze_lu`] chooses a column ordering `Q` from the pattern of `A + Aᵀ`.
//! [`SymbolicLu::factor`] then computes `P A Q = L U` one column at a time: a depth-first search
//! finds which entries of the column can be nonzero, a sparse triangular solve computes them,
//! and a pivot is chosen among the rows not yet pivoted, preferring the diagonal (which keeps the
//! ordering's sparsity) while it is large enough. [`Lu::refactor`] reuses the pivot sequence and
//! the patterns of `L` and `U` for a matrix with new values and the same pattern.

use num_complex::ComplexFloat;
use num_traits::NumCast;

use super::access::check_permutation;
use super::cholesky::inverse_permutation;
use super::solvers::reach;
use super::{CompressedMatrix, CscMatrix, Ordering};
use crate::sparse_error::SparseError;
use crate::sparse_matrix::SparseVector;
use crate::Scalar;

/// Options for an LU factorization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LuOptions {
    /// How large a pivot must be, relative to the largest candidate in its column: the
    /// diagonal entry is kept as the pivot if its magnitude is at least `pivot_threshold`
    /// times the largest, and otherwise the largest is taken. Must be in `(0, 1]`: 1.0 is true
    /// partial pivoting; smaller values keep more of the ordering's sparsity, at the cost of
    /// allowing more growth in the factors. The default, 0.1, is UMFPACK's default partial
    /// pivot tolerance.
    pub pivot_threshold: f64,
}

impl Default for LuOptions {
    fn default() -> LuOptions {
        LuOptions {
            pivot_threshold: 0.1,
        }
    }
}

/// The symbolic analysis for an LU factorization, from [`CscMatrix::analyze_lu`]: the column
/// ordering `Q`, chosen from the pattern alone, and the pattern it was chosen for.
#[derive(Clone, Debug)]
pub struct SymbolicLu {
    col_perm: Vec<usize>,
    /// Whether to keep a large enough diagonal entry as the pivot: true under a symmetric
    /// ordering, false under a column-only one (COLAMD), whose diagonal means nothing
    prefer_diagonal: bool,
    a_ptr: Vec<usize>,
    a_idx: Vec<usize>,
}

/// The outcome of a successful [`Lu::refactor`] call: whether the previous pivots were reused or
/// the matrix had to be factored afresh.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RefactorOutcome {
    /// Every pivot of the previous factorization was still large enough, so its pivot
    /// sequence and patterns were reused.
    Reused,
    /// Some reused pivot fell below the threshold, so the matrix was factored afresh, with new
    /// pivots, under the same column ordering.
    Repivoted,
}

/// An LU factorization `P A Q = L U` with threshold partial pivoting, from
/// [`SymbolicLu::factor`] or [`CscMatrix::lu`]. Keep it to solve many right-hand sides against
/// the same matrix, or to refactor a matrix with the same pattern.
///
/// `L` is unit lower triangular and `U` upper triangular, both in CSC form, and each stores
/// its whole symbolic pattern, so some entries may be zero. Row `i` of `P A Q` is row
/// `permutation()[i]` of `A`, and column `j` is column `column_permutation()[j]`.
#[derive(Clone, Debug)]
pub struct Lu<T: Scalar = f64> {
    l: CscMatrix<T>,
    u: CscMatrix<T>,
    row_perm: Vec<usize>,
    /// Row `i` of `A` is row `row_perm_inverse[i]` of `P A Q`
    row_perm_inverse: Vec<usize>,
    col_perm: Vec<usize>,
    options: LuOptions,
    /// As in `SymbolicLu`, for a fresh factorization in `refactor`
    prefer_diagonal: bool,
    a_ptr: Vec<usize>,
    a_idx: Vec<usize>,
}

impl<T: Scalar> CscMatrix<T> {
    /// Analyzes this matrix's pattern for an LU factorization whose columns, and preferred
    /// pivots, follow `ordering` applied to the pattern of `A + Aᵀ`. Only the pattern is read.
    /// Returns `SparseError::NotSquare` for a non-square matrix.
    ///
    /// Under `Ordering::Colamd`, which orders columns from the pattern of `A` itself, there is
    /// no meaningful diagonal to prefer, so each column's pivot is its largest candidate:
    /// plain partial pivoting, with the threshold unused.
    pub fn analyze_lu(&self, ordering: Ordering) -> Result<SymbolicLu, SparseError> {
        if ordering == Ordering::Colamd && !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        Ok(SymbolicLu {
            col_perm: self.ordering(ordering)?,
            prefer_diagonal: ordering != Ordering::Colamd,
            a_ptr: self.ptr.clone(),
            a_idx: self.idx.clone(),
        })
    }

    /// [`CscMatrix::analyze_lu`] with a column ordering given rather than computed: column `j`
    /// of `A Q` is column `col_perm[j]` of `A`, and the factorization prefers pivots on the
    /// diagonal of `Qᵀ A Q`, as it does under a computed ordering. For a caller that orders many
    /// related matrices from one ordering, such as the submatrices of a larger pattern. Only
    /// the pattern is read.
    ///
    /// Returns `SparseError::NotSquare` for a non-square matrix, and
    /// `SparseError::InvalidPermutation` unless `col_perm` is a permutation of its columns.
    pub fn analyze_lu_with_permutation(
        &self,
        col_perm: Vec<usize>,
    ) -> Result<SymbolicLu, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        check_permutation(&col_perm, self.shape.1)?;
        Ok(SymbolicLu {
            col_perm,
            prefer_diagonal: true,
            a_ptr: self.ptr.clone(),
            a_idx: self.idx.clone(),
        })
    }
}

impl<T: Scalar + ComplexFloat> CscMatrix<T> {
    /// Factors a square matrix as `P A Q = L U`: `analyze_lu(Ordering::Amd)` followed by
    /// `factor` with the default options (pivot threshold 0.1). Analyze once instead to factor
    /// several matrices with the same pattern, or refactor the result.
    ///
    /// Returns `SparseError::NotSquare`, or `SparseError::SingularMatrix` naming the column of
    /// `A` that has no nonzero pivot.
    pub fn lu(&self) -> Result<Lu<T>, SparseError> {
        self.analyze_lu(Ordering::Amd)?
            .factor(self, LuOptions::default())
    }

    /// Solves `A x = b` for a general square matrix by LU factorization (`lu`), failing with
    /// its errors, or `SparseError::ShapeMismatch` if `b` doesn't have one entry per row.
    /// Factor once with `lu` instead when solving several right-hand sides.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        if b.len() != self.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: (b.len(), 1),
            });
        }
        self.lu()?.solve(b)
    }
}

impl SymbolicLu {
    /// The dimension of the analyzed matrix.
    pub fn dim(&self) -> usize {
        self.col_perm.len()
    }

    /// The column ordering `Q`: column `j` of `A Q` is column `column_permutation()[j]` of `A`.
    pub fn column_permutation(&self) -> &[usize] {
        &self.col_perm
    }

    /// Factors `matrix`, which must have the analyzed pattern, as `P A Q = L U`.
    ///
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, or `SparseError::SingularMatrix` naming the column of `A` that has no
    /// nonzero pivot. Panics if `options.pivot_threshold` is not in `(0, 1]`.
    pub fn factor<T: Scalar + ComplexFloat>(
        &self,
        matrix: &CscMatrix<T>,
        options: LuOptions,
    ) -> Result<Lu<T>, SparseError> {
        let n = self.dim();
        if matrix.shape != (n, n) || matrix.ptr != self.a_ptr || matrix.idx != self.a_idx {
            return Err(SparseError::PatternMismatch);
        }
        gilbert_peierls(matrix, &self.col_perm, options, self.prefer_diagonal)
    }
}

#[cfg(feature = "rayon")]
impl<T: Scalar + ComplexFloat + Send + Sync> Lu<T> {
    /// [`Lu::solve_many`] with the right-hand sides solved in parallel (feature `rayon`).
    /// Each is solved by one thread, with its own scratch vector and exactly the operations
    /// `solve` uses, so the results are bit for bit `solve_many`'s, in the same order. Returns
    /// `SparseError::ShapeMismatch` naming the first right-hand side whose length isn't the
    /// number of rows of `A`, before solving any.
    pub fn par_solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        use rayon::prelude::*;
        for b in rhs {
            self.check_rhs(b.len())?;
        }
        let n = self.row_perm.len();
        Ok(rhs
            .par_iter()
            .map_init(
                || Vec::with_capacity(n),
                |scratch, b| {
                    let mut x = b.clone();
                    self.substitute(&mut x, scratch);
                    x
                },
            )
            .collect())
    }
}

/// The threshold as `T`'s real type, checked to be in `(0, 1]`.
fn threshold<T: ComplexFloat>(options: LuOptions) -> T::Real {
    let tau = options.pivot_threshold;
    assert!(
        tau > 0.0 && tau <= 1.0,
        "the pivot threshold must be in (0, 1], not {tau}"
    );
    <T::Real as NumCast>::from(tau).expect("a threshold in (0, 1] converts to any float")
}

/// Left-looking LU with threshold partial pivoting, after J. R. Gilbert and T. Peierls, "Sparse
/// partial pivoting in time proportional to arithmetic operations" (1988), and T. A. Davis,
/// *Direct Methods for Sparse Linear Systems* (2006), ch. 6.
///
/// Column `k` of the factors comes from solving `L x = A[:, q[k]]` with the columns of `L` found
/// so far. Until a row is pivoted it has no place in `L`'s ordering, so while factoring, `L`'s
/// row indices are rows of `A`, and `pinv[row]` is the step at which the row pivoted. A
/// depth-first search from the column's entries through the graph of `L` finds every entry of
/// `x` that can be nonzero, in an order that respects the dependencies (the reach), so the work
/// is proportional to the arithmetic. At the end, `L`'s rows are renumbered by `pinv` and each
/// column of both factors is sorted.
fn gilbert_peierls<T: Scalar + ComplexFloat>(
    a: &CscMatrix<T>,
    q: &[usize],
    options: LuOptions,
    prefer_diagonal: bool,
) -> Result<Lu<T>, SparseError> {
    let tau = threshold::<T>(options);
    let n = q.len();
    let mut pinv: Vec<Option<usize>> = vec![None; n];
    let (mut l_ptr, mut l_idx, mut l_val) = (vec![0], vec![], vec![]);
    let (mut u_ptr, mut u_idx, mut u_val) = (vec![0], vec![], vec![]);
    let mut x = vec![T::zero(); n];
    let mut visited = vec![usize::MAX; n];
    // The reach in postorder, and the depth-first search's stack of (row, next child slot)
    let mut postorder: Vec<usize> = Vec::with_capacity(n);
    let mut stack: Vec<(usize, usize)> = vec![];

    for (k, &col) in q.iter().enumerate() {
        let column = a.ptr[col]..a.ptr[col + 1];

        // The reach of A's column in the graph of L: a pivoted row's children are the rows of
        //   its column of L below the diagonal (which is stored first)
        postorder.clear();
        for &start in &a.idx[column.clone()] {
            if visited[start] == k {
                continue;
            }
            visited[start] = k;
            let first_child = |row: usize| pinv[row].map_or(0, |s| l_ptr[s] + 1);
            stack.push((start, first_child(start)));
            while let Some(top) = stack.len().checked_sub(1) {
                let (row, mut next) = stack[top];
                let end = pinv[row].map_or(0, |s| l_ptr[s + 1]);
                let mut unvisited = None;
                while next < end {
                    let child = l_idx[next];
                    next += 1;
                    if visited[child] != k {
                        unvisited = Some(child);
                        break;
                    }
                }
                stack[top].1 = next;
                match unvisited {
                    Some(child) => {
                        visited[child] = k;
                        stack.push((child, first_child(child)));
                    }
                    None => {
                        stack.pop();
                        postorder.push(row);
                    }
                }
            }
        }

        // The sparse triangular solve, in reverse postorder (a topological order)
        for p in column {
            x[a.idx[p]] = a.data[p];
        }
        for &row in postorder.iter().rev() {
            if let Some(s) = pinv[row] {
                let x_row = x[row];
                for p in l_ptr[s] + 1..l_ptr[s + 1] {
                    x[l_idx[p]] = x[l_idx[p]] - l_val[p] * x_row;
                }
            }
        }

        // Pivoted rows go to U; the largest of the rest is a candidate pivot
        let mut largest = <T::Real as num_traits::Zero>::zero();
        let mut pivot_row = None;
        for &row in postorder.iter().rev() {
            match pinv[row] {
                Some(s) => {
                    u_idx.push(s);
                    u_val.push(x[row]);
                }
                None => {
                    let magnitude = x[row].abs();
                    if magnitude > largest {
                        largest = magnitude;
                        pivot_row = Some(row);
                    }
                }
            }
        }
        let Some(mut pivot_row) = pivot_row else {
            return Err(SparseError::SingularMatrix { pivot: col });
        };
        // Prefer the diagonal of the ordered matrix while it's large enough
        if prefer_diagonal
            && pinv[col].is_none()
            && visited[col] == k
            && x[col].abs() >= tau * largest
        {
            pivot_row = col;
        }
        let pivot = x[pivot_row];
        u_idx.push(k);
        u_val.push(pivot);
        u_ptr.push(u_idx.len());
        pinv[pivot_row] = Some(k);

        l_idx.push(pivot_row);
        l_val.push(T::one());
        for &row in postorder.iter().rev() {
            if pinv[row].is_none() {
                l_idx.push(row);
                l_val.push(x[row] / pivot);
            }
            x[row] = T::zero();
        }
        l_ptr.push(l_idx.len());
    }

    let mut row_perm = vec![0; n];
    for (row, step) in pinv.iter().enumerate() {
        row_perm[step.expect("every row pivots in a nonsingular factorization")] = row;
    }
    for row in &mut l_idx {
        *row = pinv[*row].expect("every row of L has pivoted");
    }
    sort_columns(&l_ptr, &mut l_idx, &mut l_val);
    sort_columns(&u_ptr, &mut u_idx, &mut u_val);
    Ok(Lu {
        l: CompressedMatrix::from_parts_unchecked((n, n), l_ptr, l_idx, l_val),
        u: CompressedMatrix::from_parts_unchecked((n, n), u_ptr, u_idx, u_val),
        row_perm_inverse: inverse_permutation(&row_perm),
        row_perm,
        col_perm: q.to_vec(),
        options,
        prefer_diagonal,
        a_ptr: a.ptr.clone(),
        a_idx: a.idx.clone(),
    })
}

/// Sorts each column's entries by row, keeping each value with its row.
fn sort_columns<T: Copy>(ptr: &[usize], idx: &mut [usize], val: &mut [T]) {
    let mut pairs: Vec<(usize, T)> = vec![];
    for bounds in ptr.windows(2) {
        let range = bounds[0]..bounds[1];
        pairs.clear();
        pairs.extend(
            idx[range.clone()]
                .iter()
                .copied()
                .zip(val[range.clone()].iter().copied()),
        );
        pairs.sort_unstable_by_key(|&(row, _)| row);
        for (offset, &(row, value)) in pairs.iter().enumerate() {
            idx[range.start + offset] = row;
            val[range.start + offset] = value;
        }
    }
}

/// Whether a permutation is odd.
fn is_odd(perm: &[usize]) -> bool {
    let mut seen = vec![false; perm.len()];
    let mut transpositions = 0;
    for start in 0..perm.len() {
        let mut node = start;
        let mut length = 0;
        while !seen[node] {
            seen[node] = true;
            node = perm[node];
            length += 1;
        }
        transpositions += length.max(1) - 1;
    }
    transpositions % 2 == 1
}

impl<T: Scalar + ComplexFloat> Lu<T> {
    /// The unit lower triangular factor `L`, including its unit diagonal.
    pub fn l(&self) -> &CscMatrix<T> {
        &self.l
    }

    /// The upper triangular factor `U`.
    pub fn u(&self) -> &CscMatrix<T> {
        &self.u
    }

    /// The row permutation `P`: row `i` of `P A Q` is row `permutation()[i]` of `A`.
    pub fn permutation(&self) -> &[usize] {
        &self.row_perm
    }

    /// The column permutation `Q`: column `j` of `P A Q` is column `column_permutation()[j]`
    /// of `A`.
    pub fn column_permutation(&self) -> &[usize] {
        &self.col_perm
    }

    /// The determinant of `A`: the product of `U`'s diagonal, negated once for each odd
    /// permutation among `P` and `Q`.
    pub fn determinant(&self) -> T {
        let n = self.row_perm.len();
        let negate = is_odd(&self.row_perm) != is_odd(&self.col_perm);
        let sign = if negate { -T::one() } else { T::one() };
        (0..n).fold(sign, |det, j| det * self.u.data[self.u.ptr[j + 1] - 1])
    }

    /// Solves `A x = b`: permutes `b` by `P`, substitutes forward with `L` and back with `U`,
    /// and permutes the result by `Q`. Returns `SparseError::ShapeMismatch` unless `b` has one
    /// entry per row of `A`.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        let mut x = b.to_vec();
        self.solve_in_place(&mut x)?;
        Ok(x)
    }

    /// Solves `A x = b`, overwriting `rhs` (which holds `b`) with `x`. It saves returning a new
    /// vector, though permuting still uses one scratch vector of the same length. The result is
    /// bit for bit `solve`'s. Returns `SparseError::ShapeMismatch`, leaving `rhs` untouched,
    /// unless `rhs` has one entry per row of `A`.
    pub fn solve_in_place(&self, rhs: &mut [T]) -> Result<(), SparseError> {
        self.check_rhs(rhs.len())?;
        let mut scratch = Vec::with_capacity(rhs.len());
        self.substitute(rhs, &mut scratch);
        Ok(())
    }

    /// Solves `A x = b` for each right-hand side `b` in `rhs`, returning the solutions in the
    /// same order. One scratch vector serves them all, and each solution is bit for bit what
    /// `solve` gives. Returns `SparseError::ShapeMismatch` naming the first right-hand side
    /// whose length isn't the number of rows of `A`.
    pub fn solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        for b in rhs {
            self.check_rhs(b.len())?;
        }
        let mut scratch = Vec::with_capacity(self.row_perm.len());
        Ok(rhs
            .iter()
            .map(|b| {
                let mut x = b.clone();
                self.substitute(&mut x, &mut scratch);
                x
            })
            .collect())
    }

    /// Solves `Aᵀ x = b` (the plain transpose, not conjugated, for a complex matrix) with the same
    /// factors: `Aᵀ = Q Uᵀ Lᵀ P`, so it substitutes forward with `Uᵀ` and back with `Lᵀ`,
    /// reading their stored columns as rows. Returns `SparseError::ShapeMismatch` unless `b`
    /// has one entry per row of `A`.
    pub fn solve_transpose(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.check_rhs(b.len())?;
        let n = self.row_perm.len();
        // Uᵀ y = Qᵀ b: column j of U is row j of Uᵀ, its diagonal stored last
        let mut y: Vec<T> = self.col_perm.iter().map(|&col| b[col]).collect();
        for j in 0..n {
            let diagonal = self.u.ptr[j + 1] - 1;
            let mut remainder = y[j];
            for p in self.u.ptr[j]..diagonal {
                remainder = remainder - self.u.data[p] * y[self.u.idx[p]];
            }
            y[j] = remainder / self.u.data[diagonal];
        }
        // Lᵀ z = y: column j of L is row j of Lᵀ, its unit diagonal stored first
        for j in (0..n).rev() {
            let mut remainder = y[j];
            for p in self.l.ptr[j] + 1..self.l.ptr[j + 1] {
                remainder = remainder - self.l.data[p] * y[self.l.idx[p]];
            }
            y[j] = remainder;
        }
        // P x = z
        let mut x = vec![T::zero(); n];
        for (&row, &value) in self.row_perm.iter().zip(&y) {
            x[row] = value;
        }
        Ok(x)
    }

    /// Solves `L y = P b` for a sparse `b`, touching only the part of `L` the solve needs: the
    /// entries of `y` that can be nonzero are those reachable from `b`'s entries in the graph
    /// of `L` (Gilbert 1994), and only those are computed. `y`'s indices are positions in the
    /// factor's row ordering (row `i` of `P A Q`), not rows of `A`.
    ///
    /// The values are bit for bit those of `solve`'s forward substitution, which subtracts by
    /// columns of `L` in the same order, apart from the sign of a zero: the columns skipped
    /// here would only subtract zeros. Stored entries may be zero. Returns
    /// `SparseError::ShapeMismatch` unless `b` has one entry per row of `A`.
    pub fn solve_lower_sparse(&self, b: &SparseVector<T>) -> Result<SparseVector<T>, SparseError> {
        self.check_rhs(b.len())?;
        let (pattern, y) = self.forward_sparse(b);
        let values = pattern.iter().map(|&i| y[i]).collect();
        Ok(SparseVector {
            len: b.len(),
            indices: pattern,
            values,
        })
    }

    /// Solves `A x = b` for a sparse `b`: the forward substitution touches only the part of `L`
    /// that [`Lu::solve_lower_sparse`] does, and the back substitution is `solve`'s, since `x`
    /// is generally dense. The result is bit for bit `solve`'s, apart from the sign of a zero.
    /// Returns
    /// `SparseError::ShapeMismatch` unless `b` has one entry per row of `A`.
    pub fn solve_sparse(&self, b: &SparseVector<T>) -> Result<Vec<T>, SparseError> {
        self.check_rhs(b.len())?;
        let (_, mut y) = self.forward_sparse(b);
        self.back_substitute(&mut y);
        let mut x = vec![T::zero(); y.len()];
        for (&col, &value) in self.col_perm.iter().zip(y.iter()) {
            x[col] = value;
        }
        Ok(x)
    }

    /// The reach of `b` (already checked) in increasing order, and `L⁻¹ P b` as a dense vector
    /// that is zero outside it. `L` has a unit diagonal, stored first in each column.
    fn forward_sparse(&self, b: &SparseVector<T>) -> (Vec<usize>, Vec<T>) {
        let n = self.row_perm.len();
        let pattern = reach(
            n,
            &self.l.ptr,
            &self.l.idx,
            b.indices().iter().map(|&row| self.row_perm_inverse[row]),
        );
        let mut y = vec![T::zero(); n];
        for (row, value) in b.iter() {
            y[self.row_perm_inverse[row]] = value;
        }
        for &j in &pattern {
            let y_j = y[j];
            for p in self.l.ptr[j] + 1..self.l.ptr[j + 1] {
                y[self.l.idx[p]] = y[self.l.idx[p]] - self.l.data[p] * y_j;
            }
        }
        (pattern, y)
    }

    fn check_rhs(&self, len: usize) -> Result<(), SparseError> {
        let n = self.row_perm.len();
        if len == n {
            Ok(())
        } else {
            Err(SparseError::ShapeMismatch {
                left: (n, n),
                right: (len, 1),
            })
        }
    }

    /// Overwrites `x` (holding `b`, already checked) with the solution: `P b` into `scratch`,
    /// forward substitution with `L` (unit diagonal, stored first in each column) and back
    /// substitution with `U` (diagonal stored last), by columns, then through `Q` into `x`.
    fn substitute(&self, x: &mut [T], scratch: &mut Vec<T>) {
        let n = self.row_perm.len();
        scratch.clear();
        scratch.extend(self.row_perm.iter().map(|&row| x[row]));
        let y = scratch.as_mut_slice();
        for j in 0..n {
            let y_j = y[j];
            for p in self.l.ptr[j] + 1..self.l.ptr[j + 1] {
                y[self.l.idx[p]] = y[self.l.idx[p]] - self.l.data[p] * y_j;
            }
        }
        self.back_substitute(y);
        for (&col, &value) in self.col_perm.iter().zip(y.iter()) {
            x[col] = value;
        }
    }

    /// Overwrites `y` with `U⁻¹ y`: back substitution with `U` by columns, its diagonal stored
    /// last.
    fn back_substitute(&self, y: &mut [T]) {
        for j in (0..y.len()).rev() {
            let diagonal = self.u.ptr[j + 1] - 1;
            y[j] = y[j] / self.u.data[diagonal];
            let y_j = y[j];
            for p in self.u.ptr[j]..diagonal {
                y[self.u.idx[p]] = y[self.u.idx[p]] - self.u.data[p] * y_j;
            }
        }
    }

    /// Factors `matrix`, which must have the pattern of the matrix this factorization came
    /// from, reusing its column ordering, pivot sequence and the patterns of `L` and `U`, so no
    /// search is needed: each column is a sparse triangular solve over a known pattern.
    ///
    /// Each reused pivot is checked against the threshold this factorization used. If one is
    /// too small (or zero), the matrix is factored afresh with new pivots under the same column
    /// ordering, and the outcome is `RefactorOutcome::Repivoted`; otherwise it is
    /// `RefactorOutcome::Reused`.
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs, or
    /// `SparseError::SingularMatrix` if the fresh factorization finds no pivot; on an error,
    /// `self` is unchanged.
    pub fn refactor(&mut self, matrix: &CscMatrix<T>) -> Result<RefactorOutcome, SparseError> {
        let n = self.row_perm.len();
        if matrix.shape != (n, n) || matrix.ptr != self.a_ptr || matrix.idx != self.a_idx {
            return Err(SparseError::PatternMismatch);
        }
        let tau = threshold::<T>(self.options);
        let mut step_of = vec![0; n];
        for (step, &row) in self.row_perm.iter().enumerate() {
            step_of[row] = step;
        }
        let mut l_val = vec![T::zero(); self.l.data.len()];
        let mut u_val = vec![T::zero(); self.u.data.len()];
        let mut x = vec![T::zero(); n];
        for k in 0..n {
            let col = self.col_perm[k];
            for p in matrix.ptr[col]..matrix.ptr[col + 1] {
                x[step_of[matrix.idx[p]]] = matrix.data[p];
            }
            // U's column in increasing row order is a topological order, since L is lower
            //   triangular in the pivoted numbering
            let diagonal = self.u.ptr[k + 1] - 1;
            let above = self.u.ptr[k]..diagonal;
            for (&j, u_jk_slot) in self.u.idx[above.clone()].iter().zip(&mut u_val[above]) {
                let u_jk = x[j];
                x[j] = T::zero();
                *u_jk_slot = u_jk;
                let column_j = self.l.ptr[j] + 1..self.l.ptr[j + 1];
                for (&row, &l_rj) in self.l.idx[column_j.clone()].iter().zip(&l_val[column_j]) {
                    x[row] = x[row] - l_rj * u_jk;
                }
            }
            let pivot = x[k];
            let below = self.l.ptr[k] + 1..self.l.ptr[k + 1];
            let largest = self.l.idx[below.clone()]
                .iter()
                .map(|&row| x[row].abs())
                .fold(pivot.abs(), |a, b| if b > a { b } else { a });
            if pivot == T::zero() || pivot.abs() < tau * largest {
                // Too small to reuse: factor afresh, leaving self as it was if that fails
                *self =
                    gilbert_peierls(matrix, &self.col_perm, self.options, self.prefer_diagonal)?;
                return Ok(RefactorOutcome::Repivoted);
            }
            u_val[diagonal] = pivot;
            x[k] = T::zero();
            l_val[self.l.ptr[k]] = T::one();
            for p in below {
                let row = self.l.idx[p];
                l_val[p] = x[row] / pivot;
                x[row] = T::zero();
            }
        }
        self.l.data = l_val;
        self.u.data = u_val;
        Ok(RefactorOutcome::Reused)
    }
}
