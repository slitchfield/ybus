//! Sparse LU factorization: Gilbert and Peierls's left-looking algorithm with threshold partial
//! pivoting, for real and complex matrices.
//!
//! [`CscMatrix::analyze_lu`] chooses a column ordering `Q` from the pattern of `A + Aᵀ`.
//! [`SymbolicLu::factor`] then computes `P A Q = L U` one column at a time: a depth-first search
//! finds which entries of the column can be nonzero, a sparse triangular solve computes them,
//! and a pivot is chosen among the rows not yet pivoted, preferring the diagonal (which keeps the
//! ordering's sparsity) while it is large enough.

use num_complex::ComplexFloat;
use num_traits::NumCast;

use super::{CompressedMatrix, CscMatrix, Ordering};
use crate::sparse_error::SparseError;
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
    a_ptr: Vec<usize>,
    a_idx: Vec<usize>,
}

/// An LU factorization `P A Q = L U` with threshold partial pivoting, from
/// [`SymbolicLu::factor`] or [`CscMatrix::lu`]. Keep it to solve many right-hand sides against
/// the same matrix.
///
/// `L` is unit lower triangular and `U` upper triangular, both in CSC form, and each stores
/// its whole symbolic pattern, so some entries may be zero. Row `i` of `P A Q` is row
/// `permutation()[i]` of `A`, and column `j` is column `column_permutation()[j]`.
#[derive(Clone, Debug)]
pub struct Lu<T: Scalar = f64> {
    l: CscMatrix<T>,
    u: CscMatrix<T>,
    row_perm: Vec<usize>,
    col_perm: Vec<usize>,
}

impl<T: Scalar> CscMatrix<T> {
    /// Analyzes this matrix's pattern for an LU factorization whose columns, and preferred
    /// pivots, follow `ordering` applied to the pattern of `A + Aᵀ`. Only the pattern is read.
    /// Returns `SparseError::NotSquare` for a non-square matrix.
    pub fn analyze_lu(&self, ordering: Ordering) -> Result<SymbolicLu, SparseError> {
        Ok(SymbolicLu {
            col_perm: self.ordering(ordering)?,
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
        gilbert_peierls(matrix, &self.col_perm, options)
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
        if pinv[col].is_none() && visited[col] == k && x[col].abs() >= tau * largest {
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
        row_perm,
        col_perm: q.to_vec(),
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
        let n = self.row_perm.len();
        if b.len() != n {
            return Err(SparseError::ShapeMismatch {
                left: (n, n),
                right: (b.len(), 1),
            });
        }
        let mut y: Vec<T> = self.row_perm.iter().map(|&row| b[row]).collect();
        // L has a unit diagonal, stored first in each column
        for j in 0..n {
            let y_j = y[j];
            for p in self.l.ptr[j] + 1..self.l.ptr[j + 1] {
                y[self.l.idx[p]] = y[self.l.idx[p]] - self.l.data[p] * y_j;
            }
        }
        // U's diagonal is stored last in each column
        for j in (0..n).rev() {
            let diagonal = self.u.ptr[j + 1] - 1;
            y[j] = y[j] / self.u.data[diagonal];
            let y_j = y[j];
            for p in self.u.ptr[j]..diagonal {
                y[self.u.idx[p]] = y[self.u.idx[p]] - self.u.data[p] * y_j;
            }
        }
        let mut x = vec![T::zero(); n];
        for (&col, value) in self.col_perm.iter().zip(y) {
            x[col] = value;
        }
        Ok(x)
    }
}
