//! Solvers for linear systems `A x = b`: triangular solves, LU with partial pivoting,
//! Cholesky for symmetric positive definite matrices, and conjugate gradient.
//!
//! All of them need a float element type, since they divide. The factorizations work by rows,
//! so they're on `CsrMatrix`; the triangular solves and conjugate gradient work in either
//! orientation.

use std::collections::BTreeMap;

use num_traits::Float;

use super::{CompressedMatrix, CsrMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// An LU factorization with partial pivoting, `P A = L U`, from `CsrMatrix::lu`.
///
/// `L` is unit lower triangular and `U` upper triangular, both in CSR form. Keep it around to
/// solve many right-hand sides against the same matrix.
#[derive(Clone, Debug)]
pub struct Lu<T: Scalar = f64> {
    l: CsrMatrix<T>,
    u: CsrMatrix<T>,
    perm: Vec<usize>,
    swaps: usize,
}

impl<T: Scalar + Float> Lu<T> {
    /// The unit lower triangular factor `L`, including its unit diagonal.
    pub fn l(&self) -> &CsrMatrix<T> {
        &self.l
    }

    /// The upper triangular factor `U`.
    pub fn u(&self) -> &CsrMatrix<T> {
        &self.u
    }

    /// The row permutation: row `i` of `P A` is row `permutation()[i]` of `A`.
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// The determinant of the original matrix: the product of `U`'s diagonal, negated once
    /// per row swap.
    pub fn determinant(&self) -> T {
        let sign = if self.swaps.is_multiple_of(2) {
            T::one()
        } else {
            -T::one()
        };
        (0..self.perm.len()).fold(sign, |det, i| {
            det * self.u.stored_value(i, i).unwrap_or(T::zero())
        })
    }

    /// Solves `A x = b` for `x`. Returns `SparseError::ShapeMismatch` unless `b` has one
    /// entry per row of `A`.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        if b.len() != self.perm.len() {
            return Err(SparseError::ShapeMismatch {
                left: self.l.shape,
                right: (b.len(), 1),
            });
        }
        let permuted: Vec<T> = self.perm.iter().map(|&row| b[row]).collect();
        let y = self.l.solve_lower_triangular(&permuted)?;
        self.u.solve_upper_triangular(&y)
    }
}

impl<T: Scalar + Float, O: Orientation> CompressedMatrix<T, O> {
    /// Checks that `self` is square and `b_len` matches, returning the size.
    fn check_system(&self, b_len: usize) -> Result<usize, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        if b_len != self.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: (b_len, 1),
            });
        }
        Ok(self.shape.0)
    }

    /// Solves `L x = b` by forward substitution, reading only the lower triangle of `self`
    /// (entries above the diagonal are ignored). CSR substitutes row by row and CSC column by
    /// column; both subtract each `x[i]`'s terms in the same order, so they agree bit for bit.
    ///
    /// Returns `SparseError::SingularMatrix` if a diagonal entry is zero or not stored, naming
    /// the first one.
    pub fn solve_lower_triangular(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        let n = self.check_system(b.len())?;
        if O::ROW_MAJOR {
            let mut x = vec![T::zero(); n];
            for i in 0..n {
                let mut remainder = b[i];
                let mut diagonal = T::zero();
                for p in self.ptr[i]..self.ptr[i + 1] {
                    let col = self.idx[p];
                    if col < i {
                        remainder = remainder - self.data[p] * x[col];
                    } else {
                        if col == i {
                            diagonal = self.data[p];
                        }
                        break;
                    }
                }
                if diagonal == T::zero() {
                    return Err(SparseError::SingularMatrix { pivot: i });
                }
                x[i] = remainder / diagonal;
            }
            Ok(x)
        } else {
            let mut x = b.to_vec();
            for j in 0..n {
                let column = self.ptr[j]..self.ptr[j + 1];
                let below = column.start + self.idx[column.clone()].partition_point(|&row| row < j);
                if below == column.end || self.idx[below] != j || self.data[below] == T::zero() {
                    return Err(SparseError::SingularMatrix { pivot: j });
                }
                x[j] = x[j] / self.data[below];
                for p in below + 1..column.end {
                    x[self.idx[p]] = x[self.idx[p]] - self.data[p] * x[j];
                }
            }
            Ok(x)
        }
    }

    /// Solves `U x = b` by back substitution, reading only the upper triangle of `self`
    /// (entries below the diagonal are ignored). CSR substitutes row by row and CSC column by
    /// column, which subtract terms in different orders, so the two can differ in the last
    /// bits.
    ///
    /// Returns `SparseError::SingularMatrix` if a diagonal entry is zero or not stored, naming
    /// the last one.
    pub fn solve_upper_triangular(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        let n = self.check_system(b.len())?;
        if O::ROW_MAJOR {
            let mut x = vec![T::zero(); n];
            for i in (0..n).rev() {
                let mut remainder = b[i];
                let mut diagonal = T::zero();
                for p in self.ptr[i]..self.ptr[i + 1] {
                    let col = self.idx[p];
                    if col == i {
                        diagonal = self.data[p];
                    } else if col > i {
                        remainder = remainder - self.data[p] * x[col];
                    }
                }
                if diagonal == T::zero() {
                    return Err(SparseError::SingularMatrix { pivot: i });
                }
                x[i] = remainder / diagonal;
            }
            Ok(x)
        } else {
            let mut x = b.to_vec();
            for j in (0..n).rev() {
                let column = self.ptr[j]..self.ptr[j + 1];
                let at = column.start + self.idx[column.clone()].partition_point(|&row| row < j);
                if at == column.end || self.idx[at] != j || self.data[at] == T::zero() {
                    return Err(SparseError::SingularMatrix { pivot: j });
                }
                x[j] = x[j] / self.data[at];
                for p in column.start..at {
                    x[self.idx[p]] = x[self.idx[p]] - self.data[p] * x[j];
                }
            }
            Ok(x)
        }
    }

    /// Solves `A x = b` for a symmetric positive definite matrix by the conjugate gradient
    /// method, starting from zero.
    ///
    /// Stops when the residual norm is at most `tolerance * ||b||`. Each iteration costs one
    /// matrix-vector product, so this suits large systems where factoring would fill in.
    /// Symmetry is assumed, not checked. Returns `SparseError::NoConvergence` if the
    /// tolerance is not reached within `max_iterations`, or if the search direction shows the
    /// matrix is not positive definite.
    pub fn conjugate_gradient(
        &self,
        b: &[T],
        tolerance: T,
        max_iterations: usize,
    ) -> Result<Vec<T>, SparseError> {
        let n = self.check_system(b.len())?;
        let dot = |u: &[T], v: &[T]| u.iter().zip(v).fold(T::zero(), |sum, (&a, &b)| sum + a * b);

        let b_norm = dot(b, b).sqrt();
        let mut x = vec![T::zero(); n];
        if b_norm == T::zero() {
            return Ok(x);
        }
        let mut residual = b.to_vec();
        let mut direction = residual.clone();
        let mut residual_sq = dot(&residual, &residual);

        for iteration in 0..max_iterations {
            if residual_sq.sqrt() <= tolerance * b_norm {
                return Ok(x);
            }
            let a_direction = self.mul_vec(&direction)?;
            let curvature = dot(&direction, &a_direction);
            if curvature.is_nan() || curvature <= T::zero() {
                return Err(SparseError::NoConvergence {
                    iterations: iteration,
                });
            }
            let step = residual_sq / curvature;
            for i in 0..n {
                x[i] = x[i] + step * direction[i];
                residual[i] = residual[i] - step * a_direction[i];
            }
            let next_residual_sq = dot(&residual, &residual);
            let beta = next_residual_sq / residual_sq;
            for i in 0..n {
                direction[i] = residual[i] + beta * direction[i];
            }
            residual_sq = next_residual_sq;
        }

        if residual_sq.sqrt() <= tolerance * b_norm {
            Ok(x)
        } else {
            Err(SparseError::NoConvergence {
                iterations: max_iterations,
            })
        }
    }
}

impl<T: Scalar + Float> CsrMatrix<T> {
    /// Factors a square matrix as `P A = L U` with partial (row) pivoting, returning `L` and
    /// `U` as CSR matrices that store no zeros (except `L`'s unit diagonal is always stored).
    ///
    /// Rows are kept as ordered sparse maps, so only the entries that fill in are ever
    /// touched; the pivot search still looks at every remaining row for each column, so
    /// the cost is at least quadratic in the dimension.
    ///
    /// Returns `SparseError::NotSquare`, or `SparseError::SingularMatrix` if a column has no
    /// usable pivot.
    pub fn lu(&self) -> Result<Lu<T>, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        let n = self.shape.0;
        let mut rows: Vec<BTreeMap<usize, T>> = self
            .ptr
            .windows(2)
            .map(|bounds| {
                let range = bounds[0]..bounds[1];
                self.idx[range.clone()]
                    .iter()
                    .copied()
                    .zip(self.data[range].iter().copied())
                    .collect()
            })
            .collect();
        let mut perm: Vec<usize> = (0..n).collect();
        let mut swaps = 0;

        for k in 0..n {
            let magnitude = |row: &BTreeMap<usize, T>| row.get(&k).map_or(T::zero(), |v| v.abs());
            let mut best = k;
            let mut best_magnitude = magnitude(&rows[k]);
            for (i, row) in rows.iter().enumerate().skip(k + 1) {
                if magnitude(row) > best_magnitude {
                    best = i;
                    best_magnitude = magnitude(row);
                }
            }
            if best_magnitude == T::zero() {
                return Err(SparseError::SingularMatrix { pivot: k });
            }
            if best != k {
                rows.swap(k, best);
                perm.swap(k, best);
                swaps += 1;
            }

            let pivot = rows[k][&k];
            let pivot_row: Vec<(usize, T)> =
                rows[k].range(k + 1..).map(|(&j, &v)| (j, v)).collect();
            for row in rows.iter_mut().skip(k + 1) {
                let Some(&below) = row.get(&k) else {
                    continue;
                };
                let multiplier = below / pivot;
                if multiplier == T::zero() {
                    row.remove(&k);
                    continue;
                }
                row.insert(k, multiplier);
                for &(j, u_kj) in &pivot_row {
                    let current = row.get(&j).copied().unwrap_or(T::zero());
                    let updated = current - multiplier * u_kj;
                    if updated == T::zero() {
                        row.remove(&j);
                    } else {
                        row.insert(j, updated);
                    }
                }
            }
        }

        // Each row splits into L's part (left of the diagonal, plus the unit diagonal) and
        //   U's; the maps iterate in column order, so the rows come out sorted
        let (mut l_ptr, mut l_idx, mut l_data) = (vec![0], vec![], vec![]);
        let (mut u_ptr, mut u_idx, mut u_data) = (vec![0], vec![], vec![]);
        for (i, row) in rows.iter().enumerate() {
            for (&j, &value) in row.range(..i).filter(|(_, &value)| value != T::zero()) {
                l_idx.push(j);
                l_data.push(value);
            }
            l_idx.push(i);
            l_data.push(T::one());
            l_ptr.push(l_idx.len());
            for (&j, &value) in row.range(i..).filter(|(_, &value)| value != T::zero()) {
                u_idx.push(j);
                u_data.push(value);
            }
            u_ptr.push(u_idx.len());
        }
        Ok(Lu {
            l: CsrMatrix::from_parts_unchecked((n, n), l_ptr, l_idx, l_data),
            u: CsrMatrix::from_parts_unchecked((n, n), u_ptr, u_idx, u_data),
            perm,
            swaps,
        })
    }

    /// Solves `A x = b` for a general square matrix by LU factorization. Factor once with
    /// `lu` instead when solving several right-hand sides.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.check_system(b.len())?;
        self.lu()?.solve(b)
    }

    /// Computes the Cholesky factor `L` (lower triangular, in CSR form) of a symmetric positive
    /// definite matrix, so that `A = L Lᵀ`. `L` stores no zeros.
    ///
    /// Works row by row, solving a sparse triangular system for each, so cost follows the
    /// fill-in of `L`; there's no fill-reducing ordering yet. Returns `SparseError::NotSquare`,
    /// `SparseError::NotSymmetric` (judged mathematically, ignoring stored zeros), or
    /// `SparseError::NotPositiveDefinite` if a pivot is not strictly positive.
    pub fn cholesky(&self) -> Result<CsrMatrix<T>, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        if *self != self.clone().transpose() {
            return Err(SparseError::NotSymmetric);
        }
        let n = self.shape.0;
        let mut diagonal = vec![T::zero(); n];
        // Entries of L strictly below the diagonal, grouped by column: (row, value)
        let mut columns: Vec<Vec<(usize, T)>> = vec![vec![]; n];
        let (mut ptr, mut idx, mut data) = (vec![0], vec![], vec![]);

        for i in 0..n {
            // Row i of A (lower part), which is turned into row i of L one column at a time,
            //   smallest column first, since each entry depends on the smaller ones.
            let mut work: BTreeMap<usize, T> = BTreeMap::new();
            let mut a_ii = T::zero();
            for p in self.ptr[i]..self.ptr[i + 1] {
                let col = self.idx[p];
                if col < i {
                    work.insert(col, self.data[p]);
                } else {
                    if col == i {
                        a_ii = self.data[p];
                    }
                    break;
                }
            }

            let mut remaining = a_ii;
            while let Some((j, w_j)) = work.pop_first() {
                let l_ij = w_j / diagonal[j];
                if l_ij == T::zero() {
                    continue;
                }
                for &(m, l_mj) in &columns[j] {
                    let current = work.get(&m).copied().unwrap_or(T::zero());
                    work.insert(m, current - l_mj * l_ij);
                }
                columns[j].push((i, l_ij));
                // Entries leave the work map in increasing column order, so the row is sorted
                idx.push(j);
                data.push(l_ij);
                remaining = remaining - l_ij * l_ij;
            }

            if remaining.is_nan() || remaining <= T::zero() {
                return Err(SparseError::NotPositiveDefinite { pivot: i });
            }
            diagonal[i] = remaining.sqrt();
            idx.push(i);
            data.push(diagonal[i]);
            ptr.push(idx.len());
        }
        Ok(CsrMatrix::from_parts_unchecked((n, n), ptr, idx, data))
    }

    /// Solves `A x = b` for a symmetric positive definite matrix by Cholesky factorization:
    /// forward substitution with `L`, then back substitution with `Lᵀ` by rows. Fails with the
    /// errors of `cholesky`.
    pub fn solve_spd(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.check_system(b.len())?;
        let l = self.cholesky()?;
        let y = l.solve_lower_triangular(b)?;
        // Lᵀ in CSR form, reoriented in O(nnz); back substitution by rows adds each entry's
        //   terms in column order, as the solve always has
        let l_transposed: CsrMatrix<T> = l.to_csc().transpose();
        l_transposed.solve_upper_triangular(&y)
    }
}
