//! Ordered sparse LDLᵀ with static pivots, on the Cholesky factorization's symbolic
//! analysis: [`Ldlt`].

use num_complex::ComplexFloat;
use num_traits::{Float, NumCast};

use super::cholesky::{back_by_rows, forward_by_rows, SymbolicCholesky};
use super::{CompressedMatrix, CsrMatrix, Ordering};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// A symmetric matrix's inertia: how many of its eigenvalues are positive and how many negative.
/// An [`Ldlt`] factorization has no zero eigenvalue, since a zero pivot fails it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Inertia {
    /// The positive eigenvalues
    pub positive: usize,
    /// The negative eigenvalues
    pub negative: usize,
}

/// An ordered LDLᵀ factorization `P A Pᵀ = L D Lᴴ` of a symmetric (for a complex element type,
/// Hermitian) matrix with nonzero static pivots, from [`SymbolicCholesky::factor_ldlt`] or
/// [`CsrMatrix::ldlt`]. For a real matrix, `Lᴴ` is `Lᵀ`. Keep it to solve many right-hand
/// sides against the same matrix.
///
/// `L` is unit lower triangular and `D` diagonal and real. The pivots are the diagonal entries in the order the analysis chose, with no pivoting during
/// the factorization, so the factor's pattern is exactly the Cholesky factor's
/// ([`SymbolicCholesky`] serves both) and its cost is known in advance. It is the up-looking
/// LDLᵀ of Davis's LDL package (T. A. Davis, "Algorithm 849: a concise sparse Cholesky
/// factorization package," *ACM Trans. Math. Softw.* 31(4):587–591, 2005,
/// doi:10.1145/1114268.1114277).
///
/// Static pivots suit a symmetric quasi-definite matrix, `[H Bᴴ; B −G]` with `H` and `G`
/// positive definite. Such a matrix has an LDLᵀ factorization under every symmetric ordering,
/// with `D` positive on `H`'s rows and negative on `G`'s (R. J. Vanderbei, "Symmetric
/// quasidefinite matrices," *SIAM J. Optim.* 5(1):100–113, 1995, doi:10.1137/0805005). An
/// interior-point method's KKT matrix becomes one when regularized, which is the caller's
/// choice: this factorization doesn't perturb small pivots. A matrix that needs pivoting, such
/// as `[[0, 1], [1, 0]]`, fails with a zero pivot.
///
/// Static pivots lose accuracy as the regularization shrinks, since a small pivot early in the
/// order makes large entries of `L`. On random quasi-definite matrices regularized by `δ`, the
/// scaled residual grew as 1/δ: about 3e-10 at δ = 1e-8, and near machine precision from 1e-2
/// up (see the crate's `ldlt_test`). Iterative refinement against the factorization recovers
/// the lost accuracy where it matters.
///
/// By Sylvester's law of inertia, `D`'s signs are `A`'s eigenvalues' signs:
/// [`Ldlt::inertia`].
///
/// It holds `L` in CSR form, and also `Lᴴ` in CSR form for the back substitution, like
/// [`Cholesky`](super::Cholesky).
#[derive(Clone, Debug)]
pub struct Ldlt<T: Scalar = f64> {
    l: CsrMatrix<T>,
    l_adjoint: CsrMatrix<T>,
    /// `D`'s entries, real (with zero imaginary part for a complex element type)
    d: Vec<T>,
    perm: Vec<usize>,
}

impl SymbolicCholesky {
    /// Factors `matrix`, which must have the analyzed pattern and be symmetric (for a complex
    /// element type, Hermitian), as `P A Pᵀ = L D Lᴴ`, `L` unit lower triangular and `D`
    /// diagonal and real, with the pivots in the analysis's order and no pivoting.
    ///
    /// The up-looking algorithm, on the Cholesky factor's pattern: row `i` solves a sparse
    /// triangular system with the rows above it, giving `z_j = L(i, j) d_j`, then
    /// `L(i, j) = z_j / d_j`, and `d_i = a_ii − Σ_j L(i, j) conj(z_j)`. `L` stores the whole
    /// symbolic pattern, so an entry that happens to be zero stays stored.
    ///
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, `SparseError::NotSymmetric` if `A ≠ Aᴴ` (judged mathematically, ignoring
    /// stored zeros), or `SparseError::SingularMatrix` if a pivot is zero or NaN, naming its row
    /// of `A` (in the original numbering).
    pub fn factor_ldlt<T: Scalar + ComplexFloat>(
        &self,
        matrix: &CsrMatrix<T>,
    ) -> Result<Ldlt<T>, SparseError> {
        let n = self.dim();
        if matrix.shape != (n, n) || matrix.ptr != self.a_ptr || matrix.idx != self.a_idx {
            return Err(SparseError::PatternMismatch);
        }
        if !matrix.is_hermitian() {
            return Err(SparseError::NotSymmetric);
        }
        let c = matrix.permute_symmetric(&self.perm)?;

        // L's columns below the diagonal, filled as the rows are computed, as in the Cholesky's
        //   factorization: column j's entries from rows already done are
        //   col_rows/col_values[col_start[j]..col_fill[j]]
        let mut col_start = vec![0; n + 1];
        for &j in &self.l_idx {
            col_start[j + 1] += 1;
        }
        for j in 0..n {
            col_start[j + 1] += col_start[j];
        }
        for (j, start) in col_start.iter_mut().enumerate() {
            *start -= j;
        }
        let mut col_fill = col_start.clone();
        let below = self.l_idx.len() - n;
        let mut col_rows = vec![0; below];
        let mut col_values = vec![T::zero(); below];

        let mut l_values = vec![T::zero(); self.l_idx.len()];
        let mut d = vec![T::zero(); n];
        let mut work = vec![T::zero(); n];
        for i in 0..n {
            // Row i of P A Pᵀ: the part left of the diagonal into the work vector, and the
            //   (real) diagonal
            let mut remaining = <T::Real as num_traits::Zero>::zero();
            for p in c.ptr[i]..c.ptr[i + 1] {
                let col = c.idx[p];
                if col < i {
                    work[col] = c.data[p];
                } else {
                    if col == i {
                        remaining = c.data[p].re();
                    }
                    break;
                }
            }
            let diagonal_slot = self.l_ptr[i + 1] - 1;
            let row = self.l_ptr[i]..diagonal_slot;
            for (&j, l_value) in self.l_idx[row.clone()].iter().zip(&mut l_values[row]) {
                // work[j] is now a_ij less the earlier columns' parts: z_j = L(i, j) d_j
                let z_j = work[j];
                work[j] = T::zero();
                // a_im = Σ_j L(i, j) d_j conj(L(m, j)) = Σ_j z_j conj(L(m, j))
                for q in col_start[j]..col_fill[j] {
                    let m = col_rows[q];
                    work[m] = work[m] - col_values[q].conj() * z_j;
                }
                let l_ij = z_j / d[j];
                col_rows[col_fill[j]] = i;
                col_values[col_fill[j]] = l_ij;
                col_fill[j] += 1;
                *l_value = l_ij;
                remaining = remaining - (l_ij * z_j.conj()).re();
            }
            if Float::is_nan(remaining) || remaining == <T::Real as num_traits::Zero>::zero() {
                return Err(SparseError::SingularMatrix {
                    pivot: self.perm[i],
                });
            }
            d[i] = <T as NumCast>::from(remaining)
                .expect("a real value converts to its own element type");
            l_values[diagonal_slot] = T::one();
        }

        // Lᴴ in CSR form is L in CSC form, conjugated: each column's unit diagonal, then the
        //   entries below it
        let mut lt_ptr = Vec::with_capacity(n + 1);
        lt_ptr.push(0);
        let mut lt_idx = Vec::with_capacity(self.l_idx.len());
        let mut lt_values = Vec::with_capacity(self.l_idx.len());
        for j in 0..n {
            lt_idx.push(j);
            lt_values.push(T::one());
            lt_idx.extend_from_slice(&col_rows[col_start[j]..col_fill[j]]);
            lt_values.extend(
                col_values[col_start[j]..col_fill[j]]
                    .iter()
                    .map(|value| value.conj()),
            );
            lt_ptr.push(lt_idx.len());
        }
        Ok(Ldlt {
            l: CompressedMatrix::from_parts_unchecked(
                (n, n),
                self.l_ptr.clone(),
                self.l_idx.clone(),
                l_values,
            ),
            l_adjoint: CompressedMatrix::from_parts_unchecked((n, n), lt_ptr, lt_idx, lt_values),
            d,
            perm: self.perm.clone(),
        })
    }
}

impl<T: Scalar + ComplexFloat> Ldlt<T> {
    /// The unit lower triangular factor `L` of `P A Pᵀ`, with its unit diagonal stored. It
    /// stores its whole symbolic pattern, so some entries may be zero.
    pub fn l(&self) -> &CsrMatrix<T> {
        &self.l
    }

    /// `D`'s entries, one per row of `P A Pᵀ`: real, and nonzero (with zero imaginary part for a
    /// complex element type).
    pub fn d(&self) -> &[T] {
        &self.d
    }

    /// The ordering `P`: row and column `i` of `P A Pᵀ` are row and column `permutation()[i]`
    /// of `A`.
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// `A`'s inertia, from `D`'s signs (Sylvester's law of inertia).
    pub fn inertia(&self) -> Inertia {
        let positive = self
            .d
            .iter()
            .filter(|d| d.re() > <T::Real as num_traits::Zero>::zero())
            .count();
        Inertia {
            positive,
            negative: self.d.len() - positive,
        }
    }

    /// Solves `A x = b`: permutes `b`, substitutes forward with `L`, divides by `D`, substitutes
    /// back with `Lᴴ`, and permutes the result back. Returns `SparseError::ShapeMismatch` unless
    /// `b` has one entry per row of `A`.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        let mut x = b.to_vec();
        self.solve_in_place(&mut x)?;
        Ok(x)
    }

    /// Solves `A x = b`, overwriting `rhs` (which holds `b`) with `x`, bit for bit `solve`'s.
    /// Returns `SparseError::ShapeMismatch`, leaving `rhs` untouched, unless `rhs` has one entry
    /// per row of `A`.
    pub fn solve_in_place(&self, rhs: &mut [T]) -> Result<(), SparseError> {
        self.check_rhs(rhs.len())?;
        let mut scratch = Vec::with_capacity(rhs.len());
        self.substitute(rhs, &mut scratch);
        Ok(())
    }

    /// Solves `A x = b` for each right-hand side `b` in `rhs`, returning the solutions in the
    /// same order, each bit for bit what `solve` gives. Returns `SparseError::ShapeMismatch`
    /// naming the first right-hand side whose length isn't the number of rows of `A`.
    pub fn solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        for b in rhs {
            self.check_rhs(b.len())?;
        }
        let mut scratch = Vec::with_capacity(self.perm.len());
        Ok(rhs
            .iter()
            .map(|b| {
                let mut x = b.clone();
                self.substitute(&mut x, &mut scratch);
                x
            })
            .collect())
    }

    fn check_rhs(&self, len: usize) -> Result<(), SparseError> {
        if len == self.perm.len() {
            Ok(())
        } else {
            Err(SparseError::ShapeMismatch {
                left: self.l.shape,
                right: (len, 1),
            })
        }
    }

    /// Overwrites `x` (holding `b`, already checked) with the solution: `P b` into `scratch`,
    /// forward substitution with `L`, division by `D`, back substitution with `Lᴴ`, then back
    /// through `P`.
    fn substitute(&self, x: &mut [T], scratch: &mut Vec<T>) {
        scratch.clear();
        scratch.extend(self.perm.iter().map(|&row| x[row]));
        let y = scratch.as_mut_slice();
        forward_by_rows(&self.l, y);
        for (y, &d) in y.iter_mut().zip(&self.d) {
            *y = *y / d;
        }
        back_by_rows(&self.l_adjoint, y);
        for (&row, &value) in self.perm.iter().zip(y.iter()) {
            x[row] = value;
        }
    }
}

#[cfg(feature = "rayon")]
impl<T: Scalar + ComplexFloat + Send + Sync> Ldlt<T> {
    /// [`Ldlt::solve_many`] with the right-hand sides solved in parallel (feature `rayon`). Each
    /// is solved by one thread, with its own scratch vector and exactly the operations `solve`
    /// uses, so the results are bit for bit `solve_many`'s, in the same order. Returns
    /// `SparseError::ShapeMismatch` naming the first right-hand side whose length isn't the
    /// number of rows of `A`, before solving any.
    pub fn par_solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        use rayon::prelude::*;
        for b in rhs {
            self.check_rhs(b.len())?;
        }
        let n = self.perm.len();
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

impl<T: Scalar + ComplexFloat> CsrMatrix<T> {
    /// Factors a symmetric (for a complex element type, Hermitian) matrix as
    /// `P A Pᵀ = L D Lᴴ` with static pivots, ordered by approximate minimum degree:
    /// `analyze_cholesky(Ordering::Amd)` followed by
    /// [`factor_ldlt`](SymbolicCholesky::factor_ldlt). The same analysis serves the Cholesky
    /// and LDLᵀ factorizations, so analyze once to factor several matrices with the same
    /// pattern, either way.
    ///
    /// Returns `SparseError::NotSquare`, `SparseError::NotSymmetric`, or
    /// `SparseError::SingularMatrix` naming the row of `A` whose pivot was zero.
    pub fn ldlt(&self) -> Result<Ldlt<T>, SparseError> {
        self.analyze_cholesky(Ordering::Amd)?.factor_ldlt(self)
    }
}
