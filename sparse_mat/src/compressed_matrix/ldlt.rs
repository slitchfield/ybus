//! Ordered sparse LDLᵀ with static pivots, on the Cholesky factorization's symbolic
//! analysis: [`Ldlt`].

use num_complex::ComplexFloat;
use num_traits::{Float, NumCast};

use super::cholesky::{back_by_rows, forward_by_rows, SymbolicCholesky};
use super::factorization::sealed::Substitute;
use super::factorization::Factorization;
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
/// choice: [`factor_ldlt`](SymbolicCholesky::factor_ldlt) doesn't perturb small pivots, and a
/// matrix that needs pivoting, such as `[[0, 1], [1, 0]]`, fails with a zero pivot.
/// [`factor_ldlt_replacing`](SymbolicCholesky::factor_ldlt_replacing) replaces tiny pivots
/// instead, for matrices whose regularization is lost to rounding, and
/// [`CsrMatrix::symmetric_scaling`] equilibrates a badly scaled matrix first.
///
/// Static pivots lose accuracy as the regularization shrinks, since a small pivot early in the
/// order makes large entries of `L`. On random quasi-definite matrices regularized by `δ`, the
/// scaled residual grew as 1/δ: about 3e-10 at δ = 1e-8, and near machine precision from 1e-2
/// up (see the crate's `ldlt_test`). Iterative refinement against the factorization recovers
/// the lost accuracy where it matters
/// ([`Ldlt::solve_refined`], [`Ldlt::solve_refined_normwise`]).
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
    /// How many pivots [`SymbolicCholesky::factor_ldlt_replacing`] replaced
    replaced: usize,
}

/// Which pivots [`SymbolicCholesky::factor_ldlt_replacing`] replaces, and with what: a pivot
/// with `|d| ≤ threshold` becomes `value` (for a row expected positive) or `−value` (for one
/// expected negative).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PivotReplacement<'a, R> {
    /// For each row of `A`, in the original numbering, whether its pivot is expected positive
    /// (as on a quasi-definite matrix's `H` block) rather than negative (its `−G` block)
    pub positive: &'a [bool],
    /// The largest magnitude replaced, at least 0
    pub threshold: R,
    /// The replacement's magnitude, positive
    pub value: R,
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
        self.factor_ldlt_with(matrix, None)
    }

    /// [`factor_ldlt`](Self::factor_ldlt), replacing each pivot whose magnitude is at most
    /// `replacement.threshold` by `±replacement.value`, with the sign its row expects, rather
    /// than failing on it or dividing by it. [`Ldlt::replaced`] counts the replacements.
    ///
    /// **What it's for:** static pivots on a regularized KKT-shaped matrix whose entries span
    /// many orders of magnitude. There a pivot whose exact value is about the regularization
    /// can be lost to rounding, become zero and fail the factorization. Replacing it is the
    /// "dynamic" regularization of interior-point solvers for quadratic and conic programs
    /// (Altman and Gondzio, "Regularized symmetric indefinite systems in interior point
    /// methods for linear and quadratic optimization," *Optim. Methods Softw.* 11(1–4):275–302,
    /// 1999, doi:10.1080/10556789908805754, regularize the pivots as they're computed). It
    /// changes the matrix that's factored, so solve against the original with refinement
    /// ([`Ldlt::solve_refined_normwise`]).
    ///
    /// **What it leaves alone:** a pivot larger than the threshold keeps its value even if its
    /// sign isn't the expected one, so [`Ldlt::inertia`] still reports a matrix that isn't
    /// quasi-definite, as when an optimization problem's Hessian is indefinite. Replacing such
    /// a pivot by a small one of the expected sign would make `L`'s entries large, so that is
    /// for the caller to correct, by regularizing more.
    ///
    /// Returns the errors of [`factor_ldlt`](Self::factor_ldlt), and
    /// `SparseError::ShapeMismatch` unless `replacement.positive` has one entry per row.
    /// A NaN pivot is still `SparseError::SingularMatrix`, and so is a zero one with a zero
    /// threshold.
    pub fn factor_ldlt_replacing<T: Scalar + ComplexFloat>(
        &self,
        matrix: &CsrMatrix<T>,
        replacement: &PivotReplacement<'_, T::Real>,
    ) -> Result<Ldlt<T>, SparseError> {
        if replacement.positive.len() != self.dim() {
            return Err(SparseError::ShapeMismatch {
                left: (self.dim(), self.dim()),
                right: (replacement.positive.len(), 1),
            });
        }
        self.factor_ldlt_with(matrix, Some(replacement))
    }

    fn factor_ldlt_with<T: Scalar + ComplexFloat>(
        &self,
        matrix: &CsrMatrix<T>,
        replacement: Option<&PivotReplacement<'_, T::Real>>,
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
        let mut replaced = 0;
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
            if let Some(replacement) = replacement {
                if Float::abs(remaining) <= replacement.threshold {
                    remaining = if replacement.positive[self.perm[i]] {
                        replacement.value
                    } else {
                        -replacement.value
                    };
                    replaced += 1;
                }
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
            replaced,
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

    /// How many pivots [`SymbolicCholesky::factor_ldlt_replacing`] replaced: 0 from
    /// [`SymbolicCholesky::factor_ldlt`].
    pub fn replaced(&self) -> usize {
        self.replaced
    }

    /// `A`'s inertia, from `D`'s signs (Sylvester's law of inertia). With replaced pivots, it's
    /// the inertia of the matrix actually factored.
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

    /// A symmetric scaling `s`, positive, such that `S A S` (`S = diag(s)`) has every row's and
    /// column's largest magnitude near 1, by `sweeps` sweeps of Ruiz's algorithm: each divides
    /// row and column `i` by the square root of row `i`'s largest magnitude (P. A. Knight,
    /// D. Ruiz and B. Uçar, "A symmetry preserving algorithm for matrix scaling," *SIAM J.
    /// Matrix Anal. Appl.* 35(3):931–955, 2014, doi:10.1137/110825753, which shows it converges
    /// for a symmetric matrix and keeps it symmetric). An empty row keeps a factor of 1.
    ///
    /// Scale a matrix whose entries span many orders of magnitude before factoring it with
    /// static pivots: solve `(S A S) y = S b`, then `x = S y`.
    ///
    /// Returns `SparseError::NotSquare` for a non-square matrix, and
    /// `SparseError::NotSymmetric` unless `A = Aᴴ`.
    pub fn symmetric_scaling(&self, sweeps: usize) -> Result<Vec<T::Real>, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        if !self.is_hermitian() {
            return Err(SparseError::NotSymmetric);
        }
        let one = <T::Real as num_traits::One>::one();
        let zero = <T::Real as num_traits::Zero>::zero();
        let n = self.shape.0;
        let mut scale = vec![one; n];
        let mut largest = vec![zero; n];
        for _ in 0..sweeps {
            largest.fill(zero);
            for (row, largest) in largest.iter_mut().enumerate() {
                for p in self.ptr[row]..self.ptr[row + 1] {
                    let entry = self.data[p].abs() * scale[row] * scale[self.idx[p]];
                    if entry > *largest {
                        *largest = entry;
                    }
                }
            }
            for (scale, &largest) in scale.iter_mut().zip(&largest) {
                if largest > zero {
                    *scale = *scale / Float::sqrt(largest);
                }
            }
        }
        Ok(scale)
    }
}

impl<T: Scalar + ComplexFloat> Factorization<T> for Ldlt<T> {
    fn dim(&self) -> usize {
        self.perm.len()
    }
}

impl<T: Scalar + ComplexFloat> Substitute<T> for Ldlt<T> {
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
