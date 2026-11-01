//! Solvers for linear systems `A x = b`: triangular solves and conjugate gradient. The ordered
//! Cholesky factorization is in `cholesky.rs`, and LU in `lu.rs`.
//!
//! Both work in either orientation, for any `ComplexFloat` element type, real or complex.

use num_complex::ComplexFloat;
use num_traits::{Float, NumCast};

use super::properties::magnitude_squared;
use super::{CompressedMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// Checks that `self` is square and `b_len` matches, returning the size.
    pub(crate) fn check_system(&self, b_len: usize) -> Result<usize, SparseError> {
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
}

impl<T: Scalar + ComplexFloat, O: Orientation> CompressedMatrix<T, O> {
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
}

impl<T: Scalar + ComplexFloat, O: Orientation> CompressedMatrix<T, O> {
    /// Solves `A x = b` for a symmetric (for a complex element type, Hermitian) positive
    /// definite matrix by the conjugate gradient method, starting from zero.
    ///
    /// Inner products conjugate their first argument, `⟨u, v⟩ = Σ conj(uᵢ) vᵢ`, so the residual's
    /// squared norm, the curvature `⟨d, A d⟩` and the step lengths are real (`T::Real`). Stops
    /// when the residual norm is at most `tolerance * ||b||`. Each iteration costs one
    /// matrix-vector product, so this suits large systems where factoring would fill in. That
    /// `A` is Hermitian is assumed, not checked. Returns `SparseError::NoConvergence` if the
    /// tolerance is not reached within `max_iterations`, or if the search direction shows the
    /// matrix is not positive definite.
    pub fn conjugate_gradient(
        &self,
        b: &[T],
        tolerance: T::Real,
        max_iterations: usize,
    ) -> Result<Vec<T>, SparseError> {
        let n = self.check_system(b.len())?;
        let real_zero = <T::Real as num_traits::Zero>::zero();
        // For a real element type, these are the plain dot product and sum of squares
        let inner = |u: &[T], v: &[T]| {
            u.iter()
                .zip(v)
                .fold(T::zero(), |sum, (&a, &b)| sum + a.conj() * b)
        };
        let norm_squared = |u: &[T]| {
            u.iter()
                .fold(real_zero, |sum, &a| sum + magnitude_squared(a))
        };
        let as_element = |value: T::Real| {
            <T as NumCast>::from(value).expect("a real value converts to its own element type")
        };

        let b_norm = Float::sqrt(norm_squared(b));
        let mut x = vec![T::zero(); n];
        if b_norm == real_zero {
            return Ok(x);
        }
        let mut residual = b.to_vec();
        let mut direction = residual.clone();
        let mut residual_sq = norm_squared(&residual);

        for iteration in 0..max_iterations {
            if Float::sqrt(residual_sq) <= tolerance * b_norm {
                return Ok(x);
            }
            let a_direction = self.mul_vec(&direction)?;
            // Real for a Hermitian matrix; its imaginary part is rounding
            let curvature = inner(&direction, &a_direction).re();
            if Float::is_nan(curvature) || curvature <= real_zero {
                return Err(SparseError::NoConvergence {
                    iterations: iteration,
                });
            }
            let step = as_element(residual_sq / curvature);
            for i in 0..n {
                x[i] = x[i] + step * direction[i];
                residual[i] = residual[i] - step * a_direction[i];
            }
            let next_residual_sq = norm_squared(&residual);
            let beta = as_element(next_residual_sq / residual_sq);
            for i in 0..n {
                direction[i] = residual[i] + beta * direction[i];
            }
            residual_sq = next_residual_sq;
        }

        if Float::sqrt(residual_sq) <= tolerance * b_norm {
            Ok(x)
        } else {
            Err(SparseError::NoConvergence {
                iterations: max_iterations,
            })
        }
    }
}
