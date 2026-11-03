//! Accuracy tools for the factorizations: fixed-precision iterative refinement, with the
//! componentwise backward error it measures, and an estimate of the 1-norm condition number.
//!
//! **Refinement** (Skeel 1980; Higham 2002, ch. 12) improves a solution `x` of `A x = b` by
//! solving for the residual `r = b − A x` with the same factors and adding the correction. Each
//! step's componentwise backward error, `ω = maxᵢ |rᵢ| / (|A| |x| + |b|)ᵢ` (Oettli and Prager),
//! is the smallest relative change to each entry of `A` and `b` that makes `x` exact. It stops
//! when `ω` reaches the unit roundoff, stops halving, or after a number of steps, and returns the
//! best `x` seen.
//!
//! **The condition estimate** follows Hager (1984) with Higham's refinements (1988), as LAPACK's
//! `xLACON` does: a few solves with `A` and `Aᵀ` find a vector `v` with `‖v‖₁ = 1` whose
//! `‖A⁻¹ v‖₁` is close to, and never more than, `‖A⁻¹‖₁`. It is usually exact, and rarely off by
//! more than a factor of 3; times `‖A‖₁` it estimates `κ₁(A)`.

use num_complex::ComplexFloat;
use num_traits::{Float, Zero};

use super::{CompressedMatrix, CscMatrix, CsrMatrix, Orientation};
use crate::compressed_matrix::{Cholesky, Lu};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// The result of [`Lu::solve_refined`] or [`Cholesky::solve_refined`].
#[derive(Clone, Debug, PartialEq)]
pub struct Refined<T: ComplexFloat> {
    /// The solution with the smallest componentwise backward error found
    pub x: Vec<T>,
    /// How many corrections were applied to the first solution to reach `x`
    pub steps: usize,
    /// `x`'s componentwise backward error, `maxᵢ |b − A x|ᵢ / (|A| |x| + |b|)ᵢ`: 0 when every
    /// residual is zero, infinite when some residual is nonzero where `|A| |x| + |b|` is zero
    pub backward_error: T::Real,
}

/// The componentwise backward error of `x` for `A x = b`.
fn componentwise_backward_error<T: Scalar + ComplexFloat, O: Orientation>(
    a: &CompressedMatrix<T, O>,
    x: &[T],
    b: &[T],
) -> T::Real {
    let n = a.shape.0;
    let ax = a.mul_vec(x).expect("x has one entry per column");
    // |A| |x|, accumulated in whichever orientation the matrix has
    let mut abs_ax = vec![<T::Real as Zero>::zero(); n];
    for major in 0..a.ptr.len() - 1 {
        for p in a.ptr[major]..a.ptr[major + 1] {
            let (row, col) = if O::ROW_MAJOR {
                (major, a.idx[p])
            } else {
                (a.idx[p], major)
            };
            abs_ax[row] = abs_ax[row] + a.data[p].abs() * x[col].abs();
        }
    }
    let mut worst = <T::Real as Zero>::zero();
    for i in 0..n {
        let residual = (b[i] - ax[i]).abs();
        let scale = abs_ax[i] + b[i].abs();
        let ratio = if scale > <T::Real as Zero>::zero() {
            residual / scale
        } else if residual > <T::Real as Zero>::zero() {
            <T::Real as Float>::infinity()
        } else {
            <T::Real as Zero>::zero()
        };
        if ratio > worst {
            worst = ratio;
        }
    }
    worst
}

/// Fixed-precision refinement of `A x = b`, given a solve with `A`'s factors.
fn refine<T: Scalar + ComplexFloat, O: Orientation>(
    a: &CompressedMatrix<T, O>,
    b: &[T],
    max_steps: usize,
    solve: impl Fn(&[T]) -> Vec<T>,
) -> Refined<T> {
    let unit_roundoff = <T::Real as Float>::epsilon()
        / (<T::Real as num_traits::One>::one() + <T::Real as num_traits::One>::one());
    let half = <T::Real as num_traits::One>::one()
        / (<T::Real as num_traits::One>::one() + <T::Real as num_traits::One>::one());
    let mut x = solve(b);
    let mut omega = componentwise_backward_error(a, &x, b);
    let mut best = Refined {
        x: x.clone(),
        steps: 0,
        backward_error: omega,
    };
    for step in 1..=max_steps {
        if omega <= unit_roundoff {
            break;
        }
        let ax = a.mul_vec(&x).expect("x has one entry per column");
        let r: Vec<T> = b.iter().zip(&ax).map(|(&b, &ax)| b - ax).collect();
        let d = solve(&r);
        for (x_i, d_i) in x.iter_mut().zip(d) {
            *x_i = *x_i + d_i;
        }
        let previous = omega;
        omega = componentwise_backward_error(a, &x, b);
        if omega < best.backward_error {
            best = Refined {
                x: x.clone(),
                steps: step,
                backward_error: omega,
            };
        }
        // Stagnation (or a NaN): another step is unlikely to gain much
        if Float::is_nan(omega) || omega > half * previous {
            break;
        }
    }
    best
}

/// Hager's and Higham's estimate of `‖A⁻¹‖₁` for an `n × n` real `A`, given solves with `A` and
/// `Aᵀ`, as LAPACK's `xLACON`, keeping the largest estimate seen.
fn inverse_norm_1_estimate(
    n: usize,
    solve: impl Fn(&[f64]) -> Vec<f64>,
    solve_transpose: impl Fn(&[f64]) -> Vec<f64>,
) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let norm_1 = |v: &[f64]| v.iter().map(|value| value.abs()).sum::<f64>();
    let sign = |v: &[f64]| -> Vec<f64> {
        v.iter()
            .map(|&value| if value >= 0.0 { 1.0 } else { -1.0 })
            .collect()
    };
    let argmax_abs = |v: &[f64]| -> usize {
        let mut best = 0;
        for (i, value) in v.iter().enumerate() {
            if value.abs() > v[best].abs() {
                best = i;
            }
        }
        best
    };
    let mut y = solve(&vec![1.0 / n as f64; n]);
    let mut estimate = norm_1(&y);
    if n == 1 {
        return estimate;
    }
    let mut signs = sign(&y);
    let mut z = solve_transpose(&signs);
    let mut j = argmax_abs(&z);
    for _ in 2..=5 {
        let mut unit = vec![0.0; n];
        unit[j] = 1.0;
        y = solve(&unit);
        let previous = estimate;
        let current = norm_1(&y);
        estimate = estimate.max(current);
        let new_signs = sign(&y);
        // A repeated sign vector, or no growth, means the search has converged
        if new_signs == signs || current <= previous {
            break;
        }
        signs = new_signs;
        z = solve_transpose(&signs);
        let last = j;
        j = argmax_abs(&z);
        if z[last].abs() == z[j].abs() {
            break;
        }
    }
    // Higham's extra vector, of alternating signs and growing size, which catches matrices
    //   the search misses
    let alternating: Vec<f64> = (0..n)
        .map(|i| {
            let size = 1.0 + i as f64 / (n - 1) as f64;
            if i % 2 == 0 {
                size
            } else {
                -size
            }
        })
        .collect();
    let extra = 2.0 * norm_1(&solve(&alternating)) / (3.0 * n as f64);
    estimate.max(extra)
}

impl<T: Scalar + ComplexFloat> Lu<T> {
    /// Solves `A x = b`, then refines `x` by fixed-precision iterative refinement, for at most
    /// `max_steps` corrections, returning the iterate with the smallest componentwise backward
    /// error (see the module documentation). `a` must be the matrix these factors came from.
    /// Returns `SparseError::ShapeMismatch` unless `a` is n × n and `b` has n entries.
    pub fn solve_refined(
        &self,
        a: &CscMatrix<T>,
        b: &[T],
        max_steps: usize,
    ) -> Result<Refined<T>, SparseError> {
        let n = self.l().shape.0;
        check_refinement_shapes(a.shape, n, b.len())?;
        Ok(refine(a, b, max_steps, |rhs| {
            self.solve(rhs).expect("one entry per row")
        }))
    }
}

impl<T: Scalar + ComplexFloat> Cholesky<T> {
    /// [`Lu::solve_refined`] with Cholesky factors: `a` must be the matrix they came from.
    pub fn solve_refined(
        &self,
        a: &CsrMatrix<T>,
        b: &[T],
        max_steps: usize,
    ) -> Result<Refined<T>, SparseError> {
        let n = self.l().shape.0;
        check_refinement_shapes(a.shape, n, b.len())?;
        Ok(refine(a, b, max_steps, |rhs| {
            self.solve(rhs).expect("one entry per row")
        }))
    }
}

fn check_refinement_shapes(
    a_shape: (usize, usize),
    n: usize,
    b_len: usize,
) -> Result<(), SparseError> {
    if a_shape != (n, n) {
        return Err(SparseError::ShapeMismatch {
            left: (n, n),
            right: a_shape,
        });
    }
    if b_len != n {
        return Err(SparseError::ShapeMismatch {
            left: (n, n),
            right: (b_len, 1),
        });
    }
    Ok(())
}

impl Lu<f64> {
    /// An estimate of `‖A⁻¹‖₁`, from a few solves with `A` and `Aᵀ` (see the module
    /// documentation): never more than the true value, usually equal to it.
    pub fn inverse_norm_1_estimate(&self) -> f64 {
        inverse_norm_1_estimate(
            self.l().shape.0,
            |rhs| self.solve(rhs).expect("one entry per row"),
            |rhs| self.solve_transpose(rhs).expect("one entry per row"),
        )
    }

    /// An estimate of `A`'s 1-norm condition number `κ₁(A) = ‖A‖₁ ‖A⁻¹‖₁`, where `a` is the
    /// matrix these factors came from.
    pub fn condition_estimate_1(&self, a: &CscMatrix<f64>) -> f64 {
        a.norm_1() * self.inverse_norm_1_estimate()
    }
}

impl Cholesky<f64> {
    /// An estimate of `‖A⁻¹‖₁` (see [`Lu::inverse_norm_1_estimate`]). `A` is symmetric, so
    /// solves with `Aᵀ` are solves with `A`.
    pub fn inverse_norm_1_estimate(&self) -> f64 {
        let solve = |rhs: &[f64]| self.solve(rhs).expect("one entry per row");
        inverse_norm_1_estimate(self.l().shape.0, solve, solve)
    }

    /// An estimate of `κ₁(A)`, where `a` is the matrix these factors came from.
    pub fn condition_estimate_1(&self, a: &CsrMatrix<f64>) -> f64 {
        a.norm_1() * self.inverse_norm_1_estimate()
    }
}
