//! Dense vector operations the solvers share.
//!
//! They stay in this crate rather than in `ybus_sparse`, whose subject is sparse matrices:
//! [`norm_inf`] follows MATLAB's `norm(v, Inf)`, including its handling of NaN and of empty
//! vectors, which MIPS's termination and failure tests depend on.

/// The ∞-norm of a vector, as MATLAB's `norm(v, Inf)`: 0 for an empty vector, and NaN if any
/// entry is NaN.
pub(crate) fn norm_inf(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |m, &x| {
        if x.is_nan() || x.abs() > m {
            x.abs()
        } else {
            m
        }
    })
}

/// `aᵀ b`, summed in order, over the shorter length if they differ.
pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
