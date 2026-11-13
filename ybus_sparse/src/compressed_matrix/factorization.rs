//! [`Factorization`], the solves every factorization of a square matrix provides.

use crate::sparse_error::SparseError;
use crate::Scalar;

pub(crate) mod sealed {
    /// The one step each factorization supplies, so that [`Factorization`](super::Factorization)
    /// can give every one the same solves. Private, which also keeps `Factorization` to this
    /// crate's factorizations.
    pub trait Substitute<T> {
        /// Overwrites `x` (holding `b`, with `dim()` entries, already checked) with the
        /// solution of `A x = b`, using `scratch` as workspace: whatever it holds is
        /// overwritten, and its capacity is kept for the next solve
        fn substitute(&self, x: &mut [T], scratch: &mut Vec<T>);
    }
}

/// The solves of `A x = b` that every factorization of a square matrix `A` provides:
/// [`Cholesky`](super::Cholesky), [`Lu`](super::Lu) and [`Ldlt`](super::Ldlt), and the
/// [`experimental`](super::experimental) ones. Each factorization's own methods add what only it
/// has, such as [`Lu::solve_transpose`](super::Lu::solve_transpose) or the LDLᵀ factorizations'
/// inertia.
///
/// Every solve of a right-hand side runs the same operations, so a solution is bit for bit the
/// same whichever of these methods computes it. Each returns `SparseError::ShapeMismatch` unless
/// every right-hand side has one entry per row of `A`, before solving any.
///
/// The trait is sealed: only this crate's factorizations implement it.
pub trait Factorization<T: Scalar>: sealed::Substitute<T> {
    /// The number of rows (and columns) of the factored matrix.
    fn dim(&self) -> usize;

    /// Solves `A x = b`.
    fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        let mut x = b.to_vec();
        self.solve_in_place(&mut x)?;
        Ok(x)
    }

    /// Solves `A x = b`, overwriting `rhs` (which holds `b`) with `x`. It saves returning a new
    /// vector, though most factorizations still permute through one scratch vector. On an
    /// error, `rhs` is untouched.
    fn solve_in_place(&self, rhs: &mut [T]) -> Result<(), SparseError> {
        check_rhs(self, rhs.len())?;
        let mut scratch = Vec::with_capacity(rhs.len());
        self.substitute(rhs, &mut scratch);
        Ok(())
    }

    /// Solves `A x = b` for each right-hand side `b` in `rhs`, returning the solutions in the
    /// same order, with one scratch vector serving them all.
    fn solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        for b in rhs {
            check_rhs(self, b.len())?;
        }
        let mut scratch = Vec::with_capacity(self.dim());
        Ok(rhs
            .iter()
            .map(|b| {
                let mut x = b.clone();
                self.substitute(&mut x, &mut scratch);
                x
            })
            .collect())
    }

    /// [`Factorization::solve_many`] with the right-hand sides solved in parallel (feature
    /// `rayon`), each by one thread with its own scratch vector, so the results are bit for bit
    /// `solve_many`'s, in the same order.
    #[cfg(feature = "rayon")]
    fn par_solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError>
    where
        Self: Sync,
        T: Send + Sync,
    {
        use rayon::prelude::*;
        for b in rhs {
            check_rhs(self, b.len())?;
        }
        let n = self.dim();
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

/// `ShapeMismatch` unless a right-hand side of length `len` fits `factor`
fn check_rhs<T: Scalar, F: Factorization<T> + ?Sized>(
    factor: &F,
    len: usize,
) -> Result<(), SparseError> {
    let n = factor.dim();
    if len == n {
        Ok(())
    } else {
        Err(SparseError::ShapeMismatch {
            left: (n, n),
            right: (len, 1),
        })
    }
}
