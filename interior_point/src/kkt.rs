//! The reduced KKT system both solvers solve at each iteration,
//!
//! `[W Jgᵀ; Jg 0] [dx; dλ] = [r_x; r_λ]`,
//!
//! factored by LU, the default, or by a pivoted LDLᵀ ([`KktFactorization`]).
//!
//! **LU** ([`KKT_ORDERING`]) factors the matrix as assembled, as
//! MATPOWER's MIPS does with its backslash. MIPS's results match MATPOWER's with it.
//!
//! **LDLᵀ** is `sparse_mat`'s multifrontal factorization with threshold Bunch–Kaufman pivoting
//! ([`IndefiniteLdlt`]), ordered by AMD:
//!
//! - **Exactly symmetric:** the lower triangle of `W` is mirrored, so the factorization's symmetry
//!   check holds, and every diagonal is stored, so the pattern stays fixed and the analysis is
//!   kept.
//! - **No regularization:** 2×2 pivots handle the zero diagonals, of the constraints' block and of
//!   `W`'s, that defeat static pivots.
//! - **The threshold, 1e-4:** measured on the KKT matrices of DC and AC `case9241pegase` and
//!   `case_ACTIVSg25k`. MA57's 0.01 delayed so many pivots that the factor outgrew LU's, and 1e-8
//!   lost accuracy that refinement couldn't recover.
//! - **Refined:** each solve takes up to 3 steps of normwise refinement against the matrix, until
//!   the relative residual is at most 1e-12 or stops halving; it reached LU's residual on those
//!   matrices.
//! - **The inertia:** the factorization gives it, and the predictor-corrector method corrects an
//!   indefinite `W` with it (see [`crate::mehrotra`]). MIPS only swaps the factorization.

use std::cell::Cell;
use std::rc::Rc;

use sparse_mat::compressed_matrix::{
    CscMatrix, CsrMatrix, IndefiniteLdlt, IndefiniteOptions, Inertia, Lu, LuOptions, Ordering,
    SymbolicIndefinite, SymbolicLu,
};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::Factorization;

use crate::mips::KKT_ORDERING;

/// How the solvers factor the KKT matrix (see the [module documentation](self)).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum KktFactorization {
    /// LU, ordered by [`KKT_ORDERING`], as assembled: MATPOWER's
    /// choice, and the default
    #[default]
    Lu,
    /// The pivoted LDLᵀ, ordered by approximate minimum degree, with each solve refined. It
    /// reports the matrix's inertia, which the predictor-corrector method uses to correct an
    /// indefinite Hessian.
    Ldlt,
}

/// LDLᵀ's pivot threshold
const THRESHOLD: f64 = 1e-4;
/// LDLᵀ's refinement: the most steps, and the relative residual it stops at
const REFINEMENT_STEPS: usize = 3;
const REFINEMENT_TOLERANCE: f64 = 1e-12;

/// The KKT matrix's factorization, with its analysis kept across a solve's iterations while
/// the pattern holds, and counts of the analyses and refinement steps.
pub(crate) struct KktSolver {
    method: KktFactorization,
    lu: Option<SymbolicLu>,
    ldlt: Option<SymbolicIndefinite>,
    analyses: usize,
    /// Shared with each factor, which adds its solves' steps
    refinement_steps: Rc<Cell<usize>>,
}

/// One iteration's factored KKT matrix.
pub(crate) struct KktFactor {
    kind: FactorKind,
    refinement_steps: Rc<Cell<usize>>,
}

enum FactorKind {
    Lu(Box<Lu<f64>>),
    Ldlt {
        ldlt: IndefiniteLdlt<f64>,
        /// The factored matrix, which the solves refine against
        matrix: CsrMatrix<f64>,
    },
}

impl KktSolver {
    pub(crate) fn new(method: KktFactorization) -> KktSolver {
        KktSolver {
            method,
            lu: None,
            ldlt: None,
            analyses: 0,
            refinement_steps: Rc::new(Cell::new(0)),
        }
    }

    /// Factors `[W Jgᵀ; Jg 0]` (`nx + neq` square), from `W`'s triplets (duplicates summed) and
    /// `Jg` (`neq × nx`), with `delta_w` added to `W`'s diagonal for LDLᵀ (an inertia
    /// correction; LU ignores it). Returns the analysis's and factorization's errors.
    ///
    /// LU assembles the matrix as the solvers always have: `W`'s triplets, then each entry of
    /// `Jg` at `(j, nx + k)` and `(nx + k, j)`.
    pub(crate) fn factor(
        &mut self,
        nx: usize,
        w: Vec<(usize, usize, f64)>,
        jg: &CsrMatrix<f64>,
        delta_w: f64,
    ) -> Result<KktFactor, SparseError> {
        let size = nx + jg.shape().0;
        let kind = match self.method {
            KktFactorization::Lu => {
                let mut triplets = w;
                for (k, j, v) in jg.iter() {
                    triplets.push((j, nx + k, v));
                    triplets.push((nx + k, j, v));
                }
                let kkt = CscMatrix::from_triplets(size, size, triplets).expect("indices in range");
                FactorKind::Lu(Box::new(self.factor_lu(&kkt)?))
            }
            KktFactorization::Ldlt => self.factor_ldlt(nx, size, w, jg, delta_w)?,
        };
        Ok(KktFactor {
            kind,
            refinement_steps: Rc::clone(&self.refinement_steps),
        })
    }

    /// Factors `kkt` with the kept analysis, or, if there is none yet or `kkt`'s pattern
    /// differs from the one it was made for, with a fresh one, which is then kept.
    fn factor_lu(&mut self, kkt: &CscMatrix<f64>) -> Result<Lu<f64>, SparseError> {
        if let Some(symbolic) = &self.lu {
            match symbolic.factor(kkt, LuOptions::default()) {
                Err(SparseError::PatternMismatch) => {}
                factored => return factored,
            }
        }
        let symbolic = kkt.analyze_lu(KKT_ORDERING)?;
        self.analyses += 1;
        let factored = symbolic.factor(kkt, LuOptions::default());
        self.lu = Some(symbolic);
        factored
    }

    fn factor_ldlt(
        &mut self,
        nx: usize,
        size: usize,
        w: Vec<(usize, usize, f64)>,
        jg: &CsrMatrix<f64>,
        delta_w: f64,
    ) -> Result<FactorKind, SparseError> {
        // The lower triangle, every diagonal stored (with δ_w on W's), then mirrored: exactly
        //   symmetric
        let mut lower: Vec<(usize, usize, f64)> =
            w.into_iter().filter(|&(r, c, _)| r >= c).collect();
        lower.extend(jg.iter().map(|(k, j, v)| (nx + k, j, v)));
        lower.extend((0..size).map(|i| (i, i, if i < nx { delta_w } else { 0.0 })));
        let lower = CsrMatrix::from_triplets(size, size, lower).expect("indices in range");
        let mut both = Vec::with_capacity(2 * lower.values().len());
        for (r, c, v) in lower.iter() {
            both.push((r, c, v));
            if r != c {
                both.push((c, r, v));
            }
        }
        let matrix = CsrMatrix::from_triplets(size, size, both).expect("indices in range");
        let options = IndefiniteOptions {
            threshold: THRESHOLD,
        };
        if let Some(symbolic) = &self.ldlt {
            match symbolic.factor(&matrix, options) {
                Err(SparseError::PatternMismatch) => {}
                factored => {
                    return Ok(FactorKind::Ldlt {
                        ldlt: factored?,
                        matrix,
                    })
                }
            }
        }
        let symbolic = matrix.analyze_indefinite(Ordering::Amd)?;
        self.analyses += 1;
        let factored = symbolic.factor(&matrix, options);
        self.ldlt = Some(symbolic);
        Ok(FactorKind::Ldlt {
            ldlt: factored?,
            matrix,
        })
    }

    /// The analyses made so far.
    pub(crate) fn analyses(&self) -> usize {
        self.analyses
    }

    /// The refinement steps taken so far, over every solve (0 with LU).
    pub(crate) fn refinement_steps(&self) -> usize {
        self.refinement_steps.get()
    }
}

impl KktFactor {
    /// Solves the KKT system for `rhs`: with LU, directly; with LDLᵀ, refined.
    pub(crate) fn solve(&self, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
        match &self.kind {
            FactorKind::Lu(lu) => lu.solve(rhs),
            FactorKind::Ldlt { ldlt, matrix } => {
                let refined = ldlt.solve_refined_normwise(
                    matrix,
                    rhs,
                    REFINEMENT_STEPS,
                    REFINEMENT_TOLERANCE,
                )?;
                self.refinement_steps
                    .set(self.refinement_steps.get() + refined.steps);
                Ok(refined.x)
            }
        }
    }

    /// The factored matrix's inertia, with LDLᵀ (`None` with LU).
    pub(crate) fn inertia(&self) -> Option<Inertia> {
        match &self.kind {
            FactorKind::Lu(_) => None,
            FactorKind::Ldlt { ldlt, .. } => Some(ldlt.inertia()),
        }
    }
}
