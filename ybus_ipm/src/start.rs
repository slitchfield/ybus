//! Where a solve starts: from a point alone (a cold start, as MIPS always starts), or from a
//! previous solution's point and multipliers (a warm start), for a sequence of related
//! problems.
//!
//! **Why a warm start shifts.** An interior-point method started exactly at an optimum stalls.
//! There its slacks and multipliers are complementary, `z_i μ_i ≈ 0` with one of each pair at
//! zero, and the step to the boundary blocks almost every step. A warm start therefore keeps
//! the previous point `x`, its equality multipliers `λ` and its inequality multipliers `μ`. It
//! recomputes the slacks from `x` (`z = −h(x)`, and the bounds' slacks), and raises every slack
//! and inequality multiplier to at least `warm_shift` (the solvers' option). The predictor-
//! corrector method's bound slacks are the exception: they're `x − xmin` and `xmax − x` by
//! definition, so its point is pushed `warm_shift` inside its bounds instead. That moves each
//! complementary pair back off the boundary, by the same amount wherever it was active, so each
//! product `z_i μ_i` is at least `warm_shift²`. The equality multipliers, which are free, aren't
//! moved.
//!
//! That's the simplest of the shifts the literature proposes for warm starts: E. A. Yildirim and
//! S. J. Wright, "Warm-start strategies in interior-point methods for linear programming,"
//! *SIAM J. Optim.* 12(3):782–810, 2002
//! ([doi:10.1137/S1052623400369235](https://doi.org/10.1137/S1052623400369235)), and J. Gondzio
//! and A. Grothey, "A new unblocking technique to warmstart interior point methods based on
//! sensitivity analysis," *SIAM J. Optim.* 19(3):1184–1210, 2008
//! ([doi:10.1137/060678129](https://doi.org/10.1137/060678129)), whose unblocking is more
//! elaborate.
//!
//! A primal-only warm start needs nothing new: it's [`Start::Cold`] with the previous point.
//!
//! **Scale.** A [`Solution`] reports its multipliers divided by the solver's `cost_mult`, in
//! the objective's own units. A warm start multiplies them back by the `cost_mult` of the solve
//! it starts, and shifts them in those scaled units, the solver's internal ones.

use crate::solution::{Multipliers, Solution};

/// Where a solve starts (see the [module documentation](self)).
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub enum Start<'a> {
    /// From this point, with the multipliers and slacks started as the solver always starts
    /// them: exactly the solvers' original path
    Cold(&'a [f64]),
    /// From a previous solution's point and multipliers, shifted off the boundary
    Warm(&'a WarmStart),
}

/// A previous solution's point and multipliers, for [`Start::Warm`]. The multipliers are in
/// MIPS's layout, as every [`Solution`] reports them, and must fit the problem being solved:
/// one per nonlinear constraint, linear row and variable.
#[derive(Clone, Debug, PartialEq)]
pub struct WarmStart {
    /// The point
    pub x: Vec<f64>,
    /// The multipliers there
    pub multipliers: Multipliers,
}

impl WarmStart {
    /// The point and multipliers `solution` ended at.
    pub fn from_solution(solution: &Solution) -> WarmStart {
        WarmStart {
            x: solution.x.clone(),
            multipliers: solution.multipliers.clone(),
        }
    }
}

/// `multipliers` times `factor`: a solution's reported multipliers back in a solver's scaled
/// units, for `factor` its `cost_mult`
pub(crate) fn scaled(multipliers: &Multipliers, factor: f64) -> Multipliers {
    let scale = |v: &[f64]| v.iter().map(|v| v * factor).collect();
    Multipliers {
        eq_nonlin: scale(&multipliers.eq_nonlin),
        ineq_nonlin: scale(&multipliers.ineq_nonlin),
        mu_l: scale(&multipliers.mu_l),
        mu_u: scale(&multipliers.mu_u),
        lower: scale(&multipliers.lower),
        upper: scale(&multipliers.upper),
    }
}
