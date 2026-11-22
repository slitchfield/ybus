//! A primal-dual interior-point method with Mehrotra's predictor-corrector, written from the
//! literature (Mehrotra 1992; Nocedal and Wright, *Numerical Optimization*, 2nd ed., Algorithm
//! 14.3, §16.6 and ch. 19; S. J. Wright, *Primal-Dual Interior-Point Methods*, 1997) for the
//! same programs as [`mips`](crate::mips):
//!
//! minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.
//!
//! # The method
//!
//! **The constraints.** The linear rows are split as MIPS splits them: an equality `Ae x = be`
//! for `|u − l| ≤ ε`, and otherwise one or two inequalities `Ai x − bi ≤ 0`. The nonlinear
//! constraints come first, as in MIPS: the equalities are `[g(x); Ae x − be] = 0`, with Jacobian
//! `Jg = [∂g/∂x; Ae]`, and the inequalities `[h(x); Ai x − bi] + z = 0`, with Jacobian `Jh` and
//! a slack `z > 0` and multiplier `μ > 0` each. A bound side within ±1e10 is kept apart: `x`
//! stays strictly inside the bounds, whose slacks are `s_l = x − xmin` and `s_u = xmax − x`,
//! with multipliers `π_l` and `π_u`. A variable with `xmin = xmax` is an equality row, as in
//! MIPS.
//!
//! **Each iteration** takes a Newton step on the perturbed KKT conditions, reduced to
//!
//! `[W Jgᵀ; Jg 0] [dx; dλ] = [−r_d − Jhᵀ Z⁻¹ (r_z + M (h + z)) + S_l⁻¹ r_l − S_u⁻¹ r_u; −g]`
//!
//! with `W = ∇²L + Jhᵀ Z⁻¹ M Jh + S_l⁻¹ Π_l + S_u⁻¹ Π_u`. Here:
//!
//! - `∇²L` is the Lagrangian's Hessian, or the objective's without nonlinear constraints;
//! - the bounds add only diagonal terms;
//! - `r_d` is the Lagrangian's gradient;
//! - `r_z`, `r_l` and `r_u` are the right-hand sides of the linearized complementarity
//!   conditions.
//!
//! The matrix is factored once (by LU ordered by [`KKT_ORDERING`](crate::mips::KKT_ORDERING), or
//! by LDLᵀ, see [`crate::kkt`]; analyzed at the first iteration and kept while the pattern holds)
//! and solved twice:
//!
//! 1. **The predictor**, an affine step towards complementarity zero (`r_z = −Z μ`, and so on),
//!    whose longest feasible step lengths give the complementarity it would reach, `μ_aff`.
//! 2. **The corrector**, aiming at `σ μ` with `σ = (μ_aff / μ)³`, `μ` being the average
//!    complementarity, and adding the predictor's second-order term
//!    (`r_z = −Z μ − ΔZ_aff Δμ_aff + σ μ e`).
//!
//! **With nonlinear constraints,** the second-order term is damped by the affine step lengths,
//! `r_z = −Z μ − α_p,aff α_d,aff ΔZ_aff Δμ_aff + σ μ e`.
//!
//! - **Why:** the affine step can be enormous on a nonconvex problem far from its solution.
//!   With the full term, the corrector's step reached 2e12 at the AC OPF's first iteration on
//!   `case_ACTIVSg10k`, past `max_stepsize`. On `case_ACTIVSg2000`, `case3012wp`, `case3120sp`,
//!   `case3375wp` and `case8387pegase`, the steps were blocked at lengths of 1e-4 from the
//!   start, and the solve failed where MIPS converges.
//! - **The effect:** damped, all five converge. Merit-function backtracking couldn't help, since
//!   the steps were already short.
//! - **Its origin:** this form is the crate's own safeguard, chosen from those experiments (phase
//!   O3's record), not a method from the literature.
//! - **Unchanged elsewhere:** linear and quadratic programs take the full term, as Mehrotra
//!   describes it.
//!
//! **Step lengths.** The primal and dual steps are each at most `xi` of the way to the boundary.
//!
//! - **A linear program** takes them with separate lengths.
//! - **Otherwise both take the shorter of the two:** wherever the objective has curvature
//!   (Nocedal and Wright, §16.6), and wherever there are nonlinear constraints.
//! - **Why:** there the dual residual depends on `x`, through `∇f(x)` or the constraints'
//!   Jacobians. Separate lengths can then let complementarity grow, as they did on a small QP
//!   with two-sided bounds, which failed to converge with them.
//!
//! **Otherwise:**
//!
//! - A problem with no inequalities or bounds takes plain Newton steps.
//! - There is no merit function or line search by default, as in MIPS. An option,
//!   `line_search` ([`LineSearch`], `docs/LineSearchChecklist.md`), adds Wächter and Biegler's
//!   monotone barrier with their filter line search and feasibility restoration, alone or in
//!   Nocedal, Wächter and Waltz's adaptive framework, which keeps Mehrotra's steps while the
//!   KKT error falls.
//! - **The inertia correction,** with the pivoted LDLᵀ (`kkt_factorization`, see [`crate::kkt`]).
//!   On a nonconvex problem an indefinite `∇²L` can make the step a poor one, as in MIPS. LDLᵀ
//!   reveals it: the KKT matrix's inertia must be `(nx, neq)`, which holds when `W` is positive
//!   definite on `Jg`'s null space. When it isn't, `W` is shifted by `δ_w I` and refactored, as in
//!   Wächter and Biegler's Algorithm IC ("On the implementation of an interior-point filter
//!   line-search algorithm for large-scale nonlinear programming," *Math. Program.*
//!   106(1):25–57, 2006, doi:10.1007/s10107-004-0559-y, §3.1, with their defaults):
//!   - **the first shift** is 1e-4, or a third of the last iteration's (at least 1e-20);
//!   - **it grows** 100× the first time within an iteration and 8× after;
//!   - **past 1e40,** the solve fails.
//!
//!   The shift is in the iteration's record. It's off by default (`inertia_correction`), as the
//!   OPF survey favoured. With LU there's no inertia and no correction.
//! - It converges by MIPS's four measures, defined alike, with the bounds' slacks and
//!   multipliers counted as MIPS counts its bound rows, so iteration counts compare fairly.

mod filter;
mod restoration;

use ybus_sparse::compressed_matrix::CsrMatrix;

use crate::interior_point_error::InteriorPointError;
use ybus_sparse::compressed_matrix::Inertia;

use crate::kkt::{KktFactorization, KktOrdering, KktSolver};
use crate::newton::WeightedGram;
use crate::problem::{
    Constraints, NonlinearConstraints, NonlinearMultipliers, Problem, QuadraticProblem,
};
use crate::solution::{History, Multipliers, Solution, Status};
use crate::start::{self, Start};
use crate::vector::{dot, norm_inf};

/// The predictor-corrector method's options, with MIPS's termination defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct MehrotraOptions {
    /// Termination tolerance for feasibility
    pub feastol: f64,
    /// Termination tolerance for the gradient
    pub gradtol: f64,
    /// Termination tolerance for complementarity
    pub comptol: f64,
    /// Termination tolerance for the cost
    pub costtol: f64,
    /// The most iterations
    pub max_it: usize,
    /// The objective's scale factor, as MIPS's `cost_mult`
    pub cost_mult: f64,
    /// The fraction of the distance to the boundary a step may take, in `[0.5, 1)`
    pub xi: f64,
    /// The slacks' and multipliers' smallest starting value, as MIPS's `z0`
    pub z0: f64,
    /// A warm start's shift ([`Start::Warm`]): every slack and inequality multiplier, the
    /// bounds' included, starts at least this far from zero, in the solver's scaled units
    /// (the objective's times `cost_mult`). Unused
    /// by a cold start. The default, 1e-4, is provisional until the OPF sequences' survey
    /// measures it.
    pub warm_shift: f64,
    /// How far inside a bound `x0` is moved, relative to the bound's magnitude (at least 1)
    pub bound_push: f64,
    /// The most a starting point is moved inside a variable's two bounds, as a fraction of the
    /// distance between them
    pub bound_frac: f64,
    /// The smallest acceptable primal or dual step length
    pub alpha_min: f64,
    /// The threshold below which a multiplier is set to zero when its constraint isn't binding,
    /// as MIPS's `mu_threshold`
    pub mu_threshold: f64,
    /// The largest acceptable 2-norm of the reduced Newton step
    pub max_stepsize: f64,
    /// How the KKT system is factored: LU by default, or LDLᵀ (see [`crate::kkt`])
    pub kkt_factorization: KktFactorization,
    /// How the pivoted LDLᵀ orders the KKT matrix: AMD by default, or matched (see
    /// [`KktOrdering`]); LU ignores it
    pub kkt_ordering: KktOrdering,
    /// Whether to correct the KKT matrix's inertia when it's factored by LDLᵀ (see the module
    /// docs): off by default, since on the OPF survey it changed no outcome for the better and
    /// cost iterations (`docs/benchmarks/ldlt-opf-2026-11-08.md`); with LU there's no inertia,
    /// and this has no effect
    pub inertia_correction: bool,
    /// The step's globalization (see [`LineSearch`]): none by default, today's method. A filter
    /// needs `kkt_factorization` to be LDLᵀ, and turns `inertia_correction` on
    pub line_search: LineSearch,
}

/// How the predictor-corrector method accepts a step (`MehrotraOptions::line_search`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum LineSearch {
    /// Every step is taken whole, at the fraction-to-the-boundary length: the method as it was
    /// before the line search, bit for bit
    #[default]
    Off,
    /// Wächter and Biegler's monotone barrier with their filter line search (see
    /// `docs/LineSearchChecklist.md`): centered Newton steps at a fixed barrier parameter,
    /// lowered as each barrier problem is solved, each accepted or shortened by the filter, with
    /// second-order corrections and the inertia correction on. Where no step is acceptable, or
    /// none can be computed, feasibility restoration (§3.3) minimizes the constraint violation
    /// until the filter accepts a point with less of it; a solve whose restoration can't, as at
    /// a local minimizer of the violation, ends as `NumericalFailure` with its last record in
    /// [`IterationMode::Restoration`]
    MonotoneFilter,
    /// Nocedal, Wächter and Waltz's adaptive barrier (see `docs/LineSearchChecklist.md`): free
    /// mode takes `Off`'s predictor-corrector steps, with Mehrotra's barrier parameter, while
    /// the KKT error keeps falling to `0.9999` of the largest of its last five accepted values.
    /// When it doesn't, fixed mode runs `MonotoneFilter`'s method from that point, at `0.8` of
    /// its average complementarity, until the KKT error meets the same test, and free mode
    /// resumes. Wherever free mode's test always holds, the solve is `Off`'s with LDLᵀ and the
    /// inertia correction, bit for bit
    AdaptiveFilter,
}

/// Which phase of the method an iteration belongs to (`MehrotraIteration::mode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum IterationMode {
    /// Mehrotra's adaptive barrier: the method's only mode with `LineSearch::Off`, and
    /// `LineSearch::AdaptiveFilter`'s while its progress test holds
    #[default]
    Free,
    /// A fixed barrier parameter, with the filter line search
    Fixed,
    /// Feasibility restoration: one record for the whole phase, after it, with its
    /// iterations in `restoration_iterations`
    Restoration,
}

impl Default for MehrotraOptions {
    fn default() -> MehrotraOptions {
        MehrotraOptions {
            feastol: 1e-6,
            gradtol: 1e-6,
            comptol: 1e-6,
            costtol: 1e-6,
            max_it: 150,
            cost_mult: 1.0,
            xi: 0.99995,
            z0: 1.0,
            warm_shift: 1e-4,
            bound_push: 1e-2,
            bound_frac: 1e-2,
            alpha_min: 1e-8,
            mu_threshold: 1e-5,
            max_stepsize: 1e10,
            kkt_factorization: KktFactorization::Lu,
            kkt_ordering: KktOrdering::Amd,
            inertia_correction: false,
            line_search: LineSearch::Off,
        }
    }
}

/// The inertia correction's parameters, Wächter and Biegler's defaults (§3.1 of the paper the module docs cite): the first
/// shift `δ̄_w⁰`, its smallest and largest, the factor it shrinks by from one iteration to the
/// next (`κ_w⁻`), and the factors it grows by within an iteration (`κ̄_w⁺` the first time, `κ_w⁺`
/// after)
const DELTA_W_FIRST: f64 = 1e-4;
const DELTA_W_MIN: f64 = 1e-20;
const DELTA_W_MAX: f64 = 1e40;
const DELTA_W_DECREASE: f64 = 1.0 / 3.0;
const DELTA_W_FIRST_INCREASE: f64 = 100.0;
const DELTA_W_INCREASE: f64 = 8.0;

/// One iteration's record: the termination measures, the average complementarity and the
/// centering it chose, the step's 2-norm, the objective and the step lengths.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct MehrotraIteration {
    /// The feasibility measure
    pub feascond: f64,
    /// The gradient measure
    pub gradcond: f64,
    /// The complementarity measure
    pub compcond: f64,
    /// The cost measure
    pub costcond: f64,
    /// The average complementarity after the step, `μ`
    pub mu: f64,
    /// The centering the step aimed at, `σ = (μ_aff/μ)³` (0 for a plain Newton step)
    pub sigma: f64,
    /// The step's 2-norm, `‖dx‖`
    pub stepsize: f64,
    /// The objective, without `cost_mult`
    pub obj: f64,
    /// The primal step length
    pub alphap: f64,
    /// The dual step length
    pub alphad: f64,
    /// The inertia correction's shift of `W`, `δ_w` (0 when the inertia was right, and always
    /// with LU)
    pub delta_w: f64,
    /// The iteration's phase: always [`IterationMode::Free`] with `LineSearch::Off`
    pub mode: IterationMode,
    /// The constraint violation after the step, `θ = ‖[g; h + z]‖₁`, over the nonlinear and
    /// linear rows, the inequalities with their slacks
    pub theta: f64,
    /// The barrier function after the step, `φ_μ = f − μ (Σ ln z + Σ ln s_l + Σ ln s_u)`, with
    /// `f` times `cost_mult`, at the fixed barrier parameter `μ` (0 without a filter)
    pub barrier: f64,
    /// How many times the line search halved the step (0 without a filter)
    pub backtracks: usize,
    /// Whether the step was a second-order correction (false without a filter)
    pub second_order: bool,
    /// The fixed barrier parameter the step aimed at, in fixed mode (0 in free mode)
    pub barrier_parameter: f64,
    /// The filter's entries after the step (0 without a filter)
    pub filter_size: usize,
    /// In a restoration record ([`IterationMode::Restoration`]), the restoration phase's
    /// iterations, which the solve's count includes (0 otherwise)
    pub restoration_iterations: usize,
    /// With `LineSearch::AdaptiveFilter`, the original problem's KKT error `Φ` after the step,
    /// which its progress test compares (0 otherwise)
    pub kkt_error: f64,
}

/// The internal solve's hooks, which feasibility restoration sets for the solve of its own
///   problem: none for a solve the caller asked for
#[derive(Default)]
struct Hooks<'h> {
    /// Fixed mode's first barrier parameter, in place of `μ₀`
    mu_start: Option<f64>,
    /// Whether restoration may not be entered: inside restoration, so that it can't recurse
    no_restoration: bool,
    /// Ends the solve, as converged, after an iteration whose point it accepts
    stop: Option<&'h Stop<'h>>,
}

/// A solve's early end: whether to stop at a point
type Stop<'h> = dyn Fn(&[f64]) -> bool + 'h;

/// A trial point of fixed mode's line search: its evaluation and slacks, and its constraint
/// violation `θ` and barrier function `φ`
struct Trial {
    e: Evaluation,
    z: Vec<f64>,
    theta: f64,
    phi: f64,
}

/// The linear rows and bounds, split as the method uses them.
struct Split {
    /// `Ae` and `be`: the linear rows with `|u − l| ≤ ε`, then one row per fixed variable
    ae: CsrMatrix<f64>,
    be: Vec<f64>,
    /// `Ai` and `bi`, as MIPS stacks them: `[A(ilt); −A(igt); A(ibx); −A(ibx)]`
    ai: CsrMatrix<f64>,
    bi: Vec<f64>,
    /// The rows of `A` in each class, as MIPS names them
    ieq: Vec<usize>,
    igt: Vec<usize>,
    ilt: Vec<usize>,
    ibx: Vec<usize>,
    /// The fixed variables, whose equality rows follow `A(ieq)`'s in `Ae`
    fixed: Vec<usize>,
    /// The variables with a lower bound (not fixed), and the bounds
    lo: Vec<usize>,
    xmin: Vec<f64>,
    /// The variables with an upper bound (not fixed), and the bounds
    up: Vec<usize>,
    xmax: Vec<f64>,
    /// The number of rows of `A` (0 if it's dropped)
    na: usize,
}

impl Split {
    fn new(constraints: &Constraints, nx: usize) -> Result<Split, InteriorPointError> {
        // The linear constraints: none if A is absent, or if every limit is infinite (as MIPS)
        let mut a = constraints.a.as_ref();
        if a.is_some_and(|a| a.shape().1 != nx) {
            return Err(InteriorPointError::Dimension { what: "A" });
        }
        let na0 = a.map_or(0, |a| a.shape().0);
        let mut l = constraints
            .l
            .clone()
            .unwrap_or_else(|| vec![f64::NEG_INFINITY; na0]);
        let mut u = constraints
            .u
            .clone()
            .unwrap_or_else(|| vec![f64::INFINITY; na0]);
        if l.len() != na0 || u.len() != na0 {
            return Err(InteriorPointError::Dimension { what: "l or u" });
        }
        let l_free = l.iter().all(|&v| v == f64::NEG_INFINITY);
        let u_free = u.iter().all(|&v| v == f64::INFINITY);
        if na0 == 0 || (l_free && u_free) {
            a = None;
            l.clear();
            u.clear();
        }
        let na = l.len();
        let bound = |side: &Option<Vec<f64>>, default: f64| -> Result<Vec<f64>, _> {
            let values = side.clone().unwrap_or_else(|| vec![default; nx]);
            if values.len() == nx {
                Ok(values)
            } else {
                Err(InteriorPointError::Dimension {
                    what: "xmin or xmax",
                })
            }
        };
        let xmin = bound(&constraints.xmin, f64::NEG_INFINITY)?;
        let xmax = bound(&constraints.xmax, f64::INFINITY)?;

        let rows = 0..na;
        let ieq: Vec<usize> = rows
            .clone()
            .filter(|&r| (u[r] - l[r]).abs() <= f64::EPSILON)
            .collect();
        let igt: Vec<usize> = rows
            .clone()
            .filter(|&r| u[r] >= 1e10 && l[r] > -1e10)
            .collect();
        let ilt: Vec<usize> = rows
            .clone()
            .filter(|&r| l[r] <= -1e10 && u[r] < 1e10)
            .collect();
        let ibx: Vec<usize> = rows
            .filter(|&r| (u[r] - l[r]).abs() > f64::EPSILON && u[r] < 1e10 && l[r] > -1e10)
            .collect();
        let fixed: Vec<usize> = (0..nx)
            .filter(|&k| (xmax[k] - xmin[k]).abs() <= f64::EPSILON)
            .collect();
        let lo: Vec<usize> = (0..nx)
            .filter(|&k| xmin[k] > -1e10 && (xmax[k] - xmin[k]).abs() > f64::EPSILON)
            .collect();
        let up: Vec<usize> = (0..nx)
            .filter(|&k| xmax[k] < 1e10 && (xmax[k] - xmin[k]).abs() > f64::EPSILON)
            .collect();

        let columns: Vec<usize> = (0..nx).collect();
        let empty = || CsrMatrix::from_triplets(0, nx, vec![]).expect("0 × nx");
        let select = |rows: &[usize]| match a {
            Some(a) => a.select(rows, &columns).expect("rows of A"),
            None => empty(),
        };
        let fixed_rows = CsrMatrix::from_triplets(
            fixed.len(),
            nx,
            fixed
                .iter()
                .enumerate()
                .map(|(i, &k)| (i, k, 1.0))
                .collect::<Vec<_>>(),
        )
        .expect("fixed variables in range");
        let ae = CsrMatrix::vstack([select(&ieq), fixed_rows]).expect("nx columns");
        let mut be: Vec<f64> = ieq.iter().map(|&r| u[r]).collect();
        be.extend(fixed.iter().map(|&k| xmax[k]));
        let ai = CsrMatrix::vstack([select(&ilt), -select(&igt), select(&ibx), -select(&ibx)])
            .expect("nx columns");
        let mut bi: Vec<f64> = ilt.iter().map(|&r| u[r]).collect();
        bi.extend(igt.iter().map(|&r| -l[r]));
        bi.extend(ibx.iter().map(|&r| u[r]));
        bi.extend(ibx.iter().map(|&r| -l[r]));
        Ok(Split {
            ae,
            be,
            ai,
            bi,
            ieq,
            igt,
            ilt,
            ibx,
            fixed,
            xmin: lo.iter().map(|&k| xmin[k]).collect(),
            xmax: up.iter().map(|&k| xmax[k]).collect(),
            lo,
            up,
            na,
        })
    }
}

/// The longest step `α ≤ 1` keeping `values + α steps` nonnegative, times `fraction` if the
/// boundary is within reach (1 if no step heads towards it).
fn step_length(values: &[f64], steps: &[f64], fraction: f64) -> f64 {
    let limit = values
        .iter()
        .zip(steps)
        .filter(|&(_, &d)| d < 0.0)
        .map(|(v, d)| v / -d)
        .fold(f64::INFINITY, f64::min);
    if limit.is_finite() {
        (fraction * limit).min(1.0)
    } else {
        1.0
    }
}

/// The primal and dual variables of an iterate, beyond `x` and `λ`, over the inequality rows
/// and the bounds.
struct Pairs {
    z: Vec<f64>,
    mu: Vec<f64>,
    s_l: Vec<f64>,
    pi_l: Vec<f64>,
    s_u: Vec<f64>,
    pi_u: Vec<f64>,
}

impl Pairs {
    fn count(&self) -> usize {
        self.z.len() + self.s_l.len() + self.s_u.len()
    }

    /// The total complementarity, `zᵀμ + s_lᵀπ_l + s_uᵀπ_u`
    fn complementarity(&self) -> f64 {
        dot(&self.z, &self.mu) + dot(&self.s_l, &self.pi_l) + dot(&self.s_u, &self.pi_u)
    }
}

/// A step: the change in every variable.
struct Step {
    dx: Vec<f64>,
    dlam: Vec<f64>,
    dz: Vec<f64>,
    dmu: Vec<f64>,
    dpi_l: Vec<f64>,
    dpi_u: Vec<f64>,
}

impl Step {
    /// The primal and dual step lengths that keep every slack and multiplier nonnegative, at
    /// most `fraction` of the way to the boundary
    fn lengths(&self, split: &Split, p: &Pairs, fraction: f64) -> (f64, f64) {
        let ds_l: Vec<f64> = split.lo.iter().map(|&k| self.dx[k]).collect();
        let ds_u: Vec<f64> = split.up.iter().map(|&k| -self.dx[k]).collect();
        let alphap = step_length(&p.z, &self.dz, fraction)
            .min(step_length(&p.s_l, &ds_l, fraction))
            .min(step_length(&p.s_u, &ds_u, fraction));
        let alphad = step_length(&p.mu, &self.dmu, fraction)
            .min(step_length(&p.pi_l, &self.dpi_l, fraction))
            .min(step_length(&p.pi_u, &self.dpi_u, fraction));
        (alphap, alphad)
    }

    /// The total complementarity after steps of these lengths
    fn complementarity_after(&self, split: &Split, p: &Pairs, alphap: f64, alphad: f64) -> f64 {
        let after = |v: &[f64], d: &[f64], a: f64| -> Vec<f64> {
            v.iter().zip(d).map(|(v, d)| v + a * d).collect()
        };
        let ds_l: Vec<f64> = split.lo.iter().map(|&k| self.dx[k]).collect();
        let ds_u: Vec<f64> = split.up.iter().map(|&k| -self.dx[k]).collect();
        dot(
            &after(&p.z, &self.dz, alphap),
            &after(&p.mu, &self.dmu, alphad),
        ) + dot(
            &after(&p.s_l, &ds_l, alphap),
            &after(&p.pi_l, &self.dpi_l, alphad),
        ) + dot(
            &after(&p.s_u, &ds_u, alphap),
            &after(&p.pi_u, &self.dpi_u, alphad),
        )
    }
}

/// The objective (times `cost_mult`), its gradient and the constraints' values at a point, with
/// the Jacobians when there are nonlinear constraints: without them, the Jacobians are the
/// constant `Ai` and `Ae`, and `Aiᵀ` is made once.
struct Evaluation {
    f: f64,
    df: Vec<f64>,
    /// `[h(x); Ai x − bi]`
    h: Vec<f64>,
    /// `[g(x); Ae x − be]`
    g: Vec<f64>,
    jacobians: Option<Jacobians>,
}

/// The constraints' Jacobians at a point with nonlinear constraints, nonlinear rows first
struct Jacobians {
    /// `[∂h/∂x; Ai]`
    jh: CsrMatrix<f64>,
    /// `[∂g/∂x; Ae]`
    jg: CsrMatrix<f64>,
    /// `jh`'s transpose
    jh_t: CsrMatrix<f64>,
}

impl Evaluation {
    fn jh<'a>(&'a self, split: &'a Split) -> &'a CsrMatrix<f64> {
        self.jacobians.as_ref().map_or(&split.ai, |j| &j.jh)
    }

    fn jg<'a>(&'a self, split: &'a Split) -> &'a CsrMatrix<f64> {
        self.jacobians.as_ref().map_or(&split.ae, |j| &j.jg)
    }

    fn jh_t<'a>(&'a self, ai_t: &'a CsrMatrix<f64>) -> &'a CsrMatrix<f64> {
        self.jacobians.as_ref().map_or(ai_t, |j| &j.jh_t)
    }
}

/// The predictor-corrector method: minimizes `problem`'s objective from `x0` subject to its
/// nonlinear constraints and to `constraints` (see the [module](self) for the method).
///
/// `x0` is moved strictly inside its bounds first: at least `bound_push · max(1, |bound|)`
/// inside each, and at most `bound_frac` of the distance between two bounds. The slacks start
/// at `max(z0, −h(x0))` and the multipliers at `max(z0, 1/slack)`, as MIPS starts them.
///
/// It fails numerically ([`Status::NumericalFailure`]) on a singular KKT system, a NaN or too
/// long a step, a step length below `alpha_min`, or the average complementarity outside
/// `[ε, 1/ε]`.
///
/// Returns `InteriorPointError::InvalidOption` for `xi` outside `[0.5, 1)`,
/// `LineSearchNeedsLdlt` for a line search with LU, and `Dimension` for inputs of the wrong
/// length.
pub fn mehrotra<P: Problem>(
    problem: &P,
    constraints: &Constraints,
    x0: &[f64],
    options: &MehrotraOptions,
) -> Result<Solution, InteriorPointError> {
    mehrotra_from(problem, constraints, Start::Cold(x0), options)
}

/// [`mehrotra`] from `start`: a cold start from a point, exactly [`mehrotra`]'s, or a warm start
/// from a previous solution's point and multipliers (see [`crate::start`]).
///
/// A warm start moves the point only `warm_shift` inside its bounds (at most `bound_frac` of
/// the distance between two), keeps the equality multipliers, recomputes the slacks from the
/// point, and raises the inequalities' slacks and every inequality multiplier, the bounds'
/// included, to at least `warm_shift`. The bound slacks stay the point's. Returns [`mehrotra`]'s errors,
/// `InvalidOption` for a negative or non-finite `warm_shift`, and `Dimension` for a warm start
/// whose point or multipliers don't fit the problem.
pub fn mehrotra_from<P: Problem>(
    problem: &P,
    constraints: &Constraints,
    start: Start<'_>,
    options: &MehrotraOptions,
) -> Result<Solution, InteriorPointError> {
    solve_with_hooks(problem, constraints, start, options, &Hooks::default())
}

/// [`mehrotra_from`] with `hooks`. Generic over `?Sized` problems, so that restoration's solve
///   runs as `solve_with_hooks::<dyn Problem>`, one instantiation however deeply it would
//   otherwise nest
fn solve_with_hooks<P: Problem + ?Sized>(
    problem: &P,
    constraints: &Constraints,
    start: Start<'_>,
    options: &MehrotraOptions,
    hooks: &Hooks<'_>,
) -> Result<Solution, InteriorPointError> {
    // The linear rows and bounds as given, for restoration's problem, before `constraints` is
    //   the nonlinear constraints' buffer below
    let linear = constraints;
    // The warm start, if any, named apart from the cold start's own `start` below
    let (x0, warm_start): (&[f64], _) = match start {
        Start::Cold(x0) => (x0, None),
        Start::Warm(warm) => {
            if !(options.warm_shift >= 0.0 && options.warm_shift.is_finite()) {
                return Err(InteriorPointError::InvalidOption { name: "warm_shift" });
            }
            (&warm.x, Some(warm))
        }
    };
    let opt = options;
    if !(0.5..1.0).contains(&opt.xi) {
        return Err(InteriorPointError::InvalidOption { name: "xi" });
    }
    let filter = opt.line_search != LineSearch::Off;
    if filter && opt.kkt_factorization != KktFactorization::Ldlt {
        return Err(InteriorPointError::LineSearchNeedsLdlt);
    }
    // A filter's steps must be descent directions, which the inertia correction ensures
    let inertia_correction = opt.inertia_correction || filter;
    // The adaptive barrier starts in free mode; the monotone filter is always in fixed mode
    let adaptive = opt.line_search == LineSearch::AdaptiveFilter;
    let mut fixed = filter && !adaptive;
    let nx = x0.len();
    if problem.dimension() != nx {
        return Err(InteriorPointError::Dimension { what: "x0" });
    }
    let split = Split::new(constraints, nx)?;
    // The nonlinear constraints' counts, from the starting point, as MIPS takes them
    let nonlinear = problem.has_nonlinear_constraints();
    let (neqnln, niqnln) = if nonlinear {
        let c = problem.constraints(x0);
        if c.dh.shape() != (c.h.len(), nx) || c.dg.shape() != (c.g.len(), nx) {
            return Err(InteriorPointError::Dimension {
                what: "the nonlinear constraints' Jacobians",
            });
        }
        (c.g.len(), c.h.len())
    } else {
        (0, 0)
    };
    let (neq, niq) = (neqnln + split.be.len(), niqnln + split.bi.len());
    let (nlo, nup) = (split.lo.len(), split.up.len());
    // `Ai` is constant, so its transpose is made once
    let ai_t = split.ai.clone().transpose().to_csr();

    // The starting point, strictly inside the bounds: a cold start `bound_push` relative to
    //   each bound, a warm start only `warm_shift`, so that a previous solution's active bounds
    //   move as little as its other complementary pairs
    let push_from = |bound: f64| match warm_start {
        None => opt.bound_push * bound.abs().max(1.0),
        Some(_) => opt.warm_shift,
    };
    let mut x = x0.to_vec();
    for (j, &k) in split.lo.iter().enumerate() {
        let push = push_from(split.xmin[j]);
        x[k] = x[k].max(split.xmin[j] + push);
    }
    for (j, &k) in split.up.iter().enumerate() {
        let push = push_from(split.xmax[j]);
        x[k] = x[k].min(split.xmax[j] - push);
    }
    for (j, &k) in split.lo.iter().enumerate() {
        if let Ok(i) = split.up.binary_search(&k) {
            // Both bounds: within bound_frac of the distance between them, at most
            let (lower, upper) = (split.xmin[j], split.xmax[i]);
            let inset = opt.bound_frac * (upper - lower);
            let push = |bound: f64| push_from(bound).min(inset);
            x[k] = x0[k].max(lower + push(lower)).min(upper - push(upper));
        }
    }
    let slacks = |x: &[f64]| -> (Vec<f64>, Vec<f64>) {
        (
            split
                .lo
                .iter()
                .zip(&split.xmin)
                .map(|(&k, b)| x[k] - b)
                .collect(),
            split
                .up
                .iter()
                .zip(&split.xmax)
                .map(|(&k, b)| b - x[k])
                .collect(),
        )
    };
    // The nonlinear constraints' buffer, kept across evaluations (Problem::constraints_into)
    let constraints = std::cell::RefCell::new(NonlinearConstraints::empty(nx));
    let evaluate = |x: &[f64]| -> Evaluation {
        let (f, df) = problem.objective(x);
        let df = df.iter().map(|d| d * opt.cost_mult).collect();
        let linear_h = split
            .ai
            .mul_vec(x)
            .expect("Ai has nx columns")
            .iter()
            .zip(&split.bi)
            .map(|(ax, b)| ax - b)
            .collect::<Vec<f64>>();
        let linear_g = split
            .ae
            .mul_vec(x)
            .expect("Ae has nx columns")
            .iter()
            .zip(&split.be)
            .map(|(ax, b)| ax - b)
            .collect::<Vec<f64>>();
        let nonlinear_constraints = nonlinear.then(|| {
            problem.constraints_into(x, &mut constraints.borrow_mut());
            constraints.borrow()
        });
        match nonlinear_constraints {
            None => Evaluation {
                f: f * opt.cost_mult,
                df,
                h: linear_h,
                g: linear_g,
                jacobians: None,
            },
            Some(c) => {
                let jh = CsrMatrix::vstack([&c.dh, &split.ai]).expect("∂h/∂x has nx columns");
                let jg = CsrMatrix::vstack([&c.dg, &split.ae]).expect("∂g/∂x has nx columns");
                Evaluation {
                    f: f * opt.cost_mult,
                    df,
                    h: c.h.iter().copied().chain(linear_h).collect(),
                    g: c.g.iter().copied().chain(linear_g).collect(),
                    jacobians: Some(Jacobians {
                        jh_t: jh.clone().transpose().to_csr(),
                        jh,
                        jg,
                    }),
                }
            }
        }
    };

    let mut e = evaluate(&x);
    let (s_l, s_u) = slacks(&x);
    let start = |slack: f64| opt.z0.max(slack);
    let z: Vec<f64> = e.h.iter().map(|&h| start(-h)).collect();
    let multiplier = |slack: f64| opt.z0.max(1.0 / slack);
    let mut p = Pairs {
        mu: z.iter().map(|&z| multiplier(z)).collect(),
        z,
        pi_l: s_l.iter().map(|&s| multiplier(s)).collect(),
        s_l,
        pi_u: s_u.iter().map(|&s| multiplier(s)).collect(),
        s_u,
    };
    let mut lam = vec![0.0; neq];
    if let Some(warm) = warm_start {
        (lam, p) = warm_pairs(
            &start::scaled(&warm.multipliers, opt.cost_mult),
            &e.h,
            &split,
            (neqnln, niqnln),
            p,
            opt.warm_shift,
        )?;
    }
    let pairs = p.count();

    // The Lagrangian's gradient, r_d = ∇f + Jgᵀλ + Jhᵀμ − π_l + π_u
    let lagrangian_gradient = |e: &Evaluation, lam: &[f64], p: &Pairs| -> Vec<f64> {
        let mut r: Vec<f64> =
            e.df.iter()
                .zip(
                    e.jg(&split)
                        .vec_mul(lam)
                        .expect("one multiplier per equality"),
                )
                .zip(
                    e.jh(&split)
                        .vec_mul(&p.mu)
                        .expect("one multiplier per inequality"),
                )
                .map(|((d, a), b)| d + a + b)
                .collect();
        for (&k, pi) in split.lo.iter().zip(&p.pi_l) {
            r[k] -= pi;
        }
        for (&k, pi) in split.up.iter().zip(&p.pi_u) {
            r[k] += pi;
        }
        r
    };
    // MIPS's measures, with the bounds counted as MIPS's bound rows
    let conditions = |x: &[f64], g: &[f64], h: &[f64], lam: &[f64], p: &Pairs, rd: &[f64]| {
        let maxh = h
            .iter()
            .chain(&p.s_l)
            .chain(&p.s_u)
            .enumerate()
            .map(|(i, &v)| if i < h.len() { v } else { -v })
            .fold(f64::NEG_INFINITY, f64::max);
        let infeasibility = if pairs == 0 {
            norm_inf(g)
        } else {
            norm_inf(g).max(maxh)
        };
        let slack_norm = norm_inf(&p.z).max(norm_inf(&p.s_l)).max(norm_inf(&p.s_u));
        let multiplier_norm = norm_inf(lam)
            .max(norm_inf(&p.mu))
            .max(norm_inf(&p.pi_l))
            .max(norm_inf(&p.pi_u));
        let feascond = infeasibility / (1.0 + norm_inf(x).max(slack_norm));
        let gradcond = norm_inf(rd) / (1.0 + multiplier_norm);
        let compcond = p.complementarity() / (1.0 + norm_inf(x));
        (feascond, gradcond, compcond)
    };
    let average = |p: &Pairs| {
        if pairs == 0 {
            0.0
        } else {
            p.complementarity() / pairs as f64
        }
    };

    let mut rd = lagrangian_gradient(&e, &lam, &p);
    let (mut feascond, mut gradcond, mut compcond) = conditions(&x, &e.g, &e.h, &lam, &p, &rd);
    let mut f0 = e.f;
    let mut costcond = 0.0;
    // θ = ‖[g; h + z]‖₁, the constraint violation the filter measures
    let theta = |e: &Evaluation, p: &Pairs| -> f64 {
        e.g.iter().map(|g| g.abs()).sum::<f64>()
            + e.h
                .iter()
                .zip(&p.z)
                .map(|(h, z)| (h + z).abs())
                .sum::<f64>()
    };
    let mut history = vec![MehrotraIteration {
        feascond,
        gradcond,
        compcond,
        costcond,
        mu: average(&p),
        obj: e.f / opt.cost_mult,
        theta: theta(&e, &p),
        ..MehrotraIteration::default()
    }];
    let tolerances_met = |feas: f64, grad: f64, comp: f64, cost: f64| {
        feas < opt.feastol && grad < opt.gradtol && comp < opt.comptol && cost < opt.costtol
    };
    let mut converged = tolerances_met(feascond, gradcond, compcond, costcond);
    let mut failed = false;
    let mut i = 0;
    let mut kkt = KktSolver::new(opt.kkt_factorization, opt.kkt_ordering);
    let mut gram = WeightedGram::default();
    let mut hessian_buffer = CsrMatrix::from_triplets(0, 0, vec![]).expect("an empty matrix");
    // The inertia correction's last shift, and how many refactorizations it took
    let mut last_delta_w = 0.0;
    let mut corrections = 0;
    // With a filter, the fixed barrier parameter μ̄ (filter.rs). Its floor keeps the method's
    //   convergence test reachable: that test sums complementarity over every pair, as MIPS's
    //   does, where the paper's floor, ε/10, bounds each pair's
    let mu_floor = opt.comptol / (10.0 * pairs.max(1) as f64);
    let mut mu_bar = hooks.mu_start.unwrap_or(match warm_start {
        None => filter::MU_INIT,
        Some(_) => average(&p).max(mu_floor),
    });
    // The filter, and the bounds on θ it works between, from the start's (§2.3, with the
    //   paper's defaults)
    let theta_start = theta(&e, &p).max(1.0);
    let theta_min = 1e-4 * theta_start;
    let mut line_filter = filter::Filter::new(1e4 * theta_start);
    // The adaptive barrier's KKT error Φ, and the accepted values its progress test compares
    //   against (the checklist's PR 3 rules)
    let kkt_error = |rd: &[f64], e: &Evaluation, p: &Pairs| -> f64 {
        let hz: Vec<f64> = e.h.iter().zip(&p.z).map(|(h, z)| h + z).collect();
        filter::kkt_error(
            rd,
            [&e.g, &hz],
            [(&p.z, &p.mu), (&p.s_l, &p.pi_l), (&p.s_u, &p.pi_u)],
        )
    };
    let mut progress = filter::Progress::new(kkt_error(&rd, &e, &p));
    if adaptive {
        history[0].kkt_error = progress.reference();
    }

    while !converged && i < opt.max_it {
        i += 1;
        // This iteration's mode: the adaptive barrier's may change at its end
        let step_fixed = fixed;

        // The reduced KKT matrix, [W Jgᵀ; Jg 0] with W = ∇²L + Jhᵀ Z⁻¹ M Jh + Σ_l + Σ_u
        let (jh, jg, jh_t) = (e.jh(&split), e.jg(&split), e.jh_t(&ai_t));
        // Into the buffer kept across iterations; scaling in place is `scale`'s product
        if nonlinear {
            problem.lagrangian_hessian_into(
                &x,
                NonlinearMultipliers {
                    eq: &lam[..neqnln],
                    ineq: &p.mu[..niqnln],
                },
                opt.cost_mult,
                &mut hessian_buffer,
            );
        } else {
            problem.objective_hessian_into(&x, &mut hessian_buffer);
            for value in hessian_buffer.values_mut() {
                *value *= opt.cost_mult;
            }
        }
        let hessian = &hessian_buffer;
        // With curvature, the dual residual ∇f(x) + … couples x to the multipliers, so the
        //   primal and dual steps share one length (Nocedal and Wright, §16.6), and so do they
        //   with nonlinear constraints, whose Jacobians couple them too; a linear program keeps
        //   separate ones
        let curved = nonlinear || hessian.values().iter().any(|&v| v != 0.0);
        let zinv: Vec<f64> = p.z.iter().map(|z| 1.0 / z).collect();
        let weights: Vec<f64> = p.mu.iter().zip(&zinv).map(|(m, zi)| m * zi).collect();
        // W = H + Jhᵀ (D Jh), with D Jh as jh.scale_rows(&weights), through the kept plans
        let jh_weighted: Vec<f64> = (jh.values().iter())
            .zip(gram.jh_rows(hessian, jh))
            .map(|(&v, &row)| weights[row] * v)
            .collect();
        let w = gram.compute(hessian, jh, &jh_weighted);
        // The bounds' terms on W's diagonal
        let bounds: Vec<(usize, f64)> = (split.lo.iter().enumerate())
            .map(|(j, &k)| (k, p.pi_l[j] / p.s_l[j]))
            .chain(
                split
                    .up
                    .iter()
                    .enumerate()
                    .map(|(j, &k)| (k, p.pi_u[j] / p.s_u[j])),
            )
            .collect();
        // With LDLᵀ, the inertia must be (nx, neq): W positive definite on Jg's null space, so
        //   that the step is a descent direction. Otherwise W is shifted by δ_w I, more each
        //   time, as Wächter and Biegler's Algorithm IC does (see the module docs)
        let wanted = Inertia {
            positive: nx,
            negative: neq,
        };
        let mut delta_w = 0.0;
        let factor = loop {
            let Ok(factor) = kkt.factor(nx, &w, &bounds, jg, delta_w) else {
                break None;
            };
            if !inertia_correction || factor.inertia().is_none_or(|inertia| inertia == wanted) {
                break Some(factor);
            }
            kkt.recycle(factor);
            corrections += 1;
            delta_w = if delta_w == 0.0 {
                if last_delta_w == 0.0 {
                    DELTA_W_FIRST
                } else {
                    (DELTA_W_DECREASE * last_delta_w).max(DELTA_W_MIN)
                }
            } else if last_delta_w == 0.0 {
                DELTA_W_FIRST_INCREASE * delta_w
            } else {
                DELTA_W_INCREASE * delta_w
            };
            if delta_w > DELTA_W_MAX {
                break None;
            }
        };
        let Some(factor) = factor else {
            failed = true;
            break;
        };
        if delta_w > 0.0 {
            last_delta_w = delta_w;
        }

        // A step for complementarity right-hand sides r_z, r_l and r_u
        let hz: Vec<f64> = e.h.iter().zip(&p.z).map(|(h, z)| h + z).collect();
        // For constraint residuals `hz` (h + z) and `g`: the current point's, or a second-order
        //   correction's
        let solve_for =
            |r_z: &[f64], r_l: &[f64], r_u: &[f64], hz: &[f64], g: &[f64]| -> Option<Step> {
                // −r_d − Jhᵀ Z⁻¹ (r_z + M (h + z)) + S_l⁻¹ r_l − S_u⁻¹ r_u, then −g
                let inner: Vec<f64> = (0..niq)
                    .map(|k| zinv[k] * (r_z[k] + p.mu[k] * hz[k]))
                    .collect();
                let mut rhs: Vec<f64> = rd
                    .iter()
                    .zip(jh_t.mul_vec(&inner).expect("one entry per inequality"))
                    .map(|(r, a)| -r - a)
                    .collect();
                for (j, &k) in split.lo.iter().enumerate() {
                    rhs[k] += r_l[j] / p.s_l[j];
                }
                for (j, &k) in split.up.iter().enumerate() {
                    rhs[k] -= r_u[j] / p.s_u[j];
                }
                rhs.extend(g.iter().map(|g| -g));
                let solved = factor.solve(&rhs).ok()?;
                if solved.iter().any(|v| !v.is_finite()) {
                    return None;
                }
                let dx = solved[..nx].to_vec();
                let dlam = solved[nx..].to_vec();
                let jh_dx = jh.mul_vec(&dx).expect("one entry per variable");
                let dz: Vec<f64> = (0..niq).map(|k| -hz[k] - jh_dx[k]).collect();
                let dmu = (0..niq)
                    .map(|k| zinv[k] * (r_z[k] - p.mu[k] * dz[k]))
                    .collect();
                let dpi_l = (0..nlo)
                    .map(|j| (r_l[j] - p.pi_l[j] * dx[split.lo[j]]) / p.s_l[j])
                    .collect();
                let dpi_u = (0..nup)
                    .map(|j| (r_u[j] + p.pi_u[j] * dx[split.up[j]]) / p.s_u[j])
                    .collect();
                Some(Step {
                    dx,
                    dlam,
                    dz,
                    dmu,
                    dpi_l,
                    dpi_u,
                })
            };
        let solve = |r_z: &[f64], r_l: &[f64], r_u: &[f64]| solve_for(r_z, r_l, r_u, &hz, &e.g);
        // −a ∘ b
        let negated =
            |a: &[f64], b: &[f64]| -> Vec<f64> { a.iter().zip(b).map(|(a, b)| -a * b).collect() };

        // Fixed mode's complementarity right-hand sides, centered at μ̄: kept for a second-order
        //   correction, which solves with the same ones
        let centered = fixed.then(|| {
            let centered = |s: &[f64], v: &[f64]| -> Vec<f64> {
                s.iter().zip(v).map(|(s, v)| mu_bar - s * v).collect()
            };
            (
                centered(&p.z, &p.mu),
                centered(&p.s_l, &p.pi_l),
                centered(&p.s_u, &p.pi_u),
            )
        });
        let (step, sigma) = if pairs == 0 {
            // No inequalities or bounds: a plain Newton step
            (solve(&[], &[], &[]), 0.0)
        } else if let Some((r_z, r_l, r_u)) = &centered {
            // Fixed mode: the centered Newton step at the barrier parameter μ̄
            (solve(r_z, r_l, r_u), 0.0)
        } else {
            // The predictor: an affine step towards zero complementarity
            let mu_now = average(&p);
            let r_z = negated(&p.z, &p.mu);
            let r_l = negated(&p.s_l, &p.pi_l);
            let r_u = negated(&p.s_u, &p.pi_u);
            match solve(&r_z, &r_l, &r_u) {
                None => (None, 0.0),
                Some(affine) => {
                    let (ap, ad) = affine.lengths(&split, &p, 1.0);
                    let mu_aff = affine.complementarity_after(&split, &p, ap, ad) / pairs as f64;
                    let sigma = (mu_aff / mu_now).powi(3).clamp(0.0, 1.0);
                    // The corrector: centering at σμ, with the predictor's second-order term,
                    //   damped by the affine step lengths when there are nonlinear constraints
                    //   (see the module docs)
                    let target = sigma * mu_now;
                    let damping = if nonlinear { ap * ad } else { 1.0 };
                    let r_z: Vec<f64> = (0..niq)
                        .map(|k| r_z[k] - damping * affine.dz[k] * affine.dmu[k] + target)
                        .collect();
                    let r_l: Vec<f64> = (0..nlo)
                        .map(|j| {
                            r_l[j] - damping * affine.dx[split.lo[j]] * affine.dpi_l[j] + target
                        })
                        .collect();
                    // s_u's step is −dx
                    let r_u: Vec<f64> = (0..nup)
                        .map(|j| {
                            r_u[j] + damping * affine.dx[split.up[j]] * affine.dpi_u[j] + target
                        })
                        .collect();
                    (solve(&r_z, &r_l, &r_u), sigma)
                }
            }
        };
        let length = |s: &Step| (dot(&s.dx, &s.dx) + dot(&s.dlam, &s.dlam)).sqrt();
        let step = match step {
            Some(step) if length(&step) <= opt.max_stepsize => Some(step),
            // In fixed mode, a step that couldn't be computed, or that's too long, as near a
            //   point where the linearized constraints are inconsistent with the bounds, goes to
            //   feasibility restoration, as a failed line search does
            _ if fixed && !hooks.no_restoration => None,
            // In free mode, the adaptive barrier has no trial point to start fixed mode from, so
            //   it starts from the current one
            _ if adaptive && !fixed => {
                mu_bar = (filter::MONOTONE_INIT_FACTOR * average(&p)).max(mu_floor);
                line_filter.reset();
                fixed = true;
                history.push(MehrotraIteration {
                    feascond,
                    gradcond,
                    compcond,
                    costcond,
                    mu: average(&p),
                    obj: e.f / opt.cost_mult,
                    theta: theta(&e, &p),
                    kkt_error: kkt_error(&rd, &e, &p),
                    ..MehrotraIteration::default()
                });
                continue;
            }
            _ => {
                failed = true;
                break;
            }
        };

        // The step, at most xi of the way to the boundary; in fixed mode, τ of the way, with the
        //   primal and dual lengths apart, as the paper takes them (§2.2), the primal one from the
        //   filter line search
        let tau = filter::fraction(mu_bar);
        let (step, alphap, alphad, trial, backtracks, second_order) = if fixed {
            let theta_k = theta(&e, &p);
            let phi_k = filter::barrier(e.f, mu_bar, [&p.z, &p.s_l, &p.s_u]);
            let (found, backtracks) = match &step {
                None => (None, 0),
                Some(step) => {
                    // ∇φ_μᵀd over x and z: ∇fᵀdx − μ̄ (Σ dx/s_l − Σ dx/s_u + Σ dz/z)
                    let grad_phi_d = {
                        let mut d = dot(&e.df, &step.dx);
                        for (j, &k) in split.lo.iter().enumerate() {
                            d -= mu_bar * step.dx[k] / p.s_l[j];
                        }
                        for (j, &k) in split.up.iter().enumerate() {
                            d += mu_bar * step.dx[k] / p.s_u[j];
                        }
                        for (dz, z) in step.dz.iter().zip(&p.z) {
                            d -= mu_bar * dz / z;
                        }
                        d
                    };
                    let alpha_min = filter::alpha_min(theta_k, grad_phi_d, theta_min);
                    // The trial point at primal length `alpha` along a step, if it's strictly interior
                    //   and its measures are finite
                    let trial_at = |s: &Step, alpha: f64| -> Option<Trial> {
                        let x_t: Vec<f64> =
                            x.iter().zip(&s.dx).map(|(x, d)| x + alpha * d).collect();
                        let z_t: Vec<f64> =
                            p.z.iter().zip(&s.dz).map(|(z, d)| z + alpha * d).collect();
                        let (s_l, s_u) = slacks(&x_t);
                        if z_t
                            .iter()
                            .chain(&s_l)
                            .chain(&s_u)
                            .any(|v| v.is_nan() || *v <= 0.0)
                        {
                            return None;
                        }
                        let e_t = evaluate(&x_t);
                        let theta_t = e_t.g.iter().map(|g| g.abs()).sum::<f64>()
                            + e_t
                                .h
                                .iter()
                                .zip(&z_t)
                                .map(|(h, z)| (h + z).abs())
                                .sum::<f64>();
                        let phi_t = filter::barrier(e_t.f, mu_bar, [&z_t, &s_l, &s_u]);
                        (theta_t.is_finite() && phi_t.is_finite()).then_some(Trial {
                            e: e_t,
                            z: z_t,
                            theta: theta_t,
                            phi: phi_t,
                        })
                    };
                    // Whether a trial at length `alpha` is accepted (§2.3), and if so whether by the
                    //   Armijo condition, which leaves the filter as it is
                    let accept = |alpha: f64, t: &Trial| -> Option<bool> {
                        if !line_filter.acceptable(t.theta, t.phi) {
                            return None;
                        }
                        if theta_k <= theta_min && filter::switching(alpha, grad_phi_d, theta_k) {
                            filter::armijo(t.phi, phi_k, alpha, grad_phi_d).then_some(true)
                        } else {
                            filter::sufficient_decrease(t.theta, t.phi, theta_k, phi_k)
                                .then_some(false)
                        }
                    };
                    let (alpha_max, _) = step.lengths(&split, &p, tau);
                    let mut alpha = alpha_max;
                    let mut backtracks = 0;
                    let mut found: Option<(Option<Step>, f64, Trial, bool)> = None;
                    while alpha >= alpha_min {
                        if let Some(t) = trial_at(step, alpha) {
                            if let Some(armijo) = accept(alpha, &t) {
                                found = Some((None, alpha, t, armijo));
                                break;
                            }
                            // The second-order correction (§2.4), when the first trial raises θ: the
                            //   same system with the constraints' residuals accumulated
                            if backtracks == 0 && t.theta >= theta_k {
                                let (r_z, r_l, r_u) = centered.as_ref().expect("fixed mode");
                                let mut c_g: Vec<f64> =
                                    e.g.iter()
                                        .zip(&t.e.g)
                                        .map(|(g, gt)| alpha * g + gt)
                                        .collect();
                                let mut c_h: Vec<f64> = (hz.iter().zip(&t.e.h).zip(&t.z))
                                    .map(|((hz, h), z)| alpha * hz + h + z)
                                    .collect();
                                let mut theta_old = t.theta;
                                for _ in 0..filter::P_MAX {
                                    let Some(soc) = solve_for(r_z, r_l, r_u, &c_h, &c_g) else {
                                        break;
                                    };
                                    let (alpha_soc, _) = soc.lengths(&split, &p, tau);
                                    let Some(ts) = trial_at(&soc, alpha_soc) else {
                                        break;
                                    };
                                    if let Some(armijo) = accept(alpha_max, &ts) {
                                        found = Some((Some(soc), alpha_soc, ts, armijo));
                                        break;
                                    }
                                    if ts.theta > filter::KAPPA_SOC * theta_old {
                                        break;
                                    }
                                    theta_old = ts.theta;
                                    for (c, g) in c_g.iter_mut().zip(&ts.e.g) {
                                        *c = alpha_soc * *c + g;
                                    }
                                    for ((c, h), z) in c_h.iter_mut().zip(&ts.e.h).zip(&ts.z) {
                                        *c = alpha_soc * *c + h + z;
                                    }
                                }
                                if found.is_some() {
                                    break;
                                }
                            }
                        }
                        alpha *= 0.5;
                        backtracks += 1;
                    }
                    (found, backtracks)
                }
            };
            // No acceptable step: feasibility restoration (§3.3), unless this is restoration's own
            //   solve
            let Some((correction, alphap, t, armijo)) = found else {
                if hooks.no_restoration {
                    failed = true;
                    break;
                }
                // The current point's margin first, so that restoration can't return to it
                line_filter.augment(theta_k, phi_k);
                let violation = (e.g.iter().chain(&hz)).map(|v| v.abs()).fold(0.0, f64::max);
                let mu_r = mu_bar.max(violation);
                let nlin = linear.a.as_ref().map_or(0, |a| a.shape().0);
                let restoration =
                    restoration::Restoration::new(problem, &x, mu_bar, (neqnln, niqnln, nlin));
                let linear_values = linear
                    .a
                    .as_ref()
                    .map_or_else(Vec::new, |a| a.mul_vec(&x).expect("A has nx columns"));
                let (r_constraints, y0) = restoration.start(
                    linear,
                    (&e.g[..neqnln], &e.h[..niqnln], &linear_values),
                    mu_r,
                );
                // The same method in fixed mode, from restoration's own start, left where it is
                let mut r_options = *opt;
                r_options.line_search = LineSearch::MonotoneFilter;
                r_options.kkt_factorization = KktFactorization::Ldlt;
                r_options.cost_mult = 1.0;
                r_options.bound_push = 1e-10;
                r_options.bound_frac = 1e-10;
                r_options.max_it = opt.max_it.saturating_sub(i).max(1);
                // The original problem at restoration's x: each slack max(−h, μ̄), and the
                //   point's evaluation, slacks, θ and φ
                let mapped = |x_t: &[f64]| -> Option<(Evaluation, Vec<f64>, f64, f64)> {
                    let (s_l, s_u) = slacks(x_t);
                    if s_l.iter().chain(&s_u).any(|s| s.is_nan() || *s <= 0.0) {
                        return None;
                    }
                    let e_t = evaluate(x_t);
                    let z_t: Vec<f64> = e_t.h.iter().map(|h| (-h).max(mu_bar)).collect();
                    let theta_t = e_t.g.iter().map(|g| g.abs()).sum::<f64>()
                        + e_t
                            .h
                            .iter()
                            .zip(&z_t)
                            .map(|(h, z)| (h + z).abs())
                            .sum::<f64>();
                    let phi_t = filter::barrier(e_t.f, mu_bar, [&z_t, &s_l, &s_u]);
                    Some((e_t, z_t, theta_t, phi_t))
                };
                // Restoration's end (§3.3): less violation, at a point the filter accepts
                let restored = |y: &[f64]| {
                    mapped(restoration.x(y)).is_some_and(|(_, _, theta_t, phi_t)| {
                        filter::restored(&line_filter, theta_t, phi_t, theta_k)
                    })
                };
                let outcome = solve_with_hooks(
                    &restoration as &dyn Problem,
                    &r_constraints,
                    Start::Cold(&y0),
                    &r_options,
                    &Hooks {
                        mu_start: Some(mu_r),
                        no_restoration: true,
                        stop: Some(&restored),
                    },
                );
                let r_iterations = outcome.as_ref().map_or(0, |s| s.iterations);
                i += r_iterations;
                let back = outcome.ok().filter(|s| restored(&s.x)).and_then(|s| {
                    let x_t = restoration.x(&s.x).to_vec();
                    mapped(&x_t).map(|m| (x_t, m))
                });
                // Restoration ended without less violation: a local minimizer of it, so locally
                //   infeasible, or restoration failed
                let Some((x_new, (e_new, z_new, theta_new, _))) = back else {
                    history.push(MehrotraIteration {
                        mode: IterationMode::Restoration,
                        restoration_iterations: r_iterations,
                        barrier_parameter: mu_bar,
                        ..MehrotraIteration::default()
                    });
                    failed = true;
                    break;
                };
                // Resuming at the restored point: the slacks and bound multipliers centered at
                //   μ̄, and λ by least squares (§3.6), or zero if that's too large
                x = x_new;
                e = e_new;
                p.z = z_new;
                (p.s_l, p.s_u) = slacks(&x);
                let centered = |s: &[f64]| -> Vec<f64> { s.iter().map(|s| mu_bar / s).collect() };
                (p.mu, p.pi_l, p.pi_u) = (centered(&p.z), centered(&p.s_l), centered(&p.s_u));
                filter::safeguard(&mut p.mu, &p.z, mu_bar);
                filter::safeguard(&mut p.pi_l, &p.s_l, mu_bar);
                filter::safeguard(&mut p.pi_u, &p.s_u, mu_bar);
                lam = restoration::least_squares_multipliers(
                    &lagrangian_gradient(&e, &vec![0.0; neq], &p),
                    e.jg(&split),
                    opt.kkt_ordering,
                );
                rd = lagrangian_gradient(&e, &lam, &p);
                (feascond, gradcond, compcond) = conditions(&x, &e.g, &e.h, &lam, &p, &rd);
                costcond = (e.f - f0).abs() / (1.0 + f0.abs());
                f0 = e.f;
                let phi = if adaptive {
                    kkt_error(&rd, &e, &p)
                } else {
                    0.0
                };
                history.push(MehrotraIteration {
                    feascond,
                    gradcond,
                    compcond,
                    costcond,
                    mu: average(&p),
                    obj: e.f / opt.cost_mult,
                    theta: theta_new,
                    mode: IterationMode::Restoration,
                    barrier_parameter: mu_bar,
                    filter_size: line_filter.len(),
                    restoration_iterations: r_iterations,
                    kkt_error: phi,
                    ..MehrotraIteration::default()
                });
                // The adaptive barrier returns to free mode from a restored point with progress
                if adaptive && progress.accepts(phi) {
                    progress.push(phi);
                    fixed = false;
                }
                converged = tolerances_met(feascond, gradcond, compcond, costcond);
                continue;
            };
            if !armijo {
                line_filter.augment(theta_k, phi_k);
            }
            let second_order = correction.is_some();
            let chosen = correction.or(step).expect("a line search's step");
            let (_, alphad) = chosen.lengths(&split, &p, tau);
            (chosen, alphap, alphad, Some(t.e), backtracks, second_order)
        } else {
            let step = step.expect("only fixed mode goes on without a step");
            let (alphap, alphad) = if pairs == 0 {
                (1.0, 1.0)
            } else if fixed {
                step.lengths(&split, &p, tau)
            } else {
                let (alphap, alphad) = step.lengths(&split, &p, opt.xi);
                if curved {
                    let alpha = alphap.min(alphad);
                    (alpha, alpha)
                } else {
                    (alphap, alphad)
                }
            };
            (step, alphap, alphad, None, 0, false)
        };
        // The factor's storage, for the next iteration's factorization
        kkt.recycle(factor);
        let advance = |values: &mut [f64], steps: &[f64], alpha: f64| {
            for (v, d) in values.iter_mut().zip(steps) {
                *v += alpha * d;
            }
        };
        advance(&mut x, &step.dx, alphap);
        advance(&mut p.z, &step.dz, alphap);
        advance(&mut p.mu, &step.dmu, alphad);
        // The equalities' multipliers take the primal length in fixed mode, as in the paper's
        //   iteration (§2.2)
        advance(&mut lam, &step.dlam, if fixed { alphap } else { alphad });
        advance(&mut p.pi_l, &step.dpi_l, alphad);
        advance(&mut p.pi_u, &step.dpi_u, alphad);
        (p.s_l, p.s_u) = slacks(&x);
        let step_mu_bar = mu_bar;
        if fixed {
            filter::safeguard(&mut p.mu, &p.z, mu_bar);
            filter::safeguard(&mut p.pi_l, &p.s_l, mu_bar);
            filter::safeguard(&mut p.pi_u, &p.s_u, mu_bar);
        }

        // The new point, evaluated already if it was a line search's trial point
        e = match trial {
            Some(evaluated) => evaluated,
            None => evaluate(&x),
        };
        rd = lagrangian_gradient(&e, &lam, &p);
        (feascond, gradcond, compcond) = conditions(&x, &e.g, &e.h, &lam, &p, &rd);
        costcond = (e.f - f0).abs() / (1.0 + f0.abs());
        let mu_new = average(&p);
        let barrier = if fixed {
            filter::barrier(e.f, step_mu_bar, [&p.z, &p.s_l, &p.s_u])
        } else {
            0.0
        };
        // The barrier parameter falls once the barrier problem is solved to κ_ε μ̄ (§3.1):
        //   its error is the method's feasibility and gradient measures, with complementarity
        //   measured against μ̄
        if fixed {
            loop {
                let comp = filter::complementarity_error(
                    [(&p.z, &p.mu), (&p.s_l, &p.pi_l), (&p.s_u, &p.pi_u)],
                    mu_bar,
                ) / (1.0 + norm_inf(&x));
                let error = feascond.max(gradcond).max(comp);
                if error > filter::KAPPA_EPSILON * mu_bar || mu_bar <= mu_floor {
                    break;
                }
                mu_bar = filter::decrease(mu_bar, mu_floor);
                // A new barrier problem, with a new filter (§3.1)
                line_filter.reset();
            }
        }
        let phi = if adaptive {
            kkt_error(&rd, &e, &p)
        } else {
            0.0
        };
        history.push(MehrotraIteration {
            feascond,
            gradcond,
            compcond,
            costcond,
            mu: mu_new,
            sigma,
            stepsize: dot(&step.dx, &step.dx).sqrt(),
            obj: e.f / opt.cost_mult,
            alphap,
            alphad,
            delta_w,
            theta: theta(&e, &p),
            mode: if fixed {
                IterationMode::Fixed
            } else {
                IterationMode::Free
            },
            barrier,
            barrier_parameter: if fixed { step_mu_bar } else { 0.0 },
            backtracks,
            second_order,
            filter_size: line_filter.len(),
            kkt_error: phi,
            ..MehrotraIteration::default()
        });
        // The adaptive barrier's progress test: in free mode, a point without progress starts
        //   fixed mode there; in fixed mode, a point with progress returns to free mode
        if adaptive {
            if progress.accepts(phi) {
                progress.push(phi);
                fixed = false;
            } else if !fixed {
                mu_bar = (filter::MONOTONE_INIT_FACTOR * mu_new).max(mu_floor);
                line_filter.reset();
                fixed = true;
            }
        }
        if tolerances_met(feascond, gradcond, compcond, costcond) {
            converged = true;
        } else {
            // A free step of the adaptive barrier that stalls has started fixed mode instead
            if x.iter().any(|v| v.is_nan())
                || (pairs > 0
                    && !(adaptive && !step_fixed)
                    && (alphap < opt.alpha_min
                        || alphad < opt.alpha_min
                        || mu_new < f64::EPSILON
                        || mu_new > 1.0 / f64::EPSILON))
            {
                failed = true;
                break;
            }
            f0 = e.f;
        }
        // Restoration's solve ends as soon as its point would do for the original problem
        if hooks.stop.is_some_and(|stop| stop(&x)) {
            converged = true;
        }
    }
    let status = if failed {
        Status::NumericalFailure
    } else if converged {
        Status::Converged
    } else {
        Status::MaxIterations
    };

    let unscale = |v: &[f64]| -> Vec<f64> { v.iter().map(|v| v / opt.cost_mult).collect() };
    // The final iterate's multipliers, kept for a warm start before the report thresholds them
    //   (unscaled, as reported)
    let lam = unscale(&lam);
    let iterate_multipliers = report(
        &lam,
        &unscale(&p.mu),
        &unscale(&p.pi_l),
        &unscale(&p.pi_u),
        &split,
        nx,
        (neqnln, niqnln),
    );
    // Constraints well inside their limits lose a small multiplier, as in MIPS; then the scale
    //   comes off
    for (mu, &h) in p.mu.iter_mut().zip(&e.h) {
        if h < -opt.feastol && *mu < opt.mu_threshold {
            *mu = 0.0;
        }
    }
    for (pi, &s) in p
        .pi_l
        .iter_mut()
        .zip(&p.s_l)
        .chain(p.pi_u.iter_mut().zip(&p.s_u))
    {
        if s > opt.feastol && *pi < opt.mu_threshold {
            *pi = 0.0;
        }
    }
    let (mu, pi_l, pi_u) = (unscale(&p.mu), unscale(&p.pi_l), unscale(&p.pi_u));

    // The multipliers in MIPS's convention
    let multipliers = report(&lam, &mu, &pi_l, &pi_u, &split, nx, (neqnln, niqnln));
    Ok(Solution {
        x,
        f: e.f / opt.cost_mult,
        status,
        iterations: i,
        multipliers,
        iterate_multipliers,
        history: History::Mehrotra(history),
        analyses: kkt.analyses(),
        refinement_steps: kkt.refinement_steps(),
        kkt: kkt.statistics(),
        inertia_corrections: corrections,
    })
}

/// The predictor-corrector method for a quadratic or linear program: minimizes
/// `½ xᵀ H x + cᵀ x` subject to `constraints`, from `x0` (zeros if absent), as
/// [`qps_mips`](crate::mips::qps_mips) wraps MIPS.
///
/// Returns `InteriorPointError::UnboundedLinearProgram` for a linear program with no constraints
/// or bounds, and the errors of [`mehrotra`].
pub fn qps_mehrotra(
    qp: &QuadraticProblem,
    constraints: &Constraints,
    x0: Option<&[f64]>,
    options: &MehrotraOptions,
) -> Result<Solution, InteriorPointError> {
    let zeros = vec![0.0; qp.c.len()];
    qps_mehrotra_from(qp, constraints, Start::Cold(x0.unwrap_or(&zeros)), options)
}

/// [`qps_mehrotra`] from `start` (see [`mehrotra_from`]), with [`qps_mehrotra`]'s errors.
pub fn qps_mehrotra_from(
    qp: &QuadraticProblem,
    constraints: &Constraints,
    start: Start<'_>,
    options: &MehrotraOptions,
) -> Result<Solution, InteriorPointError> {
    let nx = qp.c.len();
    if qp.h.is_none()
        && constraints.a.as_ref().is_none_or(|a| a.shape().0 == 0)
        && constraints.xmin.is_none()
        && constraints.xmax.is_none()
    {
        return Err(InteriorPointError::UnboundedLinearProgram);
    }
    if let Some(h) = &qp.h {
        if h.shape() != (nx, nx) {
            return Err(InteriorPointError::Dimension { what: "H" });
        }
    }
    mehrotra_from(qp, constraints, start, options)
}

/// The multipliers in MIPS's reported layout, from the method's internal `λ`, `μ`, `π_l` and
/// `π_u`: the nonlinear constraints', then per row of A, and per variable, an equality's sign
/// deciding its side. `warm_pairs` is its inverse.
fn report(
    lam: &[f64],
    mu: &[f64],
    pi_l: &[f64],
    pi_u: &[f64],
    split: &Split,
    nx: usize,
    (neqnln, niqnln): (usize, usize),
) -> Multipliers {
    let (lam_lin, mu_lin) = (&lam[neqnln..], &mu[niqnln..]);
    let na = split.na;
    let mut mu_l = vec![0.0; na];
    let mut mu_u = vec![0.0; na];
    let (nlt, ngt, nbx) = (split.ilt.len(), split.igt.len(), split.ibx.len());
    for (k, &r) in split.ieq.iter().enumerate() {
        if lam_lin[k] < 0.0 {
            mu_l[r] = -lam_lin[k];
        } else if lam_lin[k] > 0.0 {
            mu_u[r] = lam_lin[k];
        }
    }
    for (k, &r) in split.ilt.iter().enumerate() {
        mu_u[r] = mu_lin[k];
    }
    for (k, &r) in split.igt.iter().enumerate() {
        mu_l[r] = mu_lin[nlt + k];
    }
    for (k, &r) in split.ibx.iter().enumerate() {
        mu_u[r] = mu_lin[nlt + ngt + k];
        mu_l[r] = mu_lin[nlt + ngt + nbx + k];
    }
    let mut lower = vec![0.0; nx];
    let mut upper = vec![0.0; nx];
    for (&k, &pi) in split.lo.iter().zip(pi_l) {
        lower[k] = pi;
    }
    for (&k, &pi) in split.up.iter().zip(pi_u) {
        upper[k] = pi;
    }
    for (j, &k) in split.fixed.iter().enumerate() {
        let lam = lam_lin[split.ieq.len() + j];
        if lam < 0.0 {
            lower[k] = -lam;
        } else if lam > 0.0 {
            upper[k] = lam;
        }
    }
    Multipliers {
        eq_nonlin: lam[..neqnln].to_vec(),
        ineq_nonlin: mu[..niqnln].to_vec(),
        mu_l,
        mu_u,
        lower,
        upper,
    }
}

/// A warm start's `λ` and pairs in the method's internal layout, from `multipliers` in MIPS's
/// reported layout, the inequalities' values `h` at the (already pushed) point, and the cold
/// start's pairs `cold`, whose bound slacks, the point's, it keeps: the inverse of the report
/// `mehrotra_from` builds at its end, with the inequalities' slacks and every inequality
/// multiplier, the bounds' included, raised to at least `shift`.
///
/// - `λ`: the nonlinear equalities', then each of A's equality rows (`ieq`) and each fixed
///   variable, as its upper side's multiplier less its lower side's.
/// - `μ`: the nonlinear inequalities', then A's rows as `Ai` stacks them: `ilt`, `−igt`, `ibx`,
///   `−ibx`.
/// - `π_l`, `π_u`: the bounded variables' `lower` and `upper`.
fn warm_pairs(
    multipliers: &Multipliers,
    h: &[f64],
    split: &Split,
    (neqnln, niqnln): (usize, usize),
    cold: Pairs,
    shift: f64,
) -> Result<(Vec<f64>, Pairs), InteriorPointError> {
    let m = multipliers;
    let nx = m.lower.len();
    if m.eq_nonlin.len() != neqnln
        || m.ineq_nonlin.len() != niqnln
        || m.mu_l.len() != split.na
        || m.mu_u.len() != split.na
        || m.upper.len() != nx
        || cold.s_l.len() != split.lo.len()
        || split
            .lo
            .iter()
            .chain(&split.up)
            .chain(&split.fixed)
            .any(|&k| k >= nx)
    {
        return Err(InteriorPointError::Dimension {
            what: "the warm start's multipliers",
        });
    }
    let raise = |values: &mut Vec<f64>| {
        for value in values.iter_mut() {
            *value = value.max(shift);
        }
    };
    let mut lam = m.eq_nonlin.clone();
    lam.extend(split.ieq.iter().map(|&r| m.mu_u[r] - m.mu_l[r]));
    lam.extend(split.fixed.iter().map(|&k| m.upper[k] - m.lower[k]));
    let mut mu = m.ineq_nonlin.clone();
    mu.extend(split.ilt.iter().map(|&r| m.mu_u[r]));
    mu.extend(split.igt.iter().map(|&r| m.mu_l[r]));
    mu.extend(split.ibx.iter().map(|&r| m.mu_u[r]));
    mu.extend(split.ibx.iter().map(|&r| m.mu_l[r]));
    debug_assert_eq!(mu.len(), h.len());
    raise(&mut mu);
    let z = h.iter().map(|&h| (-h).max(shift)).collect();
    let mut pi_l: Vec<f64> = split.lo.iter().map(|&k| m.lower[k]).collect();
    let mut pi_u: Vec<f64> = split.up.iter().map(|&k| m.upper[k]).collect();
    raise(&mut pi_l);
    raise(&mut pi_u);
    // The bound slacks stay the point's, x − xmin and xmax − x, as every iteration recomputes
    //   them: the push keeps them off the boundary, by warm_shift or, between two bounds closer
    //   than that allows, by bound_frac of the distance. Raised, a slack would overstate the room
    //   to the bound, and a step could cross it
    let (s_l, s_u) = (cold.s_l, cold.s_u);
    Ok((
        lam,
        Pairs {
            z,
            mu,
            s_l,
            pi_l,
            s_u,
            pi_u,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every kind of row of A (an equality, upper-only, lower-only, boxed and unbounded) and of
    ///   bound (lower-only, upper-only, both, fixed and none)
    fn split() -> Split {
        let inf = f64::INFINITY;
        let a = CsrMatrix::from_triplets(5, 5, (0..5).map(|r| (r, r, 1.0)).collect::<Vec<_>>())
            .unwrap();
        let constraints = Constraints {
            a: Some(a),
            l: Some(vec![1.0, -inf, -2.0, -3.0, -inf]),
            u: Some(vec![1.0, 2.0, inf, 3.0, inf]),
            xmin: Some(vec![0.0, -inf, -1.0, 0.5, -inf]),
            xmax: Some(vec![inf, 4.0, 1.0, 0.5, inf]),
        };
        Split::new(&constraints, 5).unwrap()
    }

    #[test]
    fn the_warm_mapping_is_the_reports_inverse() {
        let split = split();
        let (nx, nonlinear) = (5, (2, 3));
        let mut next = 0.5;
        let mut value = || {
            next += 1.25;
            next
        };
        let (mut mu_l, mut mu_u) = (vec![0.0; split.na], vec![0.0; split.na]);
        for (k, &r) in split.ieq.iter().enumerate() {
            if k % 2 == 0 {
                mu_l[r] = value()
            } else {
                mu_u[r] = value()
            }
        }
        for &r in &split.ilt {
            mu_u[r] = value();
        }
        for &r in &split.igt {
            mu_l[r] = value();
        }
        for &r in &split.ibx {
            (mu_l[r], mu_u[r]) = (value(), value());
        }
        let (mut lower, mut upper) = (vec![0.0; nx], vec![0.0; nx]);
        for &k in &split.lo {
            lower[k] = value();
        }
        for &k in &split.up {
            upper[k] = value();
        }
        for &k in &split.fixed {
            upper[k] = value();
        }
        let m = Multipliers {
            eq_nonlin: (0..nonlinear.0).map(|_| value() - 7.0).collect(),
            ineq_nonlin: (0..nonlinear.1).map(|_| value()).collect(),
            mu_l,
            mu_u,
            lower,
            upper,
        };
        assert!(!split.fixed.is_empty() && !split.ieq.is_empty() && !split.ibx.is_empty());
        let niq = nonlinear.1 + split.ilt.len() + split.igt.len() + 2 * split.ibx.len();
        let cold = Pairs {
            z: vec![1.0; niq],
            mu: vec![1.0; niq],
            s_l: vec![0.25; split.lo.len()],
            pi_l: vec![1.0; split.lo.len()],
            s_u: vec![0.75; split.up.len()],
            pi_u: vec![1.0; split.up.len()],
        };
        let (lam, p) = warm_pairs(&m, &vec![-2.0; niq], &split, nonlinear, cold, 0.0).unwrap();
        assert_eq!(p.z, vec![2.0; niq]);
        assert_eq!(
            (p.s_l[0], p.s_u[0]),
            (0.25, 0.75),
            "the bound slacks are the point's"
        );
        assert_eq!(
            report(&lam, &p.mu, &p.pi_l, &p.pi_u, &split, nx, nonlinear),
            m
        );
    }

    #[test]
    fn the_warm_start_moves_every_pair_off_the_boundary() {
        let split = split();
        let niq = split.ilt.len() + split.igt.len() + 2 * split.ibx.len();
        let m = Multipliers {
            mu_l: vec![0.0; split.na],
            mu_u: vec![0.0; split.na],
            lower: vec![0.0; 5],
            upper: vec![0.0; 5],
            ..Multipliers::default()
        };
        let cold = Pairs {
            z: vec![1.0; niq],
            mu: vec![1.0; niq],
            s_l: vec![1e-4; split.lo.len()],
            pi_l: vec![1.0; split.lo.len()],
            s_u: vec![1e-4; split.up.len()],
            pi_u: vec![1.0; split.up.len()],
        };
        let (_, p) = warm_pairs(&m, &vec![0.0; niq], &split, (0, 0), cold, 1e-3).unwrap();
        for values in [&p.z, &p.mu, &p.pi_l, &p.pi_u] {
            assert!(values.iter().all(|&v| v >= 1e-3), "{values:?}");
        }
        // The bound slacks stay the point's, even below the shift (between two close bounds):
        //   raised, they'd let a step cross its bound
        assert!(
            p.s_l.iter().chain(&p.s_u).all(|&s| s == 1e-4),
            "{:?} {:?}",
            p.s_l,
            p.s_u
        );
    }
}
