//! A primal-dual interior-point method with Mehrotra's predictor-corrector, written from the
//! literature (Mehrotra 1992; Nocedal and Wright, *Numerical Optimization*, 2nd ed., Algorithm
//! 14.3 and §16.6; S. J. Wright, *Primal-Dual Interior-Point Methods*, 1997) for the same
//! programs as [`mips`](crate::mips):
//!
//! minimize `f(x)` subject to `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.
//!
//! So far it takes problems without nonlinear constraints: unconstrained ones, linear and
//! quadratic programs, and any smooth objective under linear constraints and bounds.
//!
//! # The method
//!
//! The linear rows are split as MIPS splits them: an equality `Ae x = be` for `|u − l| ≤ ε`, and
//! otherwise one or two inequalities `Ai x − bi ≤ 0`, each with a slack `z > 0` and a multiplier
//! `μ > 0`. A bound side within ±1e10 is kept apart: `x` stays strictly inside the bounds, whose
//! slacks are `s_l = x − xmin` and `s_u = xmax − x`, with multipliers `π_l` and `π_u`. A
//! variable with `xmin = xmax` is an equality row, as in MIPS. Each iteration takes a Newton step
//! on the perturbed KKT conditions, reduced to
//!
//! `[W Jgᵀ; Jg 0] [dx; dλ] = [−r_d − Jhᵀ Z⁻¹ (r_z + M (h + z)) + S_l⁻¹ r_l − S_u⁻¹ r_u; −g]`
//!
//! with `W = ∇²f + Jhᵀ Z⁻¹ M Jh + S_l⁻¹ Π_l + S_u⁻¹ Π_u`: the bounds add only diagonal terms. Here
//! `r_d` is the Lagrangian's gradient, and `r_z`, `r_l` and `r_u` are the right-hand sides of
//! the linearized complementarity conditions. The matrix is factored once (LU ordered by
//! [`KKT_ORDERING`]) and solved twice:
//!
//! 1. **The predictor**, an affine step towards complementarity zero (`r_z = −Z μ`, and so on),
//!    whose longest feasible step lengths give the complementarity it would reach, `μ_aff`.
//! 2. **The corrector**, aiming at `σ μ` with `σ = (μ_aff / μ)³`, `μ` being the average
//!    complementarity, and adding the predictor's second-order term
//!    (`r_z = −Z μ − ΔZ_aff Δμ_aff + σ μ e`). Its primal and dual steps are each at most `xi` of the way to the boundary. A
//!    linear program takes them with separate lengths. Where the objective has curvature, both
//!    take the shorter of the two (Nocedal and Wright, §16.6): the dual residual `∇f(x) − …`
//!    depends on `x`, so separate lengths can let complementarity grow, as they did on a small
//!    QP with two-sided bounds, which failed to converge with them.
//!
//! A problem with no inequalities or bounds takes plain Newton steps. There is no merit function
//! or line search, as in MIPS. It converges by MIPS's four measures, defined alike, with the
//! bounds' slacks and multipliers counted as MIPS counts its bound rows, so iteration counts
//! compare fairly.

use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix, LuOptions};

use crate::interior_point_error::InteriorPointError;
use crate::mips::KKT_ORDERING;
use crate::problem::{Constraints, Problem, QuadraticProblem};
use crate::solution::Multipliers;
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
            bound_push: 1e-2,
            bound_frac: 1e-2,
            alpha_min: 1e-8,
            mu_threshold: 1e-5,
            max_stepsize: 1e10,
        }
    }
}

/// How a solve ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MehrotraStatus {
    /// Every termination test passed
    Converged,
    /// `max_it` iterations passed without converging
    MaxIterations,
    /// The KKT system was singular, a step was NaN or too long, a step length fell below
    /// `alpha_min`, or the average complementarity left `[ε, 1/ε]`
    NumericalFailure,
}

/// One iteration's record: the termination measures, the average complementarity and the
/// centering it chose, the step's 2-norm, the objective and the step lengths.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
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
}

/// The predictor-corrector method's result.
#[derive(Clone, Debug, PartialEq)]
pub struct MehrotraSolution {
    /// The final point
    pub x: Vec<f64>,
    /// The objective there
    pub f: f64,
    /// How the solve ended
    pub status: MehrotraStatus,
    /// The number of iterations taken
    pub iterations: usize,
    /// The multipliers, in MIPS's convention
    pub multipliers: Multipliers,
    /// One record for the starting point and one per iteration
    pub history: Vec<MehrotraIteration>,
}

impl MehrotraSolution {
    /// Whether the solve converged.
    pub fn converged(&self) -> bool {
        self.status == MehrotraStatus::Converged
    }
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

/// The predictor-corrector method: minimizes `problem`'s objective from `x0` subject to
/// `constraints` (see the [module](self) for the method).
///
/// `x0` is moved strictly inside its bounds first: at least `bound_push · max(1, |bound|)`
/// inside each, and at most `bound_frac` of the distance between two bounds. The slacks start
/// at `max(z0, −h(x0))` and the multipliers at `max(z0, 1/slack)`, as MIPS starts them.
///
/// Returns `InteriorPointError::Unsupported` for a problem with nonlinear constraints,
/// `InvalidOption` for `xi` outside `[0.5, 1)`, and `Dimension` for inputs of the wrong length.
pub fn mehrotra<P: Problem>(
    problem: &P,
    x0: &[f64],
    constraints: &Constraints,
    options: &MehrotraOptions,
) -> Result<MehrotraSolution, InteriorPointError> {
    let opt = options;
    if !(0.5..1.0).contains(&opt.xi) {
        return Err(InteriorPointError::InvalidOption { name: "xi" });
    }
    let nx = x0.len();
    if problem.dimension() != nx {
        return Err(InteriorPointError::Dimension { what: "x0" });
    }
    if problem.has_nonlinear_constraints() {
        return Err(InteriorPointError::Unsupported {
            what: "nonlinear constraints",
        });
    }
    let split = Split::new(constraints, nx)?;
    let (neq, niq) = (split.be.len(), split.bi.len());
    let (nlo, nup) = (split.lo.len(), split.up.len());
    // `Ai` and `Ae` are constant, so their transposes are made once
    let ai_t = split.ai.clone().transpose().to_csr();

    // The starting point, strictly inside the bounds
    let mut x = x0.to_vec();
    for (j, &k) in split.lo.iter().enumerate() {
        let push = opt.bound_push * split.xmin[j].abs().max(1.0);
        x[k] = x[k].max(split.xmin[j] + push);
    }
    for (j, &k) in split.up.iter().enumerate() {
        let push = opt.bound_push * split.xmax[j].abs().max(1.0);
        x[k] = x[k].min(split.xmax[j] - push);
    }
    for (j, &k) in split.lo.iter().enumerate() {
        if let Ok(i) = split.up.binary_search(&k) {
            // Both bounds: within bound_frac of the distance between them, at most
            let (lower, upper) = (split.xmin[j], split.xmax[i]);
            let inset = opt.bound_frac * (upper - lower);
            let push = |bound: f64| (opt.bound_push * bound.abs().max(1.0)).min(inset);
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
    let evaluate = |x: &[f64]| -> (f64, Vec<f64>, Vec<f64>, Vec<f64>) {
        let (f, df) = problem.objective(x);
        let df = df.iter().map(|d| d * opt.cost_mult).collect();
        let h: Vec<f64> = split
            .ai
            .mul_vec(x)
            .expect("Ai has nx columns")
            .iter()
            .zip(&split.bi)
            .map(|(ax, b)| ax - b)
            .collect();
        let g: Vec<f64> = split
            .ae
            .mul_vec(x)
            .expect("Ae has nx columns")
            .iter()
            .zip(&split.be)
            .map(|(ax, b)| ax - b)
            .collect();
        (f * opt.cost_mult, df, h, g)
    };

    let (mut f, mut df, mut h, mut g) = evaluate(&x);
    let (s_l, s_u) = slacks(&x);
    let start = |slack: f64| opt.z0.max(slack);
    let z: Vec<f64> = h.iter().map(|&h| start(-h)).collect();
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
    let pairs = p.count();

    // The Lagrangian's gradient, r_d = ∇f + Aeᵀλ + Aiᵀμ − π_l + π_u
    let lagrangian_gradient = |df: &[f64], lam: &[f64], p: &Pairs| -> Vec<f64> {
        let mut r: Vec<f64> = df
            .iter()
            .zip(split.ae.vec_mul(lam).expect("one multiplier per equality"))
            .zip(
                split
                    .ai
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

    let mut rd = lagrangian_gradient(&df, &lam, &p);
    let (mut feascond, mut gradcond, mut compcond) = conditions(&x, &g, &h, &lam, &p, &rd);
    let mut f0 = f;
    let mut costcond = 0.0;
    let mut history = vec![MehrotraIteration {
        feascond,
        gradcond,
        compcond,
        costcond,
        mu: average(&p),
        obj: f / opt.cost_mult,
        ..MehrotraIteration::default()
    }];
    let tolerances_met = |feas: f64, grad: f64, comp: f64, cost: f64| {
        feas < opt.feastol && grad < opt.gradtol && comp < opt.comptol && cost < opt.costtol
    };
    let mut converged = tolerances_met(feascond, gradcond, compcond, costcond);
    let mut failed = false;
    let mut i = 0;

    while !converged && i < opt.max_it {
        i += 1;

        // The reduced KKT matrix, [W Aeᵀ; Ae 0] with W = ∇²f + Aiᵀ Z⁻¹ M Ai + Σ_l + Σ_u
        let hessian = problem.objective_hessian(&x).scale(opt.cost_mult);
        // With curvature, the dual residual H x + c − … couples x to the multipliers, so the
        //   primal and dual steps share one length (Nocedal and Wright, §16.6); a linear
        //   program keeps separate ones
        let curved = hessian.values().iter().any(|&v| v != 0.0);
        let zinv: Vec<f64> = p.z.iter().map(|z| 1.0 / z).collect();
        let weights: Vec<f64> = p.mu.iter().zip(&zinv).map(|(m, zi)| m * zi).collect();
        let ai_weighted = split
            .ai
            .scale_rows(&weights)
            .expect("one weight per inequality");
        let w = hessian
            .try_add(&ai_t.matmul(&ai_weighted).expect("nx × niq times niq × nx"))
            .expect("both nx × nx");
        let mut triplets: Vec<(usize, usize, f64)> = w.iter().collect();
        for (j, &k) in split.lo.iter().enumerate() {
            triplets.push((k, k, p.pi_l[j] / p.s_l[j]));
        }
        for (j, &k) in split.up.iter().enumerate() {
            triplets.push((k, k, p.pi_u[j] / p.s_u[j]));
        }
        for (r, c, v) in split.ae.iter() {
            triplets.push((c, nx + r, v));
            triplets.push((nx + r, c, v));
        }
        let size = nx + neq;
        let kkt = CscMatrix::from_triplets(size, size, triplets).expect("indices in range");
        let Ok(lu) = kkt
            .analyze_lu(KKT_ORDERING)
            .and_then(|symbolic| symbolic.factor(&kkt, LuOptions::default()))
        else {
            failed = true;
            break;
        };

        // A step for complementarity right-hand sides r_z, r_l and r_u
        let hz: Vec<f64> = h.iter().zip(&p.z).map(|(h, z)| h + z).collect();
        let solve = |r_z: &[f64], r_l: &[f64], r_u: &[f64]| -> Option<Step> {
            // −r_d − Aiᵀ Z⁻¹ (r_z + M (h + z)) + S_l⁻¹ r_l − S_u⁻¹ r_u, then −g
            let inner: Vec<f64> = (0..niq)
                .map(|k| zinv[k] * (r_z[k] + p.mu[k] * hz[k]))
                .collect();
            let mut rhs: Vec<f64> = rd
                .iter()
                .zip(ai_t.mul_vec(&inner).expect("one entry per inequality"))
                .map(|(r, a)| -r - a)
                .collect();
            for (j, &k) in split.lo.iter().enumerate() {
                rhs[k] += r_l[j] / p.s_l[j];
            }
            for (j, &k) in split.up.iter().enumerate() {
                rhs[k] -= r_u[j] / p.s_u[j];
            }
            rhs.extend(g.iter().map(|g| -g));
            let solved = lu.solve(&rhs).ok()?;
            if solved.iter().any(|v| !v.is_finite()) {
                return None;
            }
            let dx = solved[..nx].to_vec();
            let dlam = solved[nx..].to_vec();
            let ai_dx = split.ai.mul_vec(&dx).expect("one entry per variable");
            let dz: Vec<f64> = (0..niq).map(|k| -hz[k] - ai_dx[k]).collect();
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
        // −a ∘ b
        let negated =
            |a: &[f64], b: &[f64]| -> Vec<f64> { a.iter().zip(b).map(|(a, b)| -a * b).collect() };

        let (step, sigma) = if pairs == 0 {
            // No inequalities or bounds: a plain Newton step
            (solve(&[], &[], &[]), 0.0)
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
                    // The corrector: centering at σμ, with the predictor's second-order term
                    let target = sigma * mu_now;
                    let r_z: Vec<f64> = (0..niq)
                        .map(|k| r_z[k] - affine.dz[k] * affine.dmu[k] + target)
                        .collect();
                    let r_l: Vec<f64> = (0..nlo)
                        .map(|j| r_l[j] - affine.dx[split.lo[j]] * affine.dpi_l[j] + target)
                        .collect();
                    // s_u's step is −dx
                    let r_u: Vec<f64> = (0..nup)
                        .map(|j| r_u[j] + affine.dx[split.up[j]] * affine.dpi_u[j] + target)
                        .collect();
                    (solve(&r_z, &r_l, &r_u), sigma)
                }
            }
        };
        let Some(step) = step else {
            failed = true;
            break;
        };
        let stepsize = (dot(&step.dx, &step.dx) + dot(&step.dlam, &step.dlam)).sqrt();
        if stepsize > opt.max_stepsize {
            failed = true;
            break;
        }

        // The step, at most xi of the way to the boundary
        let (alphap, alphad) = if pairs == 0 {
            (1.0, 1.0)
        } else {
            let (alphap, alphad) = step.lengths(&split, &p, opt.xi);
            if curved {
                let alpha = alphap.min(alphad);
                (alpha, alpha)
            } else {
                (alphap, alphad)
            }
        };
        let advance = |values: &mut [f64], steps: &[f64], alpha: f64| {
            for (v, d) in values.iter_mut().zip(steps) {
                *v += alpha * d;
            }
        };
        advance(&mut x, &step.dx, alphap);
        advance(&mut p.z, &step.dz, alphap);
        advance(&mut p.mu, &step.dmu, alphad);
        advance(&mut lam, &step.dlam, alphad);
        advance(&mut p.pi_l, &step.dpi_l, alphad);
        advance(&mut p.pi_u, &step.dpi_u, alphad);
        (p.s_l, p.s_u) = slacks(&x);

        // The new point
        (f, df, h, g) = evaluate(&x);
        rd = lagrangian_gradient(&df, &lam, &p);
        (feascond, gradcond, compcond) = conditions(&x, &g, &h, &lam, &p, &rd);
        costcond = (f - f0).abs() / (1.0 + f0.abs());
        let mu_new = average(&p);
        history.push(MehrotraIteration {
            feascond,
            gradcond,
            compcond,
            costcond,
            mu: mu_new,
            sigma,
            stepsize: dot(&step.dx, &step.dx).sqrt(),
            obj: f / opt.cost_mult,
            alphap,
            alphad,
        });
        if tolerances_met(feascond, gradcond, compcond, costcond) {
            converged = true;
        } else {
            if x.iter().any(|v| v.is_nan())
                || (pairs > 0
                    && (alphap < opt.alpha_min
                        || alphad < opt.alpha_min
                        || mu_new < f64::EPSILON
                        || mu_new > 1.0 / f64::EPSILON))
            {
                failed = true;
                break;
            }
            f0 = f;
        }
    }
    let status = if failed {
        MehrotraStatus::NumericalFailure
    } else if converged {
        MehrotraStatus::Converged
    } else {
        MehrotraStatus::MaxIterations
    };

    // Constraints well inside their limits lose a small multiplier, as in MIPS; then the scale
    //   comes off
    for (mu, &h) in p.mu.iter_mut().zip(&h) {
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
    let unscale = |v: &[f64]| -> Vec<f64> { v.iter().map(|v| v / opt.cost_mult).collect() };
    let (lam, mu, pi_l, pi_u) = (
        unscale(&lam),
        unscale(&p.mu),
        unscale(&p.pi_l),
        unscale(&p.pi_u),
    );

    // The multipliers in MIPS's convention: per row of A, and per variable
    let na = split.na;
    let mut mu_l = vec![0.0; na];
    let mut mu_u = vec![0.0; na];
    let (nlt, ngt, nbx) = (split.ilt.len(), split.igt.len(), split.ibx.len());
    for (k, &r) in split.ieq.iter().enumerate() {
        if lam[k] < 0.0 {
            mu_l[r] = -lam[k];
        } else if lam[k] > 0.0 {
            mu_u[r] = lam[k];
        }
    }
    for (k, &r) in split.ilt.iter().enumerate() {
        mu_u[r] = mu[k];
    }
    for (k, &r) in split.igt.iter().enumerate() {
        mu_l[r] = mu[nlt + k];
    }
    for (k, &r) in split.ibx.iter().enumerate() {
        mu_u[r] = mu[nlt + ngt + k];
        mu_l[r] = mu[nlt + ngt + nbx + k];
    }
    let mut lower = vec![0.0; nx];
    let mut upper = vec![0.0; nx];
    for (&k, &pi) in split.lo.iter().zip(&pi_l) {
        lower[k] = pi;
    }
    for (&k, &pi) in split.up.iter().zip(&pi_u) {
        upper[k] = pi;
    }
    for (j, &k) in split.fixed.iter().enumerate() {
        let lam = lam[split.ieq.len() + j];
        if lam < 0.0 {
            lower[k] = -lam;
        } else if lam > 0.0 {
            upper[k] = lam;
        }
    }
    Ok(MehrotraSolution {
        x,
        f: f / opt.cost_mult,
        status,
        iterations: i,
        multipliers: Multipliers {
            eq_nonlin: vec![],
            ineq_nonlin: vec![],
            mu_l,
            mu_u,
            lower,
            upper,
        },
        history,
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
) -> Result<MehrotraSolution, InteriorPointError> {
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
    let zeros = vec![0.0; nx];
    mehrotra(qp, x0.unwrap_or(&zeros), constraints, options)
}
