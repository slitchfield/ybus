// Derived from MIPS's mips.m and qps_mips.m (MIPS 1.5.2):
//   Copyright (c) 2009-2024, Power Systems Engineering Research Center (PSERC)
//   by Ray Zimmerman, PSERC Cornell
//   Covered by the 3-clause BSD License; see LICENSE-MIPS at the repository's root.

//! MIPS, the MATPOWER Interior Point Solver: a primal-dual interior-point method for nonlinear
//! programs, ported step for step from MATPOWER's BSD-licensed `mips.m` and `qps_mips.m`
//! (MIPS 1.5.2, in MATPOWER 8.1) so that its iterations, results and multipliers match
//! MATPOWER's (Wang, Murillo-Sánchez, Zimmerman and Thomas, 2007).

use ybus_sparse::compressed_matrix::{CsrMatrix, Ordering};

use crate::interior_point_error::InteriorPointError;
use crate::kkt::{KktFactorization, KktOrdering, KktSolver};
use crate::newton::WeightedGram;
use crate::problem::{
    Constraints, NonlinearConstraints, NonlinearMultipliers, Problem, QuadraticProblem,
};
use crate::solution::{History, Multipliers, Solution, Status};
use crate::start::{self, Start};
use crate::vector::{dot, norm_inf};

/// The fill-reducing ordering for the KKT system's LU factorization: COLAMD.
///
/// The system's lower-right block, `0` for the equality constraints, has zero diagonals, which
/// the LU must pivot away from by rows. AMD's symmetric ordering can't anticipate that, while
/// COLAMD bounds the fill for any row pivoting. On `case_ACTIVSg10k`'s DC OPF (a 22,949 × 22,949
/// KKT system with 103,095 entries), AMD's factors held 14 million entries and took 24 s per
/// iteration, and COLAMD's 640,000, in 0.036 s (measured 2026-11-04). On MATPOWER's cases,
/// both reproduce MATPOWER's iteration counts alike except on case2736sp's degenerate DC OPF,
/// where COLAMD takes 44 to MATPOWER's 45 and AMD 45. Neither is UMFPACK's own ordering, which
/// MATPOWER uses under Octave; a fuller study of the KKT solve is a follow-up in
/// `docs/OptimizationChecklist.md`.
///
/// The solvers analyze the KKT matrix (this ordering, and the symbolic LU) once, and keep the
/// analysis for every later iteration whose KKT matrix has the same pattern, which is all of
/// them when the problem's Jacobians and Hessians keep theirs: the ordering is chosen from the
/// pattern alone, so a kept analysis factors exactly as a fresh one would.
pub const KKT_ORDERING: Ordering = Ordering::Colamd;

/// MIPS's options, with its defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct MipsOptions {
    /// Termination tolerance for feasibility (`feastol`)
    pub feastol: f64,
    /// Termination tolerance for the gradient (`gradtol`)
    pub gradtol: f64,
    /// Termination tolerance for complementarity (`comptol`)
    pub comptol: f64,
    /// Termination tolerance for the cost (`costtol`)
    pub costtol: f64,
    /// The most iterations (`max_it`)
    pub max_it: usize,
    /// Whether to reduce a step that worsens both feasibility and the gradient
    /// (`step_control`)
    pub step_control: bool,
    /// The most step-size reductions per step, with step control (`sc.red_it`)
    pub red_it: usize,
    /// The objective's scale factor (`cost_mult`)
    pub cost_mult: f64,
    /// The fraction of the distance to the boundary a step may take (`xi`), in `[0.5, 1)`
    pub xi: f64,
    /// The centering parameter (`sigma`), in `(0, 1]`
    pub sigma: f64,
    /// The slack variables' starting value (`z0`)
    pub z0: f64,
    /// A warm start's shift ([`Start::Warm`]): every slack and inequality multiplier starts at
    /// least this far from zero, in the solver's scaled units (the objective's times
    /// `cost_mult`). Unused by a cold start. The default,
    /// 1e-4, is provisional until the OPF sequences' survey measures it.
    pub warm_shift: f64,
    /// The smallest acceptable primal or dual step length (`alpha_min`)
    pub alpha_min: f64,
    /// The bounds on step control's ratio of actual to predicted change (`rho_min`, `rho_max`)
    pub rho_min: f64,
    /// See `rho_min`
    pub rho_max: f64,
    /// The threshold below which an inequality's multiplier is set to zero when the inequality
    /// isn't binding (`mu_threshold`)
    pub mu_threshold: f64,
    /// The largest acceptable 2-norm of the reduced Newton step (`max_stepsize`)
    pub max_stepsize: f64,
    /// How the KKT system is factored: LU by default, as MATPOWER's MIPS, whose results this
    /// port matches; LDLᵀ only swaps the factorization, with no inertia correction (see
    /// [`crate::kkt`])
    pub kkt_factorization: KktFactorization,
    /// How the pivoted LDLᵀ orders the KKT matrix: AMD by default, or matched (see
    /// [`KktOrdering`]); LU ignores it
    pub kkt_ordering: KktOrdering,
}

impl Default for MipsOptions {
    fn default() -> MipsOptions {
        MipsOptions {
            feastol: 1e-6,
            gradtol: 1e-6,
            comptol: 1e-6,
            costtol: 1e-6,
            max_it: 150,
            step_control: false,
            red_it: 20,
            cost_mult: 1.0,
            xi: 0.99995,
            sigma: 0.1,
            z0: 1.0,
            warm_shift: 1e-4,
            alpha_min: 1e-8,
            rho_min: 0.95,
            rho_max: 1.05,
            mu_threshold: 1e-5,
            max_stepsize: 1e10,
            kkt_factorization: KktFactorization::Lu,
            kkt_ordering: KktOrdering::Amd,
        }
    }
}

/// One iteration's record (`output.hist`): the termination measures, the barrier parameter,
/// the step's 2-norm, the objective and the step lengths.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MipsIteration {
    /// The feasibility measure
    pub feascond: f64,
    /// The gradient measure
    pub gradcond: f64,
    /// The complementarity measure
    pub compcond: f64,
    /// The cost measure
    pub costcond: f64,
    /// The barrier parameter, `γ`
    pub gamma: f64,
    /// The Newton step's 2-norm, `‖dx‖`
    pub stepsize: f64,
    /// The objective, without `cost_mult`
    pub obj: f64,
    /// The primal step length
    pub alphap: f64,
    /// The dual step length
    pub alphad: f64,
}

/// The objective, constraints and Jacobians at a point: MIPS's `f`, `df`, `h`, `g`, and the
/// Jacobians as rows (`dh'` and `dg'` in MIPS), nonlinear ones first.
struct Evaluation {
    f: f64,
    df: Vec<f64>,
    h: Vec<f64>,
    g: Vec<f64>,
    jh: CsrMatrix<f64>,
    jg: CsrMatrix<f64>,
}

/// The fixed parts of the problem: the linear constraints split as MIPS splits them, and the
/// counts.
struct Setup<'a, P: Problem> {
    problem: &'a P,
    nx: usize,
    cost_mult: f64,
    ae: CsrMatrix<f64>,
    be: Vec<f64>,
    ai: CsrMatrix<f64>,
    bi: Vec<f64>,
    neqnln: usize,
    niqnln: usize,
    /// The nonlinear constraints' buffer, kept across evaluations (`Problem::constraints_into`)
    constraints: std::cell::RefCell<NonlinearConstraints>,
}

impl<P: Problem> Setup<'_, P> {
    /// Evaluates the problem at `x` as `mips.m` does at each new point:
    ///
    /// - `f` and `df`: the objective and its gradient, both times `cost_mult`. MIPS scales the
    ///   objective, so the multipliers it finds are scaled too, and divides `cost_mult` back out
    ///   of `f` and the multipliers when it finishes.
    /// - `h`: the nonlinear inequalities' values `h(x)`, then the linear ones' `Ai x − bi`; `g`
    ///   likewise, `g(x)` then `Ae x − be`. The nonlinear constraints come first, so the first
    ///   `niqnln` (or `neqnln`) multipliers are theirs.
    /// - `jh` and `jg`: the Jacobians with rows in the same order, `[∂h/∂x; Ai]` and
    ///   `[∂g/∂x; Ae]` (MIPS's `dh` and `dg`, transposed). Only if `full`: step control's line
    ///   search needs only the values (`mips.m` calls `gh_fcn` there with two outputs), so
    ///   otherwise they're left empty, `0 × nx`, rather than assembled.
    fn evaluate(&self, x: &[f64], full: bool) -> Evaluation {
        let (f, df) = self.problem.objective(x);
        let f = f * self.cost_mult;
        let df: Vec<f64> = df.iter().map(|d| d * self.cost_mult).collect();
        let nonlinear = self.problem.has_nonlinear_constraints().then(|| {
            self.problem
                .constraints_into(x, &mut self.constraints.borrow_mut());
            self.constraints.borrow()
        });
        let mut h = nonlinear.as_ref().map_or(vec![], |c| c.h.clone());
        let mut g = nonlinear.as_ref().map_or(vec![], |c| c.g.clone());
        let ai_x = self.ai.mul_vec(x).expect("Ai has nx columns");
        h.extend(ai_x.iter().zip(&self.bi).map(|(ax, b)| ax - b));
        let ae_x = self.ae.mul_vec(x).expect("Ae has nx columns");
        g.extend(ae_x.iter().zip(&self.be).map(|(ax, b)| ax - b));
        let (jh, jg) = match (&nonlinear, full) {
            (_, false) => {
                let empty = || CsrMatrix::from_triplets(0, self.nx, vec![]).expect("0 × nx");
                (empty(), empty())
            }
            (Some(c), true) => (
                CsrMatrix::vstack([&c.dh, &self.ai]).expect("∂h/∂x has nx columns"),
                CsrMatrix::vstack([&c.dg, &self.ae]).expect("∂g/∂x has nx columns"),
            ),
            (None, true) => (self.ai.clone(), self.ae.clone()),
        };
        Evaluation {
            f,
            df,
            h,
            g,
            jh,
            jg,
        }
    }
}

/// The termination measures: feasibility, gradient, complementarity.
fn conditions(
    e: &Evaluation,
    lx: &[f64],
    x: &[f64],
    z: &[f64],
    lam: &[f64],
    mu: &[f64],
) -> (f64, f64, f64) {
    let maxh = e.h.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let infeasibility = if e.h.is_empty() {
        norm_inf(&e.g)
    } else {
        norm_inf(&e.g).max(maxh)
    };
    let feascond = infeasibility / (1.0 + norm_inf(x).max(norm_inf(z)));
    let gradcond = norm_inf(lx) / (1.0 + norm_inf(lam).max(norm_inf(mu)));
    let compcond = dot(z, mu) / (1.0 + norm_inf(x));
    (feascond, gradcond, compcond)
}

/// `Lx = df + dg λ + dh μ`.
fn lagrangian_gradient(e: &Evaluation, lam: &[f64], mu: &[f64]) -> Vec<f64> {
    let dg_lam = e.jg.vec_mul(lam).expect("one multiplier per equality");
    let dh_mu = e.jh.vec_mul(mu).expect("one multiplier per inequality");
    e.df.iter()
        .zip(&dg_lam)
        .zip(&dh_mu)
        .map(|((df, a), b)| df + a + b)
        .collect()
}

/// MIPS: minimizes `problem`'s objective from `x0` subject to its nonlinear constraints and to
/// `constraints`, as MATPOWER's `mips.m` does.
///
/// - **Constraints:** the bounds and linear constraints are split as MIPS splits them: a row
///   with `|u − l| ≤ ε` is an equality, and a side beyond ±1e10 is unbounded, so each other
///   row becomes one or two inequalities.
/// - **Each iteration** solves the reduced KKT system `[M dg; dgᵀ 0] [dx; dλ] = [−N; −g]`,
///   with `M = Lxx + dh Z⁻¹ diag(μ) dhᵀ` and `N = Lx + dh Z⁻¹ (diag(μ) h + γ e)`, by
///   `ybus_sparse`'s LU ordered by COLAMD (see [`KKT_ORDERING`]; MATLAB uses its backslash),
///   then takes primal and
///   dual steps no longer than `xi` of the way to the boundary, and sets `γ = σ zᵀμ / niq`.
/// - **It converges** when the feasibility, gradient, complementarity and cost measures all
///   fall below their tolerances, and fails numerically on a NaN or too long a step, a step
///   length below `alpha_min`, or `γ` outside `[ε, 1/ε]`.
///
/// Returns `InteriorPointError::InvalidOption` for `xi` outside `[0.5, 1)` or `sigma` outside `(0, 1]`,
/// and `InteriorPointError::Dimension` for inputs of the wrong length.
// MATLAB's `gamma < eps || gamma > 1/eps`, which a NaN γ fails, unlike `!range.contains`
#[allow(clippy::manual_range_contains)]
pub fn mips<P: Problem>(
    problem: &P,
    constraints: &Constraints,
    x0: &[f64],
    options: &MipsOptions,
) -> Result<Solution, InteriorPointError> {
    mips_from(problem, constraints, Start::Cold(x0), options)
}

/// [`mips`] from `start`: a cold start from a point, exactly [`mips`]'s, or a warm start from a
/// previous solution's point and multipliers (see [`crate::start`]).
///
/// A warm start keeps the point and the equality multipliers, recomputes the slacks from the
/// point, and raises every slack and inequality multiplier to at least `warm_shift`; `γ` starts
/// at its update's value there, `σ zᵀμ / niq`. Returns [`mips`]'s errors, `InvalidOption` for a
/// negative or non-finite `warm_shift`, and `Dimension` for a warm start whose point or
/// multipliers don't fit the problem.
#[allow(clippy::manual_range_contains)]
pub fn mips_from<P: Problem>(
    problem: &P,
    constraints: &Constraints,
    start: Start<'_>,
    options: &MipsOptions,
) -> Result<Solution, InteriorPointError> {
    let x0: &[f64] = match start {
        Start::Cold(x0) => x0,
        Start::Warm(warm) => {
            if !(options.warm_shift >= 0.0 && options.warm_shift.is_finite()) {
                return Err(InteriorPointError::InvalidOption { name: "warm_shift" });
            }
            &warm.x
        }
    };
    let opt = options;
    if opt.xi >= 1.0 || opt.xi < 0.5 {
        return Err(InteriorPointError::InvalidOption { name: "xi" });
    }
    if opt.sigma > 1.0 || opt.sigma <= 0.0 {
        return Err(InteriorPointError::InvalidOption { name: "sigma" });
    }
    let nx = x0.len();
    if problem.dimension() != nx {
        return Err(InteriorPointError::Dimension { what: "x0" });
    }

    // The linear constraints: none if A is absent, or if every bound is infinite
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
    let l_free = constraints
        .l
        .as_ref()
        .is_none_or(|l| l.iter().all(|&v| v == f64::NEG_INFINITY));
    let u_free = constraints
        .u
        .as_ref()
        .is_none_or(|u| u.iter().all(|&v| v == f64::INFINITY));
    if na0 == 0 || (l_free && u_free) {
        a = None;
        l.clear();
        u.clear();
    }
    let na = l.len();
    let xmin = constraints
        .xmin
        .clone()
        .unwrap_or_else(|| vec![f64::NEG_INFINITY; nx]);
    let xmax = constraints
        .xmax
        .clone()
        .unwrap_or_else(|| vec![f64::INFINITY; nx]);
    if xmin.len() != nx || xmax.len() != nx {
        return Err(InteriorPointError::Dimension {
            what: "xmin or xmax",
        });
    }

    // AA = [I; A], ll = [xmin; l], uu = [xmax; u], split into equalities and inequalities
    let identity = CsrMatrix::<f64>::identity(nx);
    let aa = match a {
        Some(a) => CsrMatrix::vstack([&identity, a]).expect("A has nx columns"),
        None => identity,
    };
    let ll: Vec<f64> = xmin.iter().chain(&l).copied().collect();
    let uu: Vec<f64> = xmax.iter().chain(&u).copied().collect();
    let all = 0..nx + na;
    let ieq: Vec<usize> = all
        .clone()
        .filter(|&r| (uu[r] - ll[r]).abs() <= f64::EPSILON)
        .collect();
    let igt: Vec<usize> = all
        .clone()
        .filter(|&r| uu[r] >= 1e10 && ll[r] > -1e10)
        .collect();
    let ilt: Vec<usize> = all
        .clone()
        .filter(|&r| ll[r] <= -1e10 && uu[r] < 1e10)
        .collect();
    let ibx: Vec<usize> = all
        .filter(|&r| (uu[r] - ll[r]).abs() > f64::EPSILON && uu[r] < 1e10 && ll[r] > -1e10)
        .collect();
    // Ae = AA(ieq, :) and Ai = [AA(ilt, :); -AA(igt, :); AA(ibx, :); -AA(ibx, :)]
    let columns: Vec<usize> = (0..nx).collect();
    let rows = |rows: &[usize]| aa.select(rows, &columns).expect("rows of AA");
    let ae = rows(&ieq);
    let be: Vec<f64> = ieq.iter().map(|&r| uu[r]).collect();
    let ai = CsrMatrix::vstack([rows(&ilt), -rows(&igt), rows(&ibx), -rows(&ibx)])
        .expect("all nx columns");
    let mut bi: Vec<f64> = ilt.iter().map(|&r| uu[r]).collect();
    bi.extend(igt.iter().map(|&r| -ll[r]));
    bi.extend(ibx.iter().map(|&r| uu[r]));
    bi.extend(ibx.iter().map(|&r| -ll[r]));

    let nonlinear = problem.has_nonlinear_constraints();
    let (neqnln, niqnln) = if nonlinear {
        let c = problem.constraints(x0);
        (c.g.len(), c.h.len())
    } else {
        (0, 0)
    };
    let setup = Setup {
        problem,
        nx,
        cost_mult: opt.cost_mult,
        ae,
        be,
        ai,
        bi,
        neqnln,
        niqnln,
        constraints: std::cell::RefCell::new(NonlinearConstraints::empty(nx)),
    };

    // The starting point
    let mut x = x0.to_vec();
    let mut e = setup.evaluate(&x, true);
    let neq = e.g.len();
    let niq = e.h.len();
    let mut gamma = 1.0;
    let mut lam = vec![0.0; neq];
    let mut z = vec![opt.z0; niq];
    let mut mu = z.clone();
    for (z, &h) in z.iter_mut().zip(&e.h) {
        if h < -opt.z0 {
            *z = -h;
        }
    }
    for (mu, &z) in mu.iter_mut().zip(&z) {
        if gamma / z > opt.z0 {
            *mu = gamma / z;
        }
    }
    if let Start::Warm(warm) = start {
        (lam, z, mu) = warm_mips(
            &start::scaled(&warm.multipliers, opt.cost_mult),
            &e.h,
            (&ieq, &ilt, &igt, &ibx),
            (nx, na),
            (setup.neqnln, setup.niqnln),
            opt.warm_shift,
        )?;
        if niq > 0 {
            gamma = opt.sigma * dot(&z, &mu) / niq as f64;
        }
    }
    let mut f0 = e.f;
    let barrier = |e: &Evaluation, lam: &[f64], mu: &[f64], z: &[f64], gamma: f64| -> f64 {
        let hz: Vec<f64> = e.h.iter().zip(z).map(|(h, z)| h + z).collect();
        e.f + dot(lam, &e.g) + dot(mu, &hz) - gamma * z.iter().map(|z| z.ln()).sum::<f64>()
    };
    let mut big_l = if opt.step_control {
        barrier(&e, &lam, &mu, &z, gamma)
    } else {
        0.0
    };
    let mut lx = lagrangian_gradient(&e, &lam, &mu);
    let (mut feascond, mut gradcond, mut compcond) = conditions(&e, &lx, &x, &z, &lam, &mu);
    let mut costcond = (e.f - f0).abs() / (1.0 + f0.abs());
    let mut history = vec![MipsIteration {
        feascond,
        gradcond,
        compcond,
        costcond,
        gamma,
        stepsize: 0.0,
        obj: e.f / opt.cost_mult,
        alphap: 0.0,
        alphad: 0.0,
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

    while !converged && i < opt.max_it {
        i += 1;

        // The step: the reduced KKT system
        // Into the buffer kept across iterations; scaling in place is `scale`'s product
        if nonlinear {
            problem.lagrangian_hessian_into(
                &x,
                NonlinearMultipliers {
                    eq: &lam[..setup.neqnln],
                    ineq: &mu[..setup.niqnln],
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
        let lxx = &hessian_buffer;
        let zinv: Vec<f64> = z.iter().map(|z| 1.0 / z).collect();
        // dh Z⁻¹ (rows of dhᵀ scaled), then by μ: M = Lxx + dh Z⁻¹ diag(μ) dhᵀ
        let jh_zinv = e.jh.scale_rows(&zinv).expect("one slack per inequality");
        //   (as jh_zinv.scale_rows(&mu), μᵢ times each entry), through the kept plans
        let jh_zinv_mu: Vec<f64> = (jh_zinv.values().iter())
            .zip(gram.jh_rows(lxx, &e.jh))
            .map(|(&v, &row)| mu[row] * v)
            .collect();
        let m = gram.compute(lxx, &e.jh, &jh_zinv_mu);
        let weights: Vec<f64> = mu.iter().zip(&e.h).map(|(mu, h)| mu * h + gamma).collect();
        let n_vec: Vec<f64> = lx
            .iter()
            .zip(
                &jh_zinv
                    .vec_mul(&weights)
                    .expect("one weight per inequality"),
            )
            .map(|(lx, s)| lx + s)
            .collect();
        let rhs: Vec<f64> = n_vec.iter().chain(&e.g).map(|v| -v).collect();
        let solved = kkt
            .factor(nx, &m, &[], &e.jg, 0.0)
            .and_then(|factor| factor.solve(&rhs));
        let dxdlam = match solved {
            Ok(step) => step,
            // MATLAB's backslash returns NaN or Inf for a singular system, which MIPS treats
            //   as a numerical failure
            Err(_) => {
                failed = true;
                break;
            }
        };
        if dxdlam.iter().any(|v| v.is_nan()) || dot(&dxdlam, &dxdlam).sqrt() > opt.max_stepsize {
            failed = true;
            break;
        }
        let mut dx: Vec<f64> = dxdlam[..nx].to_vec();
        let mut dlam: Vec<f64> = dxdlam[nx..].to_vec();
        let jh_dx = e.jh.mul_vec(&dx).expect("one entry per variable");
        let mut dz: Vec<f64> = (0..niq).map(|k| -e.h[k] - z[k] - jh_dx[k]).collect();
        let mut dmu: Vec<f64> = (0..niq)
            .map(|k| -mu[k] + zinv[k] * (gamma - mu[k] * dz[k]))
            .collect();

        // Step control: a full step that worsens both feasibility and the gradient is cut back
        if opt.step_control {
            let x1: Vec<f64> = x.iter().zip(&dx).map(|(x, d)| x + d).collect();
            let e1 = setup.evaluate(&x1, true);
            let lx1 = lagrangian_gradient(&e1, &lam, &mu);
            let (feascond1, gradcond1, _) = conditions(&e1, &lx1, &x1, &z, &lam, &mu);
            if feascond1 > feascond && gradcond1 > gradcond {
                let mut alpha = 1.0;
                let lxx_dx = |d: &[f64]| lxx.mul_vec(d).expect("nx × nx");
                for _ in 0..opt.red_it {
                    let dx1: Vec<f64> = dx.iter().map(|d| alpha * d).collect();
                    let x1: Vec<f64> = x.iter().zip(&dx1).map(|(x, d)| x + d).collect();
                    let e1 = setup.evaluate(&x1, false);
                    let l1 = barrier(&e1, &lam, &mu, &z, gamma);
                    let rho = (l1 - big_l) / (dot(&lx, &dx1) + 0.5 * dot(&dx1, &lxx_dx(&dx1)));
                    if rho > opt.rho_min && rho < opt.rho_max {
                        break;
                    }
                    alpha /= 2.0;
                }
                for v in dx
                    .iter_mut()
                    .chain(&mut dz)
                    .chain(&mut dlam)
                    .chain(&mut dmu)
                {
                    *v *= alpha;
                }
            }
        }

        // Step lengths, at most xi of the way to the boundary
        let step_length = |values: &[f64], steps: &[f64]| -> f64 {
            let limit = values
                .iter()
                .zip(steps)
                .filter(|&(_, &d)| d < 0.0)
                .map(|(v, d)| v / -d)
                .fold(f64::INFINITY, f64::min);
            if limit.is_finite() {
                (opt.xi * limit).min(1.0)
            } else {
                1.0
            }
        };
        let alphap = step_length(&z, &dz);
        let alphad = step_length(&mu, &dmu);
        for k in 0..nx {
            x[k] += alphap * dx[k];
        }
        for k in 0..niq {
            z[k] += alphap * dz[k];
            mu[k] += alphad * dmu[k];
        }
        for k in 0..neq {
            lam[k] += alphad * dlam[k];
        }
        if niq > 0 {
            gamma = opt.sigma * dot(&z, &mu) / niq as f64;
        }

        // The new point
        e = setup.evaluate(&x, true);
        lx = lagrangian_gradient(&e, &lam, &mu);
        (feascond, gradcond, compcond) = conditions(&e, &lx, &x, &z, &lam, &mu);
        costcond = (e.f - f0).abs() / (1.0 + f0.abs());
        history.push(MipsIteration {
            feascond,
            gradcond,
            compcond,
            costcond,
            gamma,
            stepsize: dot(&dx, &dx).sqrt(),
            obj: e.f / opt.cost_mult,
            alphap,
            alphad,
        });
        if tolerances_met(feascond, gradcond, compcond, costcond) {
            converged = true;
        } else {
            if x.iter().any(|v| v.is_nan())
                || alphap < opt.alpha_min
                || alphad < opt.alpha_min
                || gamma < f64::EPSILON
                || gamma > 1.0 / f64::EPSILON
            {
                failed = true;
                break;
            }
            f0 = e.f;
            if opt.step_control {
                big_l = barrier(&e, &lam, &mu, &z, gamma);
            }
        }
    }
    let status = if failed {
        Status::NumericalFailure
    } else if converged {
        Status::Converged
    } else {
        Status::MaxIterations
    };

    // The final iterate's inequality multipliers, kept for a warm start before the report
    //   thresholds them (unscaled, as reported)
    let iterate_mu: Vec<f64> = mu.iter().map(|v| v / opt.cost_mult).collect();
    // Inequalities well inside their limits lose a small multiplier; then the scale comes off
    for (mu, &h) in mu.iter_mut().zip(&e.h) {
        if h < -opt.feastol && *mu < opt.mu_threshold {
            *mu = 0.0;
        }
    }
    let f = e.f / opt.cost_mult;
    for v in lam.iter_mut().chain(mu.iter_mut()) {
        *v /= opt.cost_mult;
    }

    // The linear constraints' multipliers, back on the rows of [I; A]
    let sets = (&ieq[..], &ilt[..], &igt[..], &ibx[..]);
    let nonlinear = (setup.neqnln, setup.niqnln);
    let multipliers = report(&lam, &mu, sets, (nx, na), nonlinear);
    let iterate_multipliers = report(&lam, &iterate_mu, sets, (nx, na), nonlinear);
    Ok(Solution {
        x,
        f,
        status,
        iterations: i,
        multipliers,
        iterate_multipliers,
        history: History::Mips(history),
        analyses: kkt.analyses(),
        refinement_steps: kkt.refinement_steps(),
        kkt: kkt.statistics(),
        inertia_corrections: 0,
    })
}

/// MIPS for a quadratic or linear program, as MATPOWER's `qps_mips`: minimizes
/// `½ xᵀ H x + cᵀ x` subject to `constraints`, from `x0` (zeros if absent).
///
/// Returns `InteriorPointError::UnboundedLinearProgram` for a linear program with no constraints or
/// bounds, and the errors of [`mips`].
pub fn qps_mips(
    qp: &QuadraticProblem,
    constraints: &Constraints,
    x0: Option<&[f64]>,
    options: &MipsOptions,
) -> Result<Solution, InteriorPointError> {
    let zeros = vec![0.0; qp.c.len()];
    qps_mips_from(qp, constraints, Start::Cold(x0.unwrap_or(&zeros)), options)
}

/// [`qps_mips`] from `start` (see [`mips_from`]), with [`qps_mips`]'s errors.
pub fn qps_mips_from(
    qp: &QuadraticProblem,
    constraints: &Constraints,
    start: Start<'_>,
    options: &MipsOptions,
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
    mips_from(qp, constraints, start, options)
}

/// The multipliers in MIPS's reported layout, from its internal `λ` and `μ`: each linear row's
/// back on the rows of `AA = [I; A]` (a row `r < nx` is variable `r`'s bounds), an equality's
/// sign deciding its side. `warm_mips` is its inverse.
fn report(
    lam: &[f64],
    mu: &[f64],
    (ieq, ilt, igt, ibx): (&[usize], &[usize], &[usize], &[usize]),
    (nx, na): (usize, usize),
    (neqnln, niqnln): (usize, usize),
) -> Multipliers {
    let lam_lin = &lam[neqnln..];
    let (nlt, ngt, nbx) = (ilt.len(), igt.len(), ibx.len());
    let mu_lin = &mu[niqnln..];
    let mut mu_l = vec![0.0; nx + na];
    let mut mu_u = vec![0.0; nx + na];
    for (k, &r) in ieq.iter().enumerate() {
        if lam_lin[k] < 0.0 {
            mu_l[r] = -lam_lin[k];
        } else if lam_lin[k] > 0.0 {
            mu_u[r] = lam_lin[k];
        }
    }
    for (k, &r) in igt.iter().enumerate() {
        mu_l[r] = mu_lin[nlt + k];
    }
    for (k, &r) in ibx.iter().enumerate() {
        mu_l[r] = mu_lin[nlt + ngt + nbx + k];
        mu_u[r] = mu_lin[nlt + ngt + k];
    }
    for (k, &r) in ilt.iter().enumerate() {
        mu_u[r] = mu_lin[k];
    }
    Multipliers {
        eq_nonlin: lam[..neqnln].to_vec(),
        ineq_nonlin: mu[..niqnln].to_vec(),
        mu_l: mu_l[nx..].to_vec(),
        mu_u: mu_u[nx..].to_vec(),
        lower: mu_l[..nx].to_vec(),
        upper: mu_u[..nx].to_vec(),
    }
}

/// A warm start's `(λ, z, μ)`
type WarmIterate = (Vec<f64>, Vec<f64>, Vec<f64>);

/// A warm start's `(λ, z, μ)` in MIPS's internal layout, from `multipliers` in its reported
/// layout and the inequalities' values `h` at the point: the inverse of the report `mips_from`
/// builds at its end, with the slacks `−h` and every slack and inequality multiplier raised to
/// at least `shift`. `AA = [I; A]`'s rows `ieq`, `ilt`, `igt` and `ibx` are the split's; a row
/// `r < nx` is variable `r`'s bounds (`lower`, `upper`), and row `nx + k` is `A`'s row `k`
/// (`mu_l`, `mu_u`).
fn warm_mips(
    multipliers: &Multipliers,
    h: &[f64],
    (ieq, ilt, igt, ibx): (&[usize], &[usize], &[usize], &[usize]),
    (nx, na): (usize, usize),
    (neqnln, niqnln): (usize, usize),
    shift: f64,
) -> Result<WarmIterate, InteriorPointError> {
    let m = multipliers;
    if m.eq_nonlin.len() != neqnln
        || m.ineq_nonlin.len() != niqnln
        || m.mu_l.len() != na
        || m.mu_u.len() != na
        || m.lower.len() != nx
        || m.upper.len() != nx
    {
        return Err(InteriorPointError::Dimension {
            what: "the warm start's multipliers",
        });
    }
    // Row r of AA's multipliers for its lower and upper sides
    let lower = |r: usize| if r < nx { m.lower[r] } else { m.mu_l[r - nx] };
    let upper = |r: usize| if r < nx { m.upper[r] } else { m.mu_u[r - nx] };
    // λ: the nonlinear equalities', then each equality row's, its upper side's less its lower's
    //   (the report gives one side the multiplier's magnitude and the other 0)
    let mut lam = m.eq_nonlin.clone();
    lam.extend(ieq.iter().map(|&r| upper(r) - lower(r)));
    // μ: the nonlinear inequalities', then the rows as Ai stacks them: ilt, −igt, ibx, −ibx
    let mut mu = m.ineq_nonlin.clone();
    mu.extend(ilt.iter().map(|&r| upper(r)));
    mu.extend(igt.iter().map(|&r| lower(r)));
    mu.extend(ibx.iter().map(|&r| upper(r)));
    mu.extend(ibx.iter().map(|&r| lower(r)));
    debug_assert_eq!(mu.len(), h.len());
    for value in &mut mu {
        *value = value.max(shift);
    }
    let z = h.iter().map(|&h| (-h).max(shift)).collect();
    Ok((lam, z, mu))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Multipliers with a distinct value in every slot a split reports, and 0 elsewhere: an
    ///   equality row's on one side only, as the report gives it
    fn distinct(
        (ieq, ilt, igt, ibx): (&[usize], &[usize], &[usize], &[usize]),
        (nx, na): (usize, usize),
        (neqnln, niqnln): (usize, usize),
    ) -> Multipliers {
        let mut next = 0.5;
        let mut value = || {
            next += 1.25;
            next
        };
        let (mut lower, mut upper) = (vec![0.0; nx + na], vec![0.0; nx + na]);
        for (k, &r) in ieq.iter().enumerate() {
            if k % 2 == 0 {
                lower[r] = value();
            } else {
                upper[r] = value();
            }
        }
        for &r in ilt {
            upper[r] = value();
        }
        for &r in igt {
            lower[r] = value();
        }
        for &r in ibx {
            (lower[r], upper[r]) = (value(), value());
        }
        Multipliers {
            eq_nonlin: (0..neqnln).map(|_| value() - 7.0).collect(),
            ineq_nonlin: (0..niqnln).map(|_| value()).collect(),
            mu_l: lower[nx..].to_vec(),
            mu_u: upper[nx..].to_vec(),
            lower: lower[..nx].to_vec(),
            upper: upper[..nx].to_vec(),
        }
    }

    #[test]
    fn the_warm_mapping_is_the_reports_inverse() {
        // AA = [I; A] with 4 variables and 5 rows of A, every kind of row on both, and one row
        //   unbounded both ways (8)
        let sets: (&[usize], &[usize], &[usize], &[usize]) = (&[3, 4], &[0, 5], &[1, 6], &[2, 7]);
        let (dims, nonlinear) = ((4, 5), (2, 3));
        let m = distinct(sets, dims, nonlinear);
        let niq = nonlinear.1 + sets.1.len() + sets.2.len() + 2 * sets.3.len();
        let h = vec![-1.0; niq];
        let (lam, z, mu) = warm_mips(&m, &h, sets, dims, nonlinear, 0.0).unwrap();
        assert_eq!(z, vec![1.0; niq]);
        assert_eq!(report(&lam, &mu, sets, dims, nonlinear), m);
    }

    #[test]
    fn the_warm_start_moves_every_pair_off_the_boundary() {
        // Active rows (h = 0) with zero multipliers, and an inactive one
        let sets: (&[usize], &[usize], &[usize], &[usize]) = (&[], &[0], &[1], &[]);
        let m = Multipliers {
            mu_l: vec![],
            mu_u: vec![],
            lower: vec![0.0, 0.0],
            upper: vec![0.0, 0.0],
            ..Multipliers::default()
        };
        let (_, z, mu) = warm_mips(&m, &[0.0, -3.0], sets, (2, 0), (0, 0), 1e-3).unwrap();
        assert_eq!(z, vec![1e-3, 3.0]);
        assert_eq!(mu, vec![1e-3, 1e-3]);
    }
}
