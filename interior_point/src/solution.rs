//! The solvers' results: a [`Solution`] from either solver, how its solve ended, its multipliers
//! in MIPS's convention, and a check of the KKT conditions at a point that doesn't depend on
//! which solver found it.

use crate::interior_point_error::InteriorPointError;
use crate::mehrotra::MehrotraIteration;
use crate::mips::MipsIteration;
use crate::problem::{Constraints, Problem};
use crate::vector::norm_inf;

/// How a solve ended, for either solver.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// Every termination test passed (MIPS's `exitflag` 1)
    Converged,
    /// `max_it` iterations passed without converging (0)
    MaxIterations,
    /// The solve failed numerically (−1); each solver's docs list the ways it can
    NumericalFailure,
}

/// The multipliers at a solution, as MIPS reports them. Each is nonnegative at a KKT point except
/// the equality constraints', and each is 0 for a constraint that's absent or unbounded that way.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Multipliers {
    /// One per nonlinear equality constraint (`eqnonlin`)
    pub eq_nonlin: Vec<f64>,
    /// One per nonlinear inequality constraint (`ineqnonlin`)
    pub ineq_nonlin: Vec<f64>,
    /// One per linear constraint, for its lower bound `l` (`mu_l`)
    pub mu_l: Vec<f64>,
    /// One per linear constraint, for its upper bound `u` (`mu_u`)
    pub mu_u: Vec<f64>,
    /// One per variable, for `xmin` (`lower`)
    pub lower: Vec<f64>,
    /// One per variable, for `xmax` (`upper`)
    pub upper: Vec<f64>,
}

/// A solver's result, from MIPS ([`mips`](crate::mips::mips)) or the predictor-corrector method
/// ([`mehrotra`](crate::mehrotra::mehrotra)).
#[derive(Clone, Debug, PartialEq)]
pub struct Solution {
    /// The final point
    pub x: Vec<f64>,
    /// The objective there
    pub f: f64,
    /// How the solve ended
    pub status: Status,
    /// The number of iterations taken
    pub iterations: usize,
    /// The multipliers
    pub multipliers: Multipliers,
    /// One record for the starting point and one per iteration, in the solver's own terms
    pub history: History,
    /// How many times the KKT matrix was analyzed: 1 if its pattern held from the first
    /// iteration to the last (see [`KKT_ORDERING`](crate::mips::KKT_ORDERING)), 0 if no
    /// iteration factored it
    pub analyses: usize,
    /// The refinement steps taken over every KKT solve: 0 with LU, which isn't refined
    pub refinement_steps: usize,
    /// How many times the KKT matrix was refactored with a larger `δ_w` to correct its inertia:
    /// only the predictor-corrector method corrects it, with LDLᵀ, so 0 otherwise
    pub inertia_corrections: usize,
}

impl Solution {
    /// Whether the solve converged.
    pub fn converged(&self) -> bool {
        self.status == Status::Converged
    }
}

/// Each iteration's record, from the solver that took them: one for the starting point, then
/// one per iteration.
#[derive(Clone, Debug, PartialEq)]
pub enum History {
    /// MIPS's records (`output.hist`)
    Mips(Vec<MipsIteration>),
    /// The predictor-corrector method's records
    Mehrotra(Vec<MehrotraIteration>),
}

impl History {
    /// The number of records: one more than the iterations taken.
    pub fn len(&self) -> usize {
        match self {
            History::Mips(records) => records.len(),
            History::Mehrotra(records) => records.len(),
        }
    }

    /// Whether there are no records, which a solve never returns.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// How far a point and its multipliers are from satisfying the KKT conditions, each an ∞-norm,
/// unscaled, from [`kkt_residuals`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KktResiduals {
    /// The Lagrangian's gradient: `∇f + ∂gᵀ λ + ∂hᵀ μ + Aᵀ (μ_u − μ_l) + (upper − lower)`
    pub stationarity: f64,
    /// The largest constraint violation: `|g(x)|`, `h(x)` above 0, `A x` outside `[l, u]`, and
    /// `x` outside `[xmin, xmax]`
    pub feasibility: f64,
    /// The largest product of a multiplier and its constraint's distance from its limit, over
    /// the limits that are finite
    pub complementarity: f64,
    /// The largest amount by which an inequality's multiplier is negative
    pub dual_feasibility: f64,
}

impl KktResiduals {
    /// The largest of the four residuals.
    pub fn max(&self) -> f64 {
        self.stationarity
            .max(self.feasibility)
            .max(self.complementarity)
            .max(self.dual_feasibility)
    }
}

/// The running maxima of the residuals other than stationarity.
#[derive(Default)]
struct Residuals {
    feasibility: f64,
    complementarity: f64,
    dual_feasibility: f64,
}

impl Residuals {
    /// An inequality whose `slack` is nonnegative when it's satisfied, with its multiplier
    fn inequality(&mut self, multiplier: f64, slack: f64) {
        self.feasibility = self.feasibility.max(-slack);
        self.complementarity = self.complementarity.max((multiplier * slack).abs());
        self.dual_feasibility = self.dual_feasibility.max(-multiplier);
    }

    /// The multiplier of a side with no limit, which must be zero, so only its sign counts here
    /// (any other value shows in the stationarity)
    fn unlimited(&mut self, multiplier: f64) {
        self.dual_feasibility = self.dual_feasibility.max(-multiplier);
    }
}

/// A limit counts as finite, as MIPS counts it, if it's within ±1e10.
fn finite(limit: f64) -> bool {
    limit.abs() < 1e10
}

/// The KKT residuals of `problem` under `constraints` at `x` with `multipliers`, in MIPS's
/// convention: the objective is unscaled (no `cost_mult`), the nonlinear constraints are
/// `g(x) = 0` and `h(x) ≤ 0`, and a linear constraint or a bound has one multiplier for each
/// side. Empty `mu_l` and `mu_u` count as zeros, as MIPS reports them when it drops `A`
/// (absent, or every limit infinite); so do empty nonlinear multipliers for a problem without
/// nonlinear constraints.
///
/// Returns `InteriorPointError::Dimension` if `x`, `A`, the limits or the multipliers don't
/// match the problem's sizes.
pub fn kkt_residuals<P: Problem>(
    problem: &P,
    constraints: &Constraints,
    x: &[f64],
    multipliers: &Multipliers,
) -> Result<KktResiduals, InteriorPointError> {
    let nx = problem.dimension();
    if x.len() != nx {
        return Err(InteriorPointError::Dimension { what: "x" });
    }
    let m = multipliers;
    if m.lower.len() != nx || m.upper.len() != nx {
        return Err(InteriorPointError::Dimension {
            what: "the bounds' multipliers",
        });
    }
    let (_, mut gradient) = problem.objective(x);
    let mut r = Residuals::default();

    // The nonlinear constraints
    if problem.has_nonlinear_constraints() {
        let c = problem.constraints(x);
        let eq = if m.eq_nonlin.is_empty() {
            vec![0.0; c.g.len()]
        } else {
            m.eq_nonlin.clone()
        };
        let ineq = if m.ineq_nonlin.is_empty() {
            vec![0.0; c.h.len()]
        } else {
            m.ineq_nonlin.clone()
        };
        if eq.len() != c.g.len() || ineq.len() != c.h.len() {
            return Err(InteriorPointError::Dimension {
                what: "the nonlinear constraints' multipliers",
            });
        }
        for (gr, d) in gradient
            .iter_mut()
            .zip(c.dg.vec_mul(&eq).expect("one multiplier per equality"))
        {
            *gr += d;
        }
        for (gr, d) in gradient
            .iter_mut()
            .zip(c.dh.vec_mul(&ineq).expect("one multiplier per inequality"))
        {
            *gr += d;
        }
        r.feasibility = r.feasibility.max(norm_inf(&c.g));
        for (&mu, &h) in ineq.iter().zip(&c.h) {
            r.inequality(mu, -h);
        }
    }

    // The linear constraints, l ≤ A x ≤ u
    if let Some(a) = &constraints.a {
        let (na, columns) = a.shape();
        if columns != nx {
            return Err(InteriorPointError::Dimension { what: "A" });
        }
        let limit = |side: &Option<Vec<f64>>, default: f64| -> Result<Vec<f64>, _> {
            let values = side.clone().unwrap_or_else(|| vec![default; na]);
            if values.len() == na {
                Ok(values)
            } else {
                Err(InteriorPointError::Dimension { what: "l or u" })
            }
        };
        let (l, u) = (
            limit(&constraints.l, f64::NEG_INFINITY)?,
            limit(&constraints.u, f64::INFINITY)?,
        );
        let side = |values: &[f64]| -> Result<Vec<f64>, InteriorPointError> {
            if values.is_empty() {
                Ok(vec![0.0; na])
            } else if values.len() == na {
                Ok(values.to_vec())
            } else {
                Err(InteriorPointError::Dimension {
                    what: "the linear constraints' multipliers",
                })
            }
        };
        let (mu_l, mu_u) = (side(&m.mu_l)?, side(&m.mu_u)?);
        let net: Vec<f64> = mu_u.iter().zip(&mu_l).map(|(u, l)| u - l).collect();
        for (gr, d) in gradient
            .iter_mut()
            .zip(a.vec_mul(&net).expect("one multiplier per row"))
        {
            *gr += d;
        }
        let ax = a.mul_vec(x).expect("A has nx columns");
        for k in 0..na {
            if finite(l[k]) {
                r.inequality(mu_l[k], ax[k] - l[k]);
            } else {
                r.unlimited(mu_l[k]);
            }
            if finite(u[k]) {
                r.inequality(mu_u[k], u[k] - ax[k]);
            } else {
                r.unlimited(mu_u[k]);
            }
        }
    } else if m.mu_l.iter().chain(&m.mu_u).any(|&v| v != 0.0) {
        return Err(InteriorPointError::Dimension {
            what: "the linear constraints' multipliers",
        });
    }

    // The bounds, xmin ≤ x ≤ xmax
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
    let (xmin, xmax) = (
        bound(&constraints.xmin, f64::NEG_INFINITY)?,
        bound(&constraints.xmax, f64::INFINITY)?,
    );
    for k in 0..nx {
        gradient[k] += m.upper[k] - m.lower[k];
        if finite(xmin[k]) {
            r.inequality(m.lower[k], x[k] - xmin[k]);
        } else {
            r.unlimited(m.lower[k]);
        }
        if finite(xmax[k]) {
            r.inequality(m.upper[k], xmax[k] - x[k]);
        } else {
            r.unlimited(m.upper[k]);
        }
    }

    // Adding +0.0 turns a −0.0 (the maximum of 0 and the negation of a zero multiplier) into
    //   0.0, so a satisfied condition reads as zero
    Ok(KktResiduals {
        stationarity: norm_inf(&gradient),
        feasibility: r.feasibility + 0.0,
        complementarity: r.complementarity + 0.0,
        dual_feasibility: r.dual_feasibility + 0.0,
    })
}
