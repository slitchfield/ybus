//! Generator costs: their values ([`GenCost::evaluate`], [`total_cost`]), as MATPOWER's
//! `totcost` computes them, and their form in the optimal power flow ([`OpfCosts`]), as
//! MATPOWER 8.1's `opf_setup` and `makeAy` build it.
//!
//! In the optimal power flow, a cost is either
//! - **polynomial**, at most quadratic, which becomes the quadratic term
//!   `½ Q p² + c p + k` of the objective, with `p` in per unit; or
//! - **piecewise linear** (convex), which becomes a helper variable `y` in the objective and one
//!   linear constraint `y ≥ mⱼ p − bⱼ` per segment, so that at the optimum `y` is the largest
//!   of the segments' lines, the cost.
//!
//! Before that, as `opf_setup` does, a piecewise linear cost of one segment is replaced by the
//! linear polynomial through its two points, and a polynomial's zero leading coefficients are
//! dropped while it has more than two, so that `[0, 0, c, k]` is linear and `[0, k]` stays
//! linear with a zero slope. Startup and shutdown costs play no part in either.

use crate::case::{Case, GenCost};
use crate::powerflow_error::PowerflowError;

impl GenCost {
    /// The cost ($/hr) at output `p` (MW, or MVAr for a reactive power cost), as MATPOWER's
    /// `totcost` computes it.
    ///
    /// A piecewise linear cost is extended past its first and last points along its first and
    /// last segments; one with a single point costs that point's cost everywhere, and one with
    /// none costs nothing. Startup and shutdown costs aren't included.
    pub fn evaluate(&self, p: f64) -> f64 {
        match self {
            // As polycost: the constant term, plus each coefficient times its power of p
            GenCost::Polynomial { coeffs, .. } => {
                let mut powers = coeffs.iter().rev();
                let constant = powers.next().copied().unwrap_or(0.0);
                powers
                    .zip(1..)
                    .fold(constant, |sum, (c, k)| sum + c * p.powi(k))
            }
            // As mkpp and ppval: the segment starting at the last point at or below p, or the
            //   first segment below the first point
            GenCost::PiecewiseLinear { points, .. } => match points.as_slice() {
                [] => 0.0,
                [(_, cost)] => *cost,
                _ => {
                    let segments = &points[..points.len() - 1];
                    let j = segments.iter().rposition(|&(x, _)| x <= p).unwrap_or(0);
                    let ((x0, c0), (x1, c1)) = (points[j], points[j + 1]);
                    (c1 - c0) / (x1 - x0) * (p - x0) + c0
                }
            },
        }
    }
}

/// The total cost ($/hr) of generators' outputs, as MATPOWER's `totcost` summed over them:
/// `costs[i]` at `pg_mw[i]`, plus, when `costs` holds reactive costs (twice as many rows as
/// generators) and `qg_mvar` is given, `costs[ng + i]` at `qg_mvar[i]`.
///
/// Every generator given is counted, in service or not; pass only those wanted.
///
/// # Errors
///
/// [`PowerflowError::CountMismatch`] if `costs` doesn't have one or two rows per entry of
/// `pg_mw`, or `qg_mvar` isn't as long as `pg_mw`.
pub fn total_cost(
    costs: &[GenCost],
    pg_mw: &[f64],
    qg_mvar: Option<&[f64]>,
) -> Result<f64, PowerflowError> {
    let ng = pg_mw.len();
    if costs.len() != ng && costs.len() != 2 * ng {
        return Err(PowerflowError::CountMismatch {
            field: "gencost",
            expected: format!("{ng} or {}", 2 * ng),
            found: costs.len(),
        });
    }
    let active: f64 = costs.iter().zip(pg_mw).map(|(c, &p)| c.evaluate(p)).sum();
    let reactive = match qg_mvar {
        Some(qg) if qg.len() != ng => {
            return Err(PowerflowError::CountMismatch {
                field: "qg_mvar",
                expected: ng.to_string(),
                found: qg.len(),
            })
        }
        Some(qg) if costs.len() == 2 * ng => costs[ng..]
            .iter()
            .zip(qg)
            .map(|(c, &q)| c.evaluate(q))
            .sum(),
        _ => 0.0,
    };
    Ok(active + reactive)
}

/// Polynomial costs in the optimal power flow's form, one entry per generator: the cost of
/// generator `i` is `½ quadratic[i] p² + linear[i] p + constant[i]` with `p` in per unit, as
/// `opf_setup`'s `Qpg`, `cpg` and `kpg` (or `Qqg`, `cqg` and `kqg`). A generator whose cost
/// is piecewise linear has zeros here.
#[derive(Clone, Debug, PartialEq)]
pub struct QuadraticCosts {
    /// `Q`: twice the quadratic coefficient, times `baseMVA²`
    pub quadratic: Vec<f64>,
    /// `c`: the linear coefficient times `baseMVA`
    pub linear: Vec<f64>,
    /// `k`: the constant term
    pub constant: Vec<f64>,
}

impl QuadraticCosts {
    fn zeros(n: usize) -> Self {
        QuadraticCosts {
            quadratic: vec![0.0; n],
            linear: vec![0.0; n],
            constant: vec![0.0; n],
        }
    }

    /// Whether every quadratic coefficient is zero, so that the costs are linear: what makes
    /// MATPOWER's default DC solver treat the DC optimal power flow as a linear program.
    pub fn is_linear(&self) -> bool {
        self.quadratic.iter().all(|&q| q == 0.0)
    }

    /// `Σ ½ Qᵢ pᵢ² + cᵢ pᵢ + kᵢ` at outputs `p` (p.u.)
    pub fn evaluate(&self, p: &[f64]) -> f64 {
        p.iter()
            .enumerate()
            .map(|(i, &p)| 0.5 * self.quadratic[i] * p * p + self.linear[i] * p + self.constant[i])
            .sum()
    }
}

/// The output a piecewise linear cost is a function of: a generator's real or reactive power,
/// by its position among the optimal power flow's generators.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CostVariable {
    /// The real power of the generator at this position
    Active(usize),
    /// The reactive power of the generator at this position
    Reactive(usize),
}

/// A piecewise linear cost in the optimal power flow's form, as `makeAy` builds it: its helper
/// variable `y` is constrained by `slopes[j] · p − y ≤ bounds[j]` for each segment `j`, with
/// `p` the output (p.u.) of [`PiecewiseCost::variable`].
#[derive(Clone, Debug, PartialEq)]
pub struct PiecewiseCost {
    /// The output the cost is a function of
    pub variable: CostVariable,
    /// Each segment's slope, `mⱼ = (cⱼ₊₁ − cⱼ) / (pⱼ₊₁ − pⱼ)` with `p` in per unit
    pub slopes: Vec<f64>,
    /// Each segment's bound, `bⱼ = mⱼ pⱼ − cⱼ`
    pub bounds: Vec<f64>,
    /// The cost at the last point, from which MATPOWER's OPF solvers start `y`
    pub last_cost: f64,
}

/// A case's generator costs as the optimal power flow uses them, for a set of generators (the
/// in-service ones, in MATPOWER's internal order), as `opf_setup` and `makeAy` build them.
///
/// The objective is [`OpfCosts::evaluate`]: `active`'s costs of the generators' real power,
/// `reactive`'s of their reactive power if any, and one helper variable `y` per entry of
/// `piecewise`, in that order: the real power costs' first, then the reactive ones', each in
/// generator order.
#[derive(Clone, Debug, PartialEq)]
pub struct OpfCosts {
    /// The polynomial costs of the generators' real power
    pub active: QuadraticCosts,
    /// The polynomial costs of their reactive power, when there are reactive costs and they
    /// were asked for
    pub reactive: Option<QuadraticCosts>,
    /// The piecewise linear costs, one helper variable each
    pub piecewise: Vec<PiecewiseCost>,
}

impl OpfCosts {
    /// The costs of the generators at case rows `gen_rows`, in that order, with their reactive
    /// power costs when `reactive` is set and the case has them; a DC optimal power flow
    /// ignores them, as `opf_setup` does.
    ///
    /// # Errors
    ///
    /// - [`PowerflowError::MissingField`] if the case has no costs
    /// - [`PowerflowError::RowOutOfRange`] if a row of `gen_rows` isn't a generator
    /// - [`PowerflowError::UnsupportedCost`] for a polynomial above quadratic (after zero
    ///   leading coefficients are dropped), or a piecewise linear cost with fewer than two points
    ///   or two consecutive points at the same output
    pub fn new(case: &Case, gen_rows: &[usize], reactive: bool) -> Result<Self, PowerflowError> {
        let costs = case
            .gencost
            .as_ref()
            .ok_or(PowerflowError::MissingField { field: "gencost" })?;
        let ng = case.generator.len();
        if let Some(&row) = gen_rows.iter().find(|&&row| row >= ng) {
            return Err(PowerflowError::RowOutOfRange {
                table: "gen",
                row,
                len: ng,
            });
        }
        let base_mva = case.base_mva;
        let mut piecewise = vec![];
        let active = prepare(
            costs,
            gen_rows,
            0,
            base_mva,
            CostVariable::Active,
            &mut piecewise,
        )?;
        let reactive = if reactive && costs.len() == 2 * ng {
            Some(prepare(
                costs,
                gen_rows,
                ng,
                base_mva,
                CostVariable::Reactive,
                &mut piecewise,
            )?)
        } else {
            None
        };
        Ok(OpfCosts {
            active,
            reactive,
            piecewise,
        })
    }

    /// The objective at real outputs `pg` and reactive outputs `qg` (p.u., one per generator)
    /// and helper variables `y` (one per piecewise linear cost): the polynomial costs plus the
    /// sum of `y`. `qg` is used only if there are reactive costs. At a solution, where each `y`
    /// is its cost, this is the total cost.
    pub fn evaluate(&self, pg: &[f64], qg: &[f64], y: &[f64]) -> f64 {
        let reactive = self
            .reactive
            .as_ref()
            .map_or(0.0, |costs| costs.evaluate(qg));
        self.active.evaluate(pg) + reactive + y.iter().sum::<f64>()
    }
}

/// Builds the polynomial costs of rows `offset + gen_rows[i]` of `costs`, and appends their
/// piecewise linear costs to `piecewise`.
fn prepare(
    costs: &[GenCost],
    gen_rows: &[usize],
    offset: usize,
    base_mva: f64,
    variable: fn(usize) -> CostVariable,
    piecewise: &mut Vec<PiecewiseCost>,
) -> Result<QuadraticCosts, PowerflowError> {
    let mut quadratic = QuadraticCosts::zeros(gen_rows.len());
    for (i, &gen_row) in gen_rows.iter().enumerate() {
        let row = offset + gen_row;
        let unsupported = |reason| PowerflowError::UnsupportedCost { row, reason };
        match &costs[row] {
            GenCost::PiecewiseLinear { points, .. } if points.len() > 2 => {
                // makeAy's rows, with p in per unit
                let mut segment = PiecewiseCost {
                    variable: variable(i),
                    slopes: vec![],
                    bounds: vec![],
                    last_cost: points[points.len() - 1].1,
                };
                for pair in points.windows(2) {
                    let ((x0, c0), (x1, c1)) = (pair[0], pair[1]);
                    if x0 == x1 {
                        return Err(unsupported("two consecutive points at the same output"));
                    }
                    let (p0, p1) = (x0 / base_mva, x1 / base_mva);
                    let m = (c1 - c0) / (p1 - p0);
                    segment.slopes.push(m);
                    segment.bounds.push(m * p0 - c0);
                }
                piecewise.push(segment);
            }
            GenCost::PiecewiseLinear { points, .. } => {
                let &[(x0, y0), (x1, y1)] = points.as_slice() else {
                    return Err(unsupported(
                        "a piecewise linear cost needs two points or more",
                    ));
                };
                if x0 == x1 {
                    return Err(unsupported("two consecutive points at the same output"));
                }
                // opf_setup: one segment becomes the linear polynomial m p + b, in MW
                let m = (y1 - y0) / (x1 - x0);
                quadratic.linear[i] = m * base_mva;
                quadratic.constant[i] = y0 - m * x0;
            }
            GenCost::Polynomial { coeffs, .. } => {
                // opf_setup's reduction: drop zero leading coefficients while more than two
                let mut c = coeffs.as_slice();
                while c.len() > 2 && c[0] == 0.0 {
                    c = &c[1..];
                }
                match *c {
                    [] => {}
                    [k] => quadratic.constant[i] = k,
                    [c1, k] => {
                        quadratic.linear[i] = c1 * base_mva;
                        quadratic.constant[i] = k;
                    }
                    [c2, c1, k] => {
                        quadratic.quadratic[i] = 2.0 * c2 * (base_mva * base_mva);
                        quadratic.linear[i] = c1 * base_mva;
                        quadratic.constant[i] = k;
                    }
                    _ => {
                        return Err(unsupported(
                            "polynomial costs above quadratic aren't supported",
                        ))
                    }
                }
            }
        }
    }
    Ok(quadratic)
}
