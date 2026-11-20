//! Optimal power flow, following MATPOWER's `runopf` with its default solver, MIPS (from the
//! `ybus_ipm` crate), or with the predictor-corrector method from the same crate
//! ([`OpfSolver`]).
//!
//! - [`cost`]: the generator costs, evaluated as MATPOWER's `totcost` evaluates them and
//!   prepared for the solver as `opf_setup` and `makeAy` prepare them.
//! - [`dc`]: the DC optimal power flow ([`dc::dc_opf`]), as `rundcopf` computes it with MIPS.
//! - [`ac`]: the AC optimal power flow ([`ac::ac_opf`]), as `runopf` computes it with MIPS.
//!
//! The solvers' types this module's API names are re-exported here, from their options
//! ([`MipsOptions`], [`MehrotraOptions`], [`KktFactorization`]) to their results ([`Solution`],
//! [`Status`], [`Multipliers`], [`History`], [`KktResiduals`]) and refusals ([`InteriorPointError`]), so that a caller
//! needs no direct dependency on `ybus_ipm`.
//!
//! See phase O1 of `docs/OptimizationChecklist.md` for the roadmap.

pub mod ac;
pub mod cost;
pub mod dc;
// The derivatives of AC power injections and branch flows that the AC OPF uses, as MATPOWER
//   Technical Note 2 gives them. Public only for opf_derivatives_test, and hidden from the docs:
//   it isn't part of the API.
#[doc(hidden)]
pub mod derivatives;

use std::collections::HashMap;
use std::f64::consts::PI;

// The solvers' types this module's API names, so that choosing and reading a solver needs no
//   direct dependency on ybus_ipm
pub use ybus_ipm::interior_point_error::InteriorPointError;
pub use ybus_ipm::kkt::{KktFactorization, KktOrdering, KktStatistics};
pub use ybus_ipm::mehrotra::{MehrotraIteration, MehrotraOptions};
pub use ybus_ipm::mips::{MipsIteration, MipsOptions};
pub use ybus_ipm::solution::{History, KktResiduals, Multipliers, Solution, Status};
use ybus_ipm::start::Start;
pub use ybus_ipm::start::WarmStart;
use ybus_sparse::compressed_matrix::CsrMatrix;

use crate::case::{BusType, Case};
use crate::opf::cost::{CostVariable, OpfCosts};
use crate::powerflow_error::PowerflowError;
use crate::topology::Network;

/// Where a kept model's solve starts ([`ac::AcOpfModel::solve_from`],
/// [`dc::DcOpfModel::solve_from`]), for a sequence of related cases on one network.
///
/// A previous solution is the `solver` field of an [`ac::AcOpfSolution`] or
/// [`dc::DcOpfSolution`], from the same model and from options with the same `cost_mult` (its
/// multipliers are reported unscaled, and scaled back by the new solve's). A solution that
/// doesn't fit the model's problem is refused, as `PowerflowError::Solver` with the solver's
/// `Dimension` error. See [`ybus_ipm::start`] for what a warm start carries and why it shifts.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub enum OpfStart<'a> {
    /// MATPOWER's starting point, as [`ac::AcOpfModel::solve`] and [`dc::DcOpfModel::solve`]
    /// start: exactly their path
    Cold,
    /// A previous solution's point alone, the multipliers and slacks started as a cold start's
    Primal(&'a Solution),
    /// A previous solution's point and multipliers, shifted off the boundary by the solver's
    /// `warm_shift`
    Warm(&'a Solution),
    /// A point and multipliers already in this model's layout, as `map_start` maps a solution
    /// of another model of the case onto it ([`ac::AcOpfModel::map_start`],
    /// [`dc::DcOpfModel::map_start`]), for example across an outage: shifted as `Warm` is
    Mapped(&'a WarmStart),
    /// [`OpfStart::Mapped`]'s point alone, the multipliers and slacks started as a cold start's
    MappedPrimal(&'a WarmStart),
}

impl<'a> OpfStart<'a> {
    /// The solver's start: from `x0`, MATPOWER's point, for a cold start, and from `warm`, which
    /// [`OpfStart::warm_start`] built, for a warm one
    pub(crate) fn solver_start(self, x0: &'a [f64], warm: &'a Option<WarmStart>) -> Start<'a> {
        match (self, warm) {
            (OpfStart::Primal(previous), _) => Start::Cold(&previous.x),
            (OpfStart::Warm(_), Some(warm)) => Start::Warm(warm),
            (OpfStart::Mapped(mapped), _) => Start::Warm(mapped),
            (OpfStart::MappedPrimal(mapped), _) => Start::Cold(&mapped.x),
            _ => Start::Cold(x0),
        }
    }

    /// A warm start's point and multipliers, copied from its previous solution
    pub(crate) fn warm_start(self) -> Option<WarmStart> {
        match self {
            OpfStart::Warm(previous) => Some(WarmStart::from_solution(previous)),
            _ => None,
        }
    }
}

/// What a variable or a constraint of an OPF problem belongs to, by the case's row: the key on
/// which `map_start` matches one model's problem to another's
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Key {
    /// A bus's voltage angle, and for AC its magnitude
    Va(usize),
    Vm(usize),
    /// A generator's real and (AC) reactive output
    Pg(usize),
    Qg(usize),
    /// The helper variable of a generator's piecewise cost of its real (`true`) or reactive
    /// output
    Y(usize, bool),
    /// A bus's real and (AC) reactive power balance
    P(usize),
    Q(usize),
    /// A branch's flow limit: AC's from and to ends, DC's one row
    Sf(usize),
    St(usize),
    Flow(usize),
    /// A branch's angle difference limit
    Angle(usize),
    /// Segment `j` of a generator's piecewise cost of its real (`true`) or reactive output
    Segment(usize, bool, usize),
    /// A dispatchable load's power factor
    PowerFactor(usize),
}

/// An OPF problem's variables and constraints, each by its [`Key`], in the problem's order:
/// the columns, the nonlinear equalities and inequalities, and the rows of `A`
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Layout {
    pub(crate) columns: Vec<Key>,
    pub(crate) eq: Vec<Key>,
    pub(crate) ineq: Vec<Key>,
    pub(crate) rows: Vec<Key>,
}

impl Layout {
    /// The keys of `costs`' piecewise costs, for generators at case rows `gen_rows`: their
    ///   helper columns, and their segments' rows, in `LinearRows::piecewise`'s order
    pub(crate) fn piecewise(&mut self, costs: &OpfCosts, gen_rows: &[usize]) {
        let keys: Vec<(usize, bool)> = costs
            .piecewise
            .iter()
            .map(|cost| match cost.variable {
                CostVariable::Active(k) => (gen_rows[k], true),
                CostVariable::Reactive(k) => (gen_rows[k], false),
            })
            .collect();
        self.columns
            .extend(keys.iter().map(|&(row, active)| Key::Y(row, active)));
        for (&(row, active), cost) in keys.iter().zip(&costs.piecewise) {
            self.rows
                .extend((0..cost.slopes.len()).map(|j| Key::Segment(row, active, j)));
        }
    }

    /// A warm start in this layout from `previous`, a solution in `from`'s: each variable and
    /// multiplier whose key `from` has takes its value there, from the final iterate's
    /// multipliers; a new variable takes its value in `cold`, this problem's cold start, and a
    /// new multiplier 0.
    ///
    /// Returns `PowerflowError::Solver` with the solver's `Dimension` error if `previous`
    /// doesn't fit `from`.
    pub(crate) fn map(
        &self,
        from: &Layout,
        previous: &Solution,
        cold: &[f64],
    ) -> Result<WarmStart, PowerflowError> {
        let m = &previous.iterate_multipliers;
        let fits = previous.x.len() == from.columns.len()
            && m.lower.len() == from.columns.len()
            && m.upper.len() == from.columns.len()
            && m.eq_nonlin.len() == from.eq.len()
            && m.ineq_nonlin.len() == from.ineq.len()
            && m.mu_l.len() == from.rows.len()
            && m.mu_u.len() == from.rows.len();
        if !fits {
            return Err(PowerflowError::Solver {
                error: InteriorPointError::Dimension {
                    what: "the previous solution",
                },
            });
        }
        // Each key's position in `from`, then each of ours taken from there
        let index = |keys: &[Key]| -> HashMap<Key, usize> {
            keys.iter().enumerate().map(|(i, &key)| (key, i)).collect()
        };
        let take = |ours: &[Key], theirs: &[Key], values: &[f64], new: &dyn Fn(usize) -> f64| {
            let positions = index(theirs);
            ours.iter()
                .enumerate()
                .map(|(i, key)| positions.get(key).map_or_else(|| new(i), |&j| values[j]))
                .collect::<Vec<f64>>()
        };
        let zero = |_| 0.0;
        Ok(WarmStart {
            x: take(&self.columns, &from.columns, &previous.x, &|i| cold[i]),
            multipliers: Multipliers {
                eq_nonlin: take(&self.eq, &from.eq, &m.eq_nonlin, &zero),
                ineq_nonlin: take(&self.ineq, &from.ineq, &m.ineq_nonlin, &zero),
                mu_l: take(&self.rows, &from.rows, &m.mu_l, &zero),
                mu_u: take(&self.rows, &from.rows, &m.mu_u, &zero),
                lower: take(&self.columns, &from.columns, &m.lower, &zero),
                upper: take(&self.columns, &from.columns, &m.upper, &zero),
            },
        })
    }
}

/// Which interior-point solver an optimal power flow runs, with its options.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum OpfSolver {
    /// MIPS, ported step for step from MATPOWER's `mips.m`: MATPOWER's default, and this
    /// crate's, whose results match MATPOWER's
    Mips(MipsOptions),
    /// The primal-dual method with Mehrotra's predictor-corrector
    /// ([`ybus_ipm::mehrotra`]), on the same problem, constraints and starting point
    Mehrotra(MehrotraOptions),
}

impl OpfSolver {
    /// The same solver, with its KKT matrix factored by `kkt`: LU by default, or the pivoted
    /// LDLᵀ (see [`ybus_ipm::kkt`]).
    pub fn with_kkt_factorization(self, kkt: KktFactorization) -> OpfSolver {
        match self {
            OpfSolver::Mips(mut options) => {
                options.kkt_factorization = kkt;
                OpfSolver::Mips(options)
            }
            OpfSolver::Mehrotra(mut options) => {
                options.kkt_factorization = kkt;
                OpfSolver::Mehrotra(options)
            }
        }
    }

    /// The same solver, with the pivoted LDLᵀ ordering its KKT matrix by `ordering`: AMD by
    /// default, or matched ([`KktOrdering`]). LU ignores it.
    pub fn with_kkt_ordering(self, ordering: KktOrdering) -> OpfSolver {
        match self {
            OpfSolver::Mips(mut options) => {
                options.kkt_ordering = ordering;
                OpfSolver::Mips(options)
            }
            OpfSolver::Mehrotra(mut options) => {
                options.kkt_ordering = ordering;
                OpfSolver::Mehrotra(options)
            }
        }
    }
}

/// The connected buses of type REF (internal indices), whose angles the OPF fixes, as
/// `opf_setup` finds them after `ext2int`, with no `bustypes`: a REF bus needs no generator.
/// `PowerflowError::NoReferenceBus` if there are none, where MATPOWER fails too (its OPF
/// solvers start every angle at the first reference bus's).
fn reference_buses(case: &Case, network: &Network) -> Result<Vec<usize>, PowerflowError> {
    let refs: Vec<usize> = (0..network.bus_rows.len())
        .filter(|&i| case.bus[network.bus_rows[i]].bus_type == BusType::Ref)
        .collect();
    if refs.is_empty() {
        Err(PowerflowError::NoReferenceBus)
    } else {
        Ok(refs)
    }
}

/// The OPF solvers' starting point within bounds `xmin` and `xmax`, as `dcopf_solver` and
/// `nlpopf_solver` set it: midway between the bounds (infinite ones replaced by ±1e10), or 1
/// inside a one-sided bound. The callers then set the angles, magnitudes and helper variables.
fn starting_point(xmin: &[f64], xmax: &[f64]) -> Vec<f64> {
    xmin.iter()
        .zip(xmax)
        .map(
            |(&lo, &hi)| match (lo == f64::NEG_INFINITY, hi == f64::INFINITY) {
                (true, false) => hi - 1.0,
                (false, true) => lo + 1.0,
                (true, true) => 0.0,
                (false, false) => (lo + hi) / 2.0,
            },
        )
        .collect()
}

/// Where the OPF solvers start every piecewise linear cost's helper variable: the largest of
/// the costs' last values, plus 10% of its magnitude. `None` without piecewise linear costs.
fn helper_start(costs: &OpfCosts) -> Option<f64> {
    costs
        .piecewise
        .iter()
        .map(|cost| cost.last_cost)
        .reduce(f64::max)
        .map(|largest| largest + 0.1 * largest.abs())
}

/// Linear constraints `l ≤ A x ≤ u`, built row by row.
#[derive(Default)]
struct LinearRows {
    triplets: Vec<(usize, usize, f64)>,
    l: Vec<f64>,
    u: Vec<f64>,
}

impl LinearRows {
    /// Appends the row with `entries` (column, value) and bounds `lower` and `upper`.
    fn push(&mut self, entries: impl IntoIterator<Item = (usize, f64)>, lower: f64, upper: f64) {
        let row = self.l.len();
        self.triplets
            .extend(entries.into_iter().map(|(col, value)| (row, col, value)));
        self.l.push(lower);
        self.u.push(upper);
    }

    fn len(&self) -> usize {
        self.l.len()
    }

    /// `makeAang`'s angle difference limits on the angles starting at column `va`, returning
    /// the limited branches (internal indices) in row order.
    fn angle_limits(&mut self, case: &Case, network: &Network, va: usize) -> Vec<usize> {
        let limited = dc::angle_limits(case, &network.branch_rows);
        for &(k, lower, upper) in &limited {
            let (from, to) = network.ends[k];
            self.push([(va + from, 1.0), (va + to, -1.0)], lower, upper);
        }
        limited.iter().map(|&(k, _, _)| k).collect()
    }

    /// `makeAy`'s segment constraints `mⱼ p − y ≤ bⱼ` for `costs`' piecewise linear costs, with
    /// real outputs starting at column `pg`, reactive ones at `qg`, and the helper variables at
    /// `y`.
    fn piecewise(&mut self, costs: &OpfCosts, pg: usize, qg: Option<usize>, y: usize) {
        for (j, cost) in costs.piecewise.iter().enumerate() {
            let column = match cost.variable {
                CostVariable::Active(k) => pg + k,
                CostVariable::Reactive(k) => {
                    qg.expect("reactive costs only where there are reactive outputs") + k
                }
            };
            for (&m, &b) in cost.slopes.iter().zip(&cost.bounds) {
                self.push([(column, m), (y + j, -1.0)], f64::NEG_INFINITY, b);
            }
        }
    }

    /// The matrix with `nx` columns, and the bounds.
    fn finish(self, nx: usize) -> (CsrMatrix<f64>, Vec<f64>, Vec<f64>) {
        let a =
            CsrMatrix::from_triplets(self.l.len(), nx, self.triplets).expect("columns below nx");
        (a, self.l, self.u)
    }
}

/// A multiplier of an angle difference limit, from $/rad to $/degree as MATPOWER reports it.
fn per_degree(multiplier: f64) -> f64 {
    multiplier * PI / 180.0
}
