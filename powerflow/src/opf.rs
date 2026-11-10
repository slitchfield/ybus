//! Optimal power flow, following MATPOWER's `runopf` with its default solver, MIPS (from the
//! `interior_point` crate), or with the predictor-corrector method from the same crate
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
//! needs no direct dependency on `interior_point`.
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

use std::f64::consts::PI;

// The solvers' types this module's API names, so that choosing and reading a solver needs no
//   direct dependency on interior_point
pub use interior_point::interior_point_error::InteriorPointError;
pub use interior_point::kkt::KktFactorization;
pub use interior_point::mehrotra::{MehrotraIteration, MehrotraOptions};
pub use interior_point::mips::{MipsIteration, MipsOptions};
pub use interior_point::solution::{History, KktResiduals, Multipliers, Solution, Status};
use sparse_mat::compressed_matrix::CsrMatrix;

use crate::case::{BusType, Case};
use crate::opf::cost::{CostVariable, OpfCosts};
use crate::powerflow_error::PowerflowError;
use crate::topology::Network;

/// Which interior-point solver an optimal power flow runs, with its options.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum OpfSolver {
    /// MIPS, ported step for step from MATPOWER's `mips.m`: MATPOWER's default, and this
    /// crate's, whose results match MATPOWER's
    Mips(MipsOptions),
    /// The primal-dual method with Mehrotra's predictor-corrector
    /// ([`interior_point::mehrotra`]), on the same problem, constraints and starting point
    Mehrotra(MehrotraOptions),
}

impl OpfSolver {
    /// The same solver, with its KKT matrix factored by `kkt`: LU by default, or the pivoted
    /// LDLᵀ (see [`interior_point::kkt`]).
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
