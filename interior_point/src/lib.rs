//! Interior-point solvers for nonlinear programs:
//!
//! minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.
//!
//! - [`problem`]: the programs the solvers take: a [`problem::Problem`] supplies the objective
//!   and any nonlinear constraints with their derivatives, and [`problem::Constraints`] the
//!   linear constraints and bounds.
//! - [`mips`]: MIPS, the MATPOWER Interior Point Solver, reimplemented step for step from
//!   MATPOWER's `mips.m`, so that its iterations and results can be checked against MATPOWER's.
//! - [`solution`]: the result both solvers return, a [`Solution`] with its [`Status`] and
//!   multipliers, and [`solution::kkt_residuals`], a check of the KKT conditions at a point that
//!   doesn't depend on the solver.
//! - [`mehrotra`]: a primal-dual method with Mehrotra's predictor-corrector, written from the
//!   literature, for the same programs.
//! - [`kkt`]: how the solvers factor the KKT system, by LU (the default) or by LDLᵀ.
//! - [`interior_point_error`]: why a solver couldn't start.
//!
//! The problem types, the result and the error type are re-exported here: [`Problem`],
//! [`QuadraticProblem`], [`Constraints`], [`Solution`], [`Status`] and [`InteriorPointError`].
//!
//! Every solver takes its arguments in the same order: the problem, its `constraints`, the
//! starting point `x0` and the options.
//!
//! The linear algebra is [`sparse_mat`]'s. See `docs/OptimizationChecklist.md` for the roadmap.

// No unsafe code here: the GPU path's FFI and kernel launches live in ybus_cuda, behind a safe
//   API
#![forbid(unsafe_code)]

pub mod interior_point_error;
pub mod kkt;
pub mod mehrotra;
pub mod mips;
pub mod problem;
pub mod solution;
mod vector;

pub use interior_point_error::InteriorPointError;
pub use problem::{Constraints, Problem, QuadraticProblem};
pub use solution::{Solution, Status};
