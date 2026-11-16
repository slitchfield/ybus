//! Interior-point solvers for nonlinear programs:
//!
//! minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.
//!
//! # Example
//!
//! A quadratic program, minimize `½ ‖x‖² − x₁ − 2 x₂` with `x₁ ≤ 0.5`, solved by MIPS:
//!
//! ```
//! use ybus_ipm::mips::{qps_mips, MipsOptions};
//! use ybus_ipm::{Constraints, QuadraticProblem};
//! use ybus_sparse::SparseMatrix;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let identity = SparseMatrix::from_triplets(2, 2, vec![(0, 0, 1.0), (1, 1, 1.0)])?.to_csr();
//! let qp = QuadraticProblem { h: Some(identity), c: vec![-1.0, -2.0] };
//! let constraints = Constraints {
//!     xmax: Some(vec![0.5, f64::INFINITY]),
//!     ..Constraints::default()
//! };
//! let solution = qps_mips(&qp, &constraints, None, &MipsOptions::default())?;
//! assert!(solution.converged());
//! assert!((solution.x[0] - 0.5).abs() < 1e-6 && (solution.x[1] - 2.0).abs() < 1e-6);
//! // The bound on x₁ is active, with multiplier 0.5: relaxing it by δ lowers f by about δ / 2
//! assert!((solution.multipliers.upper[0] - 0.5).abs() < 1e-6);
//! # Ok(())
//! # }
//! ```
//!
//! # Guide
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
//! # Which solver
//!
//! - **MIPS** when results must match MATPOWER's, which they do with the default LU
//!   factorization of the KKT systems.
//! - **The predictor-corrector method** for fewer iterations, written from the literature
//!   rather than to match MATPOWER.
//!
//! Either factors its KKT systems by LU or by LDLᵀ, chosen by `kkt_factorization` in its
//! options ([`kkt::KktFactorization`]). With LDLᵀ, the predictor-corrector method also
//! corrects an indefinite Hessian using the factorization's inertia.
//!
//! A general nonlinear program implements [`Problem`], with its derivatives; a quadratic or
//! linear one is a [`QuadraticProblem`], for `qps_mips` and `qps_mehrotra`.
//!
//! The linear algebra is [`ybus_sparse`]'s. See the
//! [roadmap](https://github.com/slitchfield/ybus/blob/main/docs/OptimizationChecklist.md) for
//! what's planned.

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

/// Runs the README's Rust examples as doctests, so they can't drift from the API. Compiled only
/// for `cargo test --doc`: the README is written for crates.io, and isn't the crate docs.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;
