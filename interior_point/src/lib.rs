//! Interior-point solvers for nonlinear programs:
//!
//! minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.
//!
//! - [`problem`]: the programs the solvers take: a [`problem::Problem`] supplies the objective
//!   and any nonlinear constraints with their derivatives, and [`problem::Constraints`] the
//!   linear constraints and bounds.
//! - [`mips`]: MIPS, the MATPOWER Interior Point Solver, reimplemented step for step from
//!   MATPOWER's `mips.m`, so that its iterations and results can be checked against MATPOWER's.
//! - [`interior_point_error`]: why a solver couldn't start.
//!
//! The linear algebra is [`sparse_mat`]'s. See `docs/OptimizationChecklist.md` for the roadmap.

// No unsafe code here: the GPU path's FFI and kernel launches live in ybus_cuda, behind a safe
//   API
#![forbid(unsafe_code)]

pub mod interior_point_error;
pub mod mips;
pub mod problem;
mod vector;
