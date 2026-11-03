//! Interior-point solvers for nonlinear programs:
//!
//! minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.
//!
//! - [`problem`]: the programs the solvers take: a [`problem::Problem`] supplies the objective
//!   and any nonlinear constraints with their derivatives, and [`problem::Constraints`] the
//!   linear constraints and bounds.
//! - [`mips`]: MIPS, the MATPOWER Interior Point Solver, reimplemented step for step from
//!   MATPOWER's `mips.m`, so that its iterations and results can be checked against MATPOWER's.
//!
//! The linear algebra is [`sparse_mat`]'s. See `docs/OptimizationChecklist.md` for the roadmap.

pub mod mips;
pub mod problem;
