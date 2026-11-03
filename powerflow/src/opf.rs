//! Optimal power flow, following MATPOWER's `runopf` with its default solver, MIPS (from the
//! `interior_point` crate).
//!
//! - [`cost`]: the generator costs, evaluated as MATPOWER's `totcost` evaluates them and
//!   prepared for the solver as `opf_setup` and `makeAy` prepare them.
//!
//! See phase O1 of `docs/OptimizationChecklist.md` for the roadmap.

pub mod cost;
