//! Optimal power flow, following MATPOWER's `runopf` with its default solver, MIPS (from the
//! `interior_point` crate).
//!
//! - [`cost`]: the generator costs, evaluated as MATPOWER's `totcost` evaluates them and
//!   prepared for the solver as `opf_setup` and `makeAy` prepare them.
//! - [`dc`]: the DC optimal power flow ([`dc::dc_opf`]), as `rundcopf` computes it with MIPS.
//!
//! See phase O1 of `docs/OptimizationChecklist.md` for the roadmap.

pub mod cost;
pub mod dc;
