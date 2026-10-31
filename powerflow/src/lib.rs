//! DC and AC power flow, built on [`sparse_mat`].
//!
//! So far: the [`case`] model of a power system (buses, generators, branches and costs) and
//! a reader for MATPOWER case files ([`case::Case::from_path`]). DC power flow is next; see
//! Phase 1.5 of the roadmap in `docs/PowerFlowChecklist.md`.

pub mod case;
pub mod dc;
pub mod powerflow_error;
