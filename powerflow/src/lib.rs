//! DC and AC power flow, built on [`sparse_mat`].
//!
//! - [`case`]: the model of a power system (buses, generators, branches and costs), and a reader
//!   for MATPOWER case files ([`case::Case::from_path`]).
//! - [`dc`]: the DC power flow ([`dc::dc_power_flow`]), as MATPOWER's `rundcpf` computes it.
//! - [`ac`]: the AC power flow by Newton's method ([`ac::ac_power_flow`]), as MATPOWER's
//!   `runpf` computes it with its default options.
//!
//! Both power flows are validated against MATPOWER's own results; see the roadmap in
//! `docs/PowerFlowChecklist.md`.

pub mod ac;
pub mod case;
pub mod dc;
pub mod powerflow_error;
mod topology;
