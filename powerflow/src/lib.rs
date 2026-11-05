//! DC and AC power flow, built on [`sparse_mat`].
//!
//! - [`case`]: the model of a power system (buses, generators, branches and costs), and a reader
//!   for MATPOWER case files ([`case::Case::from_path`]).
//! - [`dc`]: the DC power flow ([`dc::dc_power_flow`]), as MATPOWER's `rundcpf` computes it.
//! - [`ac`]: the AC power flow by Newton's method ([`ac::ac_power_flow`]), as MATPOWER's
//!   `runpf` computes it with its default options.
//! - [`opf`]: the optimal power flow, as MATPOWER's `runopf` computes it with MIPS; so far,
//!   the generator costs ([`opf::cost`]).
//! - [`islands`]: a network's islands and bridges ([`islands::find_islands`],
//!   [`islands::find_bridges`]), as MATPOWER's `find_islands` and `find_bridges` define them.
//!
//! Both power flows are validated against MATPOWER's own results; see the roadmap in
//! `docs/PowerFlowChecklist.md`.

// No unsafe code here: the GPU path's FFI and kernel launches live in ybus_cuda, behind a safe
//   API
#![forbid(unsafe_code)]

pub mod ac;
pub mod case;
pub mod dc;
pub mod islands;
pub mod opf;
pub mod powerflow_error;
mod topology;
