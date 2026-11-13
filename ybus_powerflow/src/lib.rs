//! DC and AC power flow, built on [`ybus_sparse`].
//!
//! - [`case`]: the model of a power system (buses, generators, branches and costs), and a reader
//!   for MATPOWER case files ([`case::Case::from_path`]).
//! - [`dc`]: the DC power flow ([`dc::dc_power_flow`]), as MATPOWER's `rundcpf` computes it.
//! - [`ac`]: the AC power flow by Newton's method ([`ac::ac_power_flow`]), as MATPOWER's
//!   `runpf` computes it with its default options.
//! - [`opf`]: the DC and AC optimal power flows ([`opf::dc::dc_opf`], [`opf::ac::ac_opf`]), as
//!   MATPOWER's `rundcopf` and `runopf` compute them with MIPS, or by the predictor-corrector
//!   method; and the generator costs ([`opf::cost`]).
//! - [`islands`]: a network's islands and bridges ([`islands::find_islands`],
//!   [`islands::find_bridges`]), as MATPOWER's `find_islands` and `find_bridges` define them.
//! - With the `cuda` feature, the DC power flow on the GPU: batched solves, single-branch
//!   outages and PTDF columns (`dc::DcStudy::gpu` and `dc::DcSensitivity::gpu`), on a `Gpu`
//!   from `ybus_cuda`. (Plain code, not links: the items exist only with the feature.)
//!
//! [`Case`] and [`PowerflowError`] are re-exported here.
//!
//! Both power flows are validated against MATPOWER's own results; see the roadmap in
//! `docs/PowerFlowChecklist.md`.

// No unsafe code here: the GPU path's FFI and kernel launches live in ybus_cuda, behind a safe
//   API
#![forbid(unsafe_code)]
// Labels each feature-gated item with the feature it needs, on docs.rs, which builds with
//   nightly and `--cfg docsrs` ([package.metadata.docs.rs] in Cargo.toml). Items gated directly
//   are labelled automatically; those inside a gated module carry their own `doc(cfg)`.
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod ac;
pub mod case;
pub mod dc;
pub mod islands;
pub mod opf;
pub mod powerflow_error;
mod topology;

pub use case::Case;
pub use powerflow_error::PowerflowError;

/// The CUDA device the GPU path runs on (feature `cuda`).
#[cfg(feature = "cuda")]
pub use ybus_cuda::Gpu;

/// Runs the README's Rust examples as doctests, so they can't drift from the API. Compiled only
/// for `cargo test --doc`: the README is written for crates.io, and isn't the crate docs.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;
