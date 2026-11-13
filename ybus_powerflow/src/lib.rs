//! DC and AC power flow and optimal power flow, built on [`ybus_sparse`] and [`ybus_ipm`],
//! each validated against MATPOWER's own results.
//!
//! # Example
//!
//! Two buses, with a generator at the reference bus and 50 MW of load at the other:
//!
//! ```
//! use ybus_powerflow::ac::ac_power_flow;
//! use ybus_powerflow::dc::dc_power_flow;
//! use ybus_powerflow::Case;
//!
//! # fn main() -> Result<(), ybus_powerflow::PowerflowError> {
//! let case = Case::parse(
//!     "mpc.version = '2';
//!      mpc.baseMVA = 100;
//!      mpc.bus = [1 3 0 0 0 0 1 1 0 135 1 1.1 0.9; 2 1 50 10 0 0 1 1 0 135 1 1.1 0.9];
//!      mpc.gen = [1 0 0 100 -100 1 100 1 200 0 0 0 0 0 0 0 0 0 0 0 0];
//!      mpc.branch = [1 2 0.01 0.1 0 0 0 0 0 0 1 -360 360];",
//! )?; // or Case::from_path("case9.m")
//! let dc = dc_power_flow(&case)?;
//! assert!((dc.pg_mw[0] - 50.0).abs() < 1e-9); // lossless
//! let ac = ac_power_flow(&case)?;
//! assert!(ac.pg_mw[0] > 50.0); // the generator also supplies the losses
//! assert!(ac.vm_pu[1] < 1.0); // and the load bus sags
//! # Ok(())
//! # }
//! ```
//!
//! # Guide
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
//! # Where to start
//!
//! - **One solve:** [`dc::dc_power_flow`], [`ac::ac_power_flow`], [`opf::dc::dc_opf`] or
//!   [`opf::ac::ac_opf`], each taking a [`Case`]; the `_with_options` variants take the
//!   options.
//! - **Many solves on one network:** build the model once ([`dc::DcModel`], [`ac::AcModel`])
//!   and solve it for each scenario. For the DC power flow, a [`dc::DcStudy`] keeps its
//!   factorization, for many injections and single-branch outages, and a
//!   [`dc::DcSensitivity`] gives PTDF and LODF columns from it.
//! - **Matching MATPOWER:** each solve follows the MATPOWER function it names, with
//!   MATPOWER's defaults, and the repository's tests compare its results with MATPOWER's.
//!
//! The [repository](https://github.com/slitchfield/ybus) holds the tests against MATPOWER, the
//! benchmarks and the finished
//! [roadmap](https://github.com/slitchfield/ybus/blob/main/docs/PowerFlowChecklist.md).

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
