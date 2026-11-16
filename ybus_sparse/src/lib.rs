//! Sparse matrices: assembly, products, factorizations and solvers, and Matrix Market I/O,
//! generic over the element type (`f64`, `f32`, integers, [`Complex`]).
//!
//! # Example
//!
//! ```
//! use ybus_sparse::{CsrMatrix, Factorization, SparseError, SparseMatrix};
//!
//! # fn main() -> Result<(), SparseError> {
//! // Assemble from (row, col, value) triplets; entries at the same position are summed
//! let a: CsrMatrix = SparseMatrix::from_triplets(
//!     3,
//!     3,
//!     vec![
//!         (0, 0, 2.0), (0, 1, -1.0),
//!         (1, 0, -1.0), (1, 1, 2.0), (1, 2, -1.0),
//!         (2, 1, -1.0), (2, 2, 2.0),
//!     ],
//! )?
//! .to_csr();
//! // Factor once, then solve
//! let cholesky = a.cholesky()?;
//! let x = cholesky.solve(&[1.0, 0.0, 1.0])?;
//! assert!(x.iter().all(|x| (x - 1.0).abs() < 1e-12));
//! # Ok(())
//! # }
//! ```
//!
//! # Guide
//!
//! - **Building:** [`SparseMatrix`] builds and edits a matrix. Its entries live in a hash map,
//!   so inserting, removing and looking up single entries is cheap. It reads and writes Matrix
//!   Market files ([`SparseMatrix::read_matrix_market`]).
//! - **Numerical work:** a [`CompressedMatrix`], as a [`CsrMatrix`] or a [`CscMatrix`], with
//!   the orientation in the type ([`compressed_matrix`]). Its pattern is fixed when it's
//!   built, so values can be refilled without changing its structure. It has the products,
//!   norms, selections and arithmetic, and the solvers:
//!   - LU with threshold partial pivoting ([`CscMatrix::lu`], or [`CscMatrix::solve`] for one
//!     solve)
//!   - Cholesky for symmetric positive definite matrices ([`CsrMatrix::cholesky`])
//!   - LDLᵀ with Bunch–Kaufman pivoting for symmetric indefinite ones ([`CsrMatrix::ldlt`])
//!   - conjugate gradient ([`CompressedMatrix::conjugate_gradient`])
//! - **Solves:** every factorization's come from one trait, [`Factorization`]: one or many
//!   right-hand sides, sparse right-hand sides, and refinement.
//! - **Errors:** every fallible operation returns a [`SparseError`].
//!
//! # Features
//!
//! - `serde`: `Serialize` and `Deserialize` for both matrix types.
//! - `rayon`: parallel products, prefixed `par_`, whose results are bit for bit the
//!   sequential versions'.
//!
//! See the [README](https://crates.io/crates/ybus_sparse) for a longer tour, and the
//! [repository](https://github.com/slitchfield/ybus) for the rest of the `ybus` workspace.

// No unsafe code here: the GPU path's FFI and kernel launches live in ybus_cuda, behind a safe
//   API
#![forbid(unsafe_code)]
// Labels each feature-gated item with the feature it needs, on docs.rs, which builds with
//   nightly and `--cfg docsrs` ([package.metadata.docs.rs] in Cargo.toml). Items gated directly
//   are labelled automatically; those inside a gated module carry their own `doc(cfg)`.
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod compressed_matrix;

/// The `SparseMatrix` type and its iterators.
pub mod sparse_matrix;

/// Element types a `SparseMatrix` can hold: any `Copy` numeric type (`f32`, `f64`, integers,
/// `num_complex::Complex`, ...) with the arithmetic operations and zero/one of `num_traits::Num`.
pub trait Scalar: Copy + PartialEq + num_traits::Num {}

impl<T: Copy + PartialEq + num_traits::Num> Scalar for T {}

/// The error type returned by `SparseMatrix` operations.
pub mod sparse_error;

pub use compressed_matrix::{CompressedMatrix, CscMatrix, CsrMatrix, Factorization};
pub use sparse_error::SparseError;
pub use sparse_matrix::SparseMatrix;

/// Re-exported so callers use the same `num-complex` version as this crate.
pub use num_complex;
/// Re-exported for convenience: `SparseMatrix<Complex<f64>>` is a first-class matrix type.
pub use num_complex::Complex;

/// Runs the README's Rust examples as doctests, so they can't drift from the API. Compiled only
/// for `cargo test --doc`: the README is written for crates.io, and isn't the crate docs.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;
