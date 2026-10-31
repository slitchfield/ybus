//! A sparse matrix library.
//!
//! See [`sparse_matrix::SparseMatrix`] for the matrix type and [`sparse_error::SparseError`]
//! for its errors.

pub mod compressed_matrix;

/// The `SparseMatrix` type and its iterators.
pub mod sparse_matrix;

/// Element types a `SparseMatrix` can hold: any `Copy` numeric type (`f32`, `f64`, integers,
/// `num_complex::Complex`, ...) with the arithmetic operations and zero/one of `num_traits::Num`.
pub trait Scalar: Copy + PartialEq + num_traits::Num {}

impl<T: Copy + PartialEq + num_traits::Num> Scalar for T {}

/// The error type returned by `SparseMatrix` operations.
pub mod sparse_error;

/// Re-exported so callers use the same `num-complex` version as this crate.
pub use num_complex;
/// Re-exported for convenience: `SparseMatrix<Complex<f64>>` is a first-class matrix type.
pub use num_complex::Complex;

/// Runs the README's Rust examples as doctests, so they can't drift from the API. Compiled only
/// for `cargo test --doc`; the README is not the crate docs, whose relative links to `docs/`
/// wouldn't resolve in rustdoc.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;
