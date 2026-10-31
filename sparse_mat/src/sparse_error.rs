use std::error::Error;
use std::fmt;

/// Errors returned by `SparseMatrix` operations.
#[derive(Debug, PartialEq)]
#[non_exhaustive]
pub enum SparseError {
    /// An index lies outside the matrix.
    OutOfBounds {
        /// The offending row index.
        row: usize,
        /// The offending column index.
        col: usize,
        /// The shape `(rows, cols)` of the matrix the index was checked against.
        shape: (usize, usize),
    },
    /// Dense input rows do not all have the same length.
    RaggedRows {
        /// The first row whose length differs from the first row's.
        row: usize,
        /// The length of the first row, which every row must match.
        expected: usize,
        /// The length of the offending row.
        found: usize,
    },
    /// A row or column index lies outside the matrix.
    OutOfRange {
        /// The offending row or column index.
        index: usize,
        /// The number of rows or columns it was checked against.
        len: usize,
    },
    /// A range is reversed or extends past the end of the matrix.
    InvalidRange {
        /// The start of the range.
        start: usize,
        /// The (exclusive) end of the range.
        end: usize,
        /// The number of rows or columns the range was checked against.
        len: usize,
    },
    /// Raw compressed (CSR or CSC) arrays are inconsistent with each other or with the
    /// requested shape.
    InvalidCompressed {
        /// What is wrong with the arrays.
        reason: &'static str,
    },
    /// A list meant to be a permutation of `0..len` is not one.
    InvalidPermutation {
        /// What is wrong with it.
        reason: &'static str,
    },
    /// The operation needs a square matrix.
    NotSquare {
        /// The shape `(rows, cols)` of the matrix.
        shape: (usize, usize),
    },
    /// A zero (or missing) pivot made the matrix singular to working precision.
    SingularMatrix {
        /// The row/column index of the zero pivot.
        pivot: usize,
    },
    /// The matrix is not symmetric, but the operation requires it.
    NotSymmetric,
    /// A non-positive pivot showed the matrix is not positive definite.
    NotPositiveDefinite {
        /// The row/column index of the offending pivot.
        pivot: usize,
    },
    /// An iterative method did not reach the requested tolerance.
    NoConvergence {
        /// How many iterations were run.
        iterations: usize,
    },
    /// Text being read is not valid Matrix Market.
    Parse {
        /// The 1-based line number the problem was found on.
        line: usize,
        /// What is wrong.
        message: String,
    },
    /// Reading or writing failed.
    Io {
        /// The underlying error's message.
        message: String,
    },
    /// The result's dimensions do not fit in a `usize`.
    ShapeOverflow,
    /// The operands' shapes are incompatible for the operation.
    ShapeMismatch {
        /// The shape `(rows, cols)` of the left operand.
        left: (usize, usize),
        /// The shape `(rows, cols)` of the right operand.
        right: (usize, usize),
    },
}

/// Formats a short, lowercase description of the error.
impl fmt::Display for SparseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SparseError::OutOfBounds { row, col, shape } => {
                write!(f, "index ({row}, {col}) out of bounds for shape {shape:?}")
            }
            SparseError::RaggedRows {
                row,
                expected,
                found,
            } => {
                write!(f, "row {row} has length {found}, expected {expected}")
            }
            SparseError::OutOfRange { index, len } => {
                write!(f, "index {index} out of range for length {len}")
            }
            SparseError::InvalidRange { start, end, len } => {
                write!(f, "range {start}..{end} invalid for length {len}")
            }
            SparseError::InvalidCompressed { reason } => {
                write!(f, "invalid compressed arrays: {reason}")
            }
            SparseError::InvalidPermutation { reason } => {
                write!(f, "invalid permutation: {reason}")
            }
            SparseError::NotSquare { shape } => {
                write!(f, "matrix of shape {shape:?} is not square")
            }
            SparseError::SingularMatrix { pivot } => {
                write!(f, "matrix is singular (zero pivot at {pivot})")
            }
            SparseError::NotSymmetric => write!(f, "matrix is not symmetric"),
            SparseError::NotPositiveDefinite { pivot } => {
                write!(f, "matrix is not positive definite (pivot {pivot})")
            }
            SparseError::NoConvergence { iterations } => {
                write!(f, "no convergence after {iterations} iterations")
            }
            SparseError::Parse { line, message } => {
                write!(f, "line {line}: {message}")
            }
            SparseError::Io { message } => write!(f, "i/o error: {message}"),
            SparseError::ShapeOverflow => {
                write!(f, "result dimensions overflow usize")
            }
            SparseError::ShapeMismatch { left, right } => {
                write!(f, "incompatible shapes {left:?} and {right:?}")
            }
        }
    }
}

/// Lets `SparseError` be used with `?` into `Box<dyn Error>` and similar.
impl Error for SparseError {}
