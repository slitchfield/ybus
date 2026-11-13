//! Reading and writing the Matrix Market coordinate format (`.mtx`).
//!
//! Supported: `matrix coordinate` files with a `real`, `integer` or `pattern` field and
//! `general`, `symmetric` or `skew-symmetric` symmetry. Dense (`array`) files, `complex`
//! values and `hermitian` symmetry are rejected with a parse error. Writing always produces
//! `real general`.

use std::fmt::Display;
use std::io::{BufRead, Write};

use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

/// The most entries reserved up front from a file's size line. The size line is untrusted
/// input, so beyond this the entry list grows as entries are actually read.
const MAX_PREALLOCATED_ENTRIES: usize = 1 << 20;

#[derive(Clone, Copy, PartialEq)]
enum Symmetry {
    General,
    Symmetric,
    SkewSymmetric,
}

fn parse_error(line: usize, message: impl Into<String>) -> SparseError {
    SparseError::Parse {
        line,
        message: message.into(),
    }
}

/// The numeric token for the negation of `token`: toggles a leading `-`, and turns a leading
/// `+` into `-`.
fn negate_token(token: &str) -> String {
    match token.strip_prefix('-') {
        Some(magnitude) => magnitude.to_string(),
        None => format!("-{}", token.strip_prefix('+').unwrap_or(token)),
    }
}

fn io_error(err: std::io::Error) -> SparseError {
    SparseError::Io {
        message: err.to_string(),
    }
}

impl<T: Scalar> SparseMatrix<T> {
    /// Parses a matrix from Matrix Market coordinate text; see the module docs for what is
    /// supported. Entries listed twice are summed. Returns `SparseError::Parse`, naming the
    /// line, for malformed input, including a skew-symmetric value whose negation the element
    /// type can't hold (any nonzero value of an unsigned type, or `i32::MIN`).
    ///
    /// For integer element types, duplicate entries whose sum overflows behave like integer
    /// `+`: a panic in debug builds and wrapping in release builds.
    pub fn from_matrix_market(text: &str) -> Result<SparseMatrix<T>, SparseError> {
        Self::read_matrix_market(text.as_bytes())
    }

    /// Reads a matrix in Matrix Market coordinate format from any buffered reader. Returns
    /// `SparseError::Parse` for malformed input and `SparseError::Io` if reading fails.
    pub fn read_matrix_market<R: BufRead>(reader: R) -> Result<SparseMatrix<T>, SparseError> {
        let mut lines = reader
            .lines()
            .enumerate()
            .map(|(index, line)| (index + 1, line));

        let (header_line, header) = match lines.next() {
            Some((number, line)) => (number, line.map_err(io_error)?),
            None => return Err(parse_error(1, "empty input")),
        };
        let (has_values, symmetry) = Self::parse_header(header_line, &header)?;

        // Size line: rows cols entries, after any comments and blank lines
        let (mut rows, mut cols, mut expected) = (0, 0, None);
        for (number, line) in lines.by_ref() {
            let line = line.map_err(io_error)?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('%') {
                continue;
            }
            let sizes: Vec<&str> = line.split_whitespace().collect();
            let parse = |token: &str| {
                token
                    .parse::<usize>()
                    .map_err(|_| parse_error(number, format!("invalid size '{token}'")))
            };
            if sizes.len() != 3 {
                return Err(parse_error(number, "size line must be: rows cols entries"));
            }
            (rows, cols, expected) = (parse(sizes[0])?, parse(sizes[1])?, Some(parse(sizes[2])?));
            break;
        }
        let expected = expected.ok_or_else(|| parse_error(header_line, "missing size line"))?;

        let mut triplets: Vec<(usize, usize, T)> =
            Vec::with_capacity(expected.min(MAX_PREALLOCATED_ENTRIES));
        let mut found = 0;
        let mut last_line = header_line;
        for (number, line) in lines {
            last_line = number;
            let line = line.map_err(io_error)?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('%') {
                continue;
            }
            if found == expected {
                return Err(parse_error(number, format!("more than {expected} entries")));
            }
            found += 1;

            let tokens: Vec<&str> = line.split_whitespace().collect();
            let wanted = if has_values { 3 } else { 2 };
            if tokens.len() != wanted {
                return Err(parse_error(
                    number,
                    format!("expected {wanted} fields, found {}", tokens.len()),
                ));
            }
            let index = |token: &str, bound: usize, what: &str| match token.parse::<usize>() {
                Ok(one_based) if (1..=bound).contains(&one_based) => Ok(one_based - 1),
                _ => Err(parse_error(
                    number,
                    format!("{what} index '{token}' out of range 1..={bound}"),
                )),
            };
            let (row, col) = (
                index(tokens[0], rows, "row")?,
                index(tokens[1], cols, "column")?,
            );
            let value = if has_values {
                T::from_str_radix(tokens[2], 10)
                    .map_err(|_| parse_error(number, format!("invalid value '{}'", tokens[2])))?
            } else {
                T::one()
            };

            triplets.push((row, col, value));
            if row != col {
                match symmetry {
                    Symmetry::General => {}
                    Symmetry::Symmetric => triplets.push((col, row, value)),
                    Symmetry::SkewSymmetric => {
                        // Parse the negated token rather than computing 0 - value, which
                        //   overflows (a panic in debug builds) for unsigned types and for the
                        //   most negative signed integer. Those values have no mirror entry in
                        //   the type, so they become a parse error instead.
                        let token = if has_values { tokens[2] } else { "1" };
                        let negated =
                            T::from_str_radix(&negate_token(token), 10).map_err(|_| {
                                parse_error(
                                    number,
                                    format!("value '{token}' has no negation in this element type"),
                                )
                            })?;
                        triplets.push((col, row, negated));
                    }
                }
            }
        }
        if found < expected {
            return Err(parse_error(
                last_line,
                format!("expected {expected} entries, found {found}"),
            ));
        }
        SparseMatrix::from_triplets(rows, cols, triplets)
    }

    /// Parses the `%%MatrixMarket` banner, returning whether entries carry values and the
    /// symmetry.
    fn parse_header(line_number: usize, header: &str) -> Result<(bool, Symmetry), SparseError> {
        let words: Vec<String> = header
            .split_whitespace()
            .map(|word| word.to_ascii_lowercase())
            .collect();
        if words.len() != 5 || words[0] != "%%matrixmarket" {
            return Err(parse_error(
                line_number,
                "expected '%%MatrixMarket matrix coordinate <field> <symmetry>'",
            ));
        }
        if words[1] != "matrix" || words[2] != "coordinate" {
            return Err(parse_error(
                line_number,
                format!(
                    "unsupported format '{} {}' (only matrix coordinate)",
                    words[1], words[2]
                ),
            ));
        }
        let has_values = match words[3].as_str() {
            "real" | "integer" => true,
            "pattern" => false,
            other => {
                return Err(parse_error(
                    line_number,
                    format!("unsupported field '{other}'"),
                ))
            }
        };
        let symmetry = match words[4].as_str() {
            "general" => Symmetry::General,
            "symmetric" => Symmetry::Symmetric,
            "skew-symmetric" => Symmetry::SkewSymmetric,
            other => {
                return Err(parse_error(
                    line_number,
                    format!("unsupported symmetry '{other}'"),
                ))
            }
        };
        Ok((has_values, symmetry))
    }
}

impl<T: Scalar + Display> SparseMatrix<T> {
    /// Writes the matrix in Matrix Market coordinate format (`real general`, 1-based
    /// indices, row-major order). Returns `SparseError::Io` if writing fails.
    pub fn write_matrix_market<W: Write>(&self, mut writer: W) -> Result<(), SparseError> {
        writeln!(writer, "%%MatrixMarket matrix coordinate real general").map_err(io_error)?;
        writeln!(
            writer,
            "{} {} {}",
            self.shape.0,
            self.shape.1,
            self.num_nonzero()
        )
        .map_err(io_error)?;
        for (row, col, val) in self.iter() {
            writeln!(writer, "{} {} {}", row + 1, col + 1, val).map_err(io_error)?;
        }
        Ok(())
    }

    /// The matrix as Matrix Market coordinate text; see `write_matrix_market`.
    pub fn to_matrix_market(&self) -> String {
        let mut bytes = Vec::new();
        self.write_matrix_market(&mut bytes)
            .expect("writing to a Vec cannot fail");
        String::from_utf8(bytes).expect("numbers and ASCII punctuation are valid UTF-8")
    }
}
