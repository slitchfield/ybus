//! The error type for reading cases and running power flows.

use std::{fmt, path::PathBuf};

/// Errors returned by `powerflow` operations.
#[derive(Debug, PartialEq)]
#[non_exhaustive]
pub enum PowerflowError {
    /// A case file couldn't be read.
    NoCaseFile {
        /// The path that was tried.
        bad_path: PathBuf,
        /// The underlying I/O error's message.
        reason: String,
    },
    /// Case file text is malformed, or uses syntax the reader doesn't support.
    Parse {
        /// The 1-based line the problem was found on.
        line: usize,
        /// What is wrong.
        message: String,
    },
    /// The case file lacks a field every case needs (`version`, `baseMVA`, `bus`, `gen` or
    /// `branch`).
    MissingField {
        /// The field's name in the file, without `mpc.`.
        field: &'static str,
    },
    /// The case file isn't in MATPOWER's version 2 format, the only one supported.
    UnsupportedVersion {
        /// The version the file states.
        found: String,
    },
    /// A bus row has a type code other than 1 (PQ), 2 (PV), 3 (reference) or 4 (isolated).
    UnknownBusType {
        /// The 1-based line of the bus row.
        line: usize,
        /// The code found.
        code: f64,
    },
    /// The system MVA base isn't a positive, finite number.
    InvalidBaseMva {
        /// The value found.
        base_mva: f64,
    },
    /// Two buses have the same label.
    DuplicateBus {
        /// The repeated label.
        label: usize,
        /// The 1-based row in the bus table of its first use.
        first_row: usize,
        /// The 1-based row in the bus table of its second use.
        second_row: usize,
    },
    /// A generator or branch refers to a bus label that no bus has.
    UnknownBus {
        /// The table the reference is in (`gen` or `branch`).
        table: &'static str,
        /// The 1-based row in that table.
        row: usize,
        /// The label referred to.
        label: usize,
    },
    /// The case has no buses.
    NoBuses,
    /// An optional table's length doesn't match the table it describes.
    CountMismatch {
        /// The field's name in the file (`bus_name` or `gencost`).
        field: &'static str,
        /// The length or lengths allowed.
        expected: String,
        /// The length found.
        found: usize,
    },
    /// An in-service branch has zero reactance, so its DC susceptance `1 / x` is infinite.
    ZeroReactance {
        /// The 1-based row in the branch table.
        row: usize,
    },
    /// An in-service branch has zero impedance (`r = x = 0`), so its AC series admittance
    /// `1 / (r + jx)` is infinite.
    ZeroImpedance {
        /// The 1-based row in the branch table.
        row: usize,
    },
    /// No bus can be the DC reference: there's no REF bus, and no PV bus with an in-service
    /// generator to take its place.
    NoReferenceBus,
    /// A connected bus can't reach any reference bus through in-service branches, so its
    /// island has no fixed angle and the DC solve is singular.
    IslandWithoutReference {
        /// The label of the first such bus, in case order.
        bus_label: usize,
    },
    /// The DC system is singular although every bus reaches a reference bus, which takes
    /// branch reactances that cancel out (for example a series capacitor in parallel with an
    /// equal and opposite line).
    SingularNetwork {
        /// The label of the bus whose angle couldn't be determined: the elimination step that
        /// failed, or the first angle that came out non-finite or over 10⁵ radians.
        bus_label: usize,
    },
    /// A DC model was asked to solve a case whose tables don't match the case it was built
    /// from.
    CaseMismatch {
        /// The model's case's table lengths: buses, generators, branches.
        expected: (usize, usize, usize),
        /// The given case's table lengths.
        found: (usize, usize, usize),
    },
}

impl fmt::Display for PowerflowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PowerflowError::NoCaseFile { bad_path, reason } => {
                write!(f, "can't read case file {}: {reason}", bad_path.display())
            }
            PowerflowError::Parse { line, message } => write!(f, "line {line}: {message}"),
            PowerflowError::MissingField { field } => {
                write!(f, "case file has no mpc.{field}")
            }
            PowerflowError::UnsupportedVersion { found } => write!(
                f,
                "case file is version {found:?}; only MATPOWER case format version 2 is supported"
            ),
            PowerflowError::UnknownBusType { line, code } => write!(
                f,
                "line {line}: bus type {code} is not 1 (PQ), 2 (PV), 3 (reference) or 4 (isolated)"
            ),
            PowerflowError::InvalidBaseMva { base_mva } => {
                write!(f, "baseMVA must be positive and finite, not {base_mva}")
            }
            PowerflowError::DuplicateBus {
                label,
                first_row,
                second_row,
            } => write!(
                f,
                "bus {label} appears twice in the bus table, at rows {first_row} and {second_row}"
            ),
            PowerflowError::UnknownBus { table, row, label } => {
                write!(
                    f,
                    "{table} row {row} refers to bus {label}, which doesn't exist"
                )
            }
            PowerflowError::NoBuses => write!(f, "case has no buses"),
            PowerflowError::CountMismatch {
                field,
                expected,
                found,
            } => write!(f, "{field} has {found} entries, but should have {expected}"),
            PowerflowError::ZeroReactance { row } => {
                write!(f, "branch row {row} is in service with zero reactance")
            }
            PowerflowError::ZeroImpedance { row } => {
                write!(f, "branch row {row} is in service with zero impedance")
            }
            PowerflowError::NoReferenceBus => write!(
                f,
                "no reference bus: no REF bus, and no PV bus with an in-service generator"
            ),
            PowerflowError::IslandWithoutReference { bus_label } => write!(
                f,
                "bus {bus_label} is in an island with no reference bus, so its angle is undefined"
            ),
            PowerflowError::SingularNetwork { bus_label } => write!(
                f,
                "the DC system is singular: the angle at bus {bus_label} can't be determined"
            ),
            PowerflowError::CaseMismatch { expected, found } => write!(
                f,
                "the DC model was built for a case with (buses, generators, branches) = \
                 {expected:?}, not {found:?}"
            ),
        }
    }
}

impl std::error::Error for PowerflowError {}
