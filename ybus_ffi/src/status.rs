//! The status codes every fallible export returns, as `#define`s in the header. A code is never
//! renumbered once released; new failures get new codes.

use powerflow::PowerflowError;

/// A status code: `YBUS_OK` (0) on success, or the code of what failed, with a message that
/// `ybus_last_error_message` returns.
#[allow(non_camel_case_types)]
pub type ybus_status = i32;

/// Success.
pub const YBUS_OK: ybus_status = 0;
/// A failure without a code of its own, such as one added to `powerflow` after this version of
/// the library: its message says what it was.
pub const YBUS_ERROR_OTHER: ybus_status = 1;
/// A bug in ybus, caught before it could cross into C. The handles passed to the call should
/// only be freed.
pub const YBUS_ERROR_PANIC: ybus_status = 2;
/// A pointer argument was null where a value is needed.
pub const YBUS_ERROR_NULL_ARGUMENT: ybus_status = 3;
/// A pointer argument wasn't aligned for its type.
pub const YBUS_ERROR_MISALIGNED_ARGUMENT: ybus_status = 4;
/// A string argument wasn't valid UTF-8.
pub const YBUS_ERROR_INVALID_UTF8: ybus_status = 5;
/// A buffer was shorter than the result, so nothing was written to it.
pub const YBUS_ERROR_BUFFER_TOO_SMALL: ybus_status = 6;
/// A value can't be held where it was given or asked for: a negative label, a bus type outside 1
/// to 4, a status other than 0 or 1, or a cost model other than 1 or 2. Nothing was changed.
pub const YBUS_ERROR_INVALID_VALUE: ybus_status = 7;
/// An array's length doesn't match the table it's for, or a cost table's offsets are
/// inconsistent. Nothing was changed.
pub const YBUS_ERROR_LENGTH_MISMATCH: ybus_status = 8;
/// A field code isn't one of the table's, or is one of its fields of the other value type.
pub const YBUS_ERROR_INVALID_FIELD: ybus_status = 9;
/// The case has no such data: no bus names, or no costs.
pub const YBUS_ERROR_MISSING_DATA: ybus_status = 10;

// `powerflow`'s errors, from 100, one per `PowerflowError` variant in its declaration order

/// A case file couldn't be read.
pub const YBUS_ERROR_NO_CASE_FILE: ybus_status = 100;
/// A case file's text is malformed, or uses syntax the reader doesn't support.
pub const YBUS_ERROR_PARSE: ybus_status = 101;
/// A case file lacks a field every case needs.
pub const YBUS_ERROR_MISSING_FIELD: ybus_status = 102;
/// A case file isn't in MATPOWER's version 2 format, the only one supported.
pub const YBUS_ERROR_UNSUPPORTED_VERSION: ybus_status = 103;
/// A bus has a type code other than 1 to 4.
pub const YBUS_ERROR_UNKNOWN_BUS_TYPE: ybus_status = 104;
/// The system MVA base isn't a positive, finite number.
pub const YBUS_ERROR_INVALID_BASE_MVA: ybus_status = 105;
/// Two buses have the same label.
pub const YBUS_ERROR_DUPLICATE_BUS: ybus_status = 106;
/// A generator or branch refers to a bus label that no bus has.
pub const YBUS_ERROR_UNKNOWN_BUS: ybus_status = 107;
/// The case has no buses.
pub const YBUS_ERROR_NO_BUSES: ybus_status = 108;
/// An optional table's length doesn't match the table it describes.
pub const YBUS_ERROR_COUNT_MISMATCH: ybus_status = 109;
/// An in-service branch has zero reactance, which the DC power flow can't use.
pub const YBUS_ERROR_ZERO_REACTANCE: ybus_status = 110;
/// An in-service branch has zero impedance, which the AC power flow can't use.
pub const YBUS_ERROR_ZERO_IMPEDANCE: ybus_status = 111;
/// No bus can be the reference: no REF bus, and no PV bus with an in-service generator.
pub const YBUS_ERROR_NO_REFERENCE_BUS: ybus_status = 112;
/// A connected bus can't reach any reference bus through in-service branches.
pub const YBUS_ERROR_ISLAND_WITHOUT_REFERENCE: ybus_status = 113;
/// The power flow's system is singular although every bus reaches a reference bus.
pub const YBUS_ERROR_SINGULAR_NETWORK: ybus_status = 114;
/// A model was given a case whose tables' lengths differ from those of its own case.
pub const YBUS_ERROR_CASE_MISMATCH: ybus_status = 115;
/// A model was given a case whose topology differs from its own case's: a changed network
/// needs a new model.
pub const YBUS_ERROR_TOPOLOGY_MISMATCH: ybus_status = 116;
/// Enforcing reactive power limits left no bus to hold a voltage.
pub const YBUS_ERROR_REACTIVE_LIMITS_INFEASIBLE: ybus_status = 117;
/// A generator's lower reactive power limit is above its upper one, so limits can't be
/// enforced.
pub const YBUS_ERROR_INVERTED_REACTIVE_LIMITS: ybus_status = 118;
/// Enforcing reactive power limits would make a reference bus PQ in a system with several.
pub const YBUS_ERROR_REFERENCE_AT_REACTIVE_LIMIT: ybus_status = 119;
/// A row index is past the end of its table.
pub const YBUS_ERROR_ROW_OUT_OF_RANGE: ybus_status = 120;
/// An outage was asked for of a branch that's already out of service.
pub const YBUS_ERROR_BRANCH_OUT_OF_SERVICE: ybus_status = 121;
/// A PTDF slack can't be used.
pub const YBUS_ERROR_INVALID_SLACK: ybus_status = 122;
/// The AC power flow didn't meet its tolerance within its iteration limit.
pub const YBUS_ERROR_NO_CONVERGENCE: ybus_status = 123;
/// A generator cost the optimal power flow can't use.
pub const YBUS_ERROR_UNSUPPORTED_COST: ybus_status = 124;
/// A generator the AC optimal power flow can't use.
pub const YBUS_ERROR_UNSUPPORTED_GENERATOR: ybus_status = 125;
/// The optimal power flow's solver refused to start, as for an option out of range.
pub const YBUS_ERROR_SOLVER: ybus_status = 126;
/// The GPU path failed. Reserved: the C library is built without the GPU path (`powerflow`'s
/// `cuda` feature), so no call returns it yet.
pub const YBUS_ERROR_GPU: ybus_status = 127;

/// The code of a `powerflow` error.
pub(crate) fn powerflow_code(error: &PowerflowError) -> ybus_status {
    match error {
        PowerflowError::NoCaseFile { .. } => YBUS_ERROR_NO_CASE_FILE,
        PowerflowError::Parse { .. } => YBUS_ERROR_PARSE,
        PowerflowError::MissingField { .. } => YBUS_ERROR_MISSING_FIELD,
        PowerflowError::UnsupportedVersion { .. } => YBUS_ERROR_UNSUPPORTED_VERSION,
        PowerflowError::UnknownBusType { .. } => YBUS_ERROR_UNKNOWN_BUS_TYPE,
        PowerflowError::InvalidBaseMva { .. } => YBUS_ERROR_INVALID_BASE_MVA,
        PowerflowError::DuplicateBus { .. } => YBUS_ERROR_DUPLICATE_BUS,
        PowerflowError::UnknownBus { .. } => YBUS_ERROR_UNKNOWN_BUS,
        PowerflowError::NoBuses => YBUS_ERROR_NO_BUSES,
        PowerflowError::CountMismatch { .. } => YBUS_ERROR_COUNT_MISMATCH,
        PowerflowError::ZeroReactance { .. } => YBUS_ERROR_ZERO_REACTANCE,
        PowerflowError::ZeroImpedance { .. } => YBUS_ERROR_ZERO_IMPEDANCE,
        PowerflowError::NoReferenceBus => YBUS_ERROR_NO_REFERENCE_BUS,
        PowerflowError::IslandWithoutReference { .. } => YBUS_ERROR_ISLAND_WITHOUT_REFERENCE,
        PowerflowError::SingularNetwork { .. } => YBUS_ERROR_SINGULAR_NETWORK,
        PowerflowError::CaseMismatch { .. } => YBUS_ERROR_CASE_MISMATCH,
        PowerflowError::TopologyMismatch { .. } => YBUS_ERROR_TOPOLOGY_MISMATCH,
        PowerflowError::ReactiveLimitsInfeasible { .. } => YBUS_ERROR_REACTIVE_LIMITS_INFEASIBLE,
        PowerflowError::InvertedReactiveLimits { .. } => YBUS_ERROR_INVERTED_REACTIVE_LIMITS,
        PowerflowError::ReferenceAtReactiveLimit { .. } => YBUS_ERROR_REFERENCE_AT_REACTIVE_LIMIT,
        PowerflowError::RowOutOfRange { .. } => YBUS_ERROR_ROW_OUT_OF_RANGE,
        PowerflowError::BranchOutOfService { .. } => YBUS_ERROR_BRANCH_OUT_OF_SERVICE,
        PowerflowError::InvalidSlack { .. } => YBUS_ERROR_INVALID_SLACK,
        PowerflowError::NoConvergence { .. } => YBUS_ERROR_NO_CONVERGENCE,
        PowerflowError::UnsupportedCost { .. } => YBUS_ERROR_UNSUPPORTED_COST,
        PowerflowError::UnsupportedGenerator { .. } => YBUS_ERROR_UNSUPPORTED_GENERATOR,
        PowerflowError::Solver { .. } => YBUS_ERROR_SOLVER,
        // `PowerflowError` is #[non_exhaustive], so a variant added later is OTHER until it gets
        //   a code. So is `Gpu`, which only exists with powerflow's `cuda` feature, which this
        //   crate can't name; a workspace build with every feature unifies it in, but the C
        //   library is built without it.
        _ => YBUS_ERROR_OTHER,
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::path::PathBuf;

    use powerflow::opf::InteriorPointError;

    use super::*;

    #[test]
    fn every_code_is_distinct() {
        let codes = [
            YBUS_OK,
            YBUS_ERROR_OTHER,
            YBUS_ERROR_PANIC,
            YBUS_ERROR_NULL_ARGUMENT,
            YBUS_ERROR_MISALIGNED_ARGUMENT,
            YBUS_ERROR_INVALID_UTF8,
            YBUS_ERROR_BUFFER_TOO_SMALL,
            YBUS_ERROR_INVALID_VALUE,
            YBUS_ERROR_LENGTH_MISMATCH,
            YBUS_ERROR_INVALID_FIELD,
            YBUS_ERROR_MISSING_DATA,
            YBUS_ERROR_GPU,
        ];
        let distinct: HashSet<_> = codes.iter().collect();
        assert_eq!(distinct.len(), codes.len());
        assert!(codes
            .iter()
            .all(|&code| !(100..YBUS_ERROR_GPU).contains(&code)));
    }

    // Every current variant has its own code, from 100 in declaration order, so none falls to
    //   OTHER. A variant added to powerflow later does, until it's listed here and given a code.
    #[test]
    fn each_powerflow_error_has_its_own_code() {
        let errors = [
            PowerflowError::NoCaseFile {
                bad_path: PathBuf::from("missing.m"),
                reason: "not found".to_string(),
            },
            PowerflowError::Parse {
                line: 1,
                message: String::new(),
            },
            PowerflowError::MissingField { field: "bus" },
            PowerflowError::UnsupportedVersion {
                found: "1".to_string(),
            },
            PowerflowError::UnknownBusType { line: 1, code: 5.0 },
            PowerflowError::InvalidBaseMva { base_mva: 0.0 },
            PowerflowError::DuplicateBus {
                label: 1,
                first_row: 1,
                second_row: 2,
            },
            PowerflowError::UnknownBus {
                table: "gen",
                row: 1,
                label: 9,
            },
            PowerflowError::NoBuses,
            PowerflowError::CountMismatch {
                field: "bus_name",
                expected: "3".to_string(),
                found: 2,
            },
            PowerflowError::ZeroReactance { row: 1 },
            PowerflowError::ZeroImpedance { row: 1 },
            PowerflowError::NoReferenceBus,
            PowerflowError::IslandWithoutReference { bus_label: 1 },
            PowerflowError::SingularNetwork { bus_label: 1 },
            PowerflowError::CaseMismatch {
                expected: (3, 1, 3),
                found: (3, 1, 2),
            },
            PowerflowError::TopologyMismatch {
                table: "bus",
                row: 1,
                field: "bus_type",
            },
            PowerflowError::ReactiveLimitsInfeasible { generators: 1 },
            PowerflowError::InvertedReactiveLimits { row: 1 },
            PowerflowError::ReferenceAtReactiveLimit { bus_label: 1 },
            PowerflowError::RowOutOfRange {
                table: "branch",
                row: 9,
                len: 3,
            },
            PowerflowError::BranchOutOfService { row: 1 },
            PowerflowError::InvalidSlack {
                reason: "no weight",
            },
            PowerflowError::NoConvergence {
                iterations: 10,
                mismatch_pu: 1.0,
            },
            PowerflowError::UnsupportedCost {
                row: 1,
                reason: "cubic",
            },
            PowerflowError::UnsupportedGenerator {
                row: 1,
                reason: "capability curve",
            },
            PowerflowError::Solver {
                error: InteriorPointError::UnboundedLinearProgram,
            },
        ];
        let codes: Vec<ybus_status> = errors.iter().map(powerflow_code).collect();
        assert_eq!(codes, (100..YBUS_ERROR_GPU).collect::<Vec<_>>());
    }
}
