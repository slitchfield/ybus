//! A case's columns, read whole by field code: one function per table and value type, with
//! `double` columns and `int64_t` columns (labels, areas, zones, bus types and statuses).
//!
//! A field code is the field's 1-based column number in MATPOWER's matrix for its table, so
//! the codes are fixed by the format itself and never renumbered. Codes repeat across tables,
//! which have their own functions.

use powerflow::case::{Branch, Bus, BusType, Gen};
use powerflow::Case;

use super::YbusCase;
use crate::boundary::{arg_mut, arg_ref, arg_slice, arg_slice_mut, guard, Failure};
use crate::status::ybus_status;

/// A column's field code: its 1-based column number in MATPOWER's matrix for the table.
#[allow(non_camel_case_types)]
pub type ybus_field = i32;

/// The bus number: a label other rows refer to, not a position (`int64_t`).
pub const YBUS_BUS_LABEL: ybus_field = 1;
/// The bus type, MATPOWER's code: 1 PQ, 2 PV, 3 reference, 4 isolated (`int64_t`).
pub const YBUS_BUS_TYPE: ybus_field = 2;
/// Real power demand (MW).
pub const YBUS_BUS_PD_MW: ybus_field = 3;
/// Reactive power demand (MVAr).
pub const YBUS_BUS_QD_MVAR: ybus_field = 4;
/// Shunt conductance, as the MW it draws at 1.0 p.u. voltage.
pub const YBUS_BUS_GS_MW: ybus_field = 5;
/// Shunt susceptance, as the MVAr it injects at 1.0 p.u. voltage.
pub const YBUS_BUS_BS_MVAR: ybus_field = 6;
/// The area number (`int64_t`).
pub const YBUS_BUS_AREA: ybus_field = 7;
/// Voltage magnitude (p.u.).
pub const YBUS_BUS_VM_PU: ybus_field = 8;
/// Voltage angle (degrees).
pub const YBUS_BUS_VA_DEG: ybus_field = 9;
/// Base voltage (kV).
pub const YBUS_BUS_BASE_KV: ybus_field = 10;
/// The loss zone (`int64_t`).
pub const YBUS_BUS_ZONE: ybus_field = 11;
/// Maximum voltage magnitude (p.u.).
pub const YBUS_BUS_VMAX_PU: ybus_field = 12;
/// Minimum voltage magnitude (p.u.).
pub const YBUS_BUS_VMIN_PU: ybus_field = 13;

/// The label of the bus the generator is connected to (`int64_t`).
pub const YBUS_GEN_BUS_LABEL: ybus_field = 1;
/// Real power output (MW).
pub const YBUS_GEN_PG_MW: ybus_field = 2;
/// Reactive power output (MVAr).
pub const YBUS_GEN_QG_MVAR: ybus_field = 3;
/// Maximum reactive power output (MVAr).
pub const YBUS_GEN_QMAX_MVAR: ybus_field = 4;
/// Minimum reactive power output (MVAr).
pub const YBUS_GEN_QMIN_MVAR: ybus_field = 5;
/// Voltage magnitude setpoint (p.u.).
pub const YBUS_GEN_VG_PU: ybus_field = 6;
/// The machine's MVA base.
pub const YBUS_GEN_MBASE_MVA: ybus_field = 7;
/// 1 if the generator is in service, 0 if not (`int64_t`).
pub const YBUS_GEN_STATUS: ybus_field = 8;
/// Maximum real power output (MW).
pub const YBUS_GEN_PMAX_MW: ybus_field = 9;
/// Minimum real power output (MW).
pub const YBUS_GEN_PMIN_MW: ybus_field = 10;
/// Lower real power output of the PQ capability curve (MW).
pub const YBUS_GEN_PC1_MW: ybus_field = 11;
/// Upper real power output of the PQ capability curve (MW).
pub const YBUS_GEN_PC2_MW: ybus_field = 12;
/// Minimum reactive power output at `PC1` (MVAr).
pub const YBUS_GEN_QC1MIN_MVAR: ybus_field = 13;
/// Maximum reactive power output at `PC1` (MVAr).
pub const YBUS_GEN_QC1MAX_MVAR: ybus_field = 14;
/// Minimum reactive power output at `PC2` (MVAr).
pub const YBUS_GEN_QC2MIN_MVAR: ybus_field = 15;
/// Maximum reactive power output at `PC2` (MVAr).
pub const YBUS_GEN_QC2MAX_MVAR: ybus_field = 16;
/// Ramp rate for load following and AGC (MW/min).
pub const YBUS_GEN_RAMP_AGC_MW_PER_MIN: ybus_field = 17;
/// Ramp rate for 10 minute reserves (MW).
pub const YBUS_GEN_RAMP_10_MW: ybus_field = 18;
/// Ramp rate for 30 minute reserves (MW).
pub const YBUS_GEN_RAMP_30_MW: ybus_field = 19;
/// Ramp rate for reactive power, 2 second timescale (MVAr/min).
pub const YBUS_GEN_RAMP_Q_MVAR_PER_MIN: ybus_field = 20;
/// Area participation factor.
pub const YBUS_GEN_APF: ybus_field = 21;

/// The label of the "from" bus, where a transformer's tap is (`int64_t`).
pub const YBUS_BRANCH_FROM_LABEL: ybus_field = 1;
/// The label of the "to" bus (`int64_t`).
pub const YBUS_BRANCH_TO_LABEL: ybus_field = 2;
/// Resistance (p.u.).
pub const YBUS_BRANCH_R: ybus_field = 3;
/// Reactance (p.u.).
pub const YBUS_BRANCH_X: ybus_field = 4;
/// Total line charging susceptance (p.u.).
pub const YBUS_BRANCH_B: ybus_field = 5;
/// Long-term MVA rating; 0 means unlimited.
pub const YBUS_BRANCH_RATE_A_MVA: ybus_field = 6;
/// Short-term MVA rating; 0 means unlimited.
pub const YBUS_BRANCH_RATE_B_MVA: ybus_field = 7;
/// Emergency MVA rating; 0 means unlimited.
pub const YBUS_BRANCH_RATE_C_MVA: ybus_field = 8;
/// Transformer off-nominal turns ratio, where 0 means a line (equivalent to 1).
pub const YBUS_BRANCH_TAP_RATIO: ybus_field = 9;
/// Transformer phase shift angle (degrees); positive is a delay.
pub const YBUS_BRANCH_SHIFT_DEG: ybus_field = 10;
/// 1 if the branch is in service, 0 if not (`int64_t`).
pub const YBUS_BRANCH_STATUS: ybus_field = 11;
/// Minimum angle difference, angle(Vf) − angle(Vt) (degrees).
pub const YBUS_BRANCH_ANGMIN_DEG: ybus_field = 12;
/// Maximum angle difference, angle(Vf) − angle(Vt) (degrees).
pub const YBUS_BRANCH_ANGMAX_DEG: ybus_field = 13;

/// How to read and write one `double` field of a row.
pub(crate) struct Float<R> {
    pub(crate) get: fn(&R) -> f64,
    pub(crate) set: fn(&mut R, f64),
}

/// How to read and write one `int64_t` field of a row.
pub(crate) struct Int<R> {
    /// The row's value, or `None` if it doesn't fit an `int64_t`
    pub(crate) get: fn(&R) -> Option<i64>,
    /// Whether the field can hold a value, and if not, why
    pub(crate) check: fn(i64) -> Result<(), &'static str>,
    /// Sets the row's value, which `check` has accepted
    pub(crate) set: fn(&mut R, i64),
}

/// One of a case's tables: its rows, and its fields by code.
pub(crate) struct Table<R: 'static> {
    pub(crate) name: &'static str,
    pub(crate) rows: fn(&Case) -> &[R],
    pub(crate) rows_mut: fn(&mut Case) -> &mut [R],
    pub(crate) float: fn(ybus_field) -> Option<Float<R>>,
    pub(crate) int: fn(ybus_field) -> Option<Int<R>>,
}

/// A `double` field, read and written as is
macro_rules! float {
    ($field:ident) => {
        Some(Float {
            get: |row| row.$field,
            set: |row, value| row.$field = value,
        })
    };
}

/// A label, area or zone: any non-negative integer
macro_rules! count {
    ($field:ident) => {
        Some(Int {
            get: |row| i64::try_from(row.$field).ok(),
            check: check_count,
            set: |row, value| row.$field = usize::try_from(value).expect("checked"),
        })
    };
}

/// A status: 1 in service, 0 not
macro_rules! status {
    ($field:ident) => {
        Some(Int {
            get: |row| Some(i64::from(row.$field)),
            check: check_status,
            set: |row, value| row.$field = value == 1,
        })
    };
}

fn check_count(value: i64) -> Result<(), &'static str> {
    usize::try_from(value)
        .map(|_| ())
        .map_err(|_| "isn't a non-negative integer that fits a size_t")
}

fn check_status(value: i64) -> Result<(), &'static str> {
    match value {
        0 | 1 => Ok(()),
        _ => Err("isn't a status, 0 or 1"),
    }
}

/// MATPOWER's code for a bus type, the inverse of `BusType::from_code`
fn bus_type_code(bus_type: BusType) -> i64 {
    match bus_type {
        BusType::PQ => 1,
        BusType::PV => 2,
        BusType::Ref => 3,
        BusType::Isolated => 4,
    }
}

fn bus_float(field: ybus_field) -> Option<Float<Bus>> {
    match field {
        YBUS_BUS_PD_MW => float!(pd_mw),
        YBUS_BUS_QD_MVAR => float!(qd_mvar),
        YBUS_BUS_GS_MW => float!(gs_mw),
        YBUS_BUS_BS_MVAR => float!(bs_mvar),
        YBUS_BUS_VM_PU => float!(vm_pu),
        YBUS_BUS_VA_DEG => float!(va_deg),
        YBUS_BUS_BASE_KV => float!(base_kv),
        YBUS_BUS_VMAX_PU => float!(vmax_pu),
        YBUS_BUS_VMIN_PU => float!(vmin_pu),
        _ => None,
    }
}

fn bus_int(field: ybus_field) -> Option<Int<Bus>> {
    match field {
        YBUS_BUS_LABEL => count!(bus_label),
        YBUS_BUS_TYPE => Some(Int {
            get: |bus| Some(bus_type_code(bus.bus_type)),
            check: |code| match code {
                1..=4 => Ok(()),
                _ => Err("isn't a bus type, 1 to 4"),
            },
            set: |bus, code| {
                bus.bus_type = BusType::from_code(code as f64).expect("checked");
            },
        }),
        YBUS_BUS_AREA => count!(area),
        YBUS_BUS_ZONE => count!(zone),
        _ => None,
    }
}

fn gen_float(field: ybus_field) -> Option<Float<Gen>> {
    match field {
        YBUS_GEN_PG_MW => float!(pg_mw),
        YBUS_GEN_QG_MVAR => float!(qg_mvar),
        YBUS_GEN_QMAX_MVAR => float!(qmax_mvar),
        YBUS_GEN_QMIN_MVAR => float!(qmin_mvar),
        YBUS_GEN_VG_PU => float!(vg_pu),
        YBUS_GEN_MBASE_MVA => float!(mbase_mva),
        YBUS_GEN_PMAX_MW => float!(pmax_mw),
        YBUS_GEN_PMIN_MW => float!(pmin_mw),
        YBUS_GEN_PC1_MW => float!(pc1_mw),
        YBUS_GEN_PC2_MW => float!(pc2_mw),
        YBUS_GEN_QC1MIN_MVAR => float!(qc1min_mvar),
        YBUS_GEN_QC1MAX_MVAR => float!(qc1max_mvar),
        YBUS_GEN_QC2MIN_MVAR => float!(qc2min_mvar),
        YBUS_GEN_QC2MAX_MVAR => float!(qc2max_mvar),
        YBUS_GEN_RAMP_AGC_MW_PER_MIN => float!(ramp_agc_mw_per_min),
        YBUS_GEN_RAMP_10_MW => float!(ramp_10_mw),
        YBUS_GEN_RAMP_30_MW => float!(ramp_30_mw),
        YBUS_GEN_RAMP_Q_MVAR_PER_MIN => float!(ramp_q_mvar_per_min),
        YBUS_GEN_APF => float!(apf),
        _ => None,
    }
}

fn gen_int(field: ybus_field) -> Option<Int<Gen>> {
    match field {
        YBUS_GEN_BUS_LABEL => count!(bus_label),
        YBUS_GEN_STATUS => status!(status),
        _ => None,
    }
}

fn branch_float(field: ybus_field) -> Option<Float<Branch>> {
    match field {
        YBUS_BRANCH_R => float!(r),
        YBUS_BRANCH_X => float!(x),
        YBUS_BRANCH_B => float!(b),
        YBUS_BRANCH_RATE_A_MVA => float!(rate_a_mva),
        YBUS_BRANCH_RATE_B_MVA => float!(rate_b_mva),
        YBUS_BRANCH_RATE_C_MVA => float!(rate_c_mva),
        YBUS_BRANCH_TAP_RATIO => float!(tap_ratio),
        YBUS_BRANCH_SHIFT_DEG => float!(shift_deg),
        YBUS_BRANCH_ANGMIN_DEG => float!(angmin_deg),
        YBUS_BRANCH_ANGMAX_DEG => float!(angmax_deg),
        _ => None,
    }
}

fn branch_int(field: ybus_field) -> Option<Int<Branch>> {
    match field {
        YBUS_BRANCH_FROM_LABEL => count!(from_label),
        YBUS_BRANCH_TO_LABEL => count!(to_label),
        YBUS_BRANCH_STATUS => status!(status),
        _ => None,
    }
}

pub(crate) fn bus_table() -> Table<Bus> {
    Table {
        name: "bus",
        rows: |case| &case.bus,
        rows_mut: |case| &mut case.bus,
        float: bus_float,
        int: bus_int,
    }
}

pub(crate) fn gen_table() -> Table<Gen> {
    Table {
        name: "gen",
        rows: |case| &case.generator,
        rows_mut: |case| &mut case.generator,
        float: gen_float,
        int: gen_int,
    }
}

pub(crate) fn branch_table() -> Table<Branch> {
    Table {
        name: "branch",
        rows: |case| &case.branch,
        rows_mut: |case| &mut case.branch,
        float: branch_float,
        int: branch_int,
    }
}

/// `InvalidField` unless the field is one of the table's of the given kind
fn field<C>(
    column: Option<C>,
    table: &'static str,
    field: ybus_field,
    kind: &'static str,
) -> Result<C, Failure> {
    column.ok_or(Failure::InvalidField { table, field, kind })
}

/// `BufferTooSmall` unless a buffer of `given` entries holds `needed`
pub(crate) fn check_room(
    argument: &'static str,
    needed: usize,
    given: usize,
) -> Result<(), Failure> {
    if given < needed {
        Err(Failure::BufferTooSmall {
            argument,
            needed,
            given,
        })
    } else {
        Ok(())
    }
}

/// Writes a `double` column into `out`, one entry per row.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `double`s.
unsafe fn get_floats<R>(
    table: &Table<R>,
    case: *const YbusCase,
    code: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: case and out are as the caller promised
        let (case, out) = unsafe { (arg_ref(case, "case")?, arg_slice_mut(out, len, "out")?) };
        let column = field((table.float)(code), table.name, code, "double")?;
        let rows = (table.rows)(&case.case);
        check_room("out", rows.len(), len)?;
        for (out, row) in out.iter_mut().zip(rows) {
            *out = (column.get)(row);
        }
        Ok(())
    })
}

/// Writes an `int64_t` column into `out`, one entry per row, or nothing if a value doesn't fit.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `int64_t`s.
unsafe fn get_ints<R>(
    table: &Table<R>,
    case: *const YbusCase,
    code: ybus_field,
    out: *mut i64,
    len: usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: case and out are as the caller promised
        let (case, out) = unsafe { (arg_ref(case, "case")?, arg_slice_mut(out, len, "out")?) };
        let column = field((table.int)(code), table.name, code, "int64_t")?;
        let rows = (table.rows)(&case.case);
        check_room("out", rows.len(), len)?;
        let values = rows
            .iter()
            .enumerate()
            .map(|(row, values)| {
                (column.get)(values).ok_or_else(|| {
                    Failure::InvalidValue(format!(
                        "{} row {row}'s field {code} doesn't fit an int64_t",
                        table.name
                    ))
                })
            })
            .collect::<Result<Vec<i64>, Failure>>()?;
        out[..values.len()].copy_from_slice(&values);
        Ok(())
    })
}

/// `LengthMismatch` unless an array of `given` entries is one per row of `rows`
pub(crate) fn check_length(
    argument: &'static str,
    rows: usize,
    given: usize,
) -> Result<(), Failure> {
    if given == rows {
        Ok(())
    } else {
        Err(Failure::LengthMismatch(format!(
            "`{argument}` has {given} entries, but its table has {rows} rows"
        )))
    }
}

/// Sets a `double` column from `values`, one entry per row.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `double`s.
unsafe fn set_floats<R>(
    table: &Table<R>,
    case: *mut YbusCase,
    code: ybus_field,
    values: *const f64,
    len: usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: case and values are as the caller promised
        let (case, values) = unsafe { (arg_mut(case, "case")?, arg_slice(values, len, "values")?) };
        let column = field((table.float)(code), table.name, code, "double")?;
        let rows = (table.rows_mut)(&mut case.case);
        check_length("values", rows.len(), len)?;
        for (row, &value) in rows.iter_mut().zip(values) {
            (column.set)(row, value);
        }
        Ok(())
    })
}

/// Sets an `int64_t` column from `values`, one entry per row, after checking every value, so a
/// column with one the field can't hold is left unchanged.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `int64_t`s.
unsafe fn set_ints<R>(
    table: &Table<R>,
    case: *mut YbusCase,
    code: ybus_field,
    values: *const i64,
    len: usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: case and values are as the caller promised
        let (case, values) = unsafe { (arg_mut(case, "case")?, arg_slice(values, len, "values")?) };
        let column = field((table.int)(code), table.name, code, "int64_t")?;
        let rows = (table.rows_mut)(&mut case.case);
        check_length("values", rows.len(), len)?;
        for (row, &value) in values.iter().enumerate() {
            (column.check)(value).map_err(|reason| {
                Failure::InvalidValue(format!(
                    "{} row {row}'s value {value} for field {code} {reason}",
                    table.name
                ))
            })?;
        }
        for (row, &value) in rows.iter_mut().zip(values) {
            (column.set)(row, value);
        }
        Ok(())
    })
}

/// Writes the bus table's `double` column `field` into `out`, which has room for `len`
/// entries: one per bus, in case order. Returns `YBUS_ERROR_INVALID_FIELD` if `field` isn't
/// one of the bus table's `double` fields, and `YBUS_ERROR_BUFFER_TOO_SMALL`, writing nothing,
/// if `len` is less than the number of buses.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_get_bus_f64(
    case: *const YbusCase,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_floats(&bus_table(), case, field, out, len) }
}

/// [`ybus_case_get_bus_f64`] for the bus table's `int64_t` fields: the label, type, area and
/// zone.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `int64_t`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_get_bus_i64(
    case: *const YbusCase,
    field: ybus_field,
    out: *mut i64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_ints(&bus_table(), case, field, out, len) }
}

/// [`ybus_case_get_bus_f64`] for the generator table: one entry per generator.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_get_gen_f64(
    case: *const YbusCase,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_floats(&gen_table(), case, field, out, len) }
}

/// [`ybus_case_get_bus_f64`] for the generator table's `int64_t` fields: the bus label and
/// status.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `int64_t`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_get_gen_i64(
    case: *const YbusCase,
    field: ybus_field,
    out: *mut i64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_ints(&gen_table(), case, field, out, len) }
}

/// [`ybus_case_get_bus_f64`] for the branch table: one entry per branch.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_get_branch_f64(
    case: *const YbusCase,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_floats(&branch_table(), case, field, out, len) }
}

/// [`ybus_case_get_bus_f64`] for the branch table's `int64_t` fields: the from and to labels
/// and the status.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or `len` writable `int64_t`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_get_branch_i64(
    case: *const YbusCase,
    field: ybus_field,
    out: *mut i64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_ints(&branch_table(), case, field, out, len) }
}

/// Sets the bus table's `double` column `field` from `values`, which has `len` entries: one
/// per bus, in case order. Nothing is checked beyond the length: data a solve can't use is
/// reported when the case is solved. Returns `YBUS_ERROR_INVALID_FIELD` if `field` isn't one
/// of the bus table's `double` fields, and `YBUS_ERROR_LENGTH_MISMATCH`, changing nothing, if
/// `len` isn't the number of buses.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_set_bus_f64(
    case: *mut YbusCase,
    field: ybus_field,
    values: *const f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_floats(&bus_table(), case, field, values, len) }
}

/// [`ybus_case_set_bus_f64`] for the bus table's `int64_t` fields: the label, type, area and
/// zone. Every value is checked first, and `YBUS_ERROR_INVALID_VALUE`, changing nothing, is
/// returned for a negative label, area or zone, or a type outside 1 to 4.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `int64_t`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_set_bus_i64(
    case: *mut YbusCase,
    field: ybus_field,
    values: *const i64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_ints(&bus_table(), case, field, values, len) }
}

/// [`ybus_case_set_bus_f64`] for the generator table: one entry per generator.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_set_gen_f64(
    case: *mut YbusCase,
    field: ybus_field,
    values: *const f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_floats(&gen_table(), case, field, values, len) }
}

/// [`ybus_case_set_bus_i64`] for the generator table's `int64_t` fields: the bus label and the
/// status, which must be 0 or 1.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `int64_t`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_set_gen_i64(
    case: *mut YbusCase,
    field: ybus_field,
    values: *const i64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_ints(&gen_table(), case, field, values, len) }
}

/// [`ybus_case_set_bus_f64`] for the branch table: one entry per branch.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_set_branch_f64(
    case: *mut YbusCase,
    field: ybus_field,
    values: *const f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_floats(&branch_table(), case, field, values, len) }
}

/// [`ybus_case_set_bus_i64`] for the branch table's `int64_t` fields: the from and to labels
/// and the status, which must be 0 or 1.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using, and `values` null or `len`
/// readable `int64_t`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_set_branch_i64(
    case: *mut YbusCase,
    field: ybus_field,
    values: *const i64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_ints(&branch_table(), case, field, values, len) }
}

/// Writes bus `row`'s name, a NUL-terminated string, into `buffer`, which has room for `len`
/// bytes, and the bytes it needs, counting the NUL, through `needed` unless it's null.
///
/// A null `buffer` with `len` 0 only asks the size. Otherwise, returns
/// `YBUS_ERROR_BUFFER_TOO_SMALL`, writing nothing into `buffer`, if the name doesn't fit.
/// Returns `YBUS_ERROR_MISSING_DATA` if the case has no bus names, `YBUS_ERROR_ROW_OUT_OF_RANGE`
/// if `row` is past the last bus, and `YBUS_ERROR_INVALID_VALUE` if the name holds a NUL.
///
/// # Safety
///
/// `case` must be null or a live case, `buffer` null or `len` writable bytes, and `needed`
/// null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_bus_name(
    case: *const YbusCase,
    row: usize,
    buffer: *mut std::ffi::c_char,
    len: usize,
    needed: *mut usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: case is as the caller promised
        let case = unsafe { arg_ref(case, "case") }?;
        let names = case
            .case
            .busnames
            .as_ref()
            .ok_or(Failure::MissingData("bus names"))?;
        let name = names
            .get(row)
            .ok_or(powerflow::PowerflowError::RowOutOfRange {
                table: "bus",
                row,
                len: names.len(),
            })?;
        if name.contains('\0') {
            return Err(Failure::InvalidValue(format!(
                "bus row {row}'s name holds a NUL"
            )));
        }
        let size = name.len() + 1;
        if !needed.is_null() {
            // SAFETY: needed is as the caller promised
            *unsafe { arg_mut(needed, "needed") }? = size;
        }
        if buffer.is_null() && len == 0 {
            return Ok(());
        }
        // SAFETY: buffer is as the caller promised
        let buffer = unsafe { arg_slice_mut(buffer.cast::<u8>(), len, "buffer") }?;
        check_room("buffer", size, len)?;
        buffer[..name.len()].copy_from_slice(name.as_bytes());
        buffer[name.len()] = 0;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use super::*;
    use crate::boundary::last_error;
    use crate::case::tests::{three_bus, THREE_BUS};
    use crate::case::ybus_case_free;
    use crate::status::{
        YBUS_ERROR_BUFFER_TOO_SMALL, YBUS_ERROR_INVALID_FIELD, YBUS_ERROR_INVALID_VALUE,
        YBUS_ERROR_LENGTH_MISMATCH, YBUS_ERROR_MISSING_DATA, YBUS_ERROR_ROW_OUT_OF_RANGE, YBUS_OK,
    };

    /// The rows of `mpc.<name>` in THREE_BUS's text, as numbers: the columns the codes name
    fn matrix(name: &str) -> Vec<Vec<f64>> {
        let start = format!("mpc.{name} = [");
        THREE_BUS
            .lines()
            .skip_while(|line| !line.starts_with(&start))
            .skip(1)
            .take_while(|line| !line.starts_with("];"))
            .map(|line| {
                line.trim_end_matches(';')
                    .split_whitespace()
                    .map(|value| value.parse().unwrap())
                    .collect()
            })
            .collect()
    }

    /// Every code from 1 to the table's width reads, through exactly one of the two getters,
    /// the column of that number in the case file: so each code is MATPOWER's column number
    fn check_table(
        name: &str,
        width: ybus_field,
        get_f64: unsafe extern "C" fn(*const YbusCase, ybus_field, *mut f64, usize) -> ybus_status,
        get_i64: unsafe extern "C" fn(*const YbusCase, ybus_field, *mut i64, usize) -> ybus_status,
    ) {
        let case = three_bus();
        let rows = matrix(name);
        for code in 1..=width {
            let expected: Vec<f64> = rows.iter().map(|row| row[code as usize - 1]).collect();
            let (mut floats, mut ints) = (vec![0.0; rows.len()], vec![0; rows.len()]);
            // SAFETY: case is live, and each buffer has room for rows.len() entries
            let (float, int) = unsafe {
                (
                    get_f64(case, code, floats.as_mut_ptr(), floats.len()),
                    get_i64(case, code, ints.as_mut_ptr(), ints.len()),
                )
            };
            match (float, int) {
                (YBUS_OK, YBUS_ERROR_INVALID_FIELD) => {
                    assert_eq!(floats, expected, "{name} {code}")
                }
                (YBUS_ERROR_INVALID_FIELD, YBUS_OK) => {
                    let ints: Vec<f64> = ints.iter().map(|&v| v as f64).collect();
                    assert_eq!(ints, expected, "{name} {code}");
                }
                other => panic!("{name} field {code}: {other:?}"),
            }
        }
        // SAFETY: case is live until this free
        unsafe { ybus_case_free(case) };
    }

    #[test]
    fn every_field_code_is_its_matpower_column() {
        check_table("bus", 13, ybus_case_get_bus_f64, ybus_case_get_bus_i64);
        check_table("gen", 21, ybus_case_get_gen_f64, ybus_case_get_gen_i64);
        check_table(
            "branch",
            13,
            ybus_case_get_branch_f64,
            ybus_case_get_branch_i64,
        );
    }

    #[test]
    fn a_column_needs_a_known_field_and_room() {
        let case = three_bus();
        let mut out = [7.0; 5];
        // SAFETY: case is live, and out has room for 5 entries
        unsafe {
            for code in [0, 14, -1, YBUS_BUS_TYPE] {
                assert_eq!(
                    ybus_case_get_bus_f64(case, code, out.as_mut_ptr(), 5),
                    YBUS_ERROR_INVALID_FIELD
                );
            }
            assert_eq!(
                last_error().unwrap().to_str().unwrap(),
                "2 isn't a double field of the bus table"
            );
            assert_eq!(
                ybus_case_get_bus_f64(case, YBUS_BUS_PD_MW, out.as_mut_ptr(), 2),
                YBUS_ERROR_BUFFER_TOO_SMALL
            );
            assert_eq!(out, [7.0; 5]);
            assert_eq!(
                ybus_case_get_bus_f64(case, YBUS_BUS_PD_MW, out.as_mut_ptr(), 5),
                YBUS_OK
            );
            assert_eq!(out, [0.0, 50.0, 80.0, 7.0, 7.0]);
            ybus_case_free(case);
        }
    }

    #[test]
    fn bus_names_are_read_one_at_a_time() {
        let case = three_bus();
        let mut buffer = [b'x'; 8];
        let mut needed = 0;
        // Bus `row`'s name into the first `len` bytes of `buffer`
        let name = move |row: usize, len: usize, buffer: &mut [u8], needed: &mut usize| {
            assert!(len <= buffer.len());
            // SAFETY: case is live, buffer has room for len bytes, and needed is writable
            unsafe { ybus_case_bus_name(case, row, buffer.as_mut_ptr().cast(), len, needed) }
        };
        // SAFETY: case is live until freed, and needed is writable
        unsafe {
            assert_eq!(name(1, 8, &mut buffer, &mut needed), YBUS_OK);
            assert_eq!((&buffer[..5], needed), (&b"Beta\0"[..], 5));
            assert_eq!(
                ybus_case_bus_name(case, 2, ptr::null_mut(), 0, &mut needed),
                YBUS_OK
            );
            assert_eq!(needed, 6);
            let mut short = [b'x'; 5];
            assert_eq!(
                name(2, 5, &mut short, &mut needed),
                YBUS_ERROR_BUFFER_TOO_SMALL
            );
            assert_eq!((short, needed), ([b'x'; 5], 6));
            assert_eq!(
                name(3, 8, &mut buffer, &mut needed),
                YBUS_ERROR_ROW_OUT_OF_RANGE
            );
            (*case).case.busnames = None;
            assert_eq!(
                name(0, 8, &mut buffer, &mut needed),
                YBUS_ERROR_MISSING_DATA
            );
            assert_eq!(
                last_error().unwrap().to_str().unwrap(),
                "the case has no bus names"
            );
            ybus_case_free(case);
        }
    }

    type GetF64 = unsafe extern "C" fn(*const YbusCase, ybus_field, *mut f64, usize) -> ybus_status;
    type GetI64 = unsafe extern "C" fn(*const YbusCase, ybus_field, *mut i64, usize) -> ybus_status;
    type SetF64 = unsafe extern "C" fn(*mut YbusCase, ybus_field, *const f64, usize) -> ybus_status;
    type SetI64 = unsafe extern "C" fn(*mut YbusCase, ybus_field, *const i64, usize) -> ybus_status;

    /// Every column of a table, read through whichever getter takes its code, as numbers
    fn columns(
        case: *const YbusCase,
        width: ybus_field,
        rows: usize,
        get: (GetF64, GetI64),
    ) -> Vec<Vec<f64>> {
        (1..=width)
            .map(|code| {
                let (mut floats, mut ints) = (vec![0.0; rows], vec![0; rows]);
                // SAFETY: case is live, and each buffer has room for rows entries
                unsafe {
                    if get.0(case, code, floats.as_mut_ptr(), rows) == YBUS_OK {
                        floats
                    } else {
                        assert_eq!(get.1(case, code, ints.as_mut_ptr(), rows), YBUS_OK);
                        ints.iter().map(|&v| v as f64).collect()
                    }
                }
            })
            .collect()
    }

    /// Setting each column to its own values rotated by a row changes that column alone, so
    /// every setter writes the field its code names
    fn check_setters(width: ybus_field, rows: usize, get: (GetF64, GetI64), set: (SetF64, SetI64)) {
        let case = three_bus();
        for code in 1..=width {
            let before = columns(case, width, rows, get);
            let mut rotated = before[code as usize - 1].clone();
            rotated.rotate_left(1);
            let ints: Vec<i64> = rotated.iter().map(|&v| v as i64).collect();
            // SAFETY: case is live and unshared, and each array has rows entries
            let status = unsafe {
                match set.0(case, code, rotated.as_ptr(), rows) {
                    YBUS_ERROR_INVALID_FIELD => set.1(case, code, ints.as_ptr(), rows),
                    status => status,
                }
            };
            assert_eq!(status, YBUS_OK, "field {code}");
            let mut expected = before;
            expected[code as usize - 1] = rotated;
            assert_eq!(columns(case, width, rows, get), expected, "field {code}");
        }
        // SAFETY: case is live until this free
        unsafe { ybus_case_free(case) };
    }

    #[test]
    fn every_setter_writes_its_own_column() {
        check_setters(
            13,
            3,
            (ybus_case_get_bus_f64, ybus_case_get_bus_i64),
            (ybus_case_set_bus_f64, ybus_case_set_bus_i64),
        );
        check_setters(
            21,
            2,
            (ybus_case_get_gen_f64, ybus_case_get_gen_i64),
            (ybus_case_set_gen_f64, ybus_case_set_gen_i64),
        );
        check_setters(
            13,
            3,
            (ybus_case_get_branch_f64, ybus_case_get_branch_i64),
            (ybus_case_set_branch_f64, ybus_case_set_branch_i64),
        );
    }

    #[test]
    fn an_integer_column_is_checked_whole_before_it_is_set() {
        let case = three_bus();
        let read = |code| {
            let mut out = [0; 3];
            // SAFETY: case is live, and out has room for 3 entries
            assert_eq!(
                unsafe { ybus_case_get_bus_i64(case, code, out.as_mut_ptr(), 3) },
                YBUS_OK
            );
            out
        };
        // SAFETY: case is live and unshared until freed, and each array has its length
        unsafe {
            for (code, values) in [
                (YBUS_BUS_LABEL, [1, -2, 3]),
                (YBUS_BUS_TYPE, [1, 2, 5]),
                (YBUS_BUS_TYPE, [0, 2, 1]),
                (YBUS_BUS_ZONE, [i64::MIN, 0, 0]),
            ] {
                let before = read(code);
                assert_eq!(
                    ybus_case_set_bus_i64(case, code, values.as_ptr(), 3),
                    YBUS_ERROR_INVALID_VALUE
                );
                assert_eq!(read(code), before);
            }
            assert_eq!(
                last_error().unwrap().to_str().unwrap(),
                format!(
                    "bus row 0's value {} for field 11 isn't a non-negative integer that fits \
                     a size_t",
                    i64::MIN
                )
            );
            let statuses = [1, 2, 0];
            assert_eq!(
                ybus_case_set_branch_i64(case, YBUS_BRANCH_STATUS, statuses.as_ptr(), 3),
                YBUS_ERROR_INVALID_VALUE
            );
            assert_eq!(
                last_error().unwrap().to_str().unwrap(),
                "branch row 1's value 2 for field 11 isn't a status, 0 or 1"
            );
            let off = [0, 1];
            assert_eq!(
                ybus_case_set_gen_i64(case, YBUS_GEN_STATUS, off.as_ptr(), 2),
                YBUS_OK
            );
            let generators = &(*case).case.generator;
            assert!(!generators[0].status && generators[1].status);

            // The length must be one per row, even for a null array
            let pd = [1.0, 2.0];
            assert_eq!(
                ybus_case_set_bus_f64(case, YBUS_BUS_PD_MW, pd.as_ptr(), 2),
                YBUS_ERROR_LENGTH_MISMATCH
            );
            assert_eq!(
                ybus_case_set_bus_f64(case, YBUS_BUS_PD_MW, ptr::null(), 0),
                YBUS_ERROR_LENGTH_MISMATCH
            );
            assert_eq!(
                last_error().unwrap().to_str().unwrap(),
                "`values` has 0 entries, but its table has 3 rows"
            );
            assert_eq!(
                ybus_case_set_bus_f64(case, YBUS_BUS_TYPE, pd.as_ptr(), 2),
                YBUS_ERROR_INVALID_FIELD
            );
            ybus_case_free(case);
        }
    }
}
