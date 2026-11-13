//! A case's generator costs, read and set whole as flat arrays: per row, a model, a startup
//! and a shutdown cost, and a run of values in one shared array, from `offsets[i]` to
//! `offsets[i + 1]`.
//!
//! A polynomial's values are its coefficients, highest order first; a piecewise linear cost's
//! are its `(p, cost)` points, interleaved as `p0, c0, p1, c1, …`. The rows are one per
//! generator for real power, optionally followed by one per generator for reactive power.

use ybus_powerflow::case::GenCost;

use super::YbusCase;
use crate::boundary::{
    arg_mut, arg_ref, arg_slice, arg_slice_mut, check_length, check_room, guard, Failure,
};
use crate::status::ybus_status;

/// A piecewise linear cost, MATPOWER's model 1: its values are `(p, cost)` points.
pub const YBUS_COST_PIECEWISE_LINEAR: i32 = 1;
/// A polynomial cost, MATPOWER's model 2: its values are coefficients, highest order first.
pub const YBUS_COST_POLYNOMIAL: i32 = 2;

/// A cost row's model, startup and shutdown costs, and values
fn flatten(cost: &GenCost) -> (i32, f64, f64, Vec<f64>) {
    match cost {
        GenCost::Polynomial {
            startup,
            shutdown,
            coeffs,
        } => (YBUS_COST_POLYNOMIAL, *startup, *shutdown, coeffs.clone()),
        GenCost::PiecewiseLinear {
            startup,
            shutdown,
            points,
        } => (
            YBUS_COST_PIECEWISE_LINEAR,
            *startup,
            *shutdown,
            points.iter().flat_map(|&(p, cost)| [p, cost]).collect(),
        ),
    }
}

/// The number of values a cost row has in the flat array
fn value_count(cost: &GenCost) -> usize {
    match cost {
        GenCost::Polynomial { coeffs, .. } => coeffs.len(),
        GenCost::PiecewiseLinear { points, .. } => 2 * points.len(),
    }
}

/// Writes the number of cost rows in `case` through `rows`, 0 if it has no costs, and the
/// number of values all of them hold through `values`: the sizes `ybus_case_get_costs` needs.
///
/// # Safety
///
/// `case` must be null or a live case, and `rows` and `values` null or writable `size_t`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_cost_sizes(
    case: *const YbusCase,
    rows: *mut usize,
    values: *mut usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: case, rows and values are as the caller promised
        let (case, rows, values) = unsafe {
            (
                arg_ref(case, "case")?,
                arg_mut(rows, "rows")?,
                arg_mut(values, "values")?,
            )
        };
        let costs = case.case.gencost.as_deref().unwrap_or_default();
        *rows = costs.len();
        *values = costs.iter().map(value_count).sum();
        Ok(())
    })
}

/// Writes the costs of `case` into flat arrays: each row's model (`YBUS_COST_POLYNOMIAL` or
/// `YBUS_COST_PIECEWISE_LINEAR`), startup and shutdown costs, and values, which are
/// `values[offsets[i]]` up to `values[offsets[i + 1]]`. `model`, `startup` and `shutdown` have
/// room for `rows` entries, `offsets` for `rows + 1`, and `values` for `values_len`.
///
/// Returns `YBUS_ERROR_MISSING_DATA` if the case has no costs, and
/// `YBUS_ERROR_BUFFER_TOO_SMALL`, writing nothing, if `rows` or `values_len` is less than
/// `ybus_case_cost_sizes` gives.
///
/// # Safety
///
/// `case` must be null or a live case; `model`, `startup` and `shutdown` null or `rows`
/// writable entries; `offsets` null or `rows + 1`; and `values` null or `values_len`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn ybus_case_get_costs(
    case: *const YbusCase,
    model: *mut i32,
    startup: *mut f64,
    shutdown: *mut f64,
    offsets: *mut usize,
    rows: usize,
    values: *mut f64,
    values_len: usize,
) -> ybus_status {
    guard(|| {
        let offsets_len = rows
            .checked_add(1)
            .ok_or_else(|| Failure::InvalidValue("`rows` is too large".to_string()))?;
        // SAFETY: every pointer is as the caller promised
        let (case, model, startup, shutdown, offsets, values) = unsafe {
            (
                arg_ref(case, "case")?,
                arg_slice_mut(model, rows, "model")?,
                arg_slice_mut(startup, rows, "startup")?,
                arg_slice_mut(shutdown, rows, "shutdown")?,
                arg_slice_mut(offsets, offsets_len, "offsets")?,
                arg_slice_mut(values, values_len, "values")?,
            )
        };
        let costs = case
            .case
            .gencost
            .as_deref()
            .ok_or(Failure::MissingData("costs"))?;
        check_room("rows", costs.len(), rows)?;
        check_room("values", costs.iter().map(value_count).sum(), values_len)?;
        let mut end = 0;
        offsets[0] = 0;
        for (row, cost) in costs.iter().enumerate() {
            let (row_model, row_startup, row_shutdown, row_values) = flatten(cost);
            (model[row], startup[row], shutdown[row]) = (row_model, row_startup, row_shutdown);
            values[end..end + row_values.len()].copy_from_slice(&row_values);
            end += row_values.len();
            offsets[row + 1] = end;
        }
        Ok(())
    })
}

/// Replaces the costs of `case` with `rows` rows given as `ybus_case_get_costs` writes them:
/// each row's model, startup and shutdown costs, and values `values[offsets[i]]` up to
/// `values[offsets[i + 1]]`. `rows` 0 removes the costs (and `values_len` must be 0).
///
/// The structure is checked before anything changes: `YBUS_ERROR_INVALID_VALUE` for a model
/// other than 1 or 2, or a piecewise linear row with an odd number of values, and
/// `YBUS_ERROR_LENGTH_MISMATCH` unless `offsets` starts at 0, never decreases and ends at
/// `values_len`. What a solve can't use, such as a row count that isn't the generator count or
/// a cubic polynomial, is reported when the case is solved.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using; `model`, `startup` and
/// `shutdown` null or `rows` readable entries; `offsets` null or `rows + 1` (or null when
/// `rows` is 0); and `values` null or `values_len`.
#[unsafe(no_mangle)]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn ybus_case_set_costs(
    case: *mut YbusCase,
    model: *const i32,
    startup: *const f64,
    shutdown: *const f64,
    offsets: *const usize,
    rows: usize,
    values: *const f64,
    values_len: usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: case is as the caller promised
        let case = unsafe { arg_mut(case, "case") }?;
        if rows == 0 {
            check_length("values", 0, values_len)?;
            case.case.gencost = None;
            return Ok(());
        }
        let offsets_len = rows
            .checked_add(1)
            .ok_or_else(|| Failure::InvalidValue("`rows` is too large".to_string()))?;
        // SAFETY: every pointer is as the caller promised
        let (model, startup, shutdown, offsets, values) = unsafe {
            (
                arg_slice(model, rows, "model")?,
                arg_slice(startup, rows, "startup")?,
                arg_slice(shutdown, rows, "shutdown")?,
                arg_slice(offsets, offsets_len, "offsets")?,
                arg_slice(values, values_len, "values")?,
            )
        };
        if offsets[0] != 0 || offsets[rows] != values_len {
            return Err(Failure::LengthMismatch(format!(
                "`offsets` runs from {} to {}, but `values` has {values_len} entries from 0",
                offsets[0], offsets[rows]
            )));
        }
        if let Some(row) = offsets.windows(2).position(|pair| pair[1] < pair[0]) {
            return Err(Failure::LengthMismatch(format!(
                "`offsets` decreases after row {row}"
            )));
        }
        let costs = (0..rows)
            .map(|row| {
                let row_values = &values[offsets[row]..offsets[row + 1]];
                let (startup, shutdown) = (startup[row], shutdown[row]);
                match model[row] {
                    YBUS_COST_POLYNOMIAL => Ok(GenCost::Polynomial {
                        startup,
                        shutdown,
                        coeffs: row_values.to_vec(),
                    }),
                    YBUS_COST_PIECEWISE_LINEAR if row_values.len() % 2 == 0 => {
                        Ok(GenCost::PiecewiseLinear {
                            startup,
                            shutdown,
                            points: row_values
                                .as_chunks::<2>()
                                .0
                                .iter()
                                .map(|&[p, cost]| (p, cost))
                                .collect(),
                        })
                    }
                    YBUS_COST_PIECEWISE_LINEAR => Err(Failure::InvalidValue(format!(
                        "cost row {row} is piecewise linear, but has an odd number of values, {}",
                        row_values.len()
                    ))),
                    other => Err(Failure::InvalidValue(format!(
                        "cost row {row}'s model {other} isn't 1 (piecewise linear) or 2 \
                         (polynomial)"
                    ))),
                }
            })
            .collect::<Result<Vec<GenCost>, Failure>>()?;
        case.case.gencost = Some(costs);
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use super::*;
    use crate::boundary::last_error;
    use crate::case::tests::three_bus;
    use crate::case::ybus_case_free;
    use crate::status::{
        YBUS_ERROR_BUFFER_TOO_SMALL, YBUS_ERROR_INVALID_VALUE, YBUS_ERROR_LENGTH_MISMATCH,
        YBUS_ERROR_MISSING_DATA, YBUS_OK,
    };

    /// A case's costs, read through the exports
    #[derive(Debug, PartialEq)]
    struct Flat {
        model: Vec<i32>,
        startup: Vec<f64>,
        shutdown: Vec<f64>,
        offsets: Vec<usize>,
        values: Vec<f64>,
    }

    fn read(case: *const YbusCase) -> Result<Flat, ybus_status> {
        let (mut rows, mut count) = (0, 0);
        // SAFETY: case is live, and every buffer has the room it's said to
        unsafe {
            assert_eq!(ybus_case_cost_sizes(case, &mut rows, &mut count), YBUS_OK);
            let mut flat = Flat {
                model: vec![0; rows],
                startup: vec![0.0; rows],
                shutdown: vec![0.0; rows],
                offsets: vec![0; rows + 1],
                values: vec![0.0; count],
            };
            match ybus_case_get_costs(
                case,
                flat.model.as_mut_ptr(),
                flat.startup.as_mut_ptr(),
                flat.shutdown.as_mut_ptr(),
                flat.offsets.as_mut_ptr(),
                rows,
                flat.values.as_mut_ptr(),
                count,
            ) {
                YBUS_OK => Ok(flat),
                status => Err(status),
            }
        }
    }

    fn write(case: *mut YbusCase, flat: &Flat) -> ybus_status {
        // SAFETY: case is live and unshared, and every array has the length it's said to
        unsafe {
            ybus_case_set_costs(
                case,
                flat.model.as_ptr(),
                flat.startup.as_ptr(),
                flat.shutdown.as_ptr(),
                flat.offsets.as_ptr(),
                flat.model.len(),
                flat.values.as_ptr(),
                flat.values.len(),
            )
        }
    }

    fn message() -> String {
        last_error().unwrap().into_string().unwrap()
    }

    #[test]
    fn the_costs_read_flat_and_round_trip() {
        let case = three_bus();
        let fixture = Flat {
            model: vec![2, 1],
            startup: vec![1500.0, 0.0],
            shutdown: vec![100.0, 0.0],
            offsets: vec![0, 3, 9],
            values: vec![0.11, 5.0, 150.0, 0.0, 0.0, 40.0, 2000.0, 80.0, 4400.0],
        };
        assert_eq!(read(case), Ok(fixture));

        // Real and reactive costs, of both models, read back exactly and as the Rust ones
        let costs = Flat {
            model: vec![1, 2, 2, 1],
            startup: vec![10.0, 20.0, 0.0, 0.0],
            shutdown: vec![1.0, 2.0, 0.0, 0.0],
            offsets: vec![0, 4, 6, 7, 11],
            values: vec![
                0.0, 0.0, 100.0, 3000.0, 25.0, 0.5, 0.0, -10.0, 1.0, 10.0, 2.0,
            ],
        };
        assert_eq!(write(case, &costs), YBUS_OK);
        assert_eq!(read(case).as_ref(), Ok(&costs));
        // SAFETY: case is live
        let rust = unsafe { (*case).case.gencost.clone() }.unwrap();
        assert_eq!(
            rust[0],
            GenCost::PiecewiseLinear {
                startup: 10.0,
                shutdown: 1.0,
                points: vec![(0.0, 0.0), (100.0, 3000.0)]
            }
        );
        assert_eq!(
            rust[2],
            GenCost::Polynomial {
                startup: 0.0,
                shutdown: 0.0,
                coeffs: vec![0.0]
            }
        );

        // No rows removes them
        // SAFETY: case is live and unshared; with no rows, no other array is read
        unsafe {
            assert_eq!(
                ybus_case_set_costs(
                    case,
                    ptr::null(),
                    ptr::null(),
                    ptr::null(),
                    ptr::null(),
                    0,
                    ptr::null(),
                    0
                ),
                YBUS_OK
            );
            assert_eq!((*case).case.gencost, None);
        }
        assert_eq!(read(case), Err(YBUS_ERROR_MISSING_DATA));
        assert_eq!(message(), "the case has no costs");
        // SAFETY: case is live until this free
        unsafe { ybus_case_free(case) };
    }

    #[test]
    fn a_malformed_cost_table_changes_nothing() {
        let case = three_bus();
        let before = read(case).unwrap();
        let bad = |change: fn(&mut Flat)| {
            let mut flat = read(case).unwrap();
            change(&mut flat);
            flat
        };
        for (flat, status, expected) in [
            (
                bad(|flat| flat.model[1] = 3),
                YBUS_ERROR_INVALID_VALUE,
                "cost row 1's model 3 isn't 1 (piecewise linear) or 2 (polynomial)",
            ),
            (
                bad(|flat| flat.offsets[1] = 4),
                YBUS_ERROR_INVALID_VALUE,
                "cost row 1 is piecewise linear, but has an odd number of values, 5",
            ),
            (
                bad(|flat| flat.offsets = vec![0, 10, 9]),
                YBUS_ERROR_LENGTH_MISMATCH,
                "`offsets` decreases after row 1",
            ),
            (
                bad(|flat| flat.offsets[2] = 8),
                YBUS_ERROR_LENGTH_MISMATCH,
                "`offsets` runs from 0 to 8, but `values` has 9 entries from 0",
            ),
            (
                bad(|flat| flat.offsets[0] = 1),
                YBUS_ERROR_LENGTH_MISMATCH,
                "`offsets` runs from 1 to 9, but `values` has 9 entries from 0",
            ),
        ] {
            assert_eq!(write(case, &flat), status, "{expected}");
            assert_eq!(message(), expected);
            assert_eq!(read(case).as_ref(), Ok(&before));
        }

        // Too little room is refused, writing nothing
        let mut model = [-1; 1];
        let mut values = [-1.0; 9];
        let (mut startup, mut shutdown, mut offsets) = ([0.0; 2], [0.0; 2], [0; 3]);
        // SAFETY: case is live, and each buffer has room for the length given
        unsafe {
            assert_eq!(
                ybus_case_get_costs(
                    case,
                    model.as_mut_ptr(),
                    startup.as_mut_ptr(),
                    shutdown.as_mut_ptr(),
                    offsets.as_mut_ptr(),
                    1,
                    values.as_mut_ptr(),
                    9
                ),
                YBUS_ERROR_BUFFER_TOO_SMALL
            );
            assert_eq!((model, values), ([-1], [-1.0; 9]));
            ybus_case_free(case);
        }
    }
}
