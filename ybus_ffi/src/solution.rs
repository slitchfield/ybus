//! What every solution handle shares: reading one of its columns whole, by the case's field
//! codes for the table the column belongs to.

use crate::boundary::{arg_ref, arg_slice_mut, check_room, guard, Failure};
use crate::case::columns::ybus_field;
use crate::status::ybus_status;

/// Writes the column of `solution` that `field` names in `table` into `out`, which has room for
/// `len` entries: one per row of the table. `column` picks the column, or `None` if the
/// solution has no such column in that table (`YBUS_ERROR_INVALID_FIELD`).
///
/// # Safety
///
/// `solution` must be null or a live handle, and `out` null or `len` writable `double`s.
pub(crate) unsafe fn get_column<S>(
    solution: *const S,
    table: &'static str,
    field: ybus_field,
    out: *mut f64,
    len: usize,
    column: fn(&S, ybus_field) -> Option<&[f64]>,
) -> ybus_status {
    guard(|| {
        // SAFETY: solution and out are as the caller promised
        let (solution, out) = unsafe {
            (
                arg_ref(solution, "solution")?,
                arg_slice_mut(out, len, "out")?,
            )
        };
        let values = column(solution, field).ok_or(Failure::InvalidField {
            table,
            field,
            kind: "solution",
        })?;
        check_room("out", values.len(), len)?;
        out[..values.len()].copy_from_slice(values);
        Ok(())
    })
}
