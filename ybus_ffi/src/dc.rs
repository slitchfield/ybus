//! The DC power flow through C: its options.

use powerflow::dc::DcOptions;
use sparse_mat::compressed_matrix::Ordering;

use crate::boundary::{arg_mut, free_handle, get_value, guard, new_handle, Failure};
use crate::status::ybus_status;

/// A fill-reducing ordering, for the DC power flow's B′.
#[allow(non_camel_case_types)]
pub type ybus_ordering = i32;

/// No reordering.
pub const YBUS_ORDERING_NATURAL: ybus_ordering = 1;
/// Reverse Cuthill–McKee, which keeps entries near the diagonal.
pub const YBUS_ORDERING_RCM: ybus_ordering = 2;
/// Approximate minimum degree, which keeps the factor's fill small: the default.
pub const YBUS_ORDERING_AMD: ybus_ordering = 3;
/// Column approximate minimum degree (COLAMD), for LU.
pub const YBUS_ORDERING_COLAMD: ybus_ordering = 4;

/// The code of an ordering
fn ordering_code(ordering: Ordering) -> Result<ybus_ordering, Failure> {
    match ordering {
        Ordering::Natural => Ok(YBUS_ORDERING_NATURAL),
        Ordering::Rcm => Ok(YBUS_ORDERING_RCM),
        Ordering::Amd => Ok(YBUS_ORDERING_AMD),
        Ordering::Colamd => Ok(YBUS_ORDERING_COLAMD),
        // Ordering is #[non_exhaustive]; only C's codes can set one here
        other => Err(Failure::InvalidValue(format!(
            "the ordering {other:?} has no code in this library"
        ))),
    }
}

/// The ordering of a code
fn ordering_of(code: ybus_ordering) -> Result<Ordering, Failure> {
    match code {
        YBUS_ORDERING_NATURAL => Ok(Ordering::Natural),
        YBUS_ORDERING_RCM => Ok(Ordering::Rcm),
        YBUS_ORDERING_AMD => Ok(Ordering::Amd),
        YBUS_ORDERING_COLAMD => Ok(Ordering::Colamd),
        other => Err(Failure::InvalidValue(format!(
            "{other} isn't an ordering, 1 to 4"
        ))),
    }
}

/// Options for the DC power flow, made by `ybus_dc_options_new` at MATPOWER's defaults and
/// freed by `ybus_dc_options_free`. A null options handle means the defaults wherever one is
/// taken.
pub struct YbusDcOptions {
    pub(crate) options: DcOptions,
}

/// Makes DC power flow options at the defaults (B′ ordered by `YBUS_ORDERING_AMD`) in `*out`,
/// which is set to `NULL` first and stays `NULL` on an error.
///
/// # Safety
///
/// `out` must be null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_options_new(out: *mut *mut YbusDcOptions) -> ybus_status {
    guard(|| {
        // SAFETY: out is as the caller promised
        unsafe {
            new_handle(out, || {
                Ok(YbusDcOptions {
                    options: DcOptions::default(),
                })
            })
        }
    })
}

/// Frees `options`, which mustn't be used again. Does nothing if `options` is null.
///
/// # Safety
///
/// `options` must be null or options made by this library and not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_options_free(options: *mut YbusDcOptions) {
    // SAFETY: as the caller promised
    unsafe { free_handle(options) }
}

/// Sets the ordering B′ is factored under. Every ordering gives the same solution, to within
/// rounding. Returns `YBUS_ERROR_INVALID_VALUE` for a code that isn't an ordering.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_options_set_ordering(
    options: *mut YbusDcOptions,
    code: ybus_ordering,
) -> ybus_status {
    guard(|| {
        // SAFETY: options is as the caller promised
        let options = unsafe { arg_mut(options, "options") }?;
        options.options.ordering = ordering_of(code)?;
        Ok(())
    })
}

/// Writes the ordering B′ is factored under through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `ybus_ordering`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_options_get_ordering(
    options: *const YbusDcOptions,
    out: *mut ybus_ordering,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(options, "options", out, |options| {
            ordering_code(options.options.ordering)
        })
    }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use super::*;
    use crate::status::{YBUS_ERROR_INVALID_VALUE, YBUS_ERROR_NULL_ARGUMENT, YBUS_OK};

    #[test]
    fn dc_options_start_at_the_defaults_and_set_each_ordering() {
        let mut options = ptr::null_mut();
        let mut code = 0;
        // SAFETY: options is live until freed, and code is writable
        unsafe {
            assert_eq!(ybus_dc_options_new(&mut options), YBUS_OK);
            assert_eq!((*options).options, DcOptions::default());
            assert_eq!(ybus_dc_options_get_ordering(options, &mut code), YBUS_OK);
            assert_eq!(code, YBUS_ORDERING_AMD);
            for (code, expected) in [
                (YBUS_ORDERING_NATURAL, Ordering::Natural),
                (YBUS_ORDERING_RCM, Ordering::Rcm),
                (YBUS_ORDERING_COLAMD, Ordering::Colamd),
                (YBUS_ORDERING_AMD, Ordering::Amd),
            ] {
                assert_eq!(ybus_dc_options_set_ordering(options, code), YBUS_OK);
                assert_eq!((*options).options.ordering, expected);
            }
            for bad in [0, 5, -1] {
                assert_eq!(
                    ybus_dc_options_set_ordering(options, bad),
                    YBUS_ERROR_INVALID_VALUE
                );
            }
            assert_eq!((*options).options.ordering, Ordering::Amd);
            assert_eq!(
                ybus_dc_options_get_ordering(ptr::null(), &mut code),
                YBUS_ERROR_NULL_ARGUMENT
            );
            ybus_dc_options_free(options);
            ybus_dc_options_free(ptr::null_mut());
        }
    }
}
