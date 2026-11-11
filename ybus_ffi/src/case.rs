//! Cases: a [`powerflow::Case`] behind an opaque `ybus_case` handle, loaded from a MATPOWER
//! file or its text, copied and freed.

use std::ffi::c_char;

use powerflow::Case;

use crate::boundary::{
    arg_mut, arg_ref, arg_slice, arg_str, free_handle, get_value, guard, new_handle,
};
use crate::status::ybus_status;

pub mod columns;
pub mod costs;

/// A power system case: its buses, generators, branches and optional costs and bus names, as
/// read from a MATPOWER case file. Made by `ybus_case_from_path`, `ybus_case_from_bytes` or
/// `ybus_case_clone`, and freed by `ybus_case_free`.
pub struct YbusCase {
    pub(crate) case: Case,
}

/// Reads and validates the MATPOWER case file at `path`, a NUL-terminated UTF-8 path, into a
/// new case in `*out`.
///
/// Returns `YBUS_ERROR_NO_CASE_FILE` if the file can't be read, `YBUS_ERROR_INVALID_UTF8` for a
/// path that isn't UTF-8, and the parser's and validation's codes for its contents. `*out` is
/// set to `NULL` first, and stays `NULL` on any error.
///
/// # Safety
///
/// `path` must be null or a NUL-terminated string, and `out` null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_from_path(
    path: *const c_char,
    out: *mut *mut YbusCase,
) -> ybus_status {
    guard(|| {
        // SAFETY: out and path are as the caller promised
        unsafe {
            new_handle(out, || {
                let path = arg_str(path, "path")?;
                Ok(YbusCase {
                    case: Case::from_path(path)?,
                })
            })
        }
    })
}

/// Parses and validates the `len` bytes of a MATPOWER case file's text at `bytes` into a new
/// case in `*out`. The text is decoded as UTF-8, or as Latin-1 if it isn't valid UTF-8.
///
/// Returns the parser's and validation's codes, such as `YBUS_ERROR_PARSE`. `*out` is set to
/// `NULL` first, and stays `NULL` on any error.
///
/// # Safety
///
/// `bytes` must be null or point to `len` readable bytes, and `out` null or a writable
/// pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_from_bytes(
    bytes: *const u8,
    len: usize,
    out: *mut *mut YbusCase,
) -> ybus_status {
    guard(|| {
        // SAFETY: out and bytes are as the caller promised
        unsafe {
            new_handle(out, || {
                let bytes = arg_slice(bytes, len, "bytes")?;
                Ok(YbusCase {
                    case: Case::from_bytes(bytes)?,
                })
            })
        }
    })
}

/// Copies `case` into a new case in `*out`, which can then be changed without changing
/// `case`. `*out` is set to `NULL` first, and stays `NULL` on any error.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_clone(
    case: *const YbusCase,
    out: *mut *mut YbusCase,
) -> ybus_status {
    guard(|| {
        // SAFETY: out and case are as the caller promised
        unsafe {
            new_handle(out, || {
                let case = arg_ref(case, "case")?;
                Ok(YbusCase {
                    case: case.case.clone(),
                })
            })
        }
    })
}

/// Frees `case`, which mustn't be used again. Does nothing if `case` is null.
///
/// # Safety
///
/// `case` must be null or a case made by this library and not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_free(case: *mut YbusCase) {
    // SAFETY: as the caller promised
    unsafe { free_handle(case) }
}

/// Writes `read(case)` through `out`: the shape of every scalar getter.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable `T`.
unsafe fn get<T>(case: *const YbusCase, out: *mut T, read: impl FnOnce(&Case) -> T) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(case, "case", out, |case| Ok(read(&case.case))) }
}

/// Writes the number of buses in `case` through `out`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_bus_count(
    case: *const YbusCase,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get(case, out, |case| case.bus.len()) }
}

/// Writes the number of generators in `case` through `out`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_gen_count(
    case: *const YbusCase,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get(case, out, |case| case.generator.len()) }
}

/// Writes the number of branches in `case` through `out`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_branch_count(
    case: *const YbusCase,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get(case, out, |case| case.branch.len()) }
}

/// Writes the system MVA base of `case`, which per-unit values are relative to, through `out`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_base_mva(case: *const YbusCase, out: *mut f64) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get(case, out, |case| case.base_mva) }
}

/// Sets the system MVA base of `case`. It isn't checked here: a base that isn't positive and
/// finite is reported when the case is solved.
///
/// # Safety
///
/// `case` must be null or a live case that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_set_base_mva(case: *mut YbusCase, base_mva: f64) -> ybus_status {
    guard(|| {
        // SAFETY: case is as the caller promised
        let case = unsafe { arg_mut(case, "case") }?;
        case.case.base_mva = base_mva;
        Ok(())
    })
}

/// Writes whether `case` has generator costs through `out`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable `bool`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_has_costs(case: *const YbusCase, out: *mut bool) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get(case, out, |case| case.gencost.is_some()) }
}

/// Writes whether `case` has bus names through `out`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable `bool`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_case_has_bus_names(
    case: *const YbusCase,
    out: *mut bool,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get(case, out, |case| case.busnames.is_some()) }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::ffi::CString;
    use std::ptr;

    use super::*;
    use crate::boundary::last_error;
    use crate::status::{
        YBUS_ERROR_INVALID_UTF8, YBUS_ERROR_NO_CASE_FILE, YBUS_ERROR_NULL_ARGUMENT,
        YBUS_ERROR_PARSE, YBUS_OK,
    };

    /// A small case with every optional table: three buses, two generators, three branches,
    /// a polynomial and a piecewise linear cost, and bus names. The C tests' `cases.h` holds
    /// the same text.
    pub(crate) const THREE_BUS: &str = include_str!("../tests/c/three_bus.m");

    /// THREE_BUS loaded through the export
    pub(crate) fn three_bus() -> *mut YbusCase {
        let mut case = ptr::null_mut();
        // SAFETY: the text is THREE_BUS's bytes, and case is writable
        let status =
            unsafe { ybus_case_from_bytes(THREE_BUS.as_ptr(), THREE_BUS.len(), &mut case) };
        assert_eq!(status, YBUS_OK);
        case
    }

    fn message() -> String {
        last_error().unwrap().into_string().unwrap()
    }

    #[test]
    fn a_case_loads_from_bytes_and_from_a_file_and_clones() {
        let from_bytes = three_bus();
        let path = std::env::temp_dir().join(format!("ybus_ffi_{}.m", std::process::id()));
        std::fs::write(&path, THREE_BUS).unwrap();
        let c_path = CString::new(path.to_str().unwrap()).unwrap();
        let (mut from_file, mut copy) = (ptr::null_mut(), ptr::null_mut());
        // SAFETY: c_path is NUL-terminated, the outs are writable, and each case is live until
        //   freed once
        unsafe {
            assert_eq!(
                ybus_case_from_path(c_path.as_ptr(), &mut from_file),
                YBUS_OK
            );
            assert_eq!(ybus_case_clone(from_bytes, &mut copy), YBUS_OK);
            assert_eq!((*from_file).case, (*from_bytes).case);
            assert_eq!((*copy).case, (*from_bytes).case);
            assert_eq!((*copy).case.bus.len(), 3);
            for case in [from_bytes, from_file, copy] {
                ybus_case_free(case);
            }
            ybus_case_free(ptr::null_mut());
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn loading_errors_have_their_codes_and_leave_out_null() {
        let missing = CString::new("/nonexistent/ybus/case.m").unwrap();
        let mut out = ptr::NonNull::<YbusCase>::dangling().as_ptr();
        // SAFETY: each pointer is NUL-terminated or null, or names len readable bytes; out is
        //   writable
        unsafe {
            assert_eq!(
                ybus_case_from_path(missing.as_ptr(), &mut out),
                YBUS_ERROR_NO_CASE_FILE
            );
            assert!(out.is_null());
            assert!(message().contains("/nonexistent/ybus/case.m"));

            let text = b"mpc.version = '2';\nmpc.bus = [1 2;\n";
            out = ptr::NonNull::dangling().as_ptr();
            assert_eq!(
                ybus_case_from_bytes(text.as_ptr(), text.len(), &mut out),
                YBUS_ERROR_PARSE
            );
            assert!(out.is_null());

            let latin1 = b"/tmp/\xff.m\0";
            assert_eq!(
                ybus_case_from_path(latin1.as_ptr().cast(), &mut out),
                YBUS_ERROR_INVALID_UTF8
            );
            assert_eq!(message(), "`path` isn't valid UTF-8");

            assert_eq!(
                ybus_case_from_path(missing.as_ptr(), ptr::null_mut()),
                YBUS_ERROR_NULL_ARGUMENT
            );
            assert_eq!(message(), "`out` is null");
            assert_eq!(
                ybus_case_from_path(ptr::null(), &mut out),
                YBUS_ERROR_NULL_ARGUMENT
            );
            assert_eq!(message(), "`path` is null");
            assert_eq!(
                ybus_case_clone(ptr::null(), &mut out),
                YBUS_ERROR_NULL_ARGUMENT
            );
            assert!(out.is_null());
        }
    }

    #[test]
    fn sizes_and_scalars_are_read_and_the_base_set() {
        let case = three_bus();
        let (mut buses, mut generators, mut branches) = (0, 0, 0);
        let (mut base, mut costs, mut names) = (0.0, false, false);
        // SAFETY: case is live until freed, and every out is writable
        unsafe {
            assert_eq!(ybus_case_bus_count(case, &mut buses), YBUS_OK);
            assert_eq!(ybus_case_gen_count(case, &mut generators), YBUS_OK);
            assert_eq!(ybus_case_branch_count(case, &mut branches), YBUS_OK);
            assert_eq!(ybus_case_base_mva(case, &mut base), YBUS_OK);
            assert_eq!(ybus_case_has_costs(case, &mut costs), YBUS_OK);
            assert_eq!(ybus_case_has_bus_names(case, &mut names), YBUS_OK);
            assert_eq!(
                (buses, generators, branches, base, costs, names),
                (3, 2, 3, 100.0, true, true)
            );

            // Not checked when set: a zero base is the solve's to report
            assert_eq!(ybus_case_set_base_mva(case, 0.0), YBUS_OK);
            assert_eq!(ybus_case_base_mva(case, &mut base), YBUS_OK);
            assert_eq!(base, 0.0);

            (*case).case.gencost = None;
            (*case).case.busnames = None;
            assert_eq!(ybus_case_has_costs(case, &mut costs), YBUS_OK);
            assert_eq!(ybus_case_has_bus_names(case, &mut names), YBUS_OK);
            assert_eq!((costs, names), (false, false));

            assert_eq!(
                ybus_case_bus_count(ptr::null(), &mut buses),
                YBUS_ERROR_NULL_ARGUMENT
            );
            assert_eq!(
                ybus_case_bus_count(case, ptr::null_mut()),
                YBUS_ERROR_NULL_ARGUMENT
            );
            assert_eq!(message(), "`out` is null");
            assert_eq!(
                ybus_case_set_base_mva(ptr::null_mut(), 1.0),
                YBUS_ERROR_NULL_ARGUMENT
            );
            ybus_case_free(case);
        }
    }
}
