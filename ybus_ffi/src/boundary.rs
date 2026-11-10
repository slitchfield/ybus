//! What every export goes through: [`guard`], which turns a failure or a panic into a status
//! code and records its message as the calling thread's last error, and the checks on the
//! pointers C passes in.

use std::any::Any;
use std::cell::RefCell;
use std::ffi::{c_char, CStr, CString};
use std::fmt;
use std::panic::{catch_unwind, AssertUnwindSafe};

use powerflow::PowerflowError;

use crate::case::columns::ybus_field;
use crate::status::{
    self, ybus_status, YBUS_ERROR_BUFFER_TOO_SMALL, YBUS_ERROR_INVALID_FIELD,
    YBUS_ERROR_INVALID_UTF8, YBUS_ERROR_INVALID_VALUE, YBUS_ERROR_LENGTH_MISMATCH,
    YBUS_ERROR_MISALIGNED_ARGUMENT, YBUS_ERROR_MISSING_DATA, YBUS_ERROR_NULL_ARGUMENT,
    YBUS_ERROR_PANIC, YBUS_OK,
};

/// Why an export failed.
#[derive(Debug)]
pub(crate) enum Failure {
    /// `powerflow` refused the call.
    Powerflow(PowerflowError),
    /// The named pointer argument was null.
    NullArgument(&'static str),
    /// The named pointer argument wasn't aligned for its type.
    MisalignedArgument(&'static str),
    /// The named string argument wasn't valid UTF-8.
    InvalidUtf8(&'static str),
    /// The named buffer was shorter than the result.
    BufferTooSmall {
        argument: &'static str,
        needed: usize,
        given: usize,
    },
    /// The named array's length doesn't match its table.
    LengthMismatch {
        argument: &'static str,
        expected: usize,
        given: usize,
    },
    /// A value can't be held where it was given or asked for, for this reason.
    InvalidValue(String),
    /// A field code isn't one of the table's fields of this kind (`double` or `int64_t`).
    InvalidField {
        table: &'static str,
        field: ybus_field,
        kind: &'static str,
    },
    /// The case has none of the named data.
    MissingData(&'static str),
    /// The call panicked, with this message.
    Panic(String),
}

impl Failure {
    /// The status code C sees.
    fn code(&self) -> ybus_status {
        match self {
            Failure::Powerflow(error) => status::powerflow_code(error),
            Failure::NullArgument(_) => YBUS_ERROR_NULL_ARGUMENT,
            Failure::MisalignedArgument(_) => YBUS_ERROR_MISALIGNED_ARGUMENT,
            Failure::InvalidUtf8(_) => YBUS_ERROR_INVALID_UTF8,
            Failure::BufferTooSmall { .. } => YBUS_ERROR_BUFFER_TOO_SMALL,
            Failure::LengthMismatch { .. } => YBUS_ERROR_LENGTH_MISMATCH,
            Failure::InvalidValue(_) => YBUS_ERROR_INVALID_VALUE,
            Failure::InvalidField { .. } => YBUS_ERROR_INVALID_FIELD,
            Failure::MissingData(_) => YBUS_ERROR_MISSING_DATA,
            Failure::Panic(_) => YBUS_ERROR_PANIC,
        }
    }
}

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Failure::Powerflow(error) => write!(f, "{error}"),
            Failure::NullArgument(argument) => write!(f, "`{argument}` is null"),
            Failure::MisalignedArgument(argument) => {
                write!(f, "`{argument}` isn't aligned for its type")
            }
            Failure::InvalidUtf8(argument) => write!(f, "`{argument}` isn't valid UTF-8"),
            Failure::BufferTooSmall {
                argument,
                needed,
                given,
            } => write!(
                f,
                "`{argument}` has room for {given} elements, but the result needs {needed}"
            ),
            Failure::LengthMismatch {
                argument,
                expected,
                given,
            } => write!(
                f,
                "`{argument}` has {given} entries, but its table has {expected} rows"
            ),
            Failure::InvalidValue(reason) => write!(f, "{reason}"),
            Failure::InvalidField { table, field, kind } => {
                write!(f, "{field} isn't a {kind} field of the {table} table")
            }
            Failure::MissingData(what) => write!(f, "the case has no {what}"),
            Failure::Panic(message) => write!(f, "a bug in ybus (a panic): {message}"),
        }
    }
}

impl From<PowerflowError> for Failure {
    fn from(error: PowerflowError) -> Failure {
        Failure::Powerflow(error)
    }
}

thread_local! {
    /// The message of this thread's last failed call. Set on failure only, so it outlives the
    /// calls that succeed after it.
    static LAST_ERROR: RefCell<Option<CString>> = const { RefCell::new(None) };
}

/// Runs an export's body: `YBUS_OK` if it succeeds, and otherwise its failure's code, with the
/// message recorded as this thread's last error. A panic becomes `YBUS_ERROR_PANIC` rather than
/// reaching C, where it would abort the process.
pub(crate) fn guard(body: impl FnOnce() -> Result<(), Failure>) -> ybus_status {
    let failure = match catch_unwind(AssertUnwindSafe(body)) {
        Ok(Ok(())) => return YBUS_OK,
        Ok(Err(failure)) => failure,
        Err(payload) => Failure::Panic(panic_message(payload.as_ref())),
    };
    // Recording formats the message, which mustn't unwind into C either; if it panics, the
    //   code is still returned
    let _ = catch_unwind(AssertUnwindSafe(|| record(&failure)));
    failure.code()
}

/// A panic's message, from the `&str` or `String` that `panic!` carries.
fn panic_message(payload: &(dyn Any + Send)) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "a panic without a message".to_string()
    }
}

/// Records `failure`'s message as this thread's last error, with any NUL replaced so it can be
/// a C string. Does nothing while the thread is being torn down.
fn record(failure: &Failure) {
    let message = failure.to_string().replace('\0', "\u{fffd}");
    let message = CString::new(message).expect("NULs are replaced");
    let _ = LAST_ERROR.try_with(|last| *last.borrow_mut() = Some(message));
}

/// This thread's last error message, if a call on it has failed.
pub(crate) fn last_error() -> Option<CString> {
    LAST_ERROR
        .try_with(|last| last.borrow().clone())
        .ok()
        .flatten()
}

/// `NullArgument` or `MisalignedArgument` unless `pointer` can be dereferenced as far as its
/// value goes: whether it points to a valid `T` is the caller's promise.
fn check_pointer<T>(pointer: *const T, argument: &'static str) -> Result<(), Failure> {
    if pointer.is_null() {
        Err(Failure::NullArgument(argument))
    } else if !pointer.is_aligned() {
        Err(Failure::MisalignedArgument(argument))
    } else {
        Ok(())
    }
}

/// The UTF-8 string C passed as a NUL-terminated pointer.
///
/// # Safety
///
/// If `pointer` is non-null, it must point to a NUL-terminated string that lives, unchanged,
/// for `'a`.
pub(crate) unsafe fn arg_str<'a>(
    pointer: *const c_char,
    argument: &'static str,
) -> Result<&'a str, Failure> {
    if pointer.is_null() {
        return Err(Failure::NullArgument(argument));
    }
    // SAFETY: non-null, and NUL-terminated and valid for 'a by the caller's promise
    let string = unsafe { CStr::from_ptr(pointer) };
    string.to_str().map_err(|_| Failure::InvalidUtf8(argument))
}

/// Makes a new handle for C: `out` is set to null first (when it can be written at all), then
/// to the boxed result of `make` if it succeeds. The caller frees it with the handle's own
/// `ybus_*_free`.
///
/// # Safety
///
/// If `out` is non-null and aligned, it must point to a writable pointer.
pub(crate) unsafe fn new_handle<T>(
    out: *mut *mut T,
    make: impl FnOnce() -> Result<T, Failure>,
) -> Result<(), Failure> {
    if !out.is_null() && out.is_aligned() {
        // SAFETY: non-null and aligned (just checked), and writable by the caller's promise
        unsafe { out.write(std::ptr::null_mut()) };
    }
    // SAFETY: as above, by the caller's promise
    let out = unsafe { arg_mut(out, "out") }?;
    *out = Box::into_raw(Box::new(make()?));
    Ok(())
}

/// The value C passed by pointer.
///
/// # Safety
///
/// If `pointer` is non-null and aligned, it must point to a valid `T` that lives, unchanged,
/// for `'a`.
pub(crate) unsafe fn arg_ref<'a, T>(
    pointer: *const T,
    argument: &'static str,
) -> Result<&'a T, Failure> {
    check_pointer(pointer, argument)?;
    // SAFETY: non-null and aligned (just checked), and valid for 'a by the caller's promise
    Ok(unsafe { &*pointer })
}

/// The value C passed by pointer, to change.
///
/// # Safety
///
/// If `pointer` is non-null and aligned, it must point to a valid `T` that nothing else reads or
/// writes for `'a`.
pub(crate) unsafe fn arg_mut<'a, T>(
    pointer: *mut T,
    argument: &'static str,
) -> Result<&'a mut T, Failure> {
    check_pointer(pointer, argument)?;
    // SAFETY: non-null and aligned (just checked), and exclusive for 'a by the caller's promise
    Ok(unsafe { &mut *pointer })
}

/// The `len` values C passed as a pointer and a length. The pointer may be null when `len` is
/// 0.
///
/// # Safety
///
/// If `len` is non-zero and `pointer` non-null and aligned, it must point to `len` valid,
/// consecutive `T`s that live, unchanged, for `'a`.
pub(crate) unsafe fn arg_slice<'a, T>(
    pointer: *const T,
    len: usize,
    argument: &'static str,
) -> Result<&'a [T], Failure> {
    if len == 0 {
        return Ok(&[]);
    }
    check_pointer(pointer, argument)?;
    // SAFETY: non-null and aligned (just checked), and len valid values for 'a by the caller's
    //   promise
    Ok(unsafe { std::slice::from_raw_parts(pointer, len) })
}

/// The buffer of `len` values C passed as a pointer and a length, to write. The pointer may be
/// null when `len` is 0.
///
/// # Safety
///
/// If `len` is non-zero and `pointer` non-null and aligned, it must point to `len` valid,
/// consecutive `T`s that nothing else reads or writes for `'a`.
pub(crate) unsafe fn arg_slice_mut<'a, T>(
    pointer: *mut T,
    len: usize,
    argument: &'static str,
) -> Result<&'a mut [T], Failure> {
    if len == 0 {
        return Ok(&mut []);
    }
    check_pointer(pointer, argument)?;
    // SAFETY: non-null and aligned (just checked), and len exclusive values for 'a by the
    //   caller's promise
    Ok(unsafe { std::slice::from_raw_parts_mut(pointer, len) })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn last_error() -> Option<String> {
        super::last_error().map(|message| message.into_string().unwrap())
    }

    #[test]
    fn a_body_that_succeeds_is_ok() {
        assert_eq!(guard(|| Ok(())), YBUS_OK);
        assert_eq!(last_error(), None);
    }

    #[test]
    fn a_failure_is_its_code_and_message() {
        assert_eq!(
            guard(|| Err(Failure::NullArgument("case"))),
            YBUS_ERROR_NULL_ARGUMENT
        );
        assert_eq!(last_error().as_deref(), Some("`case` is null"));
    }

    #[test]
    fn a_panic_is_caught() {
        // The default panic hook still prints this to stderr; the guard leaves the hook alone
        assert_eq!(
            guard(|| panic!("deliberately, by the test")),
            YBUS_ERROR_PANIC
        );
        assert_eq!(
            last_error().as_deref(),
            Some("a bug in ybus (a panic): deliberately, by the test")
        );
        let row = 7;
        assert_eq!(guard(|| panic!("row {row}")), YBUS_ERROR_PANIC);
        assert_eq!(
            last_error().as_deref(),
            Some("a bug in ybus (a panic): row 7")
        );
    }

    #[test]
    fn the_last_error_outlives_a_success() {
        guard(|| Err(Failure::InvalidUtf8("path")));
        assert_eq!(guard(|| Ok(())), YBUS_OK);
        assert_eq!(last_error().as_deref(), Some("`path` isn't valid UTF-8"));
    }

    #[test]
    fn each_thread_has_its_own_last_error() {
        guard(|| Err(Failure::NullArgument("here")));
        let there = std::thread::spawn(|| {
            let before = last_error();
            guard(|| Err(Failure::NullArgument("there")));
            (before, last_error())
        })
        .join()
        .unwrap();
        assert_eq!(there, (None, Some("`there` is null".to_string())));
        assert_eq!(last_error().as_deref(), Some("`here` is null"));
    }

    #[test]
    fn a_nul_in_a_message_is_replaced() {
        guard(|| Err(Failure::Panic("a\0b".to_string())));
        assert_eq!(
            last_error().as_deref(),
            Some("a bug in ybus (a panic): a\u{fffd}b")
        );
    }

    #[test]
    fn pointers_are_checked_before_use() {
        let value = [1u64, 2];
        // SAFETY (for each call): the pointers are either rejected before any dereference, or
        //   point into `value`
        unsafe {
            assert_eq!(*arg_ref(value.as_ptr(), "value").unwrap(), 1);
            assert!(matches!(
                arg_ref(std::ptr::null::<u64>(), "value"),
                Err(Failure::NullArgument("value"))
            ));
            let misaligned = value.as_ptr().cast::<u8>().wrapping_add(1).cast::<u64>();
            assert!(matches!(
                arg_ref(misaligned, "value"),
                Err(Failure::MisalignedArgument("value"))
            ));
            assert_eq!(arg_slice(value.as_ptr(), 2, "values").unwrap(), &[1, 2]);
            assert!(arg_slice(std::ptr::null::<u64>(), 0, "values")
                .unwrap()
                .is_empty());
            assert!(matches!(
                arg_slice(std::ptr::null::<u64>(), 2, "values"),
                Err(Failure::NullArgument("values"))
            ));
            let mut buffer = [0u64; 2];
            arg_slice_mut(buffer.as_mut_ptr(), 2, "buffer").unwrap()[1] = 5;
            assert_eq!(buffer, [0, 5]);
            assert!(arg_slice_mut(std::ptr::null_mut::<u64>(), 0, "buffer")
                .unwrap()
                .is_empty());
            *arg_mut(buffer.as_mut_ptr(), "buffer").unwrap() = 3;
            assert_eq!(buffer, [3, 5]);
        }
    }

    #[test]
    fn a_powerflow_error_is_its_code_and_message() {
        assert_eq!(
            guard(|| Err(PowerflowError::NoBuses.into())),
            crate::status::YBUS_ERROR_NO_BUSES
        );
        assert_eq!(last_error(), Some(PowerflowError::NoBuses.to_string()));
    }
}
