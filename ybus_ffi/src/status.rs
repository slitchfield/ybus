//! The status codes every fallible export returns, as `#define`s in the header. A code is never
//! renumbered once released; new failures get new codes.

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
