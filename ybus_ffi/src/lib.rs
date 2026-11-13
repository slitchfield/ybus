//! The C interface to [`ybus_powerflow`], built as the `ybus` library (`libybus.so`, `libybus.a`)
//! with the header `include/ybus.h`, which cbindgen generates from this crate.
//!
//! Every export is an `extern "C"` function prefixed `ybus_`. The ABI policy and the
//! thread-safety contract are in this crate's README; the plan is `docs/FfiChecklist.md`.
//!
//! With `ybus_cuda`, this is one of the two crates in the workspace with `unsafe` code: the
//! exports take raw pointers from C. Each `unsafe` block says why it's sound, and the rest of
//! the workspace keeps `#![forbid(unsafe_code)]`.

#![deny(unsafe_op_in_unsafe_fn)]

use std::ffi::c_char;
use std::panic::{catch_unwind, AssertUnwindSafe};

pub mod ac;
mod boundary;
pub mod case;
pub mod dc;
pub mod opf;
mod solution;
pub mod status;
#[cfg(test)]
mod testing;

/// The library's major version.
pub const YBUS_VERSION_MAJOR: u32 = 0;
/// The library's minor version.
pub const YBUS_VERSION_MINOR: u32 = 1;
/// The library's patch version.
pub const YBUS_VERSION_PATCH: u32 = 0;
/// The library's version as one number, major × 1 000 000 + minor × 1 000 + patch, so that
/// versions compare as numbers.
pub const YBUS_VERSION: u32 =
    YBUS_VERSION_MAJOR * 1_000_000 + YBUS_VERSION_MINOR * 1_000 + YBUS_VERSION_PATCH;

/// The version of the library linked at run time, encoded as `YBUS_VERSION` is. Comparing it
/// with the header's `YBUS_VERSION` shows whether the program was compiled against the same
/// version.
#[unsafe(no_mangle)]
pub extern "C" fn ybus_version() -> u32 {
    YBUS_VERSION
}

/// Copies the message of the calling thread's last failed call, a NUL-terminated UTF-8 string,
/// into `buffer`, which has room for `len` bytes, and returns the bytes the message needs,
/// counting the NUL; 0 if no call on this thread has failed.
///
/// Nothing is written unless `len` is at least that, so a message is never cut short: call it
/// with a null `buffer` and `len` 0 to learn the size, then with a buffer that size. Successful
/// calls leave the message alone, so it describes the most recent failure, not the most recent
/// call.
///
/// # Safety
///
/// Unless `buffer` is null, it must point to `len` bytes that can be written.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_last_error_message(buffer: *mut c_char, len: usize) -> usize {
    // A panic here returns 0 rather than reaching C; reading the message doesn't record one
    catch_unwind(AssertUnwindSafe(|| {
        let Some(message) = boundary::last_error() else {
            return 0;
        };
        let bytes = message.as_bytes_with_nul();
        if len >= bytes.len() {
            // SAFETY: by the caller's promise, a non-null buffer has len writable bytes; a null
            //   or misaligned one is refused, and then nothing is written
            if let Ok(out) = unsafe { boundary::arg_slice_mut(buffer.cast::<u8>(), len, "buffer") }
            {
                out[..bytes.len()].copy_from_slice(bytes);
            }
        }
        bytes.len()
    }))
    .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::{guard, Failure};

    // cbindgen can't evaluate env!, so the constants are literals, kept equal to the crate's
    //   version here
    #[test]
    fn the_version_constants_are_the_crates_version() {
        let part = |s: &str| s.parse::<u32>().unwrap();
        assert_eq!(
            (YBUS_VERSION_MAJOR, YBUS_VERSION_MINOR, YBUS_VERSION_PATCH),
            (
                part(env!("CARGO_PKG_VERSION_MAJOR")),
                part(env!("CARGO_PKG_VERSION_MINOR")),
                part(env!("CARGO_PKG_VERSION_PATCH"))
            )
        );
        assert_eq!(ybus_version(), 1_000);
    }

    /// The thread's last error message, read through the export
    fn read_message(len: usize) -> (usize, Vec<u8>) {
        let mut buffer = vec![0x7f_u8; len];
        // SAFETY: the buffer has len writable bytes
        let needed = unsafe { ybus_last_error_message(buffer.as_mut_ptr().cast(), len) };
        (needed, buffer)
    }

    #[test]
    fn a_thread_without_a_failure_has_no_message() {
        let fresh = std::thread::spawn(|| read_message(16)).join().unwrap();
        assert_eq!(fresh, (0, vec![0x7f; 16]));
    }

    #[test]
    fn the_last_error_message_is_copied_whole_or_not_at_all() {
        guard(|| Err(Failure::NullArgument("case")));
        let message = b"`case` is null\0";
        // SAFETY: a null buffer with length 0 is only a query
        let needed = unsafe { ybus_last_error_message(std::ptr::null_mut(), 0) };
        assert_eq!(needed, message.len());
        // One byte short: nothing is written
        assert_eq!(read_message(needed - 1), (needed, vec![0x7f; needed - 1]));
        // Exactly the size, and larger
        assert_eq!(read_message(needed), (needed, message.to_vec()));
        let (_, larger) = read_message(needed + 4);
        assert_eq!(&larger[..needed], message);
        assert_eq!(&larger[needed..], [0x7f; 4]);
        // A null buffer with a length is refused, without a write or a new message
        // SAFETY: a null buffer is refused before any write
        assert_eq!(
            unsafe { ybus_last_error_message(std::ptr::null_mut(), 64) },
            needed
        );
    }
}
