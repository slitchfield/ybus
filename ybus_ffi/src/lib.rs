//! The C interface to [`powerflow`], built as the `ybus` library (`libybus.so`, `libybus.a`)
//! with the header `include/ybus.h`, which cbindgen generates from this crate.
//!
//! Every export is an `extern "C"` function prefixed `ybus_`. The ABI policy and the
//! thread-safety contract are in this crate's README; the plan is `docs/FfiChecklist.md`.
//!
//! With `ybus_cuda`, this is one of the two crates in the workspace with `unsafe` code: the
//! exports take raw pointers from C. Each `unsafe` block says why it's sound, and the rest of
//! the workspace keeps `#![forbid(unsafe_code)]`.

#![deny(unsafe_op_in_unsafe_fn)]

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

#[cfg(test)]
mod tests {
    use super::*;

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
}
