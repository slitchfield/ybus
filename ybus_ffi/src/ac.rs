//! The AC power flow through C: its options.

use powerflow::ac::{AcMethod, AcOptions, ReactiveLimits};

use crate::boundary::{arg_mut, free_handle, get_value, guard, new_handle, Failure};
use crate::status::ybus_status;

/// How to solve the AC power flow.
#[allow(non_camel_case_types)]
pub type ybus_ac_method = i32;

/// Newton's method in polar form, as MATPOWER's `newtonpf` (`pf.alg = 'NR'`, the default).
pub const YBUS_AC_NEWTON: ybus_ac_method = 1;
/// Fast-decoupled, XB version, as MATPOWER's `fdpf` (`pf.alg = 'FDXB'`): B′ ignores
/// resistance.
pub const YBUS_AC_FAST_DECOUPLED_XB: ybus_ac_method = 2;
/// Fast-decoupled, BX version (`pf.alg = 'FDBX'`): B″ ignores resistance.
pub const YBUS_AC_FAST_DECOUPLED_BX: ybus_ac_method = 3;

/// Whether and how to enforce generators' reactive power limits, with MATPOWER's
/// `pf.enforce_q_lims` codes.
#[allow(non_camel_case_types)]
pub type ybus_reactive_limits = i32;

/// Ignore them, MATPOWER's default.
pub const YBUS_REACTIVE_LIMITS_IGNORE: ybus_reactive_limits = 0;
/// After each solve, hold every generator outside its limits at the limit, make its bus PQ,
/// and solve again, until none is.
pub const YBUS_REACTIVE_LIMITS_ALL_AT_ONCE: ybus_reactive_limits = 1;
/// The same, but only the largest violation each round.
pub const YBUS_REACTIVE_LIMITS_ONE_AT_A_TIME: ybus_reactive_limits = 2;

fn method_code(method: AcMethod) -> Result<ybus_ac_method, Failure> {
    match method {
        AcMethod::Newton => Ok(YBUS_AC_NEWTON),
        AcMethod::FastDecoupledXb => Ok(YBUS_AC_FAST_DECOUPLED_XB),
        AcMethod::FastDecoupledBx => Ok(YBUS_AC_FAST_DECOUPLED_BX),
        // AcMethod is #[non_exhaustive]; only C's codes can set one here
        other => Err(Failure::InvalidValue(format!(
            "the method {other:?} has no code in this library"
        ))),
    }
}

fn method_of(code: ybus_ac_method) -> Result<AcMethod, Failure> {
    match code {
        YBUS_AC_NEWTON => Ok(AcMethod::Newton),
        YBUS_AC_FAST_DECOUPLED_XB => Ok(AcMethod::FastDecoupledXb),
        YBUS_AC_FAST_DECOUPLED_BX => Ok(AcMethod::FastDecoupledBx),
        other => Err(Failure::InvalidValue(format!(
            "{other} isn't an AC method, 1 to 3"
        ))),
    }
}

fn limits_code(limits: ReactiveLimits) -> Result<ybus_reactive_limits, Failure> {
    match limits {
        ReactiveLimits::Ignore => Ok(YBUS_REACTIVE_LIMITS_IGNORE),
        ReactiveLimits::AllAtOnce => Ok(YBUS_REACTIVE_LIMITS_ALL_AT_ONCE),
        ReactiveLimits::OneAtATime => Ok(YBUS_REACTIVE_LIMITS_ONE_AT_A_TIME),
        // ReactiveLimits is #[non_exhaustive]; only C's codes can set one here
        other => Err(Failure::InvalidValue(format!(
            "the reactive limits {other:?} have no code in this library"
        ))),
    }
}

fn limits_of(code: ybus_reactive_limits) -> Result<ReactiveLimits, Failure> {
    match code {
        YBUS_REACTIVE_LIMITS_IGNORE => Ok(ReactiveLimits::Ignore),
        YBUS_REACTIVE_LIMITS_ALL_AT_ONCE => Ok(ReactiveLimits::AllAtOnce),
        YBUS_REACTIVE_LIMITS_ONE_AT_A_TIME => Ok(ReactiveLimits::OneAtATime),
        other => Err(Failure::InvalidValue(format!(
            "{other} isn't a reactive limits setting, 0 to 2"
        ))),
    }
}

/// Options for the AC power flow, made by `ybus_ac_options_new` at MATPOWER's defaults for a
/// method and freed by `ybus_ac_options_free`. A null options handle means Newton's defaults
/// wherever one is taken.
pub struct YbusAcOptions {
    pub(crate) options: AcOptions,
}

/// Makes AC power flow options at MATPOWER's defaults for `method` in `*out`: tolerance 1e-8,
/// at most 10 Newton iterations or 30 fast-decoupled ones, and reactive limits ignored. `*out`
/// is set to `NULL` first and stays `NULL` on an error, such as `YBUS_ERROR_INVALID_VALUE` for
/// a code that isn't a method.
///
/// # Safety
///
/// `out` must be null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_new(
    method: ybus_ac_method,
    out: *mut *mut YbusAcOptions,
) -> ybus_status {
    guard(|| {
        // SAFETY: out is as the caller promised
        unsafe {
            new_handle(out, || {
                Ok(YbusAcOptions {
                    options: AcOptions::for_method(method_of(method)?),
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
pub unsafe extern "C" fn ybus_ac_options_free(options: *mut YbusAcOptions) {
    // SAFETY: as the caller promised
    unsafe { free_handle(options) }
}

/// Changes options under the guard
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
unsafe fn set(
    options: *mut YbusAcOptions,
    change: impl FnOnce(&mut AcOptions) -> Result<(), Failure>,
) -> ybus_status {
    guard(|| {
        // SAFETY: options is as the caller promised
        change(&mut unsafe { arg_mut(options, "options") }?.options)
    })
}

/// Sets the solution method. Only the method changes: the iteration limit stays as it was,
/// though fast-decoupled's default is 30 to Newton's 10. Returns `YBUS_ERROR_INVALID_VALUE` for
/// a code that isn't a method.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_set_method(
    options: *mut YbusAcOptions,
    code: ybus_ac_method,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        set(options, |options| {
            options.method = method_of(code)?;
            Ok(())
        })
    }
}

/// Writes the solution method through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `ybus_ac_method`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_get_method(
    options: *const YbusAcOptions,
    out: *mut ybus_ac_method,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(options, "options", out, |o| method_code(o.options.method)) }
}

/// Sets the tolerance: convergence when the largest real or reactive power mismatch is below
/// it (p.u.), MATPOWER's `pf.tol`. Returns `YBUS_ERROR_INVALID_VALUE` unless it's positive and
/// finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_set_tolerance(
    options: *mut YbusAcOptions,
    tolerance: f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        set(options, |options| {
            if !(tolerance.is_finite() && tolerance > 0.0) {
                return Err(Failure::InvalidValue(format!(
                    "the tolerance {tolerance} isn't positive and finite"
                )));
            }
            options.tolerance = tolerance;
            Ok(())
        })
    }
}

/// Writes the tolerance (p.u.) through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_get_tolerance(
    options: *const YbusAcOptions,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(options, "options", out, |o| Ok(o.options.tolerance)) }
}

/// Sets the most iterations before giving up: Newton iterations, or for fast-decoupled, P
/// half-iterations, MATPOWER's `pf.nr.max_it` or `pf.fd.max_it`.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_set_max_iterations(
    options: *mut YbusAcOptions,
    max_iterations: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        set(options, |options| {
            options.max_iterations = max_iterations;
            Ok(())
        })
    }
}

/// Writes the iteration limit through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_get_max_iterations(
    options: *const YbusAcOptions,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(options, "options", out, |o| Ok(o.options.max_iterations)) }
}

/// Sets whether and how generators' reactive power limits are enforced. Returns
/// `YBUS_ERROR_INVALID_VALUE` for a code that isn't a setting.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_set_reactive_limits(
    options: *mut YbusAcOptions,
    code: ybus_reactive_limits,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        set(options, |options| {
            options.reactive_limits = limits_of(code)?;
            Ok(())
        })
    }
}

/// Writes the reactive limits setting through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable
/// `ybus_reactive_limits`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_options_get_reactive_limits(
    options: *const YbusAcOptions,
    out: *mut ybus_reactive_limits,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(options, "options", out, |o| {
            limits_code(o.options.reactive_limits)
        })
    }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use super::*;
    use crate::status::{YBUS_ERROR_INVALID_VALUE, YBUS_ERROR_NULL_ARGUMENT, YBUS_OK};

    fn new(method: ybus_ac_method) -> *mut YbusAcOptions {
        let mut options = ptr::null_mut();
        // SAFETY: options is writable
        assert_eq!(
            unsafe { ybus_ac_options_new(method, &mut options) },
            YBUS_OK
        );
        options
    }

    #[test]
    fn ac_options_start_at_each_methods_defaults() {
        for (code, method) in [
            (YBUS_AC_NEWTON, AcMethod::Newton),
            (YBUS_AC_FAST_DECOUPLED_XB, AcMethod::FastDecoupledXb),
            (YBUS_AC_FAST_DECOUPLED_BX, AcMethod::FastDecoupledBx),
        ] {
            let options = new(code);
            // SAFETY: options is live until freed
            unsafe {
                assert_eq!((*options).options, AcOptions::for_method(method));
                ybus_ac_options_free(options);
            }
        }
        let mut out = ptr::NonNull::dangling().as_ptr();
        // SAFETY: out is writable
        unsafe {
            assert_eq!(ybus_ac_options_new(0, &mut out), YBUS_ERROR_INVALID_VALUE);
        }
        assert!(out.is_null());
    }

    #[test]
    fn ac_options_set_and_get_each_field() {
        let options = new(YBUS_AC_NEWTON);
        let (mut method, mut tolerance, mut iterations, mut limits) = (0, 0.0, 0, -1);
        // SAFETY: options is live and unshared until freed, and every out is writable
        unsafe {
            assert_eq!(
                ybus_ac_options_set_method(options, YBUS_AC_FAST_DECOUPLED_BX),
                YBUS_OK
            );
            assert_eq!(ybus_ac_options_set_tolerance(options, 1e-6), YBUS_OK);
            assert_eq!(ybus_ac_options_set_max_iterations(options, 25), YBUS_OK);
            assert_eq!(
                ybus_ac_options_set_reactive_limits(options, YBUS_REACTIVE_LIMITS_ONE_AT_A_TIME),
                YBUS_OK
            );
            let mut expected = AcOptions::for_method(AcMethod::FastDecoupledBx);
            expected.tolerance = 1e-6;
            expected.max_iterations = 25;
            expected.reactive_limits = ReactiveLimits::OneAtATime;
            assert_eq!((*options).options, expected);

            assert_eq!(ybus_ac_options_get_method(options, &mut method), YBUS_OK);
            assert_eq!(
                ybus_ac_options_get_tolerance(options, &mut tolerance),
                YBUS_OK
            );
            assert_eq!(
                ybus_ac_options_get_max_iterations(options, &mut iterations),
                YBUS_OK
            );
            assert_eq!(
                ybus_ac_options_get_reactive_limits(options, &mut limits),
                YBUS_OK
            );
            assert_eq!(
                (method, tolerance, iterations, limits),
                (
                    YBUS_AC_FAST_DECOUPLED_BX,
                    1e-6,
                    25,
                    YBUS_REACTIVE_LIMITS_ONE_AT_A_TIME
                )
            );

            for bad in [0.0, -1e-8, f64::NAN, f64::INFINITY] {
                assert_eq!(
                    ybus_ac_options_set_tolerance(options, bad),
                    YBUS_ERROR_INVALID_VALUE
                );
            }
            assert_eq!(
                ybus_ac_options_set_method(options, 4),
                YBUS_ERROR_INVALID_VALUE
            );
            assert_eq!(
                ybus_ac_options_set_reactive_limits(options, 3),
                YBUS_ERROR_INVALID_VALUE
            );
            assert_eq!((*options).options, expected);
            assert_eq!(
                ybus_ac_options_set_tolerance(ptr::null_mut(), 1e-8),
                YBUS_ERROR_NULL_ARGUMENT
            );
            ybus_ac_options_free(options);
            ybus_ac_options_free(ptr::null_mut());
        }
    }
}
