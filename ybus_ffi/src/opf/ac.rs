//! The AC optimal power flow's options through C.

use powerflow::opf::ac::AcOpfOptions;
use powerflow::opf::OpfSolver;

use super::{
    get_float, get_kkt, get_max_iterations, get_solver, set_float, set_kkt, set_max_iterations,
    ybus_kkt_factorization, ybus_opf_solver, Float, OpfOptions, YBUS_OPF_MEHROTRA, YBUS_OPF_MIPS,
};
use crate::boundary::{free_handle, guard, new_handle, Failure};
use crate::status::ybus_status;

/// Options for the AC optimal power flow, made by `ybus_ac_opf_options_new` at the defaults for a
/// solver and freed by `ybus_ac_opf_options_free`. A null options handle means MATPOWER's defaults
/// wherever one is taken.
pub struct YbusAcOpfOptions {
    pub(crate) options: AcOpfOptions,
}

impl OpfOptions for YbusAcOpfOptions {
    fn solver(&self) -> &OpfSolver {
        &self.options.solver
    }

    fn solver_mut(&mut self) -> &mut OpfSolver {
        &mut self.options.solver
    }
}

/// Makes AC OPF options for `solver` in `*out`: MATPOWER's defaults for `YBUS_OPF_MIPS` (MIPS, with
/// `feastol` 5e-6, and `cost_mult` 1e-4), or for `YBUS_OPF_MEHROTRA` the predictor-corrector
/// method's tuned settings, with the KKT system factored by the pivoted LDLᵀ. The solver can't
/// change afterwards. `*out` is set to `NULL` first and stays `NULL` on an error, such as
/// `YBUS_ERROR_INVALID_VALUE` for a code that isn't a solver.
///
/// # Safety
///
/// `out` must be null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_new(
    solver: ybus_opf_solver,
    out: *mut *mut YbusAcOpfOptions,
) -> ybus_status {
    guard(|| {
        // SAFETY: out is as the caller promised
        unsafe {
            new_handle(out, || {
                let options = match solver {
                    YBUS_OPF_MIPS => AcOpfOptions::default(),
                    YBUS_OPF_MEHROTRA => AcOpfOptions::mehrotra(),
                    other => {
                        return Err(Failure::InvalidValue(format!(
                            "{other} isn't an OPF solver, 1 or 2"
                        )))
                    }
                };
                Ok(YbusAcOpfOptions { options })
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
pub unsafe extern "C" fn ybus_ac_opf_options_free(options: *mut YbusAcOpfOptions) {
    // SAFETY: as the caller promised
    unsafe { free_handle(options) }
}

/// Writes the solver the options are for through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `ybus_opf_solver`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_solver(
    options: *const YbusAcOpfOptions,
    out: *mut ybus_opf_solver,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_solver(options, out) }
}

/// Sets the termination tolerance for feasibility, the solvers' `feastol`: 5e-6 by default,
/// MATPOWER's `opf.violation`. Returns `YBUS_ERROR_INVALID_VALUE`, changing nothing, unless it's
/// positive and finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_set_feastol(
    options: *mut YbusAcOpfOptions,
    value: f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_float(options, Float::Feastol, value) }
}

/// Writes the `feastol` setting through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_feastol(
    options: *const YbusAcOpfOptions,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_float(options, Float::Feastol, out) }
}

/// Sets the termination tolerance for the gradient, `gradtol`. Returns `YBUS_ERROR_INVALID_VALUE`,
/// changing nothing, unless it's positive and finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_set_gradtol(
    options: *mut YbusAcOpfOptions,
    value: f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_float(options, Float::Gradtol, value) }
}

/// Writes the `gradtol` setting through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_gradtol(
    options: *const YbusAcOpfOptions,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_float(options, Float::Gradtol, out) }
}

/// Sets the termination tolerance for complementarity, `comptol`. Returns
/// `YBUS_ERROR_INVALID_VALUE`, changing nothing, unless it's positive and finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_set_comptol(
    options: *mut YbusAcOpfOptions,
    value: f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_float(options, Float::Comptol, value) }
}

/// Writes the `comptol` setting through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_comptol(
    options: *const YbusAcOpfOptions,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_float(options, Float::Comptol, out) }
}

/// Sets the termination tolerance for the cost, `costtol`. Returns `YBUS_ERROR_INVALID_VALUE`,
/// changing nothing, unless it's positive and finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_set_costtol(
    options: *mut YbusAcOpfOptions,
    value: f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_float(options, Float::Costtol, value) }
}

/// Writes the `costtol` setting through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_costtol(
    options: *const YbusAcOpfOptions,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_float(options, Float::Costtol, out) }
}

/// Sets the objective's scale factor, `cost_mult`: 1e-4, MATPOWER's for the AC OPF by default.
/// Returns `YBUS_ERROR_INVALID_VALUE`, changing nothing, unless it's positive and finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_set_cost_mult(
    options: *mut YbusAcOpfOptions,
    value: f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_float(options, Float::CostMult, value) }
}

/// Writes the `cost_mult` setting through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_cost_mult(
    options: *const YbusAcOpfOptions,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_float(options, Float::CostMult, out) }
}

/// Sets the most iterations, the solvers' `max_it`.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_set_max_iterations(
    options: *mut YbusAcOpfOptions,
    value: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_max_iterations(options, value) }
}

/// Writes the iteration limit through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_max_iterations(
    options: *const YbusAcOpfOptions,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_max_iterations(options, out) }
}

/// Sets how the solver's KKT system is factored: `YBUS_KKT_LU` or `YBUS_KKT_LDLT`. Returns
/// `YBUS_ERROR_INVALID_VALUE` for a code that isn't one.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_set_kkt_factorization(
    options: *mut YbusAcOpfOptions,
    code: ybus_kkt_factorization,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { set_kkt(options, code) }
}

/// Writes how the solver's KKT system is factored through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable
/// `ybus_kkt_factorization`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_options_get_kkt_factorization(
    options: *const YbusAcOpfOptions,
    out: *mut ybus_kkt_factorization,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_kkt(options, out) }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use powerflow::opf::KktFactorization;

    use super::*;
    use crate::opf::{YBUS_KKT_LDLT, YBUS_KKT_LU};
    use crate::status::{YBUS_ERROR_INVALID_VALUE, YBUS_ERROR_NULL_ARGUMENT, YBUS_OK};

    fn new(solver: ybus_opf_solver) -> *mut YbusAcOpfOptions {
        let mut options = ptr::null_mut();
        // SAFETY: options is writable
        assert_eq!(
            unsafe { ybus_ac_opf_options_new(solver, &mut options) },
            YBUS_OK
        );
        options
    }

    #[test]
    fn opf_options_start_at_each_solvers_defaults() {
        for (code, expected) in [
            (YBUS_OPF_MIPS, AcOpfOptions::default()),
            (YBUS_OPF_MEHROTRA, AcOpfOptions::mehrotra()),
        ] {
            let options = new(code);
            let mut solver = 0;
            // SAFETY: options is live until freed, and solver writable
            unsafe {
                assert_eq!((*options).options, expected);
                assert_eq!(
                    ybus_ac_opf_options_get_solver(options, &mut solver),
                    YBUS_OK
                );
                assert_eq!(solver, code);
                ybus_ac_opf_options_free(options);
            }
        }
        let mut out = ptr::NonNull::dangling().as_ptr();
        // SAFETY: out is writable
        assert_eq!(
            unsafe { ybus_ac_opf_options_new(3, &mut out) },
            YBUS_ERROR_INVALID_VALUE
        );
        assert!(out.is_null());
    }

    #[test]
    fn opf_options_set_each_setting_on_either_solver() {
        for code in [YBUS_OPF_MIPS, YBUS_OPF_MEHROTRA] {
            let options = new(code);
            let (mut value, mut iterations, mut kkt) = (0.0, 0, 0);
            // SAFETY: options is live and unshared until freed, and every out writable
            unsafe {
                assert_eq!(ybus_ac_opf_options_set_feastol(options, 1e-7), YBUS_OK);
                assert_eq!(ybus_ac_opf_options_set_gradtol(options, 2e-7), YBUS_OK);
                assert_eq!(ybus_ac_opf_options_set_comptol(options, 3e-7), YBUS_OK);
                assert_eq!(ybus_ac_opf_options_set_costtol(options, 4e-7), YBUS_OK);
                assert_eq!(ybus_ac_opf_options_set_cost_mult(options, 0.5), YBUS_OK);
                assert_eq!(ybus_ac_opf_options_set_max_iterations(options, 77), YBUS_OK);
                assert_eq!(
                    ybus_ac_opf_options_set_kkt_factorization(options, YBUS_KKT_LDLT),
                    YBUS_OK
                );
                let mut expected = if code == YBUS_OPF_MIPS {
                    AcOpfOptions::default()
                } else {
                    AcOpfOptions::mehrotra()
                };
                match &mut expected.solver {
                    OpfSolver::Mips(o) => {
                        (o.feastol, o.gradtol, o.comptol, o.costtol) = (1e-7, 2e-7, 3e-7, 4e-7);
                        (o.cost_mult, o.max_it) = (0.5, 77);
                        o.kkt_factorization = KktFactorization::Ldlt;
                    }
                    OpfSolver::Mehrotra(o) => {
                        (o.feastol, o.gradtol, o.comptol, o.costtol) = (1e-7, 2e-7, 3e-7, 4e-7);
                        (o.cost_mult, o.max_it) = (0.5, 77);
                        o.kkt_factorization = KktFactorization::Ldlt;
                    }
                    _ => unreachable!(),
                }
                assert_eq!((*options).options, expected);

                assert_eq!(
                    ybus_ac_opf_options_get_comptol(options, &mut value),
                    YBUS_OK
                );
                assert_eq!(value, 3e-7);
                assert_eq!(
                    ybus_ac_opf_options_get_cost_mult(options, &mut value),
                    YBUS_OK
                );
                assert_eq!(value, 0.5);
                assert_eq!(
                    ybus_ac_opf_options_get_max_iterations(options, &mut iterations),
                    YBUS_OK
                );
                assert_eq!(iterations, 77);
                assert_eq!(
                    ybus_ac_opf_options_get_kkt_factorization(options, &mut kkt),
                    YBUS_OK
                );
                assert_eq!(kkt, YBUS_KKT_LDLT);

                for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                    assert_eq!(
                        ybus_ac_opf_options_set_feastol(options, bad),
                        YBUS_ERROR_INVALID_VALUE
                    );
                    assert_eq!(
                        ybus_ac_opf_options_set_cost_mult(options, bad),
                        YBUS_ERROR_INVALID_VALUE
                    );
                }
                assert_eq!(
                    ybus_ac_opf_options_set_kkt_factorization(options, 0),
                    YBUS_ERROR_INVALID_VALUE
                );
                assert_eq!((*options).options, expected);
                assert_eq!(
                    ybus_ac_opf_options_set_kkt_factorization(options, YBUS_KKT_LU),
                    YBUS_OK
                );
                assert_eq!(
                    ybus_ac_opf_options_get_feastol(ptr::null(), &mut value),
                    YBUS_ERROR_NULL_ARGUMENT
                );
                ybus_ac_opf_options_free(options);
                ybus_ac_opf_options_free(ptr::null_mut());
            }
        }
    }
}
