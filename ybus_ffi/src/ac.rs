//! The AC power flow through C: its options, its solution, the power flow itself, and the
//! model, built once from a case and solved (or re-stamped) for many scenarios.

use ybus_powerflow::ac::{
    ac_power_flow_with_options, AcMethod, AcModel, AcOptions, AcSolution, ReactiveLimits,
};

use crate::boundary::{
    arg_mut, arg_ref, arg_slice_mut, check_room, free_handle, get_value, guard, new_handle, Failure,
};
use crate::case::columns::{
    ybus_field, YBUS_BRANCH_PF_MW, YBUS_BRANCH_PT_MW, YBUS_BRANCH_QF_MVAR, YBUS_BRANCH_QT_MVAR,
    YBUS_BUS_VA_DEG, YBUS_BUS_VM_PU, YBUS_GEN_PG_MW, YBUS_GEN_QG_MVAR,
};
use crate::case::YbusCase;
use crate::solution::get_column;
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

/// The options `options` points to, or Newton's defaults for a null pointer.
///
/// # Safety
///
/// `options` must be null or live options.
pub(crate) unsafe fn options_or_default(
    options: *const YbusAcOptions,
) -> Result<AcOptions, Failure> {
    if options.is_null() {
        Ok(AcOptions::default())
    } else {
        // SAFETY: options is as the caller promised
        Ok(unsafe { arg_ref(options, "options") }?.options)
    }
}

/// An AC power flow's results, indexed like the case's tables, made by `ybus_ac_power_flow` or
/// `ybus_ac_model_solve` and freed by `ybus_ac_solution_free`.
pub struct YbusAcSolution {
    pub(crate) solution: AcSolution,
}

/// Solves the AC power flow for `case`, as MATPOWER's `runpf` does, with `options` (`NULL`
/// for Newton's defaults), into a new solution in `*out`. `*out` is set to `NULL` first and
/// stays `NULL` on an error, such as `YBUS_ERROR_NO_CONVERGENCE` or
/// `YBUS_ERROR_SINGULAR_NETWORK`.
///
/// # Safety
///
/// `case` must be null or a live case, `options` null or live options, and `out` null or a
/// writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_power_flow(
    case: *const YbusCase,
    options: *const YbusAcOptions,
    out: *mut *mut YbusAcSolution,
) -> ybus_status {
    guard(|| {
        // SAFETY: out, case and options are as the caller promised
        unsafe {
            new_handle(out, || {
                let case = arg_ref(case, "case")?;
                let options = options_or_default(options)?;
                Ok(YbusAcSolution {
                    solution: ac_power_flow_with_options(&case.case, &options)?,
                })
            })
        }
    })
}

/// Frees `solution`, which mustn't be used again. Does nothing if `solution` is null.
///
/// # Safety
///
/// `solution` must be null or a solution made by this library and not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_free(solution: *mut YbusAcSolution) {
    // SAFETY: as the caller promised
    unsafe { free_handle(solution) }
}

/// Writes the solution's bus column `field` into `out`, which has room for `len` entries: one
/// per bus. The columns are `YBUS_BUS_VM_PU` and `YBUS_BUS_VA_DEG`, each bus's voltage
/// magnitude (p.u.) and angle (degrees); an isolated bus keeps its voltage from the case.
/// Returns `YBUS_ERROR_INVALID_FIELD` for any other code, and `YBUS_ERROR_BUFFER_TOO_SMALL`,
/// writing nothing, if `len` is less than the number of buses.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_get_bus(
    solution: *const YbusAcSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "bus", field, out, len, |s, field| match field {
            YBUS_BUS_VM_PU => Some(&s.solution.vm_pu),
            YBUS_BUS_VA_DEG => Some(&s.solution.va_deg),
            _ => None,
        })
    }
}

/// [`ybus_ac_solution_get_bus`] for the generators: `YBUS_GEN_PG_MW` and `YBUS_GEN_QG_MVAR`,
/// each generator's real (MW) and reactive (MVAr) output, as MATPOWER's `pfsoln` computes
/// them. 0 for a generator that's out of service or at an isolated bus.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_get_gen(
    solution: *const YbusAcSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "gen", field, out, len, |s, field| match field {
            YBUS_GEN_PG_MW => Some(&s.solution.pg_mw),
            YBUS_GEN_QG_MVAR => Some(&s.solution.qg_mvar),
            _ => None,
        })
    }
}

/// [`ybus_ac_solution_get_bus`] for the branches: `YBUS_BRANCH_PF_MW`, `YBUS_BRANCH_QF_MVAR`,
/// `YBUS_BRANCH_PT_MW` and `YBUS_BRANCH_QT_MVAR`, the real (MW) and reactive (MVAr) power
/// flowing into each branch at its from and to ends. 0 for a branch that's out of service or
/// touches an isolated bus.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_get_branch(
    solution: *const YbusAcSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(
            solution,
            "branch",
            field,
            out,
            len,
            |s, field| match field {
                YBUS_BRANCH_PF_MW => Some(&s.solution.pf_mw),
                YBUS_BRANCH_QF_MVAR => Some(&s.solution.qf_mvar),
                YBUS_BRANCH_PT_MW => Some(&s.solution.pt_mw),
                YBUS_BRANCH_QT_MVAR => Some(&s.solution.qt_mvar),
                _ => None,
            },
        )
    }
}

/// Writes how many iterations were taken through `out`: Newton iterations, or for
/// fast-decoupled, P half-iterations, as MATPOWER counts them, summed over every solve when
/// reactive limits are enforced.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_iterations(
    solution: *const YbusAcSolution,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(solution, "solution", out, |s| Ok(s.solution.iterations)) }
}

/// Writes the largest real or reactive power mismatch at the solution (p.u.) through `out`.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_mismatch_pu(
    solution: *const YbusAcSolution,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(solution, "solution", out, |s| Ok(s.solution.mismatch_pu)) }
}

/// Writes how many matrices were factored afresh through `out`: for Newton's method, 1 when
/// every later iteration refactored the Jacobian with the first factorization's pivots; for
/// fast-decoupled, 2, B′ and B″; 0 when no iteration was needed.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_full_factorizations(
    solution: *const YbusAcSolution,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(solution, "solution", out, |s| {
            Ok(s.solution.full_factorizations)
        })
    }
}

/// Writes the (0-based) generator rows held at a reactive power limit into `out`, which has
/// room for `len`, in the order they were found, and how many there are through `count`. A
/// null `out` with `len` 0 only asks the count. Returns `YBUS_ERROR_BUFFER_TOO_SMALL`, writing
/// nothing into `out`, if they don't fit. None are held when reactive limits are ignored.
///
/// # Safety
///
/// `solution` must be null or a live solution, `out` null or `len` writable `size_t`s, and
/// `count` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_solution_limited_generators(
    solution: *const YbusAcSolution,
    out: *mut usize,
    len: usize,
    count: *mut usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: solution and count are as the caller promised
        let (solution, count) =
            unsafe { (arg_ref(solution, "solution")?, arg_mut(count, "count")?) };
        let limited = &solution.solution.limited_generators;
        *count = limited.len();
        if out.is_null() && len == 0 {
            return Ok(());
        }
        // SAFETY: out is as the caller promised
        let out = unsafe { arg_slice_mut(out, len, "out") }?;
        check_room("out", limited.len(), len)?;
        out[..limited.len()].copy_from_slice(limited);
        Ok(())
    })
}

/// An AC model of a case's network (Ybus, Yf, Yt and the bus roles), made by
/// `ybus_ac_model_new` and freed by `ybus_ac_model_free`, for solving many scenarios of the
/// same network.
pub struct YbusAcModel {
    pub(crate) model: AcModel,
}

/// Builds the AC model of `case`'s network into a new model in `*out`, which is set to `NULL`
/// first and stays `NULL` on an error, such as `YBUS_ERROR_ZERO_IMPEDANCE`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_model_new(
    case: *const YbusCase,
    out: *mut *mut YbusAcModel,
) -> ybus_status {
    guard(|| {
        // SAFETY: out and case are as the caller promised
        unsafe {
            new_handle(out, || {
                Ok(YbusAcModel {
                    model: AcModel::new(&arg_ref(case, "case")?.case)?,
                })
            })
        }
    })
}

/// Solves the AC power flow of `case` on `model`, with `options` (`NULL` for Newton's
/// defaults), into a new solution in `*out`: the same result as `ybus_ac_power_flow` on
/// `case`, without rebuilding the network.
///
/// `case` must be the case the model was built (or last re-stamped) from, or a copy with
/// different injections, setpoints or starting voltages (`pd_mw`, `qd_mvar`, `pg_mw`,
/// `qg_mvar`, `vg_pu`, `qmax_mvar`, `qmin_mvar`, `vm_pu` and `va_deg`). Returns
/// `YBUS_ERROR_CASE_MISMATCH` if its tables' lengths differ, and
/// `YBUS_ERROR_TOPOLOGY_MISMATCH` if its topology does. Changed impedances need
/// `ybus_ac_model_restamp` first. `*out` is set to `NULL` first and stays `NULL` on an error.
///
/// # Safety
///
/// `model` must be null or a live model, `case` null or a live case, `options` null or live
/// options, and `out` null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_model_solve(
    model: *const YbusAcModel,
    case: *const YbusCase,
    options: *const YbusAcOptions,
    out: *mut *mut YbusAcSolution,
) -> ybus_status {
    guard(|| {
        // SAFETY: out, model, case and options are as the caller promised
        unsafe {
            new_handle(out, || {
                let model = arg_ref(model, "model")?;
                let case = arg_ref(case, "case")?;
                let options = options_or_default(options)?;
                Ok(YbusAcSolution {
                    solution: model.model.solve(&case.case, &options)?,
                })
            })
        }
    })
}

/// Re-stamps `model`'s Ybus, Yf and Yt from `case`'s branch impedances, charging, taps and
/// phase shifts and its bus shunts, without rebuilding the matrices. `case` must have the
/// model's topology: `YBUS_ERROR_CASE_MISMATCH` or `YBUS_ERROR_TOPOLOGY_MISMATCH` otherwise,
/// and `YBUS_ERROR_ZERO_IMPEDANCE` for an in-service branch with `r = x = 0`. On an error the
/// model is unchanged.
///
/// # Safety
///
/// `model` must be null or a live model that nothing else is using, and `case` null or a live
/// case.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_model_restamp(
    model: *mut YbusAcModel,
    case: *const YbusCase,
) -> ybus_status {
    guard(|| {
        // SAFETY: model and case are as the caller promised
        let (model, case) = unsafe { (arg_mut(model, "model")?, arg_ref(case, "case")?) };
        Ok(model.model.restamp(&case.case)?)
    })
}

/// Frees `model`, which mustn't be used again. Does nothing if `model` is null.
///
/// # Safety
///
/// `model` must be null or a model made by this library and not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_model_free(model: *mut YbusAcModel) {
    // SAFETY: as the caller promised
    unsafe { free_handle(model) }
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

    #[test]
    fn the_ac_power_flow_is_powerflows_bit_for_bit() {
        // Each method, and Newton with each kind of reactive limits
        let runs = [
            (YBUS_AC_NEWTON, YBUS_REACTIVE_LIMITS_IGNORE),
            (YBUS_AC_FAST_DECOUPLED_XB, YBUS_REACTIVE_LIMITS_IGNORE),
            (YBUS_AC_FAST_DECOUPLED_BX, YBUS_REACTIVE_LIMITS_IGNORE),
            (YBUS_AC_NEWTON, YBUS_REACTIVE_LIMITS_ALL_AT_ONCE),
            (YBUS_AC_NEWTON, YBUS_REACTIVE_LIMITS_ONE_AT_A_TIME),
        ];
        for (name, case) in crate::dc::tests::cases() {
            for (method, limits) in runs {
                let options = new(method);
                let mut solution = ptr::null_mut();
                // SAFETY: each handle is live until freed, and each out writable
                unsafe {
                    assert_eq!(
                        ybus_ac_options_set_reactive_limits(options, limits),
                        YBUS_OK
                    );
                    let status = ybus_ac_power_flow(case, options, &mut solution);
                    let what = format!("{name} method {method} limits {limits}");
                    match ac_power_flow_with_options(&(*case).case, &(*options).options) {
                        Ok(expected) => {
                            assert_eq!(status, YBUS_OK, "{what}");
                            assert_eq!((*solution).solution, expected, "{what}");
                        }
                        Err(error) => {
                            assert_eq!(status, crate::status::powerflow_code(&error), "{what}");
                            assert!(solution.is_null());
                        }
                    }
                    ybus_ac_solution_free(solution);
                    ybus_ac_options_free(options);
                }
            }
            // SAFETY: case is live until this free
            unsafe { crate::case::ybus_case_free(case) };
        }
    }

    #[test]
    fn an_ac_solution_reads_its_columns_and_scalars() {
        let case = crate::case::tests::three_bus();
        let mut solution = ptr::null_mut();
        // SAFETY: each handle is live until freed, and each buffer has room for its length
        unsafe {
            assert_eq!(
                ybus_ac_power_flow(case, ptr::null(), &mut solution),
                YBUS_OK
            );
            let expected =
                ac_power_flow_with_options(&(*case).case, &AcOptions::default()).unwrap();
            let read = |get: unsafe extern "C" fn(
                *const YbusAcSolution,
                ybus_field,
                *mut f64,
                usize,
            ) -> ybus_status,
                        field,
                        rows: usize| {
                let mut out = vec![0.0; rows];
                assert_eq!(get(solution, field, out.as_mut_ptr(), rows), YBUS_OK);
                out
            };
            assert_eq!(
                read(ybus_ac_solution_get_bus, YBUS_BUS_VM_PU, 3),
                expected.vm_pu
            );
            assert_eq!(
                read(ybus_ac_solution_get_bus, YBUS_BUS_VA_DEG, 3),
                expected.va_deg
            );
            assert_eq!(
                read(ybus_ac_solution_get_gen, YBUS_GEN_PG_MW, 2),
                expected.pg_mw
            );
            assert_eq!(
                read(ybus_ac_solution_get_gen, YBUS_GEN_QG_MVAR, 2),
                expected.qg_mvar
            );
            assert_eq!(
                read(ybus_ac_solution_get_branch, YBUS_BRANCH_PF_MW, 3),
                expected.pf_mw
            );
            assert_eq!(
                read(ybus_ac_solution_get_branch, YBUS_BRANCH_QF_MVAR, 3),
                expected.qf_mvar
            );
            assert_eq!(
                read(ybus_ac_solution_get_branch, YBUS_BRANCH_PT_MW, 3),
                expected.pt_mw
            );
            assert_eq!(
                read(ybus_ac_solution_get_branch, YBUS_BRANCH_QT_MVAR, 3),
                expected.qt_mvar
            );
            let (mut iterations, mut mismatch, mut factorizations, mut count) = (0, 0.0, 0, 9);
            assert_eq!(
                ybus_ac_solution_iterations(solution, &mut iterations),
                YBUS_OK
            );
            assert_eq!(
                ybus_ac_solution_mismatch_pu(solution, &mut mismatch),
                YBUS_OK
            );
            assert_eq!(
                ybus_ac_solution_full_factorizations(solution, &mut factorizations),
                YBUS_OK
            );
            assert_eq!(
                ybus_ac_solution_limited_generators(solution, ptr::null_mut(), 0, &mut count),
                YBUS_OK
            );
            assert_eq!(
                (iterations, mismatch, factorizations, count),
                (
                    expected.iterations,
                    expected.mismatch_pu,
                    expected.full_factorizations,
                    0
                )
            );
            assert_eq!(
                ybus_ac_solution_get_gen(
                    solution,
                    crate::case::columns::YBUS_GEN_VG_PU,
                    ptr::null_mut(),
                    0
                ),
                crate::status::YBUS_ERROR_INVALID_FIELD
            );
            ybus_ac_solution_free(solution);
            crate::case::ybus_case_free(case);
        }
    }

    #[test]
    fn an_ac_model_solves_and_restamps_as_the_power_flow_does() {
        let case = crate::case::tests::three_bus();
        let (mut model, mut solution) = (ptr::null_mut(), ptr::null_mut());
        // The model's solve of the case as it is now, against the power flow's
        let check = move |model: *const YbusAcModel, solution: &mut *mut YbusAcSolution| {
            // SAFETY: model, case and solution are live and writable
            unsafe {
                assert_eq!(
                    ybus_ac_model_solve(model, case, ptr::null(), solution),
                    YBUS_OK
                );
                let expected = ac_power_flow_with_options(&(*case).case, &AcOptions::default());
                assert_eq!(Ok(&(**solution).solution), expected.as_ref());
                ybus_ac_solution_free(*solution);
            }
        };
        // SAFETY: each handle is live until freed, and each out writable
        unsafe {
            assert_eq!(ybus_ac_model_new(case, &mut model), YBUS_OK);
            for scale in [1.0, 0.8, 1.2] {
                for bus in &mut (*case).case.bus {
                    bus.pd_mw *= scale;
                    bus.qd_mvar *= scale;
                }
                check(model, &mut solution);
            }
            // A changed reactance, re-stamped
            (&mut (*case).case.branch)[1].x = 0.3;
            assert_eq!(ybus_ac_model_restamp(model, case), YBUS_OK);
            check(model, &mut solution);

            // A refused re-stamp leaves the model as it was
            let before = (*model).model.ybus().values().to_vec();
            (&mut (*case).case.branch)[0].status = false;
            (&mut (*case).case.branch)[1].x = 0.4;
            assert_eq!(
                ybus_ac_model_restamp(model, case),
                crate::status::YBUS_ERROR_TOPOLOGY_MISMATCH
            );
            assert_eq!((*model).model.ybus().values(), &before[..]);
            ybus_ac_model_free(model);
            ybus_ac_model_free(ptr::null_mut());
            crate::case::ybus_case_free(case);
        }
    }
}
