//! The DC optimal power flow through C: its options, its solution, and the OPF itself.

use ybus_powerflow::opf::dc::{dc_opf_with_options, DcOpfOptions, DcOpfSolution};
use ybus_powerflow::opf::OpfSolver;

use super::{
    get_float, get_kkt, get_max_iterations, get_solver, set_float, set_kkt, set_max_iterations,
    status_code, ybus_kkt_factorization, ybus_opf_solver, ybus_opf_status, Float, OpfOptions,
    YBUS_OPF_MEHROTRA, YBUS_OPF_MIPS,
};
use crate::boundary::{arg_ref, free_handle, get_value, guard, new_handle, Failure};
use crate::case::columns::{
    ybus_field, YBUS_BRANCH_MU_ANGMAX, YBUS_BRANCH_MU_ANGMIN, YBUS_BRANCH_MU_SF, YBUS_BRANCH_MU_ST,
    YBUS_BRANCH_PF_MW, YBUS_BUS_LAM_P, YBUS_BUS_VA_DEG, YBUS_GEN_MU_PMAX, YBUS_GEN_MU_PMIN,
    YBUS_GEN_PG_MW,
};
use crate::case::YbusCase;
use crate::solution::get_column;
use crate::status::ybus_status;

/// Options for the DC optimal power flow, made by `ybus_dc_opf_options_new` at the defaults for a
/// solver and freed by `ybus_dc_opf_options_free`. A null options handle means MATPOWER's defaults
/// wherever one is taken.
pub struct YbusDcOpfOptions {
    pub(crate) options: DcOpfOptions,
}

impl OpfOptions for YbusDcOpfOptions {
    fn solver(&self) -> &OpfSolver {
        &self.options.solver
    }

    fn solver_mut(&mut self) -> &mut OpfSolver {
        &mut self.options.solver
    }
}

/// Makes DC OPF options for `solver` in `*out`: MATPOWER's defaults for `YBUS_OPF_MIPS` (MIPS, with
/// `feastol` 5e-6), or for `YBUS_OPF_MEHROTRA` the predictor-corrector method's tuned settings,
/// with the KKT system factored by the pivoted LDLᵀ. The solver can't change afterwards. `*out` is
/// set to `NULL` first and stays `NULL` on an error, such as `YBUS_ERROR_INVALID_VALUE` for a code
/// that isn't a solver.
///
/// # Safety
///
/// `out` must be null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_options_new(
    solver: ybus_opf_solver,
    out: *mut *mut YbusDcOpfOptions,
) -> ybus_status {
    guard(|| {
        // SAFETY: out is as the caller promised
        unsafe {
            new_handle(out, || {
                let options = match solver {
                    YBUS_OPF_MIPS => DcOpfOptions::default(),
                    YBUS_OPF_MEHROTRA => DcOpfOptions::mehrotra(),
                    other => {
                        return Err(Failure::InvalidValue(format!(
                            "{other} isn't an OPF solver, 1 or 2"
                        )))
                    }
                };
                Ok(YbusDcOpfOptions { options })
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
pub unsafe extern "C" fn ybus_dc_opf_options_free(options: *mut YbusDcOpfOptions) {
    // SAFETY: as the caller promised
    unsafe { free_handle(options) }
}

/// Writes the solver the options are for through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `ybus_opf_solver`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_options_get_solver(
    options: *const YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_set_feastol(
    options: *mut YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_get_feastol(
    options: *const YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_set_gradtol(
    options: *mut YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_get_gradtol(
    options: *const YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_set_comptol(
    options: *mut YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_get_comptol(
    options: *const YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_set_costtol(
    options: *mut YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_get_costtol(
    options: *const YbusDcOpfOptions,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_float(options, Float::Costtol, out) }
}

/// Sets the objective's scale factor, `cost_mult`: 1 by default. Returns
/// `YBUS_ERROR_INVALID_VALUE`, changing nothing, unless it's positive and finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_options_set_cost_mult(
    options: *mut YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_get_cost_mult(
    options: *const YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_set_max_iterations(
    options: *mut YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_get_max_iterations(
    options: *const YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_set_kkt_factorization(
    options: *mut YbusDcOpfOptions,
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
pub unsafe extern "C" fn ybus_dc_opf_options_get_kkt_factorization(
    options: *const YbusDcOpfOptions,
    out: *mut ybus_kkt_factorization,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_kkt(options, out) }
}

/// A DC optimal power flow's results, indexed like the case's tables, as `rundcopf` reports
/// them, made by `ybus_dc_opf` and freed by `ybus_dc_opf_solution_free`. They're reported
/// whether or not the solver converged: see `ybus_dc_opf_solution_status`.
pub struct YbusDcOpfSolution {
    pub(crate) solution: DcOpfSolution,
}

/// Solves `case`'s DC optimal power flow, as MATPOWER's `rundcopf` does, with `options` (`NULL`
/// for MATPOWER's defaults), into a new solution in `*out`.
///
/// A solve that doesn't converge still makes a solution, whose status says so. Returns
/// `YBUS_ERROR_SOLVER` if the solver refuses to start, and the case's and costs' codes, such
/// as `YBUS_ERROR_UNSUPPORTED_COST` for a cost above quadratic. `*out` is set to `NULL` first
/// and stays `NULL` on an error.
///
/// # Safety
///
/// `case` must be null or a live case, `options` null or live options, and `out` null or a
/// writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf(
    case: *const YbusCase,
    options: *const YbusDcOpfOptions,
    out: *mut *mut YbusDcOpfSolution,
) -> ybus_status {
    guard(|| {
        // SAFETY: out, case and options are as the caller promised
        unsafe {
            new_handle(out, || {
                let case = arg_ref(case, "case")?;
                let options = if options.is_null() {
                    DcOpfOptions::default()
                } else {
                    arg_ref(options, "options")?.options
                };
                Ok(YbusDcOpfSolution {
                    solution: dc_opf_with_options(&case.case, &options)?,
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
pub unsafe extern "C" fn ybus_dc_opf_solution_free(solution: *mut YbusDcOpfSolution) {
    // SAFETY: as the caller promised
    unsafe { free_handle(solution) }
}

/// Writes the solution's bus column `field` into `out`, which has room for `len` entries: one
/// per bus. The columns are `YBUS_BUS_VA_DEG`, each bus's voltage angle (degrees), and
/// `YBUS_BUS_LAM_P`, its real power price ($/MWh); an isolated bus keeps its angle from the
/// case and has no price. Returns `YBUS_ERROR_INVALID_FIELD` for any other code, and
/// `YBUS_ERROR_BUFFER_TOO_SMALL`, writing nothing, if `len` is less than the number of buses.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_solution_get_bus(
    solution: *const YbusDcOpfSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "bus", field, out, len, |s, field| match field {
            YBUS_BUS_VA_DEG => Some(&s.solution.va_deg),
            YBUS_BUS_LAM_P => Some(&s.solution.lam_p),
            _ => None,
        })
    }
}

/// [`ybus_dc_opf_solution_get_bus`] for the generators: `YBUS_GEN_PG_MW`, each generator's
/// real power output (MW), and `YBUS_GEN_MU_PMIN` and `YBUS_GEN_MU_PMAX`, the multipliers of
/// its lower and upper real power limits ($/MWh). 0 for a generator that's out of service or
/// at an isolated bus.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_solution_get_gen(
    solution: *const YbusDcOpfSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "gen", field, out, len, |s, field| match field {
            YBUS_GEN_PG_MW => Some(&s.solution.pg_mw),
            YBUS_GEN_MU_PMIN => Some(&s.solution.mu_pmin),
            YBUS_GEN_MU_PMAX => Some(&s.solution.mu_pmax),
            _ => None,
        })
    }
}

/// [`ybus_dc_opf_solution_get_bus`] for the branches: `YBUS_BRANCH_PF_MW`, each branch's real
/// power flow at the from end (MW); `YBUS_BRANCH_MU_SF` and `YBUS_BRANCH_MU_ST`, the
/// multipliers of its flow limit in each direction ($/MWh); and `YBUS_BRANCH_MU_ANGMIN` and
/// `YBUS_BRANCH_MU_ANGMAX`, those of its angle difference limits ($/degree).
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_solution_get_branch(
    solution: *const YbusDcOpfSolution,
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
                YBUS_BRANCH_MU_SF => Some(&s.solution.mu_sf),
                YBUS_BRANCH_MU_ST => Some(&s.solution.mu_st),
                YBUS_BRANCH_MU_ANGMIN => Some(&s.solution.mu_angmin),
                YBUS_BRANCH_MU_ANGMAX => Some(&s.solution.mu_angmax),
                _ => None,
            },
        )
    }
}

/// Writes the total cost ($/hr), including the costs' constant terms, through `out`.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_solution_objective(
    solution: *const YbusDcOpfSolution,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(solution, "solution", out, |s| Ok(s.solution.objective)) }
}

/// Writes how many iterations the solver took through `out`.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable `size_t`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_solution_iterations(
    solution: *const YbusDcOpfSolution,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe { get_value(solution, "solution", out, |s| Ok(s.solution.iterations)) }
}

/// Writes how the solver ended through `out`: `YBUS_OPF_CONVERGED`, MATPOWER's `success`, or
/// `YBUS_OPF_MAX_ITERATIONS` or `YBUS_OPF_NUMERICAL_FAILURE`.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable
/// `ybus_opf_status`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_opf_solution_status(
    solution: *const YbusDcOpfSolution,
    out: *mut ybus_opf_status,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(solution, "solution", out, |s| {
            Ok(status_code(s.solution.status))
        })
    }
}

#[cfg(test)]
mod tests {
    use std::ptr;

    use ybus_powerflow::opf::KktFactorization;

    use super::*;
    use crate::opf::{YBUS_KKT_LDLT, YBUS_KKT_LU};
    use crate::status::{YBUS_ERROR_INVALID_VALUE, YBUS_ERROR_NULL_ARGUMENT, YBUS_OK};

    fn new(solver: ybus_opf_solver) -> *mut YbusDcOpfOptions {
        let mut options = ptr::null_mut();
        // SAFETY: options is writable
        assert_eq!(
            unsafe { ybus_dc_opf_options_new(solver, &mut options) },
            YBUS_OK
        );
        options
    }

    #[test]
    fn opf_options_start_at_each_solvers_defaults() {
        for (code, expected) in [
            (YBUS_OPF_MIPS, DcOpfOptions::default()),
            (YBUS_OPF_MEHROTRA, DcOpfOptions::mehrotra()),
        ] {
            let options = new(code);
            let mut solver = 0;
            // SAFETY: options is live until freed, and solver writable
            unsafe {
                assert_eq!((*options).options, expected);
                assert_eq!(
                    ybus_dc_opf_options_get_solver(options, &mut solver),
                    YBUS_OK
                );
                assert_eq!(solver, code);
                ybus_dc_opf_options_free(options);
            }
        }
        let mut out = ptr::NonNull::dangling().as_ptr();
        // SAFETY: out is writable
        assert_eq!(
            unsafe { ybus_dc_opf_options_new(3, &mut out) },
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
                assert_eq!(ybus_dc_opf_options_set_feastol(options, 1e-7), YBUS_OK);
                assert_eq!(ybus_dc_opf_options_set_gradtol(options, 2e-7), YBUS_OK);
                assert_eq!(ybus_dc_opf_options_set_comptol(options, 3e-7), YBUS_OK);
                assert_eq!(ybus_dc_opf_options_set_costtol(options, 4e-7), YBUS_OK);
                assert_eq!(ybus_dc_opf_options_set_cost_mult(options, 0.5), YBUS_OK);
                assert_eq!(ybus_dc_opf_options_set_max_iterations(options, 77), YBUS_OK);
                assert_eq!(
                    ybus_dc_opf_options_set_kkt_factorization(options, YBUS_KKT_LDLT),
                    YBUS_OK
                );
                let mut expected = if code == YBUS_OPF_MIPS {
                    DcOpfOptions::default()
                } else {
                    DcOpfOptions::mehrotra()
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
                    ybus_dc_opf_options_get_comptol(options, &mut value),
                    YBUS_OK
                );
                assert_eq!(value, 3e-7);
                assert_eq!(
                    ybus_dc_opf_options_get_cost_mult(options, &mut value),
                    YBUS_OK
                );
                assert_eq!(value, 0.5);
                assert_eq!(
                    ybus_dc_opf_options_get_max_iterations(options, &mut iterations),
                    YBUS_OK
                );
                assert_eq!(iterations, 77);
                assert_eq!(
                    ybus_dc_opf_options_get_kkt_factorization(options, &mut kkt),
                    YBUS_OK
                );
                assert_eq!(kkt, YBUS_KKT_LDLT);

                for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                    assert_eq!(
                        ybus_dc_opf_options_set_feastol(options, bad),
                        YBUS_ERROR_INVALID_VALUE
                    );
                    assert_eq!(
                        ybus_dc_opf_options_set_cost_mult(options, bad),
                        YBUS_ERROR_INVALID_VALUE
                    );
                }
                assert_eq!(
                    ybus_dc_opf_options_set_kkt_factorization(options, 0),
                    YBUS_ERROR_INVALID_VALUE
                );
                assert_eq!((*options).options, expected);
                assert_eq!(
                    ybus_dc_opf_options_set_kkt_factorization(options, YBUS_KKT_LU),
                    YBUS_OK
                );
                assert_eq!(
                    ybus_dc_opf_options_get_feastol(ptr::null(), &mut value),
                    YBUS_ERROR_NULL_ARGUMENT
                );
                ybus_dc_opf_options_free(options);
                ybus_dc_opf_options_free(ptr::null_mut());
            }
        }
    }

    #[test]
    fn the_dc_opf_is_powerflows_bit_for_bit() {
        for (name, case) in crate::dc::tests::cases() {
            for solver in [YBUS_OPF_MIPS, YBUS_OPF_MEHROTRA] {
                for kkt in [YBUS_KKT_LU, YBUS_KKT_LDLT] {
                    let options = new(solver);
                    let mut solution = ptr::null_mut();
                    let what = format!("{name} solver {solver} KKT {kkt}");
                    // SAFETY: each handle is live until freed, and each out writable
                    unsafe {
                        assert_eq!(
                            ybus_dc_opf_options_set_kkt_factorization(options, kkt),
                            YBUS_OK
                        );
                        let status = ybus_dc_opf(case, options, &mut solution);
                        match dc_opf_with_options(&(*case).case, &(*options).options) {
                            Ok(expected) => {
                                assert_eq!(status, YBUS_OK, "{what}");
                                assert_eq!((*solution).solution, expected, "{what}");
                            }
                            Err(error) => {
                                assert_eq!(status, crate::status::powerflow_code(&error), "{what}");
                                assert!(solution.is_null());
                            }
                        }
                        ybus_dc_opf_solution_free(solution);
                        ybus_dc_opf_options_free(options);
                    }
                }
            }
            // SAFETY: case is live until this free
            unsafe { crate::case::ybus_case_free(case) };
        }
    }

    #[test]
    fn a_dc_opf_solution_reads_its_columns_and_scalars() {
        let case = crate::case::tests::three_bus();
        let mut solution = ptr::null_mut();
        // SAFETY: each handle is live until freed, and each buffer has room for its length
        unsafe {
            assert_eq!(ybus_dc_opf(case, ptr::null(), &mut solution), YBUS_OK);
            let expected = dc_opf_with_options(&(*case).case, &DcOpfOptions::default()).unwrap();
            type Get = unsafe extern "C" fn(
                *const YbusDcOpfSolution,
                ybus_field,
                *mut f64,
                usize,
            ) -> ybus_status;
            let read = |get: Get, field, rows: usize| {
                let mut out = vec![0.0; rows];
                assert_eq!(get(solution, field, out.as_mut_ptr(), rows), YBUS_OK);
                out
            };
            assert_eq!(
                read(ybus_dc_opf_solution_get_bus, YBUS_BUS_VA_DEG, 3),
                expected.va_deg
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_bus, YBUS_BUS_LAM_P, 3),
                expected.lam_p
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_gen, YBUS_GEN_PG_MW, 2),
                expected.pg_mw
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_gen, YBUS_GEN_MU_PMIN, 2),
                expected.mu_pmin
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_gen, YBUS_GEN_MU_PMAX, 2),
                expected.mu_pmax
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_branch, YBUS_BRANCH_PF_MW, 3),
                expected.pf_mw
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_branch, YBUS_BRANCH_MU_SF, 3),
                expected.mu_sf
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_branch, YBUS_BRANCH_MU_ST, 3),
                expected.mu_st
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_branch, YBUS_BRANCH_MU_ANGMIN, 3),
                expected.mu_angmin
            );
            assert_eq!(
                read(ybus_dc_opf_solution_get_branch, YBUS_BRANCH_MU_ANGMAX, 3),
                expected.mu_angmax
            );
            let (mut objective, mut iterations, mut status) = (0.0, 0, 9);
            assert_eq!(
                ybus_dc_opf_solution_objective(solution, &mut objective),
                YBUS_OK
            );
            assert_eq!(
                ybus_dc_opf_solution_iterations(solution, &mut iterations),
                YBUS_OK
            );
            assert_eq!(ybus_dc_opf_solution_status(solution, &mut status), YBUS_OK);
            assert_eq!(
                (objective, iterations),
                (expected.objective, expected.iterations)
            );
            assert_eq!(status, crate::opf::YBUS_OPF_CONVERGED);
            assert_eq!(
                ybus_dc_opf_solution_get_bus(
                    solution,
                    crate::case::columns::YBUS_BUS_VM_PU,
                    ptr::null_mut(),
                    0
                ),
                crate::status::YBUS_ERROR_INVALID_FIELD
            );
            ybus_dc_opf_solution_free(solution);
            crate::case::ybus_case_free(case);
        }
    }
}
