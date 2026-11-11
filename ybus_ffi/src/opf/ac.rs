//! The AC optimal power flow through C: its options, its solution, and the OPF itself.

use powerflow::opf::ac::{ac_opf_with_options, AcOpfOptions, AcOpfSolution};
use powerflow::opf::OpfSolver;

use super::{
    get_float, get_kkt, get_max_iterations, get_solver, set_float, set_kkt, set_max_iterations,
    status_code, ybus_kkt_factorization, ybus_opf_solver, ybus_opf_status, Float, OpfOptions,
    YBUS_OPF_MEHROTRA, YBUS_OPF_MIPS,
};
use crate::boundary::{arg_ref, free_handle, get_value, guard, new_handle, Failure};
use crate::case::columns::{
    ybus_field, YBUS_BRANCH_MU_ANGMAX, YBUS_BRANCH_MU_ANGMIN, YBUS_BRANCH_MU_SF, YBUS_BRANCH_MU_ST,
    YBUS_BRANCH_PF_MW, YBUS_BRANCH_PT_MW, YBUS_BRANCH_QF_MVAR, YBUS_BRANCH_QT_MVAR, YBUS_BUS_LAM_P,
    YBUS_BUS_LAM_Q, YBUS_BUS_MU_VMAX, YBUS_BUS_MU_VMIN, YBUS_BUS_VA_DEG, YBUS_BUS_VM_PU,
    YBUS_GEN_MU_PMAX, YBUS_GEN_MU_PMIN, YBUS_GEN_MU_QMAX, YBUS_GEN_MU_QMIN, YBUS_GEN_PG_MW,
    YBUS_GEN_QG_MVAR, YBUS_GEN_VG_PU,
};
use crate::case::YbusCase;
use crate::solution::get_column;
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

/// An AC optimal power flow's results, indexed like the case's tables, as `runopf` reports
/// them, made by `ybus_ac_opf` and freed by `ybus_ac_opf_solution_free`. They're reported
/// whether or not the solver converged: see `ybus_ac_opf_solution_status`.
pub struct YbusAcOpfSolution {
    pub(crate) solution: AcOpfSolution,
}

/// Solves `case`'s AC optimal power flow, as MATPOWER's `runopf` does, with `options` (`NULL`
/// for MATPOWER's defaults), into a new solution in `*out`.
///
/// A solve that doesn't converge still makes a solution, whose status says so. Returns
/// `YBUS_ERROR_SOLVER` if the solver refuses to start, `YBUS_ERROR_UNSUPPORTED_GENERATOR` for
/// a generator with a PQ capability curve, and the case's and costs' codes. `*out` is set to
/// `NULL` first and stays `NULL` on an error.
///
/// # Safety
///
/// `case` must be null or a live case, `options` null or live options, and `out` null or a
/// writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf(
    case: *const YbusCase,
    options: *const YbusAcOpfOptions,
    out: *mut *mut YbusAcOpfSolution,
) -> ybus_status {
    guard(|| {
        // SAFETY: out, case and options are as the caller promised
        unsafe {
            new_handle(out, || {
                let case = arg_ref(case, "case")?;
                let options = if options.is_null() {
                    AcOpfOptions::default()
                } else {
                    arg_ref(options, "options")?.options
                };
                Ok(YbusAcOpfSolution {
                    solution: ac_opf_with_options(&case.case, &options)?,
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
pub unsafe extern "C" fn ybus_ac_opf_solution_free(solution: *mut YbusAcOpfSolution) {
    // SAFETY: as the caller promised
    unsafe { free_handle(solution) }
}

/// Writes the solution's bus column `field` into `out`, which has room for `len` entries: one
/// per bus. The columns are `YBUS_BUS_VM_PU` and `YBUS_BUS_VA_DEG`, each bus's voltage; its
/// prices `YBUS_BUS_LAM_P` ($/MWh) and `YBUS_BUS_LAM_Q` ($/MVAr-hr); and the multipliers of its
/// voltage magnitude limits, `YBUS_BUS_MU_VMIN` and `YBUS_BUS_MU_VMAX` ($/p.u.-hr). An isolated
/// bus keeps its voltage from the case and has no prices. Returns `YBUS_ERROR_INVALID_FIELD`
/// for any other code, and `YBUS_ERROR_BUFFER_TOO_SMALL`, writing nothing, if `len` is less
/// than the number of buses.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_solution_get_bus(
    solution: *const YbusAcOpfSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "bus", field, out, len, |s, field| match field {
            YBUS_BUS_VM_PU => Some(&s.solution.vm_pu),
            YBUS_BUS_VA_DEG => Some(&s.solution.va_deg),
            YBUS_BUS_LAM_P => Some(&s.solution.lam_p),
            YBUS_BUS_LAM_Q => Some(&s.solution.lam_q),
            YBUS_BUS_MU_VMAX => Some(&s.solution.mu_vmax),
            YBUS_BUS_MU_VMIN => Some(&s.solution.mu_vmin),
            _ => None,
        })
    }
}

/// [`ybus_ac_opf_solution_get_bus`] for the generators: `YBUS_GEN_PG_MW` and
/// `YBUS_GEN_QG_MVAR`, each generator's output; `YBUS_GEN_VG_PU`, its voltage setpoint, its
/// bus's voltage magnitude; and the multipliers of its real and reactive power limits,
/// `YBUS_GEN_MU_PMIN` and `YBUS_GEN_MU_PMAX` ($/MWh), and `YBUS_GEN_MU_QMIN` and
/// `YBUS_GEN_MU_QMAX` ($/MVAr-hr). Out-of-service generators, and those at an isolated bus,
/// report 0 for each but `VG_PU`, which keeps the case's.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_solution_get_gen(
    solution: *const YbusAcOpfSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "gen", field, out, len, |s, field| match field {
            YBUS_GEN_PG_MW => Some(&s.solution.pg_mw),
            YBUS_GEN_QG_MVAR => Some(&s.solution.qg_mvar),
            YBUS_GEN_VG_PU => Some(&s.solution.vg_pu),
            YBUS_GEN_MU_PMAX => Some(&s.solution.mu_pmax),
            YBUS_GEN_MU_PMIN => Some(&s.solution.mu_pmin),
            YBUS_GEN_MU_QMAX => Some(&s.solution.mu_qmax),
            YBUS_GEN_MU_QMIN => Some(&s.solution.mu_qmin),
            _ => None,
        })
    }
}

/// [`ybus_ac_opf_solution_get_bus`] for the branches: the real and reactive power flowing into
/// each branch at its ends, `YBUS_BRANCH_PF_MW`, `YBUS_BRANCH_QF_MVAR`, `YBUS_BRANCH_PT_MW` and
/// `YBUS_BRANCH_QT_MVAR`; the multipliers of its apparent power limit at each end,
/// `YBUS_BRANCH_MU_SF` and `YBUS_BRANCH_MU_ST` ($/MVA-hr); and those of its angle difference
/// limits, `YBUS_BRANCH_MU_ANGMIN` and `YBUS_BRANCH_MU_ANGMAX` ($/degree).
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_solution_get_branch(
    solution: *const YbusAcOpfSolution,
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
                YBUS_BRANCH_MU_SF => Some(&s.solution.mu_sf),
                YBUS_BRANCH_MU_ST => Some(&s.solution.mu_st),
                YBUS_BRANCH_MU_ANGMIN => Some(&s.solution.mu_angmin),
                YBUS_BRANCH_MU_ANGMAX => Some(&s.solution.mu_angmax),
                _ => None,
            },
        )
    }
}

/// Writes the total cost ($/hr) through `out`.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable `double`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_ac_opf_solution_objective(
    solution: *const YbusAcOpfSolution,
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
pub unsafe extern "C" fn ybus_ac_opf_solution_iterations(
    solution: *const YbusAcOpfSolution,
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
pub unsafe extern "C" fn ybus_ac_opf_solution_status(
    solution: *const YbusAcOpfSolution,
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

    #[test]
    fn the_ac_opf_is_powerflows_bit_for_bit() {
        for (name, case) in crate::dc::tests::cases() {
            for solver in [YBUS_OPF_MIPS, YBUS_OPF_MEHROTRA] {
                for kkt in [YBUS_KKT_LU, YBUS_KKT_LDLT] {
                    let options = new(solver);
                    let mut solution = ptr::null_mut();
                    let what = format!("{name} solver {solver} KKT {kkt}");
                    // SAFETY: each handle is live until freed, and each out writable
                    unsafe {
                        assert_eq!(
                            ybus_ac_opf_options_set_kkt_factorization(options, kkt),
                            YBUS_OK
                        );
                        let status = ybus_ac_opf(case, options, &mut solution);
                        match ac_opf_with_options(&(*case).case, &(*options).options) {
                            Ok(expected) => {
                                assert_eq!(status, YBUS_OK, "{what}");
                                assert_eq!((*solution).solution, expected, "{what}");
                            }
                            Err(error) => {
                                assert_eq!(status, crate::status::powerflow_code(&error), "{what}");
                                assert!(solution.is_null());
                            }
                        }
                        ybus_ac_opf_solution_free(solution);
                        ybus_ac_opf_options_free(options);
                    }
                }
            }
            // SAFETY: case is live until this free
            unsafe { crate::case::ybus_case_free(case) };
        }
    }

    #[test]
    fn three_bus_has_a_capability_curve_the_ac_opf_refuses_until_it_is_zeroed() {
        let case = crate::case::tests::three_bus();
        let mut solution = ptr::null_mut();
        // SAFETY: each handle is live until freed, and each buffer has room for its length
        unsafe {
            assert_eq!(
                ybus_ac_opf(case, ptr::null(), &mut solution),
                crate::status::YBUS_ERROR_UNSUPPORTED_GENERATOR
            );
            for generator in &mut (*case).case.generator {
                generator.pc1_mw = 0.0;
                generator.pc2_mw = 0.0;
            }
            assert_eq!(ybus_ac_opf(case, ptr::null(), &mut solution), YBUS_OK);
            let expected = ac_opf_with_options(&(*case).case, &AcOpfOptions::default()).unwrap();
            type Get = unsafe extern "C" fn(
                *const YbusAcOpfSolution,
                ybus_field,
                *mut f64,
                usize,
            ) -> ybus_status;
            let read = |get: Get, field, rows: usize| {
                let mut out = vec![0.0; rows];
                assert_eq!(get(solution, field, out.as_mut_ptr(), rows), YBUS_OK);
                out
            };
            let s = &(*solution).solution;
            assert_eq!(
                read(ybus_ac_opf_solution_get_bus, YBUS_BUS_VM_PU, 3),
                s.vm_pu
            );
            assert_eq!(
                read(ybus_ac_opf_solution_get_bus, YBUS_BUS_LAM_Q, 3),
                s.lam_q
            );
            assert_eq!(
                read(ybus_ac_opf_solution_get_bus, YBUS_BUS_MU_VMAX, 3),
                s.mu_vmax
            );
            assert_eq!(
                read(ybus_ac_opf_solution_get_gen, YBUS_GEN_VG_PU, 2),
                s.vg_pu
            );
            assert_eq!(
                read(ybus_ac_opf_solution_get_gen, YBUS_GEN_MU_QMIN, 2),
                s.mu_qmin
            );
            assert_eq!(
                read(ybus_ac_opf_solution_get_branch, YBUS_BRANCH_QT_MVAR, 3),
                s.qt_mvar
            );
            assert_eq!(
                read(ybus_ac_opf_solution_get_branch, YBUS_BRANCH_MU_ST, 3),
                s.mu_st
            );
            assert_eq!(*s, expected);
            let (mut objective, mut status) = (0.0, 9);
            assert_eq!(
                ybus_ac_opf_solution_objective(solution, &mut objective),
                YBUS_OK
            );
            assert_eq!(ybus_ac_opf_solution_status(solution, &mut status), YBUS_OK);
            assert_eq!(
                (objective, status),
                (expected.objective, crate::opf::status_code(expected.status))
            );
            ybus_ac_opf_solution_free(solution);
            crate::case::ybus_case_free(case);
        }
    }
}
