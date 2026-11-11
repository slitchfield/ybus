//! The DC power flow through C: its options, its solution, the power flow itself, and the
//! model, built once from a case and solved for many scenarios.

use powerflow::dc::{dc_power_flow_with_options, DcModel, DcOptions, DcSolution, DcSolver};
use sparse_mat::compressed_matrix::Ordering;

use crate::boundary::{arg_mut, arg_ref, free_handle, get_value, guard, new_handle, Failure};
use crate::case::columns::{ybus_field, YBUS_BRANCH_PF_MW, YBUS_BUS_VA_DEG, YBUS_GEN_PG_MW};
use crate::case::YbusCase;
use crate::solution::get_column;
use crate::status::ybus_status;

/// A fill-reducing ordering, for the DC power flow's B′.
#[allow(non_camel_case_types)]
pub type ybus_ordering = i32;

/// No reordering.
pub const YBUS_ORDERING_NATURAL: ybus_ordering = 1;
/// Reverse Cuthill–McKee, which keeps entries near the diagonal.
pub const YBUS_ORDERING_RCM: ybus_ordering = 2;
/// Approximate minimum degree, which keeps the factor's fill small: the default.
pub const YBUS_ORDERING_AMD: ybus_ordering = 3;
/// Column approximate minimum degree (COLAMD), for LU.
pub const YBUS_ORDERING_COLAMD: ybus_ordering = 4;

/// The code of an ordering
fn ordering_code(ordering: Ordering) -> Result<ybus_ordering, Failure> {
    match ordering {
        Ordering::Natural => Ok(YBUS_ORDERING_NATURAL),
        Ordering::Rcm => Ok(YBUS_ORDERING_RCM),
        Ordering::Amd => Ok(YBUS_ORDERING_AMD),
        Ordering::Colamd => Ok(YBUS_ORDERING_COLAMD),
        // Ordering is #[non_exhaustive]; only C's codes can set one here
        other => Err(Failure::InvalidValue(format!(
            "the ordering {other:?} has no code in this library"
        ))),
    }
}

/// The ordering of a code
fn ordering_of(code: ybus_ordering) -> Result<Ordering, Failure> {
    match code {
        YBUS_ORDERING_NATURAL => Ok(Ordering::Natural),
        YBUS_ORDERING_RCM => Ok(Ordering::Rcm),
        YBUS_ORDERING_AMD => Ok(Ordering::Amd),
        YBUS_ORDERING_COLAMD => Ok(Ordering::Colamd),
        other => Err(Failure::InvalidValue(format!(
            "{other} isn't an ordering, 1 to 4"
        ))),
    }
}

/// Options for the DC power flow, made by `ybus_dc_options_new` at MATPOWER's defaults and
/// freed by `ybus_dc_options_free`. A null options handle means the defaults wherever one is
/// taken.
pub struct YbusDcOptions {
    pub(crate) options: DcOptions,
}

/// Makes DC power flow options at the defaults (B′ ordered by `YBUS_ORDERING_AMD`) in `*out`,
/// which is set to `NULL` first and stays `NULL` on an error.
///
/// # Safety
///
/// `out` must be null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_options_new(out: *mut *mut YbusDcOptions) -> ybus_status {
    guard(|| {
        // SAFETY: out is as the caller promised
        unsafe {
            new_handle(out, || {
                Ok(YbusDcOptions {
                    options: DcOptions::default(),
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
pub unsafe extern "C" fn ybus_dc_options_free(options: *mut YbusDcOptions) {
    // SAFETY: as the caller promised
    unsafe { free_handle(options) }
}

/// Sets the ordering B′ is factored under. Every ordering gives the same solution, to within
/// rounding. Returns `YBUS_ERROR_INVALID_VALUE` for a code that isn't an ordering.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_options_set_ordering(
    options: *mut YbusDcOptions,
    code: ybus_ordering,
) -> ybus_status {
    guard(|| {
        // SAFETY: options is as the caller promised
        let options = unsafe { arg_mut(options, "options") }?;
        options.options.ordering = ordering_of(code)?;
        Ok(())
    })
}

/// Writes the ordering B′ is factored under through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `ybus_ordering`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_options_get_ordering(
    options: *const YbusDcOptions,
    out: *mut ybus_ordering,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(options, "options", out, |options| {
            ordering_code(options.options.ordering)
        })
    }
}

/// The options `options` points to, or the defaults for a null pointer.
///
/// # Safety
///
/// `options` must be null or live options.
pub(crate) unsafe fn options_or_default(
    options: *const YbusDcOptions,
) -> Result<DcOptions, Failure> {
    if options.is_null() {
        Ok(DcOptions::default())
    } else {
        // SAFETY: options is as the caller promised
        Ok(unsafe { arg_ref(options, "options") }?.options)
    }
}

/// The factorization a DC power flow solved B′ with.
#[allow(non_camel_case_types)]
pub type ybus_dc_solver = i32;

/// Cholesky, which applies unless some branch has a negative reactance.
pub const YBUS_DC_SOLVER_CHOLESKY: ybus_dc_solver = 1;
/// LU with threshold partial pivoting, the fallback when B′ isn't positive definite.
pub const YBUS_DC_SOLVER_LU: ybus_dc_solver = 2;
/// No solve was needed: every connected bus is a reference bus.
pub const YBUS_DC_SOLVER_NO_SOLVE: ybus_dc_solver = 3;

/// A DC power flow's results, indexed like the case's tables, made by `ybus_dc_power_flow` or
/// `ybus_dc_model_solve` and freed by `ybus_dc_solution_free`.
pub struct YbusDcSolution {
    pub(crate) solution: DcSolution,
}

/// Solves the DC power flow for `case`, as MATPOWER's `rundcpf` does, with `options` (`NULL`
/// for the defaults), into a new solution in `*out`. `*out` is set to `NULL` first and stays
/// `NULL` on an error, such as `YBUS_ERROR_ISLAND_WITHOUT_REFERENCE` or
/// `YBUS_ERROR_SINGULAR_NETWORK`.
///
/// # Safety
///
/// `case` must be null or a live case, `options` null or live options, and `out` null or a
/// writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_power_flow(
    case: *const YbusCase,
    options: *const YbusDcOptions,
    out: *mut *mut YbusDcSolution,
) -> ybus_status {
    guard(|| {
        // SAFETY: out, case and options are as the caller promised
        unsafe {
            new_handle(out, || {
                let case = arg_ref(case, "case")?;
                let options = options_or_default(options)?;
                Ok(YbusDcSolution {
                    solution: dc_power_flow_with_options(&case.case, &options)?,
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
pub unsafe extern "C" fn ybus_dc_solution_free(solution: *mut YbusDcSolution) {
    // SAFETY: as the caller promised
    unsafe { free_handle(solution) }
}

/// Writes the solution's bus column `field` into `out`, which has room for `len` entries: one
/// per bus. The only one is `YBUS_BUS_VA_DEG`, each bus's voltage angle (degrees); an isolated
/// bus keeps its angle from the case. Returns `YBUS_ERROR_INVALID_FIELD` for any other code,
/// and `YBUS_ERROR_BUFFER_TOO_SMALL`, writing nothing, if `len` is less than the number of
/// buses.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_solution_get_bus(
    solution: *const YbusDcSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "bus", field, out, len, |s, field| match field {
            YBUS_BUS_VA_DEG => Some(&s.solution.va_deg),
            _ => None,
        })
    }
}

/// [`ybus_dc_solution_get_bus`] for the generators: `YBUS_GEN_PG_MW`, each generator's real
/// power output (MW), where the first generator at each reference bus picks up the slack. 0
/// for a generator that's out of service or at an isolated bus.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_solution_get_gen(
    solution: *const YbusDcSolution,
    field: ybus_field,
    out: *mut f64,
    len: usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_column(solution, "gen", field, out, len, |s, field| match field {
            YBUS_GEN_PG_MW => Some(&s.solution.pg_mw),
            _ => None,
        })
    }
}

/// [`ybus_dc_solution_get_bus`] for the branches: `YBUS_BRANCH_PF_MW`, each branch's real power
/// flow at the from end (MW); the flow at the to end is its negation, since DC power flow is
/// lossless. 0 for a branch that's out of service or touches an isolated bus.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or `len` writable `double`s.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_solution_get_branch(
    solution: *const YbusDcSolution,
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
                _ => None,
            },
        )
    }
}

/// Writes which factorization solved B′ through `out`.
///
/// # Safety
///
/// `solution` must be null or a live solution, and `out` null or a writable `ybus_dc_solver`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_solution_solver(
    solution: *const YbusDcSolution,
    out: *mut ybus_dc_solver,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(solution, "solution", out, |s| match s.solution.solver {
            DcSolver::Cholesky => Ok(YBUS_DC_SOLVER_CHOLESKY),
            DcSolver::Lu => Ok(YBUS_DC_SOLVER_LU),
            DcSolver::NoSolve => Ok(YBUS_DC_SOLVER_NO_SOLVE),
        })
    }
}

/// A DC model of a case's network (B′, Bf and the bus roles), made by `ybus_dc_model_new` and
/// freed by `ybus_dc_model_free`, for solving many scenarios of the same network.
pub struct YbusDcModel {
    pub(crate) model: DcModel,
}

/// Builds the DC model of `case`'s network into a new model in `*out`, which is set to `NULL`
/// first and stays `NULL` on an error, such as `YBUS_ERROR_ZERO_REACTANCE`.
///
/// # Safety
///
/// `case` must be null or a live case, and `out` null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_model_new(
    case: *const YbusCase,
    out: *mut *mut YbusDcModel,
) -> ybus_status {
    guard(|| {
        // SAFETY: out and case are as the caller promised
        unsafe {
            new_handle(out, || {
                Ok(YbusDcModel {
                    model: DcModel::new(&arg_ref(case, "case")?.case)?,
                })
            })
        }
    })
}

/// Solves the DC power flow of `case`'s injections on `model`, with `options` (`NULL` for the
/// defaults), into a new solution in `*out`: the same result as `ybus_dc_power_flow` on
/// `case`, without rebuilding the network. B′ is still factored on each call.
///
/// `case` must be the case the model was built from, or a copy with different injections or
/// reference angles (`pd_mw`, `gs_mw`, `pg_mw` and `va_deg`). Returns
/// `YBUS_ERROR_CASE_MISMATCH` if its tables' lengths differ, and
/// `YBUS_ERROR_TOPOLOGY_MISMATCH` if its topology does: a bus's label or type, a generator's
/// bus or status, or a branch's ends or status. Other changes, such as a reactance, are
/// ignored; either needs a new model. `*out` is set to `NULL` first and stays `NULL` on an
/// error.
///
/// # Safety
///
/// `model` must be null or a live model, `case` null or a live case, `options` null or live
/// options, and `out` null or a writable pointer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_model_solve(
    model: *const YbusDcModel,
    case: *const YbusCase,
    options: *const YbusDcOptions,
    out: *mut *mut YbusDcSolution,
) -> ybus_status {
    guard(|| {
        // SAFETY: out, model, case and options are as the caller promised
        unsafe {
            new_handle(out, || {
                let model = arg_ref(model, "model")?;
                let case = arg_ref(case, "case")?;
                let options = options_or_default(options)?;
                Ok(YbusDcSolution {
                    solution: model.model.solve(&case.case, &options)?,
                })
            })
        }
    })
}

/// Frees `model`, which mustn't be used again. Does nothing if `model` is null.
///
/// # Safety
///
/// `model` must be null or a model made by this library and not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ybus_dc_model_free(model: *mut YbusDcModel) {
    // SAFETY: as the caller promised
    unsafe { free_handle(model) }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::ptr;

    use super::*;
    use crate::status::{YBUS_ERROR_INVALID_VALUE, YBUS_ERROR_NULL_ARGUMENT, YBUS_OK};

    #[test]
    fn dc_options_start_at_the_defaults_and_set_each_ordering() {
        let mut options = ptr::null_mut();
        let mut code = 0;
        // SAFETY: options is live until freed, and code is writable
        unsafe {
            assert_eq!(ybus_dc_options_new(&mut options), YBUS_OK);
            assert_eq!((*options).options, DcOptions::default());
            assert_eq!(ybus_dc_options_get_ordering(options, &mut code), YBUS_OK);
            assert_eq!(code, YBUS_ORDERING_AMD);
            for (code, expected) in [
                (YBUS_ORDERING_NATURAL, Ordering::Natural),
                (YBUS_ORDERING_RCM, Ordering::Rcm),
                (YBUS_ORDERING_COLAMD, Ordering::Colamd),
                (YBUS_ORDERING_AMD, Ordering::Amd),
            ] {
                assert_eq!(ybus_dc_options_set_ordering(options, code), YBUS_OK);
                assert_eq!((*options).options.ordering, expected);
            }
            for bad in [0, 5, -1] {
                assert_eq!(
                    ybus_dc_options_set_ordering(options, bad),
                    YBUS_ERROR_INVALID_VALUE
                );
            }
            assert_eq!((*options).options.ordering, Ordering::Amd);
            assert_eq!(
                ybus_dc_options_get_ordering(ptr::null(), &mut code),
                YBUS_ERROR_NULL_ARGUMENT
            );
            ybus_dc_options_free(options);
            ybus_dc_options_free(ptr::null_mut());
        }
    }

    /// The case at `path`, through the export
    pub(crate) fn load(path: &std::path::Path) -> *mut YbusCase {
        let path = std::ffi::CString::new(path.to_str().unwrap()).unwrap();
        let mut case = ptr::null_mut();
        // SAFETY: path is NUL-terminated, and case writable
        let status = unsafe { crate::case::ybus_case_from_path(path.as_ptr(), &mut case) };
        assert_eq!(status, YBUS_OK);
        case
    }

    /// The test cases: three_bus.m, and the MATPOWER cases with committed references, when
    /// they're fetched
    pub(crate) fn cases() -> Vec<(String, *mut YbusCase)> {
        let mut cases = vec![("three_bus".to_string(), crate::case::tests::three_bus())];
        for name in crate::testing::REFERENCE_CASES {
            if let Some(path) = crate::testing::matpower_case(name) {
                cases.push((name.to_string(), load(&path)));
            }
        }
        cases
    }

    #[test]
    fn the_dc_power_flow_is_powerflows_bit_for_bit() {
        for (name, case) in cases() {
            for code in [
                YBUS_ORDERING_AMD,
                YBUS_ORDERING_NATURAL,
                YBUS_ORDERING_RCM,
                YBUS_ORDERING_COLAMD,
            ] {
                let mut options = ptr::null_mut();
                let mut solution = ptr::null_mut();
                // SAFETY: each handle is live until freed, and each out writable
                unsafe {
                    assert_eq!(ybus_dc_options_new(&mut options), YBUS_OK);
                    assert_eq!(ybus_dc_options_set_ordering(options, code), YBUS_OK);
                    let status = ybus_dc_power_flow(case, options, &mut solution);
                    let expected = dc_power_flow_with_options(&(*case).case, &(*options).options);
                    match expected {
                        Ok(expected) => {
                            assert_eq!(status, YBUS_OK, "{name} {code}");
                            assert_eq!((*solution).solution, expected, "{name} {code}");
                        }
                        Err(error) => {
                            assert_eq!(status, crate::status::powerflow_code(&error));
                            assert!(solution.is_null());
                        }
                    }
                    ybus_dc_solution_free(solution);
                    ybus_dc_options_free(options);
                }
            }
            // SAFETY: case is live until this free
            unsafe { crate::case::ybus_case_free(case) };
        }
    }

    #[test]
    fn a_dc_solution_reads_its_columns_by_the_cases_codes() {
        let case = crate::case::tests::three_bus();
        let mut solution = ptr::null_mut();
        let (mut va, mut pg, mut pf) = ([0.0; 3], [0.0; 2], [0.0; 3]);
        let mut solver = 0;
        // SAFETY: each handle is live until freed, and each buffer has room for its length
        unsafe {
            assert_eq!(
                ybus_dc_power_flow(case, ptr::null(), &mut solution),
                YBUS_OK
            );
            let expected =
                dc_power_flow_with_options(&(*case).case, &DcOptions::default()).unwrap();
            assert_eq!(
                ybus_dc_solution_get_bus(solution, YBUS_BUS_VA_DEG, va.as_mut_ptr(), 3),
                YBUS_OK
            );
            assert_eq!(
                ybus_dc_solution_get_gen(solution, YBUS_GEN_PG_MW, pg.as_mut_ptr(), 2),
                YBUS_OK
            );
            assert_eq!(
                ybus_dc_solution_get_branch(solution, YBUS_BRANCH_PF_MW, pf.as_mut_ptr(), 3),
                YBUS_OK
            );
            assert_eq!(
                (&va[..], &pg[..], &pf[..]),
                (
                    &expected.va_deg[..],
                    &expected.pg_mw[..],
                    &expected.pf_mw[..]
                )
            );
            assert_eq!(ybus_dc_solution_solver(solution, &mut solver), YBUS_OK);
            assert_eq!(solver, YBUS_DC_SOLVER_CHOLESKY);
            assert_eq!(
                ybus_dc_solution_get_bus(
                    solution,
                    crate::case::columns::YBUS_BUS_VM_PU,
                    va.as_mut_ptr(),
                    3
                ),
                crate::status::YBUS_ERROR_INVALID_FIELD
            );
            assert_eq!(
                ybus_dc_solution_get_branch(solution, YBUS_BRANCH_PF_MW, pf.as_mut_ptr(), 2),
                crate::status::YBUS_ERROR_BUFFER_TOO_SMALL
            );
            ybus_dc_solution_free(solution);
            crate::case::ybus_case_free(case);
        }
    }

    #[test]
    fn a_dc_model_solves_scenarios_as_the_power_flow_does() {
        let case = crate::case::tests::three_bus();
        let (mut model, mut solution) = (ptr::null_mut(), ptr::null_mut());
        // SAFETY: each handle is live until freed, and each out writable
        unsafe {
            assert_eq!(ybus_dc_model_new(case, &mut model), YBUS_OK);
            for scale in [1.0, 0.5, 1.5] {
                for bus in &mut (*case).case.bus {
                    bus.pd_mw *= scale;
                }
                assert_eq!(
                    ybus_dc_model_solve(model, case, ptr::null(), &mut solution),
                    YBUS_OK
                );
                let expected = dc_power_flow_with_options(&(*case).case, &DcOptions::default());
                assert_eq!(Ok(&(*solution).solution), expected.as_ref());
                ybus_dc_solution_free(solution);
            }
            // Another topology is refused
            (&mut (*case).case.branch)[0].status = false;
            assert_eq!(
                ybus_dc_model_solve(model, case, ptr::null(), &mut solution),
                crate::status::YBUS_ERROR_TOPOLOGY_MISMATCH
            );
            assert!(solution.is_null());
            ybus_dc_model_free(model);
            ybus_dc_model_free(ptr::null_mut());
            crate::case::ybus_case_free(case);
        }
    }
}
