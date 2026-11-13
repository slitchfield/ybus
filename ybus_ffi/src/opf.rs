//! The optimal power flows through C: what the DC and AC OPFs share, their solver codes and the
//! settings both solvers have, read and written on either OPF's options.

use ybus_powerflow::opf::{KktFactorization, OpfSolver};

use crate::boundary::{arg_mut, get_value, guard, Failure};
use crate::status::ybus_status;

pub mod ac;
pub mod dc;

/// Which interior-point solver an optimal power flow runs.
#[allow(non_camel_case_types)]
pub type ybus_opf_solver = i32;

/// MIPS, ported step for step from MATPOWER's `mips.m`: MATPOWER's default, whose results match
/// MATPOWER's.
pub const YBUS_OPF_MIPS: ybus_opf_solver = 1;
/// The primal-dual method with Mehrotra's predictor-corrector, on the same problem.
pub const YBUS_OPF_MEHROTRA: ybus_opf_solver = 2;

/// How the solver's KKT system is factored.
#[allow(non_camel_case_types)]
pub type ybus_kkt_factorization = i32;

/// LU: MATPOWER's choice, and MIPS's default.
pub const YBUS_KKT_LU: ybus_kkt_factorization = 1;
/// The pivoted LDLᵀ: the predictor-corrector method's default.
pub const YBUS_KKT_LDLT: ybus_kkt_factorization = 2;

/// How an optimal power flow's solver ended, with MIPS's `exitflag` codes.
#[allow(non_camel_case_types)]
pub type ybus_opf_status = i32;

/// Every termination test passed (`exitflag` 1): MATPOWER's `success`.
pub const YBUS_OPF_CONVERGED: ybus_opf_status = 1;
/// The iteration limit passed without converging (`exitflag` 0).
pub const YBUS_OPF_MAX_ITERATIONS: ybus_opf_status = 0;
/// The solve failed numerically (`exitflag` -1).
pub const YBUS_OPF_NUMERICAL_FAILURE: ybus_opf_status = -1;

/// The code of an OPF's status
pub(crate) fn status_code(status: ybus_powerflow::opf::Status) -> ybus_opf_status {
    match status {
        ybus_powerflow::opf::Status::Converged => YBUS_OPF_CONVERGED,
        ybus_powerflow::opf::Status::MaxIterations => YBUS_OPF_MAX_ITERATIONS,
        ybus_powerflow::opf::Status::NumericalFailure => YBUS_OPF_NUMERICAL_FAILURE,
    }
}

/// An OPF options handle: either OPF's, each holding an [`OpfSolver`].
pub(crate) trait OpfOptions {
    fn solver(&self) -> &OpfSolver;
    fn solver_mut(&mut self) -> &mut OpfSolver;
}

/// The settings both solvers have, by reference into either's options.
struct Common<'a> {
    feastol: &'a mut f64,
    gradtol: &'a mut f64,
    comptol: &'a mut f64,
    costtol: &'a mut f64,
    cost_mult: &'a mut f64,
    max_it: &'a mut usize,
    kkt_factorization: &'a mut KktFactorization,
}

/// The common settings of `solver`'s options
fn common(solver: &mut OpfSolver) -> Result<Common<'_>, Failure> {
    match solver {
        OpfSolver::Mips(o) => Ok(Common {
            feastol: &mut o.feastol,
            gradtol: &mut o.gradtol,
            comptol: &mut o.comptol,
            costtol: &mut o.costtol,
            cost_mult: &mut o.cost_mult,
            max_it: &mut o.max_it,
            kkt_factorization: &mut o.kkt_factorization,
        }),
        OpfSolver::Mehrotra(o) => Ok(Common {
            feastol: &mut o.feastol,
            gradtol: &mut o.gradtol,
            comptol: &mut o.comptol,
            costtol: &mut o.costtol,
            cost_mult: &mut o.cost_mult,
            max_it: &mut o.max_it,
            kkt_factorization: &mut o.kkt_factorization,
        }),
        // OpfSolver is #[non_exhaustive]; only C's codes can choose one here
        other => Err(Failure::InvalidValue(format!(
            "the solver {other:?} has no code in this library"
        ))),
    }
}

/// One of the `double` settings both solvers have
#[derive(Clone, Copy)]
pub(crate) enum Float {
    Feastol,
    Gradtol,
    Comptol,
    Costtol,
    CostMult,
}

impl Float {
    fn name(self) -> &'static str {
        match self {
            Float::Feastol => "feastol",
            Float::Gradtol => "gradtol",
            Float::Comptol => "comptol",
            Float::Costtol => "costtol",
            Float::CostMult => "cost_mult",
        }
    }

    fn field(self, common: Common<'_>) -> &mut f64 {
        match self {
            Float::Feastol => common.feastol,
            Float::Gradtol => common.gradtol,
            Float::Comptol => common.comptol,
            Float::Costtol => common.costtol,
            Float::CostMult => common.cost_mult,
        }
    }
}

/// Sets a `double` setting, which must be positive and finite.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
pub(crate) unsafe fn set_float<O: OpfOptions>(
    options: *mut O,
    setting: Float,
    value: f64,
) -> ybus_status {
    guard(|| {
        // SAFETY: options is as the caller promised
        let options = unsafe { arg_mut(options, "options") }?;
        if !(value.is_finite() && value > 0.0) {
            return Err(Failure::InvalidValue(format!(
                "the {} {value} isn't positive and finite",
                setting.name()
            )));
        }
        *setting.field(common(options.solver_mut())?) = value;
        Ok(())
    })
}

/// Writes a `double` setting through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `double`.
pub(crate) unsafe fn get_float<O: OpfOptions>(
    options: *const O,
    setting: Float,
    out: *mut f64,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(options, "options", out, |options| {
            let mut solver = *options.solver();
            Ok(*setting.field(common(&mut solver)?))
        })
    }
}

/// Sets the iteration limit, the solvers' `max_it`.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
pub(crate) unsafe fn set_max_iterations<O: OpfOptions>(
    options: *mut O,
    value: usize,
) -> ybus_status {
    guard(|| {
        // SAFETY: options is as the caller promised
        let options = unsafe { arg_mut(options, "options") }?;
        *common(options.solver_mut())?.max_it = value;
        Ok(())
    })
}

/// Writes the iteration limit through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `size_t`.
pub(crate) unsafe fn get_max_iterations<O: OpfOptions>(
    options: *const O,
    out: *mut usize,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(options, "options", out, |options| {
            let mut solver = *options.solver();
            Ok(*common(&mut solver)?.max_it)
        })
    }
}

/// Sets the KKT factorization from its code.
///
/// # Safety
///
/// `options` must be null or live options that nothing else is using.
pub(crate) unsafe fn set_kkt<O: OpfOptions>(
    options: *mut O,
    code: ybus_kkt_factorization,
) -> ybus_status {
    guard(|| {
        // SAFETY: options is as the caller promised
        let options = unsafe { arg_mut(options, "options") }?;
        let kkt = match code {
            YBUS_KKT_LU => KktFactorization::Lu,
            YBUS_KKT_LDLT => KktFactorization::Ldlt,
            other => {
                return Err(Failure::InvalidValue(format!(
                    "{other} isn't a KKT factorization, 1 or 2"
                )))
            }
        };
        *common(options.solver_mut())?.kkt_factorization = kkt;
        Ok(())
    })
}

/// Writes the KKT factorization's code through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable
/// `ybus_kkt_factorization`.
pub(crate) unsafe fn get_kkt<O: OpfOptions>(
    options: *const O,
    out: *mut ybus_kkt_factorization,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(options, "options", out, |options| {
            let mut solver = *options.solver();
            match *common(&mut solver)?.kkt_factorization {
                KktFactorization::Lu => Ok(YBUS_KKT_LU),
                KktFactorization::Ldlt => Ok(YBUS_KKT_LDLT),
                // KktFactorization is #[non_exhaustive]; only C's codes can set one here
                other => Err(Failure::InvalidValue(format!(
                    "the KKT factorization {other:?} has no code in this library"
                ))),
            }
        })
    }
}

/// Writes the solver's code through `out`.
///
/// # Safety
///
/// `options` must be null or live options, and `out` null or a writable `ybus_opf_solver`.
pub(crate) unsafe fn get_solver<O: OpfOptions>(
    options: *const O,
    out: *mut ybus_opf_solver,
) -> ybus_status {
    // SAFETY: as the caller promised
    unsafe {
        get_value(options, "options", out, |options| match options.solver() {
            OpfSolver::Mips(_) => Ok(YBUS_OPF_MIPS),
            OpfSolver::Mehrotra(_) => Ok(YBUS_OPF_MEHROTRA),
            other => Err(Failure::InvalidValue(format!(
                "the solver {other:?} has no code in this library"
            ))),
        })
    }
}
