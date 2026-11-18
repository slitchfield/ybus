// Warm starts of the kept OPF models (OpfStart, solve_from): a cold start is solve's path bit for
//   bit; a warm step from one load level's solution to the next meets the warm-start survey's
//   bar (O3's without LAM_Q) against the cold solve of the next, and a primal step runs; and a
//   previous solution from another model is refused.

mod common;

use common::reference::relative_difference;
use ybus_powerflow::case::Case;
use ybus_powerflow::opf::ac::{AcOpfModel, AcOpfOptions, AcOpfSolution};
use ybus_powerflow::opf::dc::{DcOpfModel, DcOpfOptions, DcOpfSolution};
use ybus_powerflow::opf::{
    InteriorPointError, KktResiduals, OpfSolver, OpfStart, Solution, Status,
};
use ybus_powerflow::powerflow_error::PowerflowError;

/// O3's bar: the objective's relative difference, the prices', and the residuals' factor over
///   the reference's, or their floor
const OBJECTIVE: f64 = 1e-6;
const PRICE: f64 = 1e-3;
const RESIDUAL_FACTOR: f64 = 100.0;
const RESIDUAL_FLOOR: f64 = 1e-6;

fn scaled(case: &Case, factor: f64) -> Case {
    let mut case = case.clone();
    for bus in &mut case.bus {
        bus.pd_mw *= factor;
        bus.qd_mvar *= factor;
    }
    case
}

fn case(name: &str) -> Option<Case> {
    common::matpower_case(name).map(|path| Case::from_path(&path).unwrap())
}

/// A solve's results for the bar: converged, the objective, LAM_P, the residuals, and the
///   iterations and status, reported
struct Outcome<'a> {
    converged: bool,
    objective: f64,
    lam_p: &'a [f64],
    residuals: KktResiduals,
    iterations: usize,
    status: Status,
}

fn ac(s: &AcOpfSolution) -> Outcome<'_> {
    Outcome {
        converged: s.converged(),
        objective: s.objective,
        lam_p: &s.lam_p,
        residuals: s.residuals,
        iterations: s.iterations,
        status: s.status,
    }
}

fn dc(s: &DcOpfSolution) -> Outcome<'_> {
    Outcome {
        converged: s.converged(),
        objective: s.objective,
        lam_p: &s.lam_p,
        residuals: s.residuals,
        iterations: s.iterations,
        status: s.status,
    }
}

/// What the survey's bar (O3's without LAM_Q) finds wrong with `step` against the cold solve `cold` of the same case: a
///   different outcome, and, where the cold solve converged, its values and residuals
fn misses(step: &Outcome, cold: &Outcome) -> Vec<String> {
    let mut found = vec![];
    if step.converged != cold.converged {
        found.push(format!(
            "converged {} (cold {})",
            step.converged, cold.converged
        ));
    }
    if !cold.converged {
        return found;
    }
    let objective = (step.objective - cold.objective).abs() / cold.objective.abs().max(1.0);
    if objective > OBJECTIVE || objective.is_nan() {
        found.push(format!("objective {objective:e}"));
    }
    // LAM_P only: LAM_Q is the survey's to report, not bar, since two cold solvers already
    //   differ in it by up to 9.4e-3 (docs/OptimizationChecklist.md, the warm-start rule)
    let lam_p = relative_difference(step.lam_p, cold.lam_p);
    if lam_p > PRICE || lam_p.is_nan() {
        found.push(format!("LAM_P {lam_p:e}"));
    }
    let (a, b) = (step.residuals, cold.residuals);
    for (name, a, b) in [
        ("stationarity", a.stationarity, b.stationarity),
        ("feasibility", a.feasibility, b.feasibility),
        ("complementarity", a.complementarity, b.complementarity),
        ("dual feasibility", a.dual_feasibility, b.dual_feasibility),
    ] {
        if a > RESIDUAL_FACTOR * b && a > RESIDUAL_FLOOR {
            found.push(format!("{name} {a:e} (cold {b:e})"));
        }
    }
    found
}

fn bits(s: &Solution) -> String {
    format!(
        "{:?} {:?} {} {} {:?}",
        s.x,
        s.multipliers,
        s.f.to_bits(),
        s.iterations,
        s.status
    )
}

const CASES: [&str; 2] = ["case118", "case2383wp"];

/// The warm starts' shift in these tests. The default is the sequences survey's to choose; on
///   case2383wp's DC OPF the predictor-corrector method fails from 1e-4 but not from 1e-2
const SHIFT: f64 = 1e-2;

fn shifted(solver: OpfSolver) -> OpfSolver {
    match solver {
        OpfSolver::Mips(mut options) => {
            options.warm_shift = SHIFT;
            OpfSolver::Mips(options)
        }
        OpfSolver::Mehrotra(mut options) => {
            options.warm_shift = SHIFT;
            OpfSolver::Mehrotra(options)
        }
        other => other,
    }
}

#[test]
fn a_cold_start_is_solves_path_bit_for_bit() {
    for name in CASES {
        let Some(case) = case(name) else { continue };
        let changed = scaled(&case, 1.02);
        let model = DcOpfModel::new(&case).unwrap();
        for options in [DcOpfOptions::default(), DcOpfOptions::mehrotra()] {
            let a = model.solve(&changed, &options).unwrap();
            let b = model
                .solve_from(&changed, &options, OpfStart::Cold)
                .unwrap();
            assert_eq!(bits(&a.solver), bits(&b.solver), "{name} DC");
        }
        let model = AcOpfModel::new(&case).unwrap();
        for options in [AcOpfOptions::default(), AcOpfOptions::mehrotra()] {
            let a = model.solve(&changed, &options).unwrap();
            let b = model
                .solve_from(&changed, &options, OpfStart::Cold)
                .unwrap();
            assert_eq!(bits(&a.solver), bits(&b.solver), "{name} AC");
        }
    }
}

/// Checks the warm step against `cold`, the same case solved cold, by O3's bar, and reports
///   the primal step's outcome without requiring it: a primal-only start keeps the point on its
///   active bounds with cold multipliers, and on case2383wp's DC OPF the predictor-corrector
///   method fails from there. How often is the sequences survey's to measure.
fn check_steps(label: &str, cold: &Outcome, warm: &Outcome, primal: &Outcome) {
    let found = misses(warm, cold);
    assert!(
        found.is_empty(),
        "{label}, warm: {found:?}; {:?} {} in {} (cold {:?} {} in {})",
        warm.status,
        warm.objective,
        warm.iterations,
        cold.status,
        cold.objective,
        cold.iterations
    );
    eprintln!(
        "{label}: cold {}, warm {}, primal {} ({:?}{})",
        cold.iterations,
        warm.iterations,
        primal.iterations,
        primal.status,
        if misses(primal, cold).is_empty() {
            ""
        } else {
            ", outside the bar"
        }
    );
}

#[test]
fn a_dc_warm_step_from_the_last_solution_meets_o3s_bar() {
    for name in CASES {
        let Some(case) = case(name) else { continue };
        let (first, next) = (case.clone(), scaled(&case, 1.01));
        let model = DcOpfModel::new(&case).unwrap();
        for (solver, options) in [
            ("MIPS", DcOpfOptions::default()),
            ("predictor-corrector", DcOpfOptions::mehrotra()),
        ] {
            let mut options = options;
            options.solver = shifted(options.solver);
            let previous = model.solve(&first, &options).unwrap();
            if !previous.converged() {
                // No solution to start from (case2383wp's DC on MIPS, a numerical failure)
                eprintln!("{name}, {solver}: the first solve didn't converge");
                continue;
            }
            let cold = model.solve(&next, &options).unwrap();
            let warm = model
                .solve_from(&next, &options, OpfStart::Warm(&previous.solver))
                .unwrap();
            let primal = model
                .solve_from(&next, &options, OpfStart::Primal(&previous.solver))
                .unwrap();
            check_steps(
                &format!("{name} DC, {solver}"),
                &dc(&cold),
                &dc(&warm),
                &dc(&primal),
            );
            // A warm start that fell back to a cold one would meet the bar trivially, and a
            //   primal one that carried the multipliers would be the warm one
            if cold.converged() && warm.converged() {
                assert!(
                    warm.iterations < cold.iterations,
                    "{name}, {solver}: no fewer iterations warm"
                );
            }
            assert_ne!(
                bits(&warm.solver),
                bits(&primal.solver),
                "{name}, {solver}: the primal step is the warm one"
            );
        }
    }
}

#[test]
fn an_ac_warm_step_from_the_last_solution_meets_o3s_bar() {
    for name in CASES {
        let Some(case) = case(name) else { continue };
        let (first, next) = (case.clone(), scaled(&case, 1.01));
        let model = AcOpfModel::new(&case).unwrap();
        for (solver, options) in [
            ("MIPS", AcOpfOptions::default()),
            ("predictor-corrector", AcOpfOptions::mehrotra()),
        ] {
            let mut options = options;
            options.solver = shifted(options.solver);
            let previous = model.solve(&first, &options).unwrap();
            if !previous.converged() {
                // No solution to start from (case2383wp's DC on MIPS, a numerical failure)
                eprintln!("{name}, {solver}: the first solve didn't converge");
                continue;
            }
            let cold = model.solve(&next, &options).unwrap();
            let warm = model
                .solve_from(&next, &options, OpfStart::Warm(&previous.solver))
                .unwrap();
            let primal = model
                .solve_from(&next, &options, OpfStart::Primal(&previous.solver))
                .unwrap();
            check_steps(
                &format!("{name} AC, {solver}"),
                &ac(&cold),
                &ac(&warm),
                &ac(&primal),
            );
            // A warm start that fell back to a cold one would meet the bar trivially, and a
            //   primal one that carried the multipliers would be the warm one
            if cold.converged() && warm.converged() {
                assert!(
                    warm.iterations < cold.iterations,
                    "{name}, {solver}: no fewer iterations warm"
                );
            }
            assert_ne!(
                bits(&warm.solver),
                bits(&primal.solver),
                "{name}, {solver}: the primal step is the warm one"
            );
        }
    }
}

#[test]
fn a_solution_from_another_model_is_refused() {
    let (Some(small), Some(other)) = (case("case9"), case("case30")) else {
        return;
    };
    let dc_options = DcOpfOptions::default();
    let foreign = DcOpfModel::new(&other)
        .unwrap()
        .solve(&other, &dc_options)
        .unwrap();
    let model = DcOpfModel::new(&small).unwrap();
    for start in [
        OpfStart::Warm(&foreign.solver),
        OpfStart::Primal(&foreign.solver),
    ] {
        let result = model.solve_from(&small, &dc_options, start);
        assert!(
            matches!(
                result,
                Err(PowerflowError::Solver {
                    error: InteriorPointError::Dimension { .. }
                })
            ),
            "DC: {:?}",
            result.err()
        );
    }
    let ac_options = AcOpfOptions::mehrotra();
    let foreign = AcOpfModel::new(&other)
        .unwrap()
        .solve(&other, &ac_options)
        .unwrap();
    let model = AcOpfModel::new(&small).unwrap();
    for start in [
        OpfStart::Warm(&foreign.solver),
        OpfStart::Primal(&foreign.solver),
    ] {
        let result = model.solve_from(&small, &ac_options, start);
        assert!(
            matches!(
                result,
                Err(PowerflowError::Solver {
                    error: InteriorPointError::Dimension { .. }
                })
            ),
            "AC: {:?}",
            result.err()
        );
    }
}
