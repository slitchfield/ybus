// Checks that the DC solve is backward stable on MATPOWER cases: its angles solve B′θ = P′ to
//   within rounding, which unlike the differences from MATPOWER doesn't grow with a case's
//   size or conditioning. The survey records the same measure for every case
//   (docs/benchmarks/), where it's tracked against fill-in. The cases are fetched, not
//   committed (see tests/common), so these tests skip without them, unless
//   YBUS_REQUIRE_CASES=1.

mod common;
use common::accuracy::backward_error;
use powerflow::case::Case;
use powerflow::dc::{DcModel, DcSolution, DcSolver};

/// The largest backward error allowed: 100 ε, where ε is `f64::EPSILON` (2.2e-16). These cases
/// measured 0.07 ε (case300) to 0.38 ε (case118), and the worst of every MATPOWER case in the
/// survey is 10.7 ε (case2383wp, whose factor fills in heavily). A backward-stable solve stays
/// within a small multiple of ε, while a 1e-9° error in a single angle of case118 already
/// measures 3,700 ε.
const LIMIT: f64 = 100.0 * f64::EPSILON;

/// Parses, models and solves a fetched MATPOWER case, or `None` if MATPOWER isn't fetched.
fn solve(name: &str) -> Option<(Case, DcModel, DcSolution)> {
    let path = common::matpower_case(name)?;
    let case = Case::from_path(&path).unwrap();
    let model = DcModel::new(&case).unwrap();
    let solution = model.solve(&case).unwrap();
    Some((case, model, solution))
}

fn assert_backward_stable(name: &str, solver: DcSolver) {
    let Some((case, model, solution)) = solve(name) else {
        return; // not fetched; the skip was reported on stderr
    };
    assert_eq!(solution.solver, solver, "{name}");
    let eta = backward_error(&case, &model, &solution.va_deg);
    assert!(
        eta <= LIMIT,
        "{name}: backward error {eta:e} exceeds {LIMIT:e}"
    );
}

#[test]
fn cholesky_solves_are_backward_stable() {
    for name in ["case14", "case30", "case118", "case_ACTIVSg2000"] {
        assert_backward_stable(name, DcSolver::Cholesky);
    }
}

#[test]
fn lu_solves_are_backward_stable() {
    // Negative series susceptances make these B′ indefinite
    for name in ["case60nordic", "case300"] {
        assert_backward_stable(name, DcSolver::Lu);
    }
}

#[test]
fn the_check_is_not_vacuous() {
    let Some((case, model, solution)) = solve("case118") else {
        return;
    };
    // A 1e-9° error in one PQ bus's angle, far below anything that matters, fails it
    let mut wrong = solution.va_deg.clone();
    wrong[1] += 1e-9;
    let eta = backward_error(&case, &model, &wrong);
    assert!(eta > LIMIT, "{eta:e}");
    // And so does a NaN
    wrong[1] = f64::NAN;
    assert!(backward_error(&case, &model, &wrong).is_nan());
}
