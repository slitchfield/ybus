// Compares the DC power flow with MATPOWER's rundcpf on case14, case30 and case118, whose
//   results are committed in tests/data/matpower-reference/. The cases themselves are fetched,
//   not committed (see tests/common), so these tests skip without them, unless
//   YBUS_REQUIRE_CASES=1.

mod common;
use common::reference::{self, Deviation};
use powerflow::case::Case;
use powerflow::dc::{dc_power_flow, DcSolution, DcSolver};

/// The largest differences allowed from MATPOWER's results. The files hold MATPOWER's exact
/// doubles, so what differs is rounding: Octave's sparse backslash reorders B′ to limit fill-in
/// and this crate doesn't, so they add in different orders. The largest differences measured
/// were 1.4e-14° and 2.8e-13 MW on case14, 1.2e-14° and 5.5e-14 MW on case30, and 2.3e-13° and
/// 3.6e-12 MW on case118. 1e-9 leaves room for other platforms and summation orders, and any
/// real error would be orders of magnitude larger.
const TOLERANCE: Deviation = Deviation {
    va_deg: 1e-9,
    pf_mw: 1e-9,
    pg_mw: 1e-9,
};

/// Parses and solves a fetched MATPOWER case, or `None` if MATPOWER isn't fetched.
fn solve(name: &str) -> Option<(Case, DcSolution)> {
    let path = common::matpower_case(name)?;
    let case = Case::from_path(&path).unwrap();
    let solution = dc_power_flow(&case).unwrap();
    Some((case, solution))
}

fn compare_with_matpower(name: &str) {
    let Some((case, solution)) = solve(name) else {
        return; // not fetched; the skip was reported on stderr
    };
    let reference = reference::read(&reference::committed_dir(), name).unwrap();
    let deviation = reference::deviation(&case, &solution, &reference).unwrap();
    assert!(
        deviation.within(&TOLERANCE),
        "{name}: {deviation:?} exceeds {TOLERANCE:?}"
    );
    // None of these cases has a negative reactance, so B′ is positive definite
    assert_eq!(solution.solver, DcSolver::Cholesky);
}

#[test]
fn case14_matches_matpower() {
    compare_with_matpower("case14");
}

#[test]
fn case30_matches_matpower() {
    compare_with_matpower("case30");
}

#[test]
fn case118_matches_matpower() {
    compare_with_matpower("case118");
}

#[test]
fn the_comparison_is_not_vacuous() {
    let Some((case, solution)) = solve("case14") else {
        return;
    };
    let committed = reference::committed_dir();

    // Another case's results don't describe this case
    let other = reference::read(&committed, "case30").unwrap();
    assert!(reference::deviation(&case, &solution, &other).is_err());

    // An error far smaller than any that matters, or a NaN, fails the tolerance
    let reference = reference::read(&committed, "case14").unwrap();
    let mut wrong = solution.clone();
    wrong.pf_mw[0] += 1e-6;
    let deviation = reference::deviation(&case, &wrong, &reference).unwrap();
    assert!(!deviation.within(&TOLERANCE), "{deviation:?}");
    let mut wrong = solution;
    wrong.va_deg[1] = f64::NAN;
    let deviation = reference::deviation(&case, &wrong, &reference).unwrap();
    assert!(deviation.va_deg.is_nan() && !deviation.within(&TOLERANCE));
}
