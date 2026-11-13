//! The DC optimal power flow against MATPOWER's `rundcopf` (MIPS forced): on the committed
//! reference cases, and on any case named in `YBUS_DCOPF_CASES` whose results
//! `scripts/generate-reference.sh` generated locally (never committed).

mod common;

use std::path::Path;

use common::reference::{committed_dir, dcopf_deviation, read_dcopf};
use ybus_powerflow::case::Case;
use ybus_powerflow::opf::dc::dc_opf;

/// How far each quantity may be from MATPOWER's (MW, degrees, $/MWh, $/degree), as for the
/// power flows. The largest difference measured on the committed cases is 6.9e-12 (case118's
/// dispatch), at MATPOWER's iteration counts.
const TOLERANCE: f64 = 1e-9;

/// The objective's relative difference: 1.7e-15 measured (case14).
const OBJECTIVE: f64 = 1e-12;

fn check(name: &str, case_path: &Path, reference_dir: &Path) {
    let case = Case::from_path(case_path).unwrap();
    let reference = read_dcopf(reference_dir, name, "dcopf").unwrap();
    let solution = dc_opf(&case).unwrap();
    assert_eq!(
        (solution.iterations, solution.converged()),
        (reference.iterations, reference.success),
        "{name}: iterations and convergence"
    );
    for (quantity, deviation) in dcopf_deviation(&solution, &reference) {
        let tolerance = if quantity.starts_with("objective") {
            OBJECTIVE
        } else {
            TOLERANCE
        };
        assert!(
            deviation <= tolerance,
            "{name}: {quantity} differs by {deviation:e}"
        );
    }
    // The default solver's optimum, where it found one, confirms MIPS's
    if solution.converged() && reference.default_objective.is_finite() {
        let relative = (solution.objective - reference.default_objective).abs()
            / reference.default_objective.abs().max(1.0);
        assert!(
            relative < 1e-6,
            "{name}: {} vs {} by {}",
            solution.objective,
            reference.default_objective,
            reference.default_solver
        );
    }
}

#[test]
fn committed_cases_match_rundcopf() {
    for name in ["case14", "case30", "case118"] {
        if let Some(path) = common::matpower_case(name) {
            check(name, &path, &committed_dir());
        }
    }
}

/// The local comparisons' bars, for cases whose optimum may not be unique. Where MATPOWER and
/// `ybus_powerflow` both converge at the same iteration count, the objectives and prices agree, but
/// on degenerate problems (linear costs, several generators at the same marginal cost) the
/// dispatch, angles and flows are one point of a set of optima, and differ: by up to 61.5 MW on
/// case8387pegase, as MATPOWER's MP-Core and legacy implementations differ from each other. So
/// those aren't compared here; the committed cases, whose optima are unique, compare them.
const LOCAL_OBJECTIVE: f64 = 1e-10;
const LOCAL_PRICE: f64 = 1e-8;
/// The objective's bar against MATPOWER's default solver's optimum: MIPS's tolerance
const DEFAULT_OBJECTIVE: f64 = 1e-6;
/// The objective's bar where MATPOWER's run took Octave's least-squares fallback (see
/// `generate_reference.m`) and still converged: MIPS's own tolerance
const FALLBACK_OBJECTIVE: f64 = 1e-6;

/// Compares one local case, returning what differs.
fn check_local(name: &str, case_path: &Path, reference_dir: &Path) -> Vec<String> {
    let case = Case::from_path(case_path).unwrap();
    let reference = read_dcopf(reference_dir, name, "dcopf").unwrap();
    let solution = dc_opf(&case).unwrap();
    let mut problems = vec![];
    if solution.converged() != reference.success {
        problems.push(format!(
            "{name}: {:?} after {} iterations, MATPOWER {} after {}",
            solution.status,
            solution.iterations,
            if reference.success {
                "converged"
            } else {
                "failed"
            },
            reference.iterations
        ));
        return problems;
    }
    // Iterations, unless MATPOWER failed (both fail, after however long) or Octave's fallback
    //   changed its steps
    if reference.success && !reference.singular && solution.iterations != reference.iterations {
        problems.push(format!(
            "{name}: {} iterations, MATPOWER {}",
            solution.iterations, reference.iterations
        ));
    }
    if !reference.success {
        return problems;
    }
    // The default solver's optimum (GLPK for linear costs), independently of MIPS
    if reference.default_objective.is_finite() {
        let relative = (solution.objective - reference.default_objective).abs()
            / reference.default_objective.abs().max(1.0);
        if relative.is_nan() || relative > DEFAULT_OBJECTIVE {
            problems.push(format!(
                "{name}: objective {} vs {} by {}",
                solution.objective, reference.default_objective, reference.default_solver
            ));
        }
    }
    for (quantity, deviation) in dcopf_deviation(&solution, &reference) {
        let tolerance = match quantity {
            "objective (relative)" if reference.singular => FALLBACK_OBJECTIVE,
            "objective (relative)" => LOCAL_OBJECTIVE,
            "va_deg" | "pg_mw" | "pf_mw" => continue,
            _ if reference.singular => continue,
            _ => LOCAL_PRICE,
        };
        if deviation.is_nan() || deviation > tolerance {
            problems.push(format!("{name}: {quantity} differs by {deviation:e}"));
        }
    }
    problems
}

/// Local comparisons: `YBUS_DCOPF_CASES` names cases, separated by spaces, whose MATPOWER
/// results are in the local reference directory next to the normalized cases (made by
/// `scripts/generate-reference.sh`). Every case is checked, and the differences reported
/// together.
#[test]
fn local_cases_match_rundcopf() {
    let Ok(names) = std::env::var("YBUS_DCOPF_CASES") else {
        return;
    };
    let mut problems = vec![];
    for name in names.split_whitespace() {
        let path =
            common::normalized_case(name).unwrap_or_else(|| panic!("{name}: no normalized case"));
        let dir = path
            .parent()
            .and_then(|dir| dir.to_str())
            .map(|dir| dir.replace("-normalized", "-reference"))
            .unwrap();
        eprintln!("{name}: starting");
        let start = std::time::Instant::now();
        let found = check_local(name, &path, Path::new(&dir));
        eprintln!(
            "{name}: {} in {:.1} s",
            if found.is_empty() { "ok" } else { "differs" },
            start.elapsed().as_secs_f64()
        );
        problems.extend(found);
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
