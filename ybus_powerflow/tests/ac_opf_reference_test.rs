//! The AC optimal power flow against MATPOWER's `runopf` (MIPS): on the committed reference
//! cases, and on any case named in `YBUS_ACOPF_CASES` whose results
//! `scripts/generate-reference.sh` generated locally (never committed).

mod common;

use std::path::Path;

use common::reference::{acopf_deviation, committed_dir, read_acopf, relative_difference};
use ybus_powerflow::case::Case;
use ybus_powerflow::opf::ac::ac_opf;

/// How far each quantity may be from MATPOWER's (p.u., degrees, MW, MVAr and the prices' units)
/// on the committed cases. The largest difference measured is 5.4e-9, in case30's voltage-limit
/// multipliers (values in the hundreds of $/p.u.-hr), at MATPOWER's iteration counts; the power
/// flows' 1e-9 is a little tight for the AC OPF's multipliers.
const TOLERANCE: f64 = 1e-8;

/// The objective's relative difference: 1.1e-15 measured (case14).
const OBJECTIVE: f64 = 1e-12;

#[test]
fn committed_cases_match_runopf() {
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();
        let reference = read_acopf(&committed_dir(), name, "opf").unwrap();
        let solution = ac_opf(&case).unwrap();
        assert_eq!(
            (solution.iterations, solution.converged()),
            (reference.iterations, reference.success),
            "{name}: iterations and convergence"
        );
        for (quantity, deviation) in acopf_deviation(&solution, &reference) {
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
    }
}

/// The local comparisons' bars. Where both converge: the objective, and the bus prices relative
/// to their magnitude (at least 1). Not the primal solution, which degenerate problems leave
/// non-unique (as for the DC OPF), and not the individual bound multipliers: where branches of
/// near-zero impedance lock neighbouring buses together at a voltage limit, the split of the
/// limit's price between their multipliers is ill-determined. On case3012wp `mu_vmax` differs by
/// up to 0.43 at buses whose differences cancel (99% of their magnitude), and `LAM_Q` shifts by
/// 1.16e-4 over the region, while the objective agrees within 1.3e-12. Measured over all local
/// cases (except case8387pegase, below): objectives within 1.1e-11 (case141), `LAM_P` within
/// 2.4e-7 and `LAM_Q` within 9.4e-5 (case3012wp).
///
/// case8387pegase takes 82 iterations to MATPOWER's 80, with its reactive dispatch 246 MVAr and
/// its objective 1.5e-6 away: its optimum isn't unique, and MATPOWER's legacy code differs from
/// MP-Core there by as much (246 MVAr, 9.3e-7). It is reported as a difference.
const LOCAL_OBJECTIVE: f64 = 1e-10;
const LOCAL_LAM_P: f64 = 1e-6;
const LOCAL_LAM_Q: f64 = 2e-4;
/// The objective's bar where MATPOWER's run took Octave's least-squares fallback
const FALLBACK_OBJECTIVE: f64 = 1e-6;

/// Compares one local case, returning what differs.
fn check_local(name: &str, case_path: &Path, reference_dir: &Path) -> Vec<String> {
    let case = Case::from_path(case_path).unwrap();
    let reference = read_acopf(reference_dir, name, "opf").unwrap();
    let solution = match ac_opf(&case) {
        Ok(solution) => solution,
        Err(err) => return vec![format!("{name}: {err}")],
    };
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
    if reference.success && !reference.singular && solution.iterations != reference.iterations {
        problems.push(format!(
            "{name}: {} iterations, MATPOWER {}",
            solution.iterations, reference.iterations
        ));
    }
    if !reference.success {
        return problems;
    }
    let objective =
        (solution.objective - reference.objective).abs() / reference.objective.abs().max(1.0);
    let objective_bar = if reference.singular {
        FALLBACK_OBJECTIVE
    } else {
        LOCAL_OBJECTIVE
    };
    let mut checks = vec![("objective", objective, objective_bar)];
    if !reference.singular {
        checks.push((
            "lam_p",
            relative_difference(&solution.lam_p, &reference.lam_p),
            LOCAL_LAM_P,
        ));
        checks.push((
            "lam_q",
            relative_difference(&solution.lam_q, &reference.lam_q),
            LOCAL_LAM_Q,
        ));
    }
    for (quantity, difference, bar) in checks {
        if difference.is_nan() || difference > bar {
            problems.push(format!(
                "{name}: {quantity} differs by {difference:e} (relative)"
            ));
        }
    }
    problems
}

/// Local comparisons: `YBUS_ACOPF_CASES` names cases, separated by spaces, whose MATPOWER
/// results are in the local reference directory next to the normalized cases (made by
/// `scripts/generate-reference.sh`). Every case is checked, and the differences reported
/// together; run with `--nocapture` to follow each case.
#[test]
fn local_cases_match_runopf() {
    let Ok(names) = std::env::var("YBUS_ACOPF_CASES") else {
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
        for problem in &found {
            eprintln!("  {problem}");
        }
        problems.extend(found);
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
