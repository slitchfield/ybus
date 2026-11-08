//! The DC and AC optimal power flows on the predictor-corrector solver (`OpfSolver::Mehrotra`),
//! against MATPOWER's `rundcopf` and `runopf` (MIPS): on the committed reference cases, and on
//! any case named in `YBUS_DCOPF_CASES` or `YBUS_ACOPF_CASES` whose results
//! `scripts/generate-reference.sh` generated locally. The bar is phase O3's:
//!
//! - converge wherever MATPOWER's MIPS converges;
//! - the objective within 1e-6 (relative) of MATPOWER's;
//! - the bus prices close to MATPOWER's (`PRICE`, below);
//! - the KKT residuals no worse than MIPS's on the same case (`RESIDUAL_FACTOR`, below).
//!
//! Iterations are reported, not required, and pinned on the committed cases. With the KKT
//! matrix factored by the pivoted LDLᵀ (`KktFactorization::Ldlt`), both solvers are held to the
//! same bar on the committed cases, with their iterations reported only. The primal solution
//! and the bound multipliers aren't compared, since degenerate optima leave them non-unique, as
//! in the reference tests.

mod common;

use std::path::Path;

use common::reference::{committed_dir, read_acopf, read_dcopf, relative_difference};
use interior_point::solution::KktResiduals;
use powerflow::case::Case;
use powerflow::opf::ac::{ac_opf_with_options, AcOpfOptions, AcOpfSolution};
use powerflow::opf::dc::{dc_opf_with_options, DcOpfOptions, DcOpfSolution};
use powerflow::opf::{KktFactorization, OpfSolver};

/// The objective's bar, relative to its magnitude (at least 1)
const OBJECTIVE: f64 = 1e-6;

/// The bus prices' bar, relative to each price's magnitude (at least 1). Both solvers stop once
/// their measures are below 1e-6, by different paths, and prices are multipliers, determined
/// less tightly than the objective: the largest differences measured on the committed cases
/// are 9.7e-5 (`LAM_P`, case30 AC) and 3.5e-4 (`LAM_Q`, case30 AC), 2026-11-07.
const PRICE: f64 = 1e-3;

/// How many times MIPS's own KKT residual on the same case each of this solver's may be, above
/// the floors below (decided with the maintainer, 2026-11-07). The residuals are unscaled, so on
/// the AC OPF (`cost_mult` 1e-4) they're in $/hr and reach 0.5 for MIPS itself; MIPS's
/// `mu_threshold` zeroing of small multipliers, which both solvers apply, perturbs
/// stationarity too. So the bar is MIPS's own residual, not an absolute one. The largest
/// ratio measured on the committed cases is 51 (complementarity, case14 AC: 4.3e-3 against
/// MIPS's 8.4e-5, which is 4.3e-7 in the solver's own scale).
const RESIDUAL_FACTOR: f64 = 100.0;

/// The floors under which a residual always passes: the solvers' own tolerances, expressed in
/// the residuals' unscaled units. The termination measures are relative and scaled (feasibility
/// over `1 + ‖x‖∞`, complementarity over `1 + ‖x‖∞` and times `cost_mult`, the gradient times
/// `cost_mult` over `1 + cost_mult ‖multipliers‖∞`), so a converged solve's absolute residuals
/// can exceed the tolerances themselves: case1354pegase's AC feasibility is 7.3e-6 against a
/// `feastol` of 5e-6, with `‖x‖∞` above 1.
struct Floors {
    feasibility: f64,
    complementarity: f64,
    stationarity: f64,
}

impl Floors {
    fn new(x: &[f64], multipliers: f64, cost_mult: f64) -> Floors {
        let x_norm = x.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        Floors {
            feasibility: 5e-6 * (1.0 + x_norm),
            complementarity: 1e-6 * (1.0 + x_norm) / cost_mult,
            stationarity: 1e-6 * (1.0 + cost_mult * multipliers) / cost_mult,
        }
    }
}

/// The largest magnitude among a result's multipliers
fn largest_multiplier(m: &interior_point::solution::Multipliers) -> f64 {
    m.eq_nonlin
        .iter()
        .chain(&m.ineq_nonlin)
        .chain(&m.mu_l)
        .chain(&m.mu_u)
        .chain(&m.lower)
        .chain(&m.upper)
        .fold(0.0_f64, |largest, v| largest.max(v.abs()))
}

/// This solver's iterations on the committed cases (2026-11-07), pinned as in
/// `interior_point`'s tests: DC, then AC. The AC counts follow the damped corrector (phase O3's
/// PR 4); with the full one they were 7, 11 and 9.
fn expected_iterations(name: &str) -> (usize, usize) {
    match name {
        "case14" => (6, 9),
        "case30" => (3, 11),
        "case118" => (5, 11),
        _ => panic!("no expected iterations for {name}"),
    }
}

/// What one case's comparison found
struct Comparison {
    iterations: usize,
    matpower_iterations: usize,
    objective: f64,
    lam_p: f64,
    lam_q: f64,
    /// This solver's KKT residuals and MIPS's on the same case
    residuals: KktResiduals,
    mips_residuals: KktResiduals,
    /// The residuals that always pass, from this solver's point
    floors: Floors,
}

fn report(name: &str, model: &str, c: &Comparison) {
    let (r, m) = (&c.residuals, &c.mips_residuals);
    eprintln!(
        "{name} {model}: {} iterations (MATPOWER {}), objective {:.1e}, LAM_P {:.1e}, \
         LAM_Q {:.1e}; residuals (MIPS's) stationarity {:.1e} ({:.1e}), feasibility {:.1e} \
         ({:.1e}), complementarity {:.1e} ({:.1e})",
        c.iterations,
        c.matpower_iterations,
        c.objective,
        c.lam_p,
        c.lam_q,
        r.stationarity,
        m.stationarity,
        r.feasibility,
        m.feasibility,
        r.complementarity,
        m.complementarity,
    );
}

fn compare_dc(
    name: &str,
    case_path: &Path,
    dir: &Path,
    options: &DcOpfOptions,
) -> (DcOpfSolution, Comparison, bool) {
    let case = Case::from_path(case_path).unwrap();
    let reference = read_dcopf(dir, name, "dcopf").unwrap();
    let s = dc_opf_with_options(&case, options).unwrap();
    let c = Comparison {
        iterations: s.iterations,
        matpower_iterations: reference.iterations,
        objective: (s.objective - reference.objective).abs() / reference.objective.abs().max(1.0),
        lam_p: relative_difference(&s.lam_p, &reference.lam_p),
        lam_q: 0.0,
        residuals: s.residuals,
        mips_residuals: dc_opf_with_options(&case, &DcOpfOptions::default())
            .unwrap()
            .residuals,
        floors: Floors::new(
            s.solver.x(),
            largest_multiplier(s.solver.multipliers()),
            1.0,
        ),
    };
    (s, c, reference.success)
}

fn compare_ac(
    name: &str,
    case_path: &Path,
    dir: &Path,
    options: &AcOpfOptions,
) -> (AcOpfSolution, Comparison, bool) {
    let case = Case::from_path(case_path).unwrap();
    let reference = read_acopf(dir, name, "opf").unwrap();
    let s = ac_opf_with_options(&case, options).unwrap();
    let c = Comparison {
        iterations: s.iterations,
        matpower_iterations: reference.iterations,
        objective: (s.objective - reference.objective).abs() / reference.objective.abs().max(1.0),
        lam_p: relative_difference(&s.lam_p, &reference.lam_p),
        lam_q: relative_difference(&s.lam_q, &reference.lam_q),
        residuals: s.residuals,
        mips_residuals: ac_opf_with_options(&case, &AcOpfOptions::default())
            .unwrap()
            .residuals,
        floors: Floors::new(
            s.solver.x(),
            largest_multiplier(s.solver.multipliers()),
            1e-4,
        ),
    };
    (s, c, reference.success)
}

/// What the bar finds wrong with a converged-or-not comparison
fn problems(
    name: &str,
    model: &str,
    converged: bool,
    success: bool,
    c: &Comparison,
) -> Vec<String> {
    let mut found = vec![];
    if success && !converged {
        found.push(format!(
            "{name} {model}: didn't converge where MATPOWER did"
        ));
    }
    if success && converged {
        let (r, m) = (&c.residuals, &c.mips_residuals);
        let against_mips = |mips: f64, floor: f64| (RESIDUAL_FACTOR * mips).max(floor);
        for (what, value, bar) in [
            ("objective", c.objective, OBJECTIVE),
            ("LAM_P", c.lam_p, PRICE),
            ("LAM_Q", c.lam_q, PRICE),
            (
                "stationarity",
                r.stationarity,
                against_mips(m.stationarity, c.floors.stationarity),
            ),
            (
                "feasibility",
                r.feasibility,
                against_mips(m.feasibility, c.floors.feasibility),
            ),
            (
                "complementarity",
                r.complementarity,
                against_mips(m.complementarity, c.floors.complementarity),
            ),
            (
                "dual feasibility",
                r.dual_feasibility,
                against_mips(m.dual_feasibility, 1e-6),
            ),
        ] {
            if value.is_nan() || value > bar {
                found.push(format!("{name} {model}: {what} {value:e} over {bar:e}"));
            }
        }
    }
    found
}

#[test]
fn the_defaults_run_mips() {
    assert!(matches!(DcOpfOptions::default().solver, OpfSolver::Mips(_)));
    assert!(matches!(AcOpfOptions::default().solver, OpfSolver::Mips(_)));
    let OpfSolver::Mehrotra(dc) = DcOpfOptions::mehrotra().solver else {
        panic!("DcOpfOptions::mehrotra chooses the predictor-corrector method");
    };
    assert_eq!(dc.feastol, 5e-6);
    // The DC OPF keeps LU and the AC OPF factors by LDLᵀ, without the inertia correction, as
    //   the survey decided (docs/benchmarks/ldlt-opf-2026-11-08.md)
    assert_eq!(dc.kkt_factorization, KktFactorization::Lu);
    let OpfSolver::Mehrotra(ac) = AcOpfOptions::mehrotra().solver else {
        panic!("AcOpfOptions::mehrotra chooses the predictor-corrector method");
    };
    assert_eq!((ac.feastol, ac.cost_mult), (5e-6, 1e-4));
    assert_eq!(ac.kkt_factorization, KktFactorization::Ldlt);
    assert!(!ac.inertia_correction);
}

#[test]
fn mips_analyzes_the_opf_once() {
    // MIPS's OPF keeps its KKT pattern too, so its analysis is kept for the whole solve
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();
        let dc = dc_opf_with_options(&case, &DcOpfOptions::default()).unwrap();
        let ac = ac_opf_with_options(&case, &AcOpfOptions::default()).unwrap();
        assert_eq!(dc.solver.analyses(), 1, "{name} DC");
        assert_eq!(ac.solver.analyses(), 1, "{name} AC");
    }
}

#[test]
fn dc_residuals_are_small_on_both_solvers() {
    // Both solvers' residuals come from the OPF's own kkt_residuals, so the bar relative to
    //   MIPS can't catch a fault there; the DC OPF (cost_mult 1) has small absolute residuals
    //   on both, at most 1.9e-6 measured (MIPS on case30), which anchors it
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();
        for options in [DcOpfOptions::default(), DcOpfOptions::mehrotra()] {
            let s = dc_opf_with_options(&case, &options).unwrap();
            assert!(s.residuals.max() < 1e-5, "{name}: {:?}", s.residuals);
        }
    }
}

#[test]
fn committed_cases_meet_the_bar() {
    let mut found = vec![];
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let (dc, ac) = expected_iterations(name);
        let (s, c, success) = compare_dc(name, &path, &committed_dir(), &DcOpfOptions::mehrotra());
        report(name, "DC", &c);
        found.extend(problems(name, "DC", s.converged(), success, &c));
        assert_eq!(c.iterations, dc, "{name} DC: iterations");
        // The OPF's matrices keep their patterns, so one analysis serves the whole solve
        assert_eq!(s.solver.analyses(), 1, "{name} DC: analyses");
        let (s, c, success) = compare_ac(name, &path, &committed_dir(), &AcOpfOptions::mehrotra());
        report(name, "AC", &c);
        found.extend(problems(name, "AC", s.converged(), success, &c));
        assert_eq!(c.iterations, ac, "{name} AC: iterations");
        assert_eq!(s.solver.analyses(), 1, "{name} AC: analyses");
    }
    assert!(found.is_empty(), "{}", found.join("\n"));
}

/// The local reference directory beside a normalized case
fn local_dir(path: &Path) -> String {
    path.parent()
        .and_then(|dir| dir.to_str())
        .map(|dir| dir.replace("-normalized", "-reference"))
        .unwrap()
}

/// Local comparisons: `YBUS_DCOPF_CASES` and `YBUS_ACOPF_CASES` name cases, separated by
/// spaces, as for the reference tests. Every case is checked and the differences reported
/// together; run with `--nocapture` to follow each case.
///
/// Known (2026-11-07): case300's AC `LAM_Q` differs by 6.4e-3, over `PRICE`, at bus 9533, whose
/// voltage both solvers leave just below its upper limit (1.0583 and 1.0593 against 1.06)
/// within their tolerance; the price is 0.07 $/MVAr-hr, so the difference is 0.0065. Its
/// objective agrees to 1.2e-8. Whether `LAM_Q` needs a bar of its own is left to the survey of
/// every case (phase O3's PR 4), as decided with the maintainer.
#[test]
fn local_cases_meet_the_bar() {
    let mut found = vec![];
    for (variable, model) in [("YBUS_DCOPF_CASES", "DC"), ("YBUS_ACOPF_CASES", "AC")] {
        let Ok(names) = std::env::var(variable) else {
            continue;
        };
        for name in names.split_whitespace() {
            let path = common::normalized_case(name)
                .unwrap_or_else(|| panic!("{name}: no normalized case"));
            let dir = local_dir(&path);
            let start = std::time::Instant::now();
            let (converged, success, c) = if model == "DC" {
                let (s, c, success) =
                    compare_dc(name, &path, Path::new(&dir), &DcOpfOptions::mehrotra());
                (s.converged(), success, c)
            } else {
                let (s, c, success) =
                    compare_ac(name, &path, Path::new(&dir), &AcOpfOptions::mehrotra());
                (s.converged(), success, c)
            };
            report(name, model, &c);
            eprintln!(
                "  converged {converged} (MATPOWER {success}) in {:.1} s",
                start.elapsed().as_secs_f64()
            );
            found.extend(problems(name, model, converged, success, &c));
        }
    }
    assert!(found.is_empty(), "{}", found.join("\n"));
}

#[test]
fn committed_cases_meet_the_bar_with_ldlt() {
    let mut found = vec![];
    // Only LDLᵀ refines, so its solves show it was used: (MIPS, predictor-corrector)
    let mut refinement_steps = (0, 0);
    let mut solved = 0;
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        solved += 1;
        for (solver, mut dc, mut ac) in [
            ("MIPS", DcOpfOptions::default(), AcOpfOptions::default()),
            (
                "predictor-corrector",
                DcOpfOptions::mehrotra(),
                AcOpfOptions::mehrotra(),
            ),
        ] {
            dc.solver = dc.solver.with_kkt_factorization(KktFactorization::Ldlt);
            ac.solver = ac.solver.with_kkt_factorization(KktFactorization::Ldlt);
            let (s, c, success) = compare_dc(name, &path, &committed_dir(), &dc);
            report(name, &format!("DC, {solver}, LDLᵀ"), &c);
            found.extend(problems(name, "DC", s.converged(), success, &c));
            assert_eq!(s.solver.analyses(), 1, "{name} DC, {solver}: analyses");
            let (s, c, success) = compare_ac(name, &path, &committed_dir(), &ac);
            report(name, &format!("AC, {solver}, LDLᵀ"), &c);
            found.extend(problems(name, "AC", s.converged(), success, &c));
            assert_eq!(s.solver.analyses(), 1, "{name} AC, {solver}: analyses");
            let steps = s.solver.refinement_steps();
            if solver == "MIPS" {
                refinement_steps.0 += steps;
            } else {
                refinement_steps.1 += steps;
            }
        }
    }
    eprintln!("refinement steps over the AC solves: {refinement_steps:?}");
    // Without MATPOWER's cases (an unfetched checkout) nothing was solved, as the other tests
    //   skip; YBUS_REQUIRE_CASES makes that a failure in matpower_case itself
    assert!(solved == 0 || (refinement_steps.0 > 0 && refinement_steps.1 > 0));
    assert!(found.is_empty(), "{}", found.join("\n"));
}
