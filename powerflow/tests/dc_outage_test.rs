// Tests for DC outages by rank-k updates (DcStudy::solve_outage): each against refactoring,
//   a new model of the case with the branches out of service, solved afresh, including
//   which error it reports; and single outages against the LODF.

mod common;

use common::dc_fixture::{fixture, B, G, L};
use powerflow::case::Case;
use powerflow::dc::{dc_power_flow, DcModel, DcSolution, DcSolver, PtdfSlack};
use powerflow::powerflow_error::PowerflowError;

/// The case with `rows` out of service, solved by a new model
fn refactored(case: &Case, rows: &[usize]) -> Result<DcSolution, PowerflowError> {
    let mut outaged = case.clone();
    for &row in rows {
        outaged.branch[row].status = false;
    }
    dc_power_flow(&outaged)
}

/// The largest difference between two solutions, relative to the largest magnitude of each
/// quantity in `expected` (angles, flows, generator outputs), at least 1.
fn relative_deviation(actual: &DcSolution, expected: &DcSolution) -> f64 {
    let worst = |a: &[f64], e: &[f64]| {
        let scale = e.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        a.iter()
            .zip(e)
            .fold(0.0_f64, |m, (a, e)| m.max((a - e).abs()))
            / scale
    };
    worst(&actual.va_deg, &expected.va_deg)
        .max(worst(&actual.pf_mw, &expected.pf_mw))
        .max(worst(&actual.pg_mw, &expected.pg_mw))
}

/// Checks the outage of `rows` against refactoring: the same error, or results within 1e-9
/// relative. Returns the deviation (0 for matching errors).
fn check(case: &Case, model: &DcModel, rows: &[usize]) -> f64 {
    let study = model.study().unwrap();
    let ours = study.solve_outage(case, rows);
    match (ours, refactored(case, rows)) {
        (Ok(ours), Ok(expected)) => {
            let deviation = relative_deviation(&ours, &expected);
            assert!(deviation <= 1e-9, "outage {rows:?}: {deviation:e}");
            assert_eq!(ours.solver, study.solver());
            for &row in rows {
                assert_eq!(ours.pf_mw[row], 0.0);
            }
            deviation
        }
        (Err(ours), Err(expected)) => {
            assert_eq!(ours, expected, "outage {rows:?}");
            0.0
        }
        (ours, expected) => panic!("outage {rows:?}: {ours:?} against {expected:?}"),
    }
}

/// Every set of up to three of `case`'s branch rows, checked against refactoring. Returns how
/// many were solved and how many were errors in both.
fn check_every_small_set(case: &Case) -> (usize, usize) {
    let model = DcModel::new(case).unwrap();
    let nl = case.branch.len();
    let mut sets = vec![];
    for i in 0..nl {
        sets.push(vec![i]);
        for j in i + 1..nl {
            sets.push(vec![i, j]);
            for k in j + 1..nl {
                sets.push(vec![i, j, k]);
            }
        }
    }
    let (mut solved, mut errors) = (0, 0);
    for rows in sets {
        // Out-of-service rows can't be outaged; refactoring would just ignore them
        if rows.iter().any(|&row| !case.branch[row].status) {
            let study = model.study().unwrap();
            let row = *rows.iter().find(|&&row| !case.branch[row].status).unwrap();
            assert_eq!(
                study.solve_outage(case, &rows).unwrap_err(),
                PowerflowError::BranchOutOfService { row }
            );
            continue;
        }
        match model.study().unwrap().solve_outage(case, &rows) {
            Ok(_) => solved += 1,
            Err(_) => errors += 1,
        }
        check(case, &model, &rows);
    }
    (solved, errors)
}

/// The ring of the sensitivity tests, plus a pendant bus: buses 1–4 in a ring with a chord
/// 1–3, a phase shifter on 2–3, a tap on 3–4, a parallel pair on 1–2, a pendant bus 6 off bus
/// 4 (a bridge), and bus 5 isolated with an out-of-service branch to bus 4. Bus 1 is the
/// reference; bus 3 has a second generator.
fn ring() -> (Vec<B>, Vec<G>, Vec<L>) {
    (
        vec![
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 40.0, 0.0, 0.0),
            B(3, 2, 30.0, 0.0, 0.0),
            B(4, 1, 50.0, 0.0, 0.0),
            B(5, 4, 10.0, 0.0, 0.0),
            B(6, 1, 15.0, 0.0, 0.0),
        ],
        vec![G(1, 0.0, true), G(3, 20.0, true)],
        vec![
            L(1, 2, 0.1, 0.0, 0.0, true),
            L(1, 2, 0.3, 0.0, 0.0, true),
            L(2, 3, 0.2, 0.0, 8.0, true),
            L(3, 4, 0.15, 1.05, 0.0, true),
            L(4, 1, 0.25, 0.0, 0.0, true),
            L(1, 3, 0.4, 0.0, 0.0, true),
            L(4, 5, 0.1, 0.0, 0.0, false),
            L(4, 6, 0.2, 0.0, -5.0, true),
        ],
    )
}

#[test]
fn every_small_outage_matches_refactoring() {
    let (buses, gens, branches) = ring();
    let case = fixture(&buses, &gens, &branches);
    let (solved, errors) = check_every_small_set(&case);
    // Taking out the pendant branch (row 7) strands bus 6, and some pairs and triples strand
    //   more; the rest solve
    assert!(solved > 0 && errors > 0, "{solved} solved, {errors} errors");
    eprintln!("ring: {solved} outages solved, {errors} islanding in both");
}

/// A radial feeder off the reference bus whose far end has the smallest label, so a cut near
/// the root strands buses whose smallest (in case order) is deep in the feeder; a cycle hanging
/// off the feeder, whose branches are no bridges; and a second reference bus behind a bridge,
/// which keeps its side reachable when that bridge's neighbours go.
#[test]
fn radial_feeders_and_several_references_match_refactoring() {
    let buses = vec![
        // Listed so the feeder's far end (bus 2) comes before its root side in case order
        B(2, 1, 5.0, 0.0, 0.0),
        B(1, 3, 0.0, 0.0, 0.0),
        B(5, 1, 10.0, 0.0, 0.0),
        B(4, 1, 8.0, 0.0, 0.0),
        B(3, 1, 6.0, 0.0, 0.0),
        B(9, 1, 7.0, 0.0, 0.0),
        B(10, 1, 4.0, 0.0, 0.0),
        B(7, 1, 3.0, 0.0, 0.0),
        B(8, 3, 0.0, 0.0, 0.0),
    ];
    let gens = vec![G(1, 0.0, true), G(8, 10.0, true)];
    let branches = vec![
        // The feeder 1 → 5 → 4 → 3 → 2
        L(1, 5, 0.1, 0.0, 0.0, true),
        L(5, 4, 0.1, 0.0, 0.0, true),
        L(4, 3, 0.1, 0.0, 0.0, true),
        L(3, 2, 0.1, 0.0, 0.0, true),
        // A cycle 4 → 9 → 10 → 4 off the feeder
        L(4, 9, 0.2, 0.0, 0.0, true),
        L(9, 10, 0.2, 0.0, 0.0, true),
        L(10, 4, 0.2, 0.0, 0.0, true),
        // Bus 7 off bus 5, and the second reference bus 8 behind the bridge 7–8
        L(5, 7, 0.15, 0.0, 0.0, true),
        L(7, 8, 0.15, 0.0, 0.0, true),
    ];
    let case = fixture(&buses, &gens, &branches);
    let (solved, errors) = check_every_small_set(&case);
    assert!(solved > 0 && errors > 0, "{solved} solved, {errors} errors");
    eprintln!("feeder: {solved} outages solved, {errors} islanding in both");
}

#[test]
fn outages_match_refactoring_on_the_lu_path() {
    // A series capacitor in place of the chord makes B′ indefinite
    let (buses, gens, mut branches) = ring();
    branches[5] = L(1, 3, -0.05, 0.0, 0.0, true);
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    assert_eq!(model.study().unwrap().solver(), DcSolver::Lu);
    check_every_small_set(&case);
}

#[test]
fn an_island_with_its_own_reference_is_solved() {
    // Bus 3 as a second reference: taking out 2–3, 3–4 and 1–3 leaves bus 3 on its own, with
    //   its own reference, which the DC power flow solves (as refactoring does)
    let (mut buses, gens, branches) = ring();
    buses[2].1 = 3;
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    assert_eq!(model.ref_buses(), [0, 2]);
    let study = model.study().unwrap();
    let alone = study.solve_outage(&case, &[2, 3, 5]).unwrap();
    // Alone, bus 3's generator is its island's slack, so it serves exactly its bus's 30 MW
    assert!((alone.pg_mw[1] - 30.0).abs() <= 1e-9, "{}", alone.pg_mw[1]);
    check(&case, &model, &[2, 3, 5]);
    check_every_small_set(&case);
}

#[test]
fn outage_rows_are_checked() {
    let (buses, gens, branches) = ring();
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    assert_eq!(
        study.solve_outage(&case, &[8]).unwrap_err(),
        PowerflowError::RowOutOfRange {
            table: "branch",
            row: 8,
            len: 8
        }
    );
    assert_eq!(
        study.solve_outage(&case, &[6]).unwrap_err(),
        PowerflowError::BranchOutOfService { row: 6 }
    );
    // A row given twice is one outage; no rows is the base case
    assert_eq!(
        study.solve_outage(&case, &[0, 0]).unwrap(),
        study.solve_outage(&case, &[0]).unwrap()
    );
    let base = study.solve(&case).unwrap();
    assert!(relative_deviation(&study.solve_outage(&case, &[]).unwrap(), &base) <= 1e-15);
    // A different case is refused before any outage
    let mut other = case.clone();
    other.bus.pop();
    assert!(matches!(
        study.solve_outages(&other, &[vec![0]]),
        Err(PowerflowError::CaseMismatch { .. })
    ));
}

#[test]
fn a_batch_matches_one_at_a_time() {
    let (buses, gens, branches) = ring();
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    let outages: Vec<Vec<usize>> = vec![vec![0], vec![7], vec![2, 4], vec![6], vec![1, 3, 5]];
    let batch = study.solve_outages(&case, &outages).unwrap();
    for (rows, result) in outages.iter().zip(batch) {
        assert_eq!(result, study.solve_outage(&case, rows), "{rows:?}");
    }
}

#[test]
fn single_outages_agree_with_the_lodf() {
    let Some(path) = common::matpower_case("case30") else {
        return;
    };
    let case = Case::from_path(&path).unwrap();
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    let base = study.solve(&case).unwrap();
    let sensitivity = study.sensitivity(PtdfSlack::Reference).unwrap();
    let mut checked = 0;
    for k in 0..case.branch.len() {
        let column = &sensitivity.lodf_columns(&[k]).unwrap()[0];
        if column[0].is_nan() {
            continue; // a bridge
        }
        let after = study.solve_outage(&case, &[k]).unwrap();
        for (l, lodf) in column.iter().enumerate() {
            let predicted = base.pf_mw[l] + lodf * base.pf_mw[k];
            assert!(
                (after.pf_mw[l] - predicted).abs() <= 1e-9,
                "outage {k}, branch {l}: {} vs {predicted}",
                after.pf_mw[l]
            );
        }
        checked += 1;
    }
    assert_eq!(checked, case.branch.len() - 3); // case30's three bridges
}

#[test]
fn outages_match_refactoring_on_the_committed_cases() {
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            return;
        };
        let case = Case::from_path(&path).unwrap();
        let model = DcModel::new(&case).unwrap();
        let nl = case.branch.len();
        let mut worst: f64 = 0.0;
        let (mut solved, mut islanding) = (0, 0);
        // Every single outage, and fixed pairs and triples spread over the branches
        let mut sets: Vec<Vec<usize>> = (0..nl).map(|k| vec![k]).collect();
        sets.extend((0..nl).map(|k| vec![k, (k + 7) % nl]));
        sets.extend(
            (0..nl)
                .step_by(3)
                .map(|k| vec![k, (k + 5) % nl, (k + 11) % nl]),
        );
        for rows in sets {
            let mut rows = rows;
            rows.sort_unstable();
            rows.dedup();
            match model.study().unwrap().solve_outage(&case, &rows) {
                Ok(_) => solved += 1,
                Err(_) => islanding += 1,
            }
            worst = worst.max(check(&case, &model, &rows));
        }
        eprintln!("{name}: {solved} outages within {worst:e} of refactoring, {islanding} islanding in both");
    }
}
