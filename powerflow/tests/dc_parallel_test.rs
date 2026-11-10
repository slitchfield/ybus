// Tests for the parallel batches (feature `rayon`): par_solve_outages, par_ptdf_columns and
//   par_lodf_columns, each bit for bit its sequential version.
#![cfg(feature = "rayon")]

mod common;

use powerflow::case::Case;
use powerflow::dc::{DcModel, DcOptions, DcSolution, PtdfSlack};
use powerflow::powerflow_error::PowerflowError;

fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

fn same_solution(a: &DcSolution, b: &DcSolution) -> bool {
    same_bits(&a.va_deg, &b.va_deg)
        && same_bits(&a.pf_mw, &b.pf_mw)
        && same_bits(&a.pg_mw, &b.pg_mw)
        && a.solver == b.solver
}

#[test]
fn parallel_batches_match_sequential_bit_for_bit() {
    let Some(path) = common::matpower_case("case118") else {
        return; // not fetched; the skip was reported on stderr
    };
    let case = Case::from_path(&path).unwrap();
    let model = DcModel::new(&case).unwrap();
    let study = model.study(&DcOptions::default()).unwrap();
    let nl = case.branch.len();

    // Every single outage, and some pairs (islanding ones included)
    let mut outages: Vec<Vec<usize>> = (0..nl).map(|k| vec![k]).collect();
    outages.extend((0..nl).step_by(5).map(|k| vec![k, (k + 9) % nl]));
    let parallel = study.par_solve_outages(&case, &outages).unwrap();
    let sequential = study.solve_outages(&case, &outages).unwrap();
    let mut islanding = 0;
    for (rows, (p, s)) in outages.iter().zip(parallel.iter().zip(&sequential)) {
        match (p, s) {
            (Ok(p), Ok(s)) => assert!(same_solution(p, s), "{rows:?}"),
            (Err(p), Err(s)) => {
                assert_eq!(p, s, "{rows:?}");
                islanding += 1;
            }
            _ => panic!("{rows:?}: {p:?} against {s:?}"),
        }
    }
    assert!(islanding > 0);

    // PTDF and LODF columns, at a single and a distributed slack
    let buses: Vec<usize> = (0..case.bus.len()).collect();
    let branches: Vec<usize> = (0..nl).collect();
    let weights: Vec<f64> = case.bus.iter().map(|bus| bus.pd_mw).collect();
    for slack in [PtdfSlack::Reference, PtdfSlack::Weights(weights)] {
        let sensitivity = study.sensitivity(slack).unwrap();
        let (p, s) = (
            sensitivity.par_ptdf_columns(&buses).unwrap(),
            sensitivity.ptdf_columns(&buses).unwrap(),
        );
        assert!(p.iter().zip(&s).all(|(p, s)| same_bits(p, s)));
        let (p, s) = (
            sensitivity.par_lodf_columns(&branches).unwrap(),
            sensitivity.lodf_columns(&branches).unwrap(),
        );
        assert!(p.iter().zip(&s).all(|(p, s)| same_bits(p, s)));
    }
}

#[test]
fn parallel_errors_are_the_sequential_ones() {
    let Some(path) = common::matpower_case("case30") else {
        return;
    };
    let mut case = Case::from_path(&path).unwrap();
    // Two branches out of service: asking for both, the first is the error, as sequentially
    case.branch[3].status = false;
    case.branch[20].status = false;
    let model = DcModel::new(&case).unwrap();
    let study = model.study(&DcOptions::default()).unwrap();
    let sensitivity = study.sensitivity(PtdfSlack::Reference).unwrap();
    let rows = [0, 20, 5, 3, 7];
    assert_eq!(
        sensitivity.par_lodf_columns(&rows).unwrap_err(),
        PowerflowError::BranchOutOfService { row: 20 }
    );
    assert_eq!(
        sensitivity.lodf_columns(&rows).unwrap_err(),
        PowerflowError::BranchOutOfService { row: 20 }
    );
    // A different case fails the whole batch, as sequentially
    let mut other = case.clone();
    other.bus.pop();
    assert!(matches!(
        study.par_solve_outages(&other, &[vec![0]]),
        Err(PowerflowError::CaseMismatch { .. })
    ));
}
