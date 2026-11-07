//! The DC power flow's GPU path (feature `cuda`) against the CPU's: batched solves and PTDF
//! columns and transfers. The results
//! must agree to within `BAR`, relative to the largest magnitude of each quantity, and errors
//! must be the same. The tests skip, saying so, on a machine without CUDA unless
//! `YBUS_REQUIRE_CUDA=1`, and without a case unless `YBUS_REQUIRE_CASES=1`.

#![cfg(feature = "cuda")]

mod common;

use powerflow::case::{BusType, Case};
use powerflow::dc::{DcModel, DcSolution, DcSolver, PtdfSlack};
use powerflow::powerflow_error::PowerflowError;
use powerflow::Gpu;

/// Cholesky cases, and case300, whose negative reactances make B′ indefinite, so it's LU
const CASES: &[(&str, DcSolver)] = &[
    ("case14", DcSolver::Cholesky),
    ("case118", DcSolver::Cholesky),
    ("case300", DcSolver::Lu),
    ("case2383wp", DcSolver::Cholesky),
];

const BAR: f64 = 1e-10;

fn gpu(test: &str) -> Option<Gpu> {
    let gpu = Gpu::new().expect("CUDA is present but failed to start");
    if gpu.is_none() {
        if std::env::var("YBUS_REQUIRE_CUDA").is_ok_and(|v| v == "1") {
            panic!("{test}: no CUDA device, and YBUS_REQUIRE_CUDA=1");
        }
        eprintln!("{test}: skipped, no CUDA device");
    }
    gpu
}

/// The largest difference between `a` and `b`, relative to `b`'s largest magnitude
fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let scale = b
        .iter()
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    a.iter()
        .zip(b)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
        / scale
}

/// Panics unless the GPU's solution agrees with the CPU's
fn assert_agree(what: &str, gpu: &DcSolution, cpu: &DcSolution) {
    assert_eq!(gpu.solver, cpu.solver, "{what}: solver");
    for (name, g, c) in [
        ("va_deg", &gpu.va_deg, &cpu.va_deg),
        ("pf_mw", &gpu.pf_mw, &cpu.pf_mw),
        ("pg_mw", &gpu.pg_mw, &cpu.pg_mw),
    ] {
        let difference = relative_difference(g, c);
        assert!(difference < BAR, "{what}: {name} differs by {difference:e}");
    }
}

/// `case` with every load scaled by `factor`
fn scaled(case: &Case, factor: f64) -> Case {
    let mut scaled = case.clone();
    for bus in &mut scaled.bus {
        bus.pd_mw *= factor;
    }
    scaled
}

#[test]
fn batched_solves_match_the_cpu() {
    let Some(gpu) = gpu("batched_solves_match_the_cpu") else {
        return;
    };
    for &(name, solver) in CASES {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();
        let model = DcModel::new(&case).unwrap();
        let study = model.study().unwrap();
        assert_eq!(study.solver(), solver, "{name}");
        let variants: Vec<Case> = (0..7)
            .map(|i| scaled(&case, 0.7 + 0.1 * i as f64))
            .collect();
        // A case with a branch fewer doesn't match the model
        let mut mismatched = case.clone();
        mismatched.branch.pop();
        // Blocks of 3, so the batch spans blocks and repeats a case within one
        let mut cases: Vec<&Case> = variants.iter().collect();
        cases.insert(2, &mismatched);
        cases.push(&variants[0]);
        let mut on_gpu = study.gpu(&gpu).unwrap();
        on_gpu.set_block_width(3);
        assert_eq!(on_gpu.block_width(), 3);
        let results = on_gpu.solve_many(&cases).unwrap();
        assert_eq!(results.len(), cases.len());
        for (i, (result, case)) in results.iter().zip(&cases).enumerate() {
            let what = format!("{name} case {i}");
            match (result, study.solve(case)) {
                (Ok(gpu), Ok(cpu)) => assert_agree(&what, gpu, &cpu),
                (Err(gpu), Err(cpu)) => assert_eq!(gpu, &cpu, "{what}"),
                (gpu, cpu) => panic!("{what}: GPU {gpu:?}, CPU {cpu:?}"),
            }
        }
        assert!(matches!(
            results[2],
            Err(PowerflowError::CaseMismatch { .. })
        ));
        eprintln!("{name}: {} batched solves agree", cases.len());
    }
}

#[test]
fn a_network_of_reference_buses_needs_no_gpu_solve() {
    let Some(gpu) = gpu("a_network_of_reference_buses_needs_no_gpu_solve") else {
        return;
    };
    let Some(path) = common::matpower_case("case14") else {
        return;
    };
    // Every bus with a generator and of type REF: nothing is left to solve
    let mut case = Case::from_path(&path).unwrap();
    let template = case.generator[0].clone();
    for bus in &mut case.bus {
        bus.bus_type = powerflow::case::BusType::Ref;
    }
    let labels: Vec<usize> = case.bus.iter().map(|bus| bus.bus_label).collect();
    case.generator = labels
        .into_iter()
        .map(|label| {
            let mut generator = template.clone();
            generator.bus_label = label;
            generator
        })
        .collect();
    case.gencost = None;
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    assert_eq!(study.solver(), DcSolver::NoSolve);
    let mut on_gpu = study.gpu(&gpu).unwrap();
    assert_eq!(on_gpu.block_width(), 0);
    let results = on_gpu.solve_many(&[&case]).unwrap();
    assert_eq!(results[0].as_ref().unwrap(), &study.solve(&case).unwrap());
}

/// Panics unless every GPU column agrees with the CPU's, relative to the CPU column's largest
/// magnitude (a column of zeros must be zeros)
fn assert_columns_agree(what: &str, gpu: &[Vec<f64>], cpu: &[Vec<f64>]) {
    assert_eq!(gpu.len(), cpu.len(), "{what}: columns");
    for (j, (g, c)) in gpu.iter().zip(cpu).enumerate() {
        let difference = relative_difference(g, c);
        assert!(
            difference < BAR,
            "{what}: column {j} differs by {difference:e}"
        );
    }
}

/// The slacks to test on `case`: MATPOWER's, another bus (so the sensitivity factors its own
/// reduced Bbus), and one distributed in proportion to the loads
fn slacks(case: &Case, model: &DcModel) -> Vec<(&'static str, PtdfSlack)> {
    let other = model.bus_rows()[model.bus_rows().len() / 2];
    let loads: Vec<f64> = case.bus.iter().map(|bus| bus.pd_mw.max(0.0)).collect();
    vec![
        ("reference", PtdfSlack::Reference),
        ("another bus", PtdfSlack::Bus(other)),
        ("distributed", PtdfSlack::Weights(loads)),
    ]
}

#[test]
fn ptdf_columns_and_transfers_match_the_cpu() {
    let Some(gpu) = gpu("ptdf_columns_and_transfers_match_the_cpu") else {
        return;
    };
    for &(name, _) in CASES {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();
        let model = DcModel::new(&case).unwrap();
        let study = model.study().unwrap();
        let nbus = case.bus.len();
        let all: Vec<usize> = (0..nbus).collect();
        for (slack_name, slack) in slacks(&case, &model) {
            let what = format!("{name}, {slack_name} slack");
            let sensitivity = study.sensitivity(slack).unwrap();
            let mut on_gpu = sensitivity.gpu(&gpu).unwrap();
            // Several blocks, the last one partial
            on_gpu.set_block_width(64);
            assert_columns_agree(
                &format!("{what} columns"),
                &on_gpu.ptdf_columns(&all).unwrap(),
                &sensitivity.ptdf_columns(&all).unwrap(),
            );
            let slack_bus = sensitivity.slack_bus();
            let mut transfers: Vec<(usize, usize)> =
                (0..nbus).map(|i| (i, (i * 7 + 3) % nbus)).collect();
            transfers.extend([
                (0, 0),
                (slack_bus, 1),
                (1, slack_bus),
                (slack_bus, slack_bus),
            ]);
            // Transfers, then columns again: the offset goes back and forth
            assert_columns_agree(
                &format!("{what} transfers"),
                &on_gpu.ptdf_transfers(&transfers).unwrap(),
                &sensitivity.ptdf_transfers(&transfers).unwrap(),
            );
            assert_columns_agree(
                &format!("{what} columns again"),
                &on_gpu.ptdf_columns(&all[..10]).unwrap(),
                &sensitivity.ptdf_columns(&all[..10]).unwrap(),
            );
            assert_eq!(
                on_gpu.ptdf_columns(&[nbus]).unwrap_err(),
                sensitivity.ptdf_columns(&[nbus]).unwrap_err(),
            );
        }
        eprintln!("{name}: PTDF columns and transfers agree");
    }
}

#[test]
fn an_isolated_bus_gets_the_cpus_ptdf() {
    let Some(gpu) = gpu("an_isolated_bus_gets_the_cpus_ptdf") else {
        return;
    };
    let Some(path) = common::matpower_case("case14") else {
        return;
    };
    // Bus 8 hangs off bus 7 alone, so isolating it strands nothing
    let mut case = Case::from_path(&path).unwrap();
    let row = case.bus.iter().position(|bus| bus.bus_label == 8).unwrap();
    case.bus[row].bus_type = BusType::Isolated;
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    let all: Vec<usize> = (0..case.bus.len()).collect();
    let transfers = [(row, 0), (0, row), (row, row)];
    for (slack_name, slack) in slacks(&case, &model) {
        let sensitivity = study.sensitivity(slack).unwrap();
        let mut on_gpu = sensitivity.gpu(&gpu).unwrap();
        let columns = on_gpu.ptdf_columns(&all).unwrap();
        assert!(columns[row].iter().all(|&v| v == 0.0), "{slack_name}");
        assert_columns_agree(
            slack_name,
            &columns,
            &sensitivity.ptdf_columns(&all).unwrap(),
        );
        assert_columns_agree(
            slack_name,
            &on_gpu.ptdf_transfers(&transfers).unwrap(),
            &sensitivity.ptdf_transfers(&transfers).unwrap(),
        );
    }
}
