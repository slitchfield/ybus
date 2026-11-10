//! The DC power flow's GPU path (feature `cuda`) against the CPU's: batched solves, PTDF
//! columns and transfers, and single-branch outages, solved whole and screened. The results
//! must agree to within `BAR`, relative to the largest magnitude of each quantity, and errors
//! must be the same. The tests skip, saying so, on a machine without CUDA unless
//! `YBUS_REQUIRE_CUDA=1`, and without a case unless `YBUS_REQUIRE_CASES=1`.

#![cfg(feature = "cuda")]

mod common;

use common::dc_fixture::{fixture, B, G, L};

use powerflow::case::{BusType, Case};
use powerflow::dc::{DcModel, DcOptions, DcSolution, DcSolver, OutageLoading, PtdfSlack};
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
        let study = model.study(&DcOptions::default()).unwrap();
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
    let study = model.study(&DcOptions::default()).unwrap();
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
        let study = model.study(&DcOptions::default()).unwrap();
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
    let study = model.study(&DcOptions::default()).unwrap();
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

/// The most loaded rated branch of a CPU solution, as `OutageLoading` defines it, computed here
/// independently: (case branch row, |Pf| / rate_a)
fn cpu_loading(case: &Case, solution: &DcSolution) -> (Option<usize>, f64) {
    let mut best = (None, 0.0);
    for (row, branch) in case.branch.iter().enumerate() {
        if branch.rate_a_mva > 0.0 {
            let loading = solution.pf_mw[row].abs() / branch.rate_a_mva;
            if loading > best.1 {
                best = (Some(row), loading);
            }
        }
    }
    best
}

/// Panics unless the GPU's screening of an outage agrees with the CPU's solution: the same
/// loading to within rounding, and the same branch unless another's loading is as close
fn assert_loading_agrees(what: &str, case: &Case, gpu: &OutageLoading, cpu: &DcSolution) {
    let (branch, loading) = cpu_loading(case, cpu);
    let tolerance = BAR * loading.max(1.0);
    assert!(
        (gpu.loading - loading).abs() <= tolerance,
        "{what}: loading {} against {loading}",
        gpu.loading
    );
    if gpu.branch != branch {
        let row = gpu
            .branch
            .unwrap_or_else(|| panic!("{what}: no branch, CPU {branch:?}"));
        let theirs = cpu.pf_mw[row].abs() / case.branch[row].rate_a_mva;
        assert!(
            (theirs - loading).abs() <= tolerance,
            "{what}: branch {row} ({theirs}) against {branch:?} ({loading})"
        );
    }
}

/// Every single-branch outage of `case` (and a row past the table) on the GPU, solved whole
/// and screened, against the CPU's `solve_outages`. Returns how many solved and how many
/// were errors in both.
fn check_outages(what: &str, gpu: &Gpu, case: &Case, width: usize) -> (usize, usize) {
    let model = DcModel::new(case).unwrap();
    let study = model.study(&DcOptions::default()).unwrap();
    let mut rows: Vec<usize> = (0..case.branch.len()).collect();
    rows.insert(rows.len() / 2, case.branch.len());
    let sets: Vec<Vec<usize>> = rows.iter().map(|&row| vec![row]).collect();
    let cpu = study.solve_outages(case, &sets).unwrap();
    let mut on_gpu = study.gpu(gpu).unwrap();
    on_gpu.set_block_width(width);
    let solved = on_gpu.solve_branch_outages(case, &rows).unwrap();
    let screened = on_gpu.outage_loadings(case, &rows).unwrap();
    assert_eq!((solved.len(), screened.len()), (rows.len(), rows.len()));
    let (mut ok, mut errors) = (0, 0);
    for (((row, cpu), solved), screened) in rows.iter().zip(&cpu).zip(&solved).zip(&screened) {
        let what = format!("{what}, outage of row {row}");
        match (cpu, solved, screened) {
            (Ok(cpu), Ok(solved), Ok(screened)) => {
                assert_agree(&what, solved, cpu);
                assert_eq!(solved.pf_mw[*row], 0.0, "{what}");
                assert_loading_agrees(&what, case, screened, cpu);
                ok += 1;
            }
            (Err(cpu), Err(solved), Err(screened)) => {
                assert_eq!((solved, screened), (cpu, cpu), "{what}");
                errors += 1;
            }
            (cpu, solved, screened) => {
                panic!("{what}: CPU {cpu:?}, GPU {solved:?}, screened {screened:?}")
            }
        }
    }
    (ok, errors)
}

#[test]
fn outages_match_the_cpu() {
    let Some(gpu) = gpu("outages_match_the_cpu") else {
        return;
    };
    for &(name, _) in CASES {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();
        let shifters = case.branch.iter().filter(|b| b.shift_deg != 0.0).count();
        let (ok, errors) = check_outages(name, &gpu, &case, 16);
        assert!(ok > 0, "{name}");
        eprintln!("{name}: {ok} outages solved and {errors} errors agree ({shifters} shifters)");
    }
}

#[test]
fn outages_between_two_references_match_the_cpu() {
    let Some(gpu) = gpu("outages_between_two_references_match_the_cpu") else {
        return;
    };
    let Some(path) = common::matpower_case("case14") else {
        return;
    };
    // Bus 2 (with a generator) as a second reference: branch 1–2 then touches no unknown, and
    //   the branches from bus 2 couple to a reference end
    let mut case = Case::from_path(&path).unwrap();
    let row = case.bus.iter().position(|bus| bus.bus_label == 2).unwrap();
    case.bus[row].bus_type = BusType::Ref;
    case.bus[row].va_deg = -3.0;
    // Width 1, the smallest block
    let (ok, errors) = check_outages("case14 with two references", &gpu, &case, 1);
    eprintln!("case14 with two references: {ok} solved and {errors} errors agree");
}

#[test]
fn a_singular_outage_is_left_to_the_cpu() {
    let Some(gpu) = gpu("a_singular_outage_is_left_to_the_cpu") else {
        return;
    };
    // Lines 1–2 and 2–3, a series capacitor on 1–3 that cancels them exactly, and a last line
    //   1–3: without the last, B′ is singular though every bus still reaches the reference.
    //   With b = 10, 10 and −5 the GPU's update has w exactly zero, which it flags; with
    //   x = 0.1, 0.11 and −(0.1 + 0.11) the reciprocals round, so w is tiny instead and the
    //   angles come out near 1e14 radians
    for (x12, x23, x13) in [(0.1, 0.1, -0.2), (0.1, 0.11, -(0.1 + 0.11))] {
        let buses = [
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 50.0, 0.0, 0.0),
            B(3, 1, 50.0, 0.0, 0.0),
        ];
        let branches = [
            L(1, 2, x12, 0.0, 0.0, true),
            L(2, 3, x23, 0.0, 0.0, true),
            L(1, 3, x13, 0.0, 0.0, true),
            L(1, 3, 0.1, 0.0, 0.0, true),
        ];
        let case = fixture(&buses, &[G(1, 100.0, true)], &branches);
        let model = DcModel::new(&case).unwrap();
        let study = model.study(&DcOptions::default()).unwrap();
        assert!(matches!(
            study.solve_outage(&case, &[3]),
            Err(PowerflowError::SingularNetwork { .. })
        ));
        let what = format!("singular fixture {x12}, {x23}, {x13}");
        let (ok, errors) = check_outages(&what, &gpu, &case, 4);
        assert!(ok > 0 && errors > 0, "{what}: {ok} solved, {errors} errors");
    }
}
