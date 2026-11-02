// Tests for the DC sensitivities: PTDF against MATPOWER's makePTDF on case14 and case30
//   (committed) and on every local case whose results exist, and against the DC power flow
//   itself on fixtures: a PTDF column is the change in flows per unit of power injected at a
//   bus and withdrawn at the slack.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::dc_fixture::{fixture, B, G, L};
use powerflow::case::{BusType, Case};
use powerflow::dc::{dc_power_flow, DcModel, DcSolver, PtdfSlack};
use powerflow::powerflow_error::PowerflowError;

/// A labelled dense matrix as `generate_reference.m` writes it: the column labels, and each
/// row's label with its values.
struct Labelled {
    columns: Vec<usize>,
    rows: Vec<(usize, Vec<f64>)>,
}

fn read_matrix(path: &Path) -> Option<Labelled> {
    let text = fs::read_to_string(path).ok()?;
    let mut lines = text.lines();
    let header = lines.next().unwrap();
    let columns = header
        .split(',')
        .skip(1)
        .map(|label| label.parse().unwrap())
        .collect();
    let rows = lines
        .map(|line| {
            let mut fields = line.split(',');
            let row = fields.next().unwrap().parse().unwrap();
            (row, fields.map(|value| value.parse().unwrap()).collect())
        })
        .collect();
    Some(Labelled { columns, rows })
}

/// The case bus row with each label
fn bus_row(case: &Case, label: usize) -> usize {
    case.bus
        .iter()
        .position(|bus| bus.bus_label == label)
        .unwrap()
}

/// Each bus's load as distributed-slack weights, as generate_reference.m gives makePTDF
fn load_weights(case: &Case) -> Vec<f64> {
    case.bus
        .iter()
        .map(|bus| {
            if bus.bus_type == BusType::Isolated {
                0.0
            } else {
                bus.pd_mw
            }
        })
        .collect()
}

/// The largest absolute difference between our PTDF columns and MATPOWER's labelled
/// columns (rows labelled by 1-based branch row, columns by bus number).
fn ptdf_deviation(case: &Case, reference: &Labelled, slack: PtdfSlack) -> f64 {
    let model = DcModel::new(case).unwrap();
    let study = model.study().unwrap();
    let sensitivity = study.sensitivity(slack).unwrap();
    let buses: Vec<usize> = reference
        .columns
        .iter()
        .map(|&label| bus_row(case, label))
        .collect();
    let columns = sensitivity.ptdf_columns(&buses).unwrap();
    let mut worst: f64 = 0.0;
    for (row, values) in &reference.rows {
        for (column, value) in columns.iter().zip(values) {
            worst = worst.max((column[row - 1] - value).abs());
        }
    }
    worst
}

#[test]
fn ptdf_matches_makeptdf_on_the_committed_cases() {
    let dir = common::reference::committed_dir();
    for name in ["case14", "case30"] {
        let Some(path) = common::matpower_case(name) else {
            return; // not fetched; the skip was reported on stderr
        };
        let case = Case::from_path(&path).unwrap();
        let reference = read_matrix(&dir.join(format!("{name}_ptdf.csv"))).unwrap();
        let single = ptdf_deviation(&case, &reference, PtdfSlack::Reference);
        let reference = read_matrix(&dir.join(format!("{name}_ptdf_distributed.csv"))).unwrap();
        let distributed =
            ptdf_deviation(&case, &reference, PtdfSlack::Weights(load_weights(&case)));
        eprintln!("{name}: PTDF within {single:e}, with distributed slack {distributed:e}");
        assert!(single <= 1e-12 && distributed <= 1e-12, "{name}");
    }
}

#[test]
fn rows_columns_transfers_and_the_full_matrix_agree() {
    let Some(path) = common::matpower_case("case30") else {
        return;
    };
    let case = Case::from_path(&path).unwrap();
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    for slack in [
        PtdfSlack::Reference,
        PtdfSlack::Bus(7),
        PtdfSlack::Weights(load_weights(&case)),
    ] {
        let sensitivity = study.sensitivity(slack.clone()).unwrap();
        let full = sensitivity.ptdf();
        assert_eq!(
            (full.len(), full[0].len()),
            (case.branch.len(), case.bus.len())
        );
        // Rows come from one solve per branch, the full matrix from one per bus
        let branches: Vec<usize> = (0..case.branch.len()).collect();
        let rows = sensitivity.ptdf_rows(&branches).unwrap();
        for (row, full_row) in rows.iter().zip(&full) {
            for (a, b) in row.iter().zip(full_row) {
                assert!((a - b).abs() <= 1e-13, "{slack:?}: {a} vs {b}");
            }
        }
        // A transfer is the difference of two columns, whatever the slack
        let transfers = sensitivity.ptdf_transfers(&[(3, 20), (0, 29)]).unwrap();
        for (transfer, (from, to)) in transfers.iter().zip([(3, 20), (0, 29)]) {
            for (k, value) in transfer.iter().enumerate() {
                let difference = full[k][from] - full[k][to];
                assert!((value - difference).abs() <= 1e-13, "{slack:?}");
            }
        }
    }
}

/// A PTDF column by the DC power flow: the change in each branch's flow (p.u.) when bus
/// `bus` (a case row) serves 1 p.u. less load, the slack generator making it up.
fn by_power_flow(case: &Case, bus: usize) -> Vec<f64> {
    let base = dc_power_flow(case).unwrap();
    let mut injected = case.clone();
    injected.bus[bus].pd_mw -= case.base_mva;
    let changed = dc_power_flow(&injected).unwrap();
    changed
        .pf_mw
        .iter()
        .zip(&base.pf_mw)
        .map(|(after, before)| (after - before) / case.base_mva)
        .collect()
}

fn assert_columns_close(actual: &[f64], expected: &[f64], what: &str) {
    for (k, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!((a - e).abs() <= 1e-12, "{what}, branch {k}: {a} vs {e}");
    }
}

/// Four buses: a ring 1–2–3–4 with a chord 1–3, a phase shifter on 2–3, a tap on 3–4, a
/// parallel pair on 1–2, and bus 5 isolated with an out-of-service branch to bus 4. Bus 1 is
/// the reference.
fn ring() -> (Vec<B>, Vec<G>, Vec<L>) {
    (
        vec![
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 40.0, 0.0, 0.0),
            B(3, 2, 30.0, 0.0, 0.0),
            B(4, 1, 50.0, 0.0, 0.0),
            B(5, 4, 10.0, 0.0, 0.0),
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
        ],
    )
}

#[test]
fn a_ptdf_column_is_the_change_in_flows() {
    let (buses, gens, mut branches) = ring();
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    assert_eq!(study.solver(), DcSolver::Cholesky);
    let sensitivity = study.sensitivity(PtdfSlack::Reference).unwrap();
    assert_eq!(sensitivity.slack_bus(), 0);
    let columns = sensitivity.ptdf_columns(&[1, 2, 3, 4, 0]).unwrap();
    for (bus, column) in [1, 2, 3].iter().zip(&columns) {
        assert_columns_close(column, &by_power_flow(&case, *bus), &format!("bus {bus}"));
    }
    // The isolated bus and the slack itself move nothing; the out-of-service branch never
    //   carries flow
    assert!(columns[3]
        .iter()
        .chain(&columns[4])
        .all(|&value| value == 0.0));
    assert!(columns.iter().all(|column| column[6] == 0.0));

    // A series capacitor in place of the chord makes B′ indefinite, so LU serves
    branches[5] = L(1, 3, -0.05, 0.0, 0.0, true);
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    assert_eq!(study.solver(), DcSolver::Lu);
    let sensitivity = study.sensitivity(PtdfSlack::Reference).unwrap();
    for bus in 1..4 {
        let column = sensitivity.ptdf_columns(&[bus]).unwrap();
        assert_columns_close(
            &column[0],
            &by_power_flow(&case, bus),
            &format!("LU, bus {bus}"),
        );
    }
}

#[test]
fn another_slack_bus_and_several_references_follow_makeptdf() {
    // Bus 3 as slack: the injection is withdrawn there. By the power flow, that's the case
    //   with bus 3 as the only reference
    let (mut buses, gens, branches) = ring();
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    let at_three = study.sensitivity(PtdfSlack::Bus(2)).unwrap();
    assert_eq!(at_three.slack_bus(), 2);
    buses[0].1 = 2;
    buses[2].1 = 3;
    let moved = fixture(&buses, &gens, &branches);
    for bus in [0, 1, 3] {
        let column = at_three.ptdf_columns(&[bus]).unwrap();
        assert_columns_close(
            &column[0],
            &by_power_flow(&moved, bus),
            &format!("bus {bus}"),
        );
    }

    // Two references (buses 1 and 3): makePTDF's slack is the first, and it treats bus 3 as
    //   an ordinary bus, so the PTDF is the one-reference case's
    let (mut buses, gens, branches) = ring();
    buses[2].1 = 3;
    let two = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&two).unwrap();
    assert_eq!(model.ref_buses(), [0, 2]);
    let study = model.study().unwrap();
    let sensitivity = study.sensitivity(PtdfSlack::Reference).unwrap();
    assert_eq!(sensitivity.slack_bus(), 0);
    for bus in [1, 2, 3] {
        let column = sensitivity.ptdf_columns(&[bus]).unwrap();
        assert_columns_close(
            &column[0],
            &by_power_flow(&case, bus),
            &format!("bus {bus}"),
        );
    }
}

#[test]
fn a_distributed_slack_withdraws_in_proportion() {
    // Weights 1 and 3 on buses 2 and 4: injecting at bus j and withdrawing a quarter at bus 2
    //   and three quarters at bus 4 is the single-slack column for j, less a quarter of bus
    //   2's and three quarters of bus 4's
    let (buses, gens, branches) = ring();
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    let single = study.sensitivity(PtdfSlack::Reference).unwrap();
    let weights = vec![0.0, 1.0, 0.0, 3.0, 0.0];
    let distributed = study.sensitivity(PtdfSlack::Weights(weights)).unwrap();
    let h = single.ptdf_columns(&[0, 1, 2, 3]).unwrap();
    for bus in 0..4 {
        let expected: Vec<f64> = (0..case.branch.len())
            .map(|k| h[bus][k] - 0.25 * h[1][k] - 0.75 * h[3][k])
            .collect();
        let column = distributed.ptdf_columns(&[bus]).unwrap();
        assert_columns_close(&column[0], &expected, &format!("bus {bus}"));
        // Rows agree with columns under a distributed slack too
        let rows = distributed.ptdf_rows(&[0, 2, 5]).unwrap();
        for (row, k) in rows.iter().zip([0, 2, 5]) {
            assert!((row[bus] - expected[k]).abs() <= 1e-12);
        }
    }
}

#[test]
fn unusable_slacks_and_rows_are_errors() {
    let (buses, gens, branches) = ring();
    let case = fixture(&buses, &gens, &branches);
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    let invalid = |slack| match study.sensitivity(slack) {
        Err(PowerflowError::InvalidSlack { .. }) => true,
        other => panic!("{other:?}"),
    };
    assert!(invalid(PtdfSlack::Bus(4))); // isolated
    assert!(invalid(PtdfSlack::Weights(vec![1.0; 4]))); // one short
    assert!(invalid(PtdfSlack::Weights(vec![0.0, 0.0, 0.0, 0.0, 1.0]))); // isolated bus
    assert!(invalid(PtdfSlack::Weights(vec![1.0, -1.0, 0.0, 0.0, 0.0]))); // sums to zero
    assert_eq!(
        study.sensitivity(PtdfSlack::Bus(5)).unwrap_err(),
        PowerflowError::RowOutOfRange {
            table: "bus",
            row: 5,
            len: 5
        }
    );
    let sensitivity = study.sensitivity(PtdfSlack::Reference).unwrap();
    assert_eq!(
        sensitivity.ptdf_rows(&[7]).unwrap_err(),
        PowerflowError::RowOutOfRange {
            table: "branch",
            row: 7,
            len: 7
        }
    );
    assert!(sensitivity.ptdf_columns(&[9]).is_err());
    assert!(sensitivity.ptdf_transfers(&[(0, 9)]).is_err());
}

#[test]
fn two_islands_need_a_slack_each() {
    // Buses 1–2 and 3–4, each with a reference: the DC power flow solves them, but a single
    //   slack can't balance both
    let case = fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 10.0, 0.0, 0.0),
            B(3, 3, 0.0, 0.0, 0.0),
            B(4, 1, 10.0, 0.0, 0.0),
        ],
        &[G(1, 0.0, true), G(3, 0.0, true)],
        &[L(1, 2, 0.1, 0.0, 0.0, true), L(3, 4, 0.1, 0.0, 0.0, true)],
    );
    let model = DcModel::new(&case).unwrap();
    let study = model.study().unwrap();
    assert_eq!(
        study.sensitivity(PtdfSlack::Reference).unwrap_err(),
        PowerflowError::IslandWithoutReference { bus_label: 3 }
    );
}

/// The local reference results' directory, beside the normalized cases, or `None` to skip.
fn local_reference_dir() -> Option<PathBuf> {
    let normalized = common::normalized_case("case14")?.parent()?.to_path_buf();
    let name = normalized
        .file_name()?
        .to_str()?
        .strip_suffix("-normalized")?;
    let dir = normalized.with_file_name(format!("{name}-reference"));
    dir.is_dir().then_some(dir)
}

#[test]
fn ptdf_matches_makeptdf_on_every_local_case() {
    let Some(dir) = local_reference_dir() else {
        eprintln!("skipping: no local reference results (scripts/generate-reference.sh)");
        return;
    };
    let normalized = common::normalized_case("case14").unwrap();
    let normalized = normalized.parent().unwrap();
    let mut names: Vec<String> = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|entry| {
            let name = entry.unwrap().file_name().into_string().unwrap();
            name.strip_suffix("_ptdf.csv").map(str::to_string)
        })
        .collect();
    names.sort();
    let mut worst = (0.0, String::new());
    for name in &names {
        let case = Case::from_path(normalized.join(format!("{name}.m"))).unwrap();
        let reference = read_matrix(&dir.join(format!("{name}_ptdf.csv"))).unwrap();
        let mut deviation = ptdf_deviation(&case, &reference, PtdfSlack::Reference);
        if let Some(reference) = read_matrix(&dir.join(format!("{name}_ptdf_distributed.csv"))) {
            let weights = PtdfSlack::Weights(load_weights(&case));
            deviation = deviation.max(ptdf_deviation(&case, &reference, weights));
        }
        assert!(
            deviation <= 1e-9,
            "{name}: PTDF differs from makePTDF's by {deviation:e}"
        );
        if deviation > worst.0 {
            worst = (deviation, name.clone());
        }
    }
    eprintln!(
        "{} cases' PTDF match makePTDF's; the largest difference is {:e} ({})",
        names.len(),
        worst.0,
        worst.1
    );
}
