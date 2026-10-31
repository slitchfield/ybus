// Tests for the DC model (Bbus assembly): matrices, phase-shifter injections, which elements
//   are dropped, bus roles and the errors. Expected values are worked out by hand.

mod common;
use common::dc_fixture::{fixture, B, G, L};
use powerflow::dc::DcModel;
use powerflow::powerflow_error::PowerflowError;

use std::f64::consts::PI;

/// Three buses in a triangle: 1–2 (x = 0.1), 2–3 (x = 0.2) and 1–3 (x = 0.25), so the
/// susceptances are 10, 5 and 4. Bus 1 is the reference, with the only generator.
fn triangle() -> (Vec<B>, Vec<G>, Vec<L>) {
    (
        vec![
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 50.0, 0.0, 0.0),
            B(3, 1, 30.0, 0.0, 0.0),
        ],
        vec![G(1, 80.0, true)],
        vec![
            L(1, 2, 0.1, 0.0, 0.0, true),
            L(2, 3, 0.2, 0.0, 0.0, true),
            L(1, 3, 0.25, 0.0, 0.0, true),
        ],
    )
}

fn model(tables: &(Vec<B>, Vec<G>, Vec<L>)) -> Result<DcModel, PowerflowError> {
    DcModel::new(&fixture(&tables.0, &tables.1, &tables.2))
}

#[test]
fn triangle_matrices() {
    let model = model(&triangle()).unwrap();
    assert_eq!(
        model.bbus().to_dense(),
        vec![
            vec![14.0, -10.0, -4.0],
            vec![-10.0, 15.0, -5.0],
            vec![-4.0, -5.0, 9.0]
        ]
    );
    assert_eq!(
        model.bf().to_dense(),
        vec![
            vec![10.0, -10.0, 0.0],
            vec![0.0, 5.0, -5.0],
            vec![4.0, 0.0, -4.0]
        ]
    );
    assert_eq!(model.pbusinj_pu(), [0.0; 3]);
    assert_eq!(model.pfinj_pu(), [0.0; 3]);
    assert_eq!(model.base_mva(), 100.0);
    assert_eq!(model.bus_rows(), [0, 1, 2]);
    assert_eq!(model.branch_rows(), [0, 1, 2]);
    assert_eq!(model.gen_rows(), [0]);
    assert_eq!(model.ref_buses(), [0]);
    assert!(model.pv_buses().is_empty());
    assert_eq!(model.pq_buses(), [1, 2]);
}

#[test]
fn parallel_branches_sum() {
    // Two x = 0.2 lines are one x = 0.1 line. 1 / 0.2 = 5 exactly, so the sums are exact.
    let mut doubled = triangle();
    doubled.2[0] = L(1, 2, 0.2, 0.0, 0.0, true);
    doubled.2.push(L(2, 1, 0.2, 0.0, 0.0, true));
    let doubled = model(&doubled).unwrap();
    let single = model(&triangle()).unwrap();
    assert_eq!(doubled.bbus().to_dense(), single.bbus().to_dense());
    assert_eq!(doubled.bf().shape(), (4, 3));
    assert_eq!(doubled.bf().to_dense()[3], vec![-5.0, 5.0, 0.0]);
}

#[test]
fn bbus_is_exactly_symmetric_with_many_parallel_branches() {
    // cholesky requires exact symmetry. Values whose sums round differently by order catch
    //   a stamp order that differs between (i, j) and (j, i).
    let mut tables = triangle();
    for x in [0.3, 0.7, 0.11, 0.13, 0.17, 0.19] {
        tables.2.push(L(1, 2, x, 0.0, 0.0, true));
        tables.2.push(L(3, 2, x, 0.0, 0.0, true));
    }
    assert!(model(&tables).unwrap().bbus().is_symmetric());
}

#[test]
fn tap_ratio_divides_susceptance_and_zero_means_one() {
    let mut tables = triangle();
    tables.2[0] = L(1, 2, 0.1, 2.0, 0.0, true);
    let model_tap = model(&tables).unwrap();
    assert_eq!(model_tap.bf().to_dense()[0], vec![5.0, -5.0, 0.0]);
    tables.2[0] = L(1, 2, 0.1, 1.0, 0.0, true);
    let model_one = model(&tables).unwrap();
    assert_eq!(
        model_one.bbus().to_dense(),
        model(&triangle()).unwrap().bbus().to_dense()
    );
}

#[test]
fn phase_shifter_injections() {
    // pfinj = b · (−shift · π/180) at the from end; pbusinj = Cftᵀ · pfinj
    let mut tables = triangle();
    tables.2[0] = L(1, 2, 0.1, 0.0, 10.0, true);
    let model = model(&tables).unwrap();
    let expected = 10.0 * (-10.0 * PI / 180.0);
    assert_eq!(model.pfinj_pu(), [expected, 0.0, 0.0]);
    assert_eq!(model.pbusinj_pu(), [expected, -expected, 0.0]);
    // The shift doesn't change the susceptances
    assert_eq!(
        model.bbus().to_dense(),
        self::model(&triangle()).unwrap().bbus().to_dense()
    );
}

#[test]
fn out_of_service_branches_are_dropped() {
    let mut tables = triangle();
    tables.2[1].5 = false;
    let model = model(&tables).unwrap();
    assert_eq!(model.branch_rows(), [0, 2]);
    assert_eq!(
        model.bbus().to_dense(),
        vec![
            vec![14.0, -10.0, -4.0],
            vec![-10.0, 10.0, 0.0],
            vec![-4.0, 0.0, 4.0]
        ]
    );
}

#[test]
fn isolated_buses_and_their_elements_are_dropped() {
    // Bus 9 is isolated, with an in-service generator and an in-service branch to bus 2. Both
    //   count as out of service, as in MATPOWER's ext2int.
    let mut tables = triangle();
    tables.0.insert(1, B(9, 4, 20.0, 0.0, 7.5));
    tables.1.push(G(9, 20.0, true));
    tables.2.push(L(9, 2, 0.1, 0.0, 0.0, true));
    let model = model(&tables).unwrap();
    assert_eq!(model.bus_rows(), [0, 2, 3]);
    assert_eq!(model.gen_rows(), [0]);
    assert_eq!(model.branch_rows(), [0, 1, 2]);
    assert_eq!(
        model.bbus().to_dense(),
        self::model(&triangle()).unwrap().bbus().to_dense()
    );
}

#[test]
fn bus_labels_need_not_be_consecutive() {
    let relabel = |label: usize| [0, 40, 7, 1000][label];
    let (buses, gens, branches) = triangle();
    let relabeled = (
        buses
            .iter()
            .map(|b| B(relabel(b.0), b.1, b.2, b.3, b.4))
            .collect(),
        gens.iter().map(|g| G(relabel(g.0), g.1, g.2)).collect(),
        branches
            .iter()
            .map(|l| L(relabel(l.0), relabel(l.1), l.2, l.3, l.4, l.5))
            .collect(),
    );
    let model = model(&relabeled).unwrap();
    assert_eq!(
        model.bbus().to_dense(),
        self::model(&triangle()).unwrap().bbus().to_dense()
    );
    assert_eq!(model.ref_buses(), [0]);
}

#[test]
fn bus_roles_follow_bustypes() {
    // A PV bus with a generator is PV; one whose only generator is off is PQ
    let mut tables = triangle();
    tables.0[1].1 = 2;
    tables.0[2].1 = 2;
    tables.1.push(G(2, 10.0, true));
    tables.1.push(G(3, 10.0, false));
    let roles = model(&tables).unwrap();
    assert_eq!(roles.ref_buses(), [0]);
    assert_eq!(roles.pv_buses(), [1]);
    assert_eq!(roles.pq_buses(), [2]);

    // A REF bus without an in-service generator is PQ, and the first PV bus becomes the
    //   reference
    tables.1[0].2 = false;
    let promoted = model(&tables).unwrap();
    assert_eq!(promoted.ref_buses(), [1]);
    assert!(promoted.pv_buses().is_empty());
    assert_eq!(promoted.pq_buses(), [0, 2]);

    // Several REF buses with generators are all references
    let mut tables = triangle();
    tables.0[2].1 = 3;
    tables.1.push(G(3, 10.0, true));
    assert_eq!(model(&tables).unwrap().ref_buses(), [0, 2]);
}

#[test]
fn no_reference_bus_is_an_error() {
    let mut tables = triangle();
    tables.1[0].2 = false;
    assert_eq!(model(&tables).unwrap_err(), PowerflowError::NoReferenceBus);
}

#[test]
fn zero_reactance_is_an_error_only_in_service() {
    let mut tables = triangle();
    tables.2[1] = L(2, 3, 0.0, 0.0, 0.0, true);
    assert_eq!(
        model(&tables).unwrap_err(),
        PowerflowError::ZeroReactance { row: 2 }
    );
    tables.2[1].5 = false;
    assert!(model(&tables).is_ok());
}

#[test]
fn island_without_reference_is_an_error() {
    // Buses 4 and 5 are connected to each other but not to the reference
    let mut tables = triangle();
    tables.0.push(B(4, 1, 10.0, 0.0, 0.0));
    tables.0.push(B(5, 1, 10.0, 0.0, 0.0));
    tables.2.push(L(4, 5, 0.1, 0.0, 0.0, true));
    assert_eq!(
        model(&tables).unwrap_err(),
        PowerflowError::IslandWithoutReference { bus_label: 4 }
    );

    // With its own reference bus and generator, the island is fine
    tables.0[4].1 = 3;
    tables.1.push(G(5, 20.0, true));
    assert_eq!(model(&tables).unwrap().ref_buses(), [0, 4]);
}

#[test]
fn an_edited_case_is_validated() {
    // Only parsing validates, so these edits reach the model unchecked. Without the model's
    //   own check, an unknown label silently dropped its generator (case14 then solved with
    //   one generator's 40 MW missing), and a duplicate one gave a misleading island error.
    let tables = triangle();
    let valid = fixture(&tables.0, &tables.1, &tables.2);

    let mut case = valid.clone();
    case.generator[0].bus_label = 99;
    assert_eq!(
        DcModel::new(&case).unwrap_err(),
        PowerflowError::UnknownBus {
            table: "gen",
            row: 1,
            label: 99
        }
    );

    let mut case = valid.clone();
    case.bus[2].bus_label = 1;
    assert_eq!(
        DcModel::new(&case).unwrap_err(),
        PowerflowError::DuplicateBus {
            label: 1,
            first_row: 1,
            second_row: 3
        }
    );

    let mut case = valid;
    case.base_mva = 0.0;
    assert_eq!(
        DcModel::new(&case).unwrap_err(),
        PowerflowError::InvalidBaseMva { base_mva: 0.0 }
    );
}
