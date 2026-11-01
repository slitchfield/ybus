// Tests for the AC Newton Jacobian, built through an assembly map: it must equal, bit for bit,
//   the matrix-operation build it replaced, kept here as `reference_jacobian`.

mod common;
use common::ac_fixture::{ac_fixture, AcBranch as L, AcBus as B, AcGen as G};
use powerflow::ac::{AcModel, AcOptions};
use powerflow::case::Case;
use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix};
use sparse_mat::Complex;

use std::f64::consts::PI;

/// The Jacobian as Phase 4.5 built it: `dS/dVa` and `dS/dVm` from Ybus with matrix operations
/// (MATPOWER's `dSbus_dV`), then the four blocks selected over the PV and PQ buses and stacked.
fn reference_jacobian(model: &AcModel, v: &[Complex<f64>]) -> CscMatrix<f64> {
    let ybus = model.ybus();
    let current = ybus.mul_vec(v).unwrap();
    let v_unit: Vec<Complex<f64>> = v.iter().map(|value| value / value.norm()).collect();
    let j_v: Vec<Complex<f64>> = v.iter().map(|value| value * Complex::i()).collect();
    let ybus_v = ybus.scale_cols(v).unwrap();
    let ds_dva = (&CsrMatrix::from_diagonal(current.iter().copied()) - &ybus_v)
        .conj()
        .scale_rows(&j_v)
        .unwrap();
    let ds_dvm = &ybus
        .scale_cols(&v_unit)
        .unwrap()
        .conj()
        .scale_rows(v)
        .unwrap()
        + &CsrMatrix::from_diagonal(current.iter().zip(&v_unit).map(|(i, u)| i.conj() * u));
    let pvpq: Vec<usize> = model
        .pv_buses()
        .iter()
        .chain(model.pq_buses())
        .copied()
        .collect();
    let pq = model.pq_buses();
    let block = |m: &CsrMatrix<f64>, rows: &[usize], cols: &[usize]| m.select(rows, cols).unwrap();
    let (re_va, im_va) = (ds_dva.real(), ds_dva.imag());
    let (re_vm, im_vm) = (ds_dvm.real(), ds_dvm.imag());
    let top = CsrMatrix::hstack([block(&re_va, &pvpq, &pvpq), block(&re_vm, &pvpq, pq)]).unwrap();
    let bottom = CsrMatrix::hstack([block(&im_va, pq, &pvpq), block(&im_vm, pq, pq)]).unwrap();
    CsrMatrix::vstack([top, bottom]).unwrap().to_csc()
}

fn assert_same_bits(got: &CscMatrix<f64>, want: &CscMatrix<f64>, what: &str) {
    assert_eq!(got.shape(), want.shape(), "{what}");
    assert_eq!(got.pointers(), want.pointers(), "{what}");
    assert_eq!(got.indices(), want.indices(), "{what}");
    let bits = |m: &CscMatrix<f64>| m.values().iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(got), bits(want), "{what}");
}

/// Voltage states to compare at: the case's start, a flat start, the solution, and the solution
/// with every angle and magnitude nudged differently.
fn states(case: &Case, model: &AcModel) -> Vec<Vec<Complex<f64>>> {
    let rows = model.bus_rows();
    let polar = |vm: &dyn Fn(usize, usize) -> f64, va: &dyn Fn(usize, usize) -> f64| {
        rows.iter()
            .enumerate()
            .map(|(k, &row)| Complex::from_polar(vm(k, row), va(k, row) * PI / 180.0))
            .collect::<Vec<_>>()
    };
    let mut states = vec![
        polar(&|_, row| case.bus[row].vm_pu, &|_, row| {
            case.bus[row].va_deg
        }),
        polar(&|_, _| 1.0, &|_, _| 0.0),
    ];
    if let Ok(solution) = model.solve(case, AcOptions::default()) {
        states.push(polar(&|_, row| solution.vm_pu[row], &|_, row| {
            solution.va_deg[row]
        }));
        states.push(polar(
            &|k, row| solution.vm_pu[row] * (1.0 + 0.003 * ((k % 7) as f64 - 3.0)),
            &|k, row| solution.va_deg[row] + 0.7 * ((k % 5) as f64 - 2.0),
        ));
    }
    states
}

fn check(case: &Case, name: &str) {
    let model = AcModel::new(case).unwrap();
    for (k, v) in states(case, &model).iter().enumerate() {
        assert_same_bits(
            &model.jacobian(v),
            &reference_jacobian(&model, v),
            &format!("{name}, state {k}"),
        );
    }
}

/// A four-bus network with a PV bus, PQ buses, a shunt, line charging and a tapped transformer.
fn network(shift: f64) -> Case {
    ac_fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            B(2, 2, 20.0, 5.0, 0.0, 0.0, 1.0, 0.0),
            B(3, 1, 45.0, 15.0, 0.0, 10.0, 1.0, 0.0),
            B(4, 1, 40.0, 12.0, 0.0, 0.0, 1.0, 0.0),
        ],
        &[
            G(1, 0.0, 0.0, 100.0, -100.0, 1.04, true),
            G(2, 40.0, 0.0, 60.0, -40.0, 1.02, true),
        ],
        &[
            L(1, 2, 0.02, 0.06, 0.03, 0.0, 0.0, true),
            L(2, 3, 0.03, 0.09, 0.02, 0.0, 0.0, true),
            L(1, 3, 0.01, 0.05, 0.02, 0.0, 0.0, true),
            L(3, 4, 0.005, 0.08, 0.0, 0.97, shift, true),
        ],
    )
}

#[test]
fn jacobian_equals_the_matrix_operation_build_on_fixtures() {
    check(&network(0.0), "fixture");
    check(&network(-4.0), "fixture with a phase shifter");
}

#[test]
fn jacobian_equals_the_matrix_operation_build_on_matpower_cases() {
    for name in ["case14", "case30", "case118", "case2383wp"] {
        let Some(path) = common::matpower_case(name) else {
            continue; // not fetched; the skip was reported on stderr
        };
        check(&Case::from_path(&path).unwrap(), name);
    }
}
