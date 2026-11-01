// Tests for the AC model (Ybus, Yf, Yt and Sbus): the branch π-model, shunts, which elements are
//   dropped, re-stamping and the errors. Expected values are computed independently here from
//   the π-model's formulas.

mod common;
use common::ac_fixture::{ac_fixture, AcBranch as L, AcBus as B, AcGen as G};
use powerflow::ac::AcModel;
use powerflow::case::Case;
use powerflow::powerflow_error::PowerflowError;
use sparse_mat::Complex;

use std::f64::consts::PI;

fn c(re: f64, im: f64) -> Complex<f64> {
    Complex::new(re, im)
}

fn close(a: Complex<f64>, b: Complex<f64>) -> bool {
    (a - b).norm() <= 1e-12 * b.norm().max(1.0)
}

/// Two buses joined by a transformer with resistance, charging, an off-nominal tap and a phase
/// shift, and shunts at both buses.
fn transformer_case() -> Case {
    ac_fixture(
        &[
            B(1, 3, 0.0, 0.0, 10.0, -20.0, 1.0, 0.0),
            B(2, 1, 50.0, 20.0, 0.0, 15.0, 1.0, 0.0),
        ],
        &[G(1, 50.0, 10.0, 100.0, -100.0, 1.02, true)],
        &[L(1, 2, 0.01, 0.1, 0.02, 0.95, 10.0, true)],
    )
}

#[test]
fn transformer_stamps_match_the_pi_model() {
    let model = AcModel::new(&transformer_case()).unwrap();
    // The π-model by hand
    let ys = c(1.0, 0.0) / c(0.01, 0.1);
    let tau = Complex::from_polar(0.95, 10.0 * PI / 180.0);
    let ytt = ys + c(0.0, 0.01);
    let yff = ytt / (0.95 * 0.95);
    let yft = -ys / tau.conj();
    let ytf = -ys / tau;
    let (shunt1, shunt2) = (c(0.1, -0.2), c(0.0, 0.15));

    let ybus = model.ybus().to_dense();
    for (got, want) in [
        (ybus[0][0], yff + shunt1),
        (ybus[0][1], yft),
        (ybus[1][0], ytf),
        (ybus[1][1], ytt + shunt2),
    ] {
        assert!(close(got, want), "{got} vs {want}");
    }
    let (yf, yt) = (model.yf().to_dense(), model.yt().to_dense());
    assert!(close(yf[0][0], yff) && close(yf[0][1], yft));
    assert!(close(yt[0][0], ytf) && close(yt[0][1], ytt));
    // A phase shift makes Ybus unsymmetric
    assert!(!model.ybus().is_symmetric());
    assert_eq!(model.base_mva(), 100.0);
}

/// A triangle of lines with resistance, and a fourth bus on a transformer from bus 3.
fn network(charging: f64, tap: f64, shift: f64) -> Case {
    ac_fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            B(2, 2, 20.0, 5.0, 0.0, 0.0, 1.0, 0.0),
            B(3, 1, 40.0, 15.0, 3.0, 12.0, 1.0, 0.0),
            B(4, 1, 30.0, 10.0, 0.0, -5.0, 1.0, 0.0),
        ],
        &[
            G(1, 60.0, 0.0, 100.0, -100.0, 1.0, true),
            G(2, 30.0, 0.0, 50.0, -50.0, 1.01, true),
        ],
        &[
            L(1, 2, 0.02, 0.06, charging, 0.0, 0.0, true),
            L(2, 3, 0.03, 0.09, charging, 0.0, 0.0, true),
            L(1, 3, 0.01, 0.05, charging, 0.0, 0.0, true),
            L(3, 4, 0.0, 0.08, 0.0, tap, shift, true),
        ],
    )
}

#[test]
fn rows_sum_to_the_shunts_without_charging_or_taps() {
    // With no charging and unit taps, each branch's stamps cancel across a row, leaving the
    //   bus shunt
    let model = AcModel::new(&network(0.0, 0.0, 0.0)).unwrap();
    let shunts = [c(0.0, 0.0), c(0.0, 0.0), c(0.03, 0.12), c(0.0, -0.05)];
    for (row, shunt) in model.ybus().to_dense().iter().zip(shunts) {
        let sum: Complex<f64> = row.iter().sum();
        assert!(close(sum, shunt), "{sum} vs {shunt}");
    }
}

#[test]
fn ybus_is_symmetric_unless_phase_shifted() {
    // Off-nominal taps and charging keep Ybus symmetric (Yft = Ytf for a real tap); a phase
    //   shift doesn't. It's never Hermitian, since series admittances aren't real.
    let tapped = AcModel::new(&network(0.04, 1.05, 0.0)).unwrap();
    assert!(tapped.ybus().is_symmetric());
    assert!(!tapped.ybus().is_hermitian());
    let shifted = AcModel::new(&network(0.04, 1.05, -3.0)).unwrap();
    assert!(!shifted.ybus().is_symmetric());
    // Every stamped position is stored, including the shunt diagonal of a bus without a shunt
    assert_eq!(tapped.ybus().num_stored(), 4 + 2 * 4);
}

#[test]
fn sbus_is_generation_less_load() {
    let case = network(0.0, 0.0, 0.0);
    let model = AcModel::new(&case).unwrap();
    let sbus = model.sbus_pu(&case);
    let expected = [c(0.6, 0.0), c(0.1, -0.05), c(-0.4, -0.15), c(-0.3, -0.1)];
    for (got, want) in sbus.iter().zip(expected) {
        assert!(close(*got, want), "{got} vs {want}");
    }
}

#[test]
fn restamping_equals_a_fresh_build() {
    let case = network(0.04, 1.05, -3.0);
    let mut model = AcModel::new(&case).unwrap();
    let mut changed = case.clone();
    changed.branch[0].r = 0.025;
    changed.branch[1].x = 0.12;
    changed.branch[2].b = 0.1;
    changed.branch[3].tap_ratio = 0.98;
    changed.branch[3].shift_deg = 5.0;
    changed.bus[2].bs_mvar = -7.0;
    model.restamp(&changed).unwrap();
    let fresh = AcModel::new(&changed).unwrap();
    for (got, want) in [
        (model.ybus(), fresh.ybus()),
        (model.yf(), fresh.yf()),
        (model.yt(), fresh.yt()),
    ] {
        assert_eq!(got.pointers(), want.pointers());
        assert_eq!(got.indices(), want.indices());
        assert_eq!(got.values(), want.values(), "bit for bit");
    }

    // Errors leave the model as it was
    let before = model.ybus().clone();
    let mut zero = changed.clone();
    zero.branch[1].r = 0.0;
    zero.branch[1].x = 0.0;
    assert_eq!(
        model.restamp(&zero),
        Err(PowerflowError::ZeroImpedance { row: 2 })
    );
    let mut shorter = changed.clone();
    shorter.branch.pop();
    assert!(matches!(
        model.restamp(&shorter),
        Err(PowerflowError::CaseMismatch { .. })
    ));
    assert_eq!(model.ybus().values(), before.values());
}

#[test]
fn zero_impedance_is_refused_only_in_service() {
    let mut case = network(0.0, 0.0, 0.0);
    case.branch[2].r = 0.0;
    case.branch[2].x = 0.0;
    assert_eq!(
        AcModel::new(&case).err(),
        Some(PowerflowError::ZeroImpedance { row: 3 })
    );
    // Out of service, it's dropped before the check
    case.branch[2].status = false;
    let model = AcModel::new(&case).unwrap();
    assert_eq!(model.branch_rows(), [0, 1, 3]);
    // Zero resistance alone is a lossless line, which is fine
    assert!(AcModel::new(&network(0.0, 0.0, 0.0)).is_ok());
}

#[test]
fn isolated_buses_and_out_of_service_elements_are_dropped() {
    let mut case = network(0.0, 0.0, 0.0);
    case.bus[3].bus_type = powerflow::case::BusType::Isolated;
    case.generator[1].status = false;
    let model = AcModel::new(&case).unwrap();
    assert_eq!(model.bus_rows(), [0, 1, 2]);
    assert_eq!(model.branch_rows(), [0, 1, 2]);
    assert_eq!(model.gen_rows(), [0]);
    // Bus 2 has lost its generator, so it's PQ
    assert_eq!(model.ref_buses(), [0]);
    assert!(model.pv_buses().is_empty());
    assert_eq!(model.pq_buses(), [1, 2]);
    assert_eq!(model.ybus().shape(), (3, 3));
}
