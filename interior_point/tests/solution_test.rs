// Tests for the solver-independent KKT check (interior_point::solution::kkt_residuals): zero at
//   a known exact optimum, each residual catching its own kind of error, and small at MIPS's
//   solutions of MATPOWER's MIPS test problems, nonlinear ones included.

mod common;

use common::{dense, qp4_constraints, qp4_hessian, Hs71, Nonlinear3};
use interior_point::mips::{mips, qps_mips, MipsOptions};
use interior_point::problem::{Constraints, QuadraticProblem};
use interior_point::solution::{kkt_residuals, Multipliers};

/// t_qps_mips's qp2: min ½ xᵀ [[1, −1], [−1, 2]] x − 2x₁ − 6x₂ subject to three rows ≤ (2, 2, 3)
/// and x ≥ 0, whose optimum is x = (2/3, 4/3) with the first two rows binding
fn qp2() -> (QuadraticProblem, Constraints) {
    (
        QuadraticProblem {
            h: Some(dense(&[&[1.0, -1.0], &[-1.0, 2.0]])),
            c: vec![-2.0, -6.0],
        },
        Constraints {
            a: Some(dense(&[&[1.0, 1.0], &[-1.0, 2.0], &[2.0, 1.0]])),
            l: None,
            u: Some(vec![2.0, 2.0, 3.0]),
            xmin: Some(vec![0.0; 2]),
            xmax: None,
        },
    )
}

/// qp2's exact optimum and multipliers: ∇f = H x + c = (−8/3, −10/3) = −Aᵀ μ_u with
/// μ_u = (28/9, 4/9, 0)
fn qp2_optimum() -> (Vec<f64>, Multipliers) {
    (
        vec![2.0 / 3.0, 4.0 / 3.0],
        Multipliers {
            mu_l: vec![0.0; 3],
            mu_u: vec![28.0 / 9.0, 4.0 / 9.0, 0.0],
            lower: vec![0.0; 2],
            upper: vec![0.0; 2],
            ..Multipliers::default()
        },
    )
}

#[test]
fn an_exact_optimum_has_no_residual() {
    let (qp, constraints) = qp2();
    let (x, multipliers) = qp2_optimum();
    let r = kkt_residuals(&qp, &constraints, &x, &multipliers).unwrap();
    assert!(r.max() < 1e-14, "{r:?}");
}

#[test]
fn each_residual_catches_its_own_error() {
    let (qp, constraints) = qp2();
    let (x, multipliers) = qp2_optimum();

    // A wrong multiplier breaks stationarity, by its row's coefficients
    let mut m = multipliers.clone();
    m.mu_u[0] += 0.5;
    let r = kkt_residuals(&qp, &constraints, &x, &m).unwrap();
    assert!((r.stationarity - 0.5).abs() < 1e-14, "{r:?}");

    // A point past the second row's limit is infeasible by that much
    let r = kkt_residuals(
        &qp,
        &constraints,
        &[2.0 / 3.0, 4.0 / 3.0 + 0.1],
        &multipliers,
    )
    .unwrap();
    assert!((r.feasibility - 0.2).abs() < 1e-14, "{r:?}");

    // A multiplier on a row that isn't binding breaks complementarity: the third row is 1/3 below
    //   its limit
    let mut m = multipliers.clone();
    m.mu_u[2] = 0.3;
    let r = kkt_residuals(&qp, &constraints, &x, &m).unwrap();
    assert!((r.complementarity - 0.1).abs() < 1e-14, "{r:?}");

    // A negative multiplier is dual infeasible
    let mut m = multipliers.clone();
    m.lower[1] = -0.25;
    let r = kkt_residuals(&qp, &constraints, &x, &m).unwrap();
    assert!((r.dual_feasibility - 0.25).abs() < 1e-14, "{r:?}");
}

#[test]
fn mips_solutions_have_small_residuals() {
    let options = MipsOptions::default();

    let qp = QuadraticProblem {
        h: Some(qp4_hessian()),
        c: vec![0.0; 4],
    };
    let s = qps_mips(
        &qp,
        &qp4_constraints(),
        Some(&[1.0, 0.0, 0.0, 1.0]),
        &options,
    )
    .unwrap();
    let r = kkt_residuals(&qp, &qp4_constraints(), &s.x, &s.multipliers).unwrap();
    assert!(r.max() < 1e-5, "qp4: {r:?}");

    let s = mips(
        &Nonlinear3,
        &Constraints::default(),
        &[1.0, 1.0, 0.0],
        &options,
    )
    .unwrap();
    let r = kkt_residuals(&Nonlinear3, &Constraints::default(), &s.x, &s.multipliers).unwrap();
    assert!(r.max() < 1e-5, "nonlinear3: {r:?}");

    let bounds = Constraints {
        xmin: Some(vec![1.0; 4]),
        xmax: Some(vec![5.0; 4]),
        ..Constraints::default()
    };
    let s = mips(&Hs71, &bounds, &[1.0, 5.0, 5.0, 1.0], &options).unwrap();
    let r = kkt_residuals(&Hs71, &bounds, &s.x, &s.multipliers).unwrap();
    assert!(r.max() < 1e-5, "hs71: {r:?}");
}

#[test]
fn mismatched_multipliers_are_errors() {
    let (qp, constraints) = qp2();
    let (x, mut multipliers) = qp2_optimum();
    multipliers.mu_u.pop();
    assert!(kkt_residuals(&qp, &constraints, &x, &multipliers).is_err());
    let (_, multipliers) = qp2_optimum();
    assert!(kkt_residuals(&qp, &constraints, &[0.0], &multipliers).is_err());
}
