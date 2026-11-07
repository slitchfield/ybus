// Tests the predictor-corrector method (interior_point::mehrotra) on the test problems of
//   MATPOWER's MIPS tests that have no nonlinear constraints (tests/common): it must converge
//   wherever MATPOWER's MIPS does, with the objective within 1e-6 (relative) of MIPS's, the
//   solution and multipliers close to MIPS's and to each problem's known answer, and small KKT
//   residuals; and fail where MIPS fails. Its iteration counts are printed beside MIPS's.

mod common;

use common::{
    dense, known, qp4_constraints, qp4_hessian, references, Banana, Nonlinear2, Quadratic3,
};
use interior_point::interior_point_error::InteriorPointError;
use interior_point::mehrotra::{mehrotra, qps_mehrotra, MehrotraOptions, MehrotraSolution};
use interior_point::problem::{Constraints, Problem, QuadraticProblem};
use interior_point::solution::kkt_residuals;

/// The bar for the objective, relative to the larger magnitude (at least 1)
const OBJECTIVE: f64 = 1e-6;

/// How far the solution and the multipliers may be from MIPS's, relative to the larger magnitude
/// (at least 1). Both stop once their measures are below 1e-6, from different paths, so they
/// agree only to about that.
const AGREEMENT: f64 = 1e-4;

/// The largest KKT residual accepted, unscaled
const RESIDUAL: f64 = 1e-5;

/// Each problem's iteration count, as this method takes it (2026-11-07). The method is
/// deterministic, so a change to it that changes a count shows here, as MIPS's tests pin
/// MATPOWER's counts; a change that's meant to must update them, with the reason.
fn expected_iterations(name: &str) -> usize {
    match name {
        "banana" => 6,
        "quadratic3" | "qps_quadratic3" => 2,
        "qp4" | "qps_qp4" => 13,
        "lp3" => 10,
        "qp2" => 8,
        _ => panic!("no expected iteration count for {name}"),
    }
}

fn close_within(ours: &[f64], theirs: &[f64], tolerance: f64, what: &str) -> f64 {
    assert_eq!(ours.len(), theirs.len(), "{what}: lengths");
    let mut worst = 0.0_f64;
    for (k, (a, b)) in ours.iter().zip(theirs).enumerate() {
        let scale = a.abs().max(b.abs()).max(1.0);
        let difference = (a - b).abs() / scale;
        worst = worst.max(difference);
        assert!(
            difference <= tolerance,
            "{what}[{k}]: ours {a:e}, MIPS's {b:e}"
        );
    }
    worst
}

/// Checks a converged solution against MATPOWER's MIPS run of the same problem (without step
/// control), and its KKT residuals; prints the iterations against MIPS's
fn check<P: Problem>(
    name: &str,
    problem: &P,
    constraints: &Constraints,
    solution: &MehrotraSolution,
) {
    let runs = references();
    let reference = &runs[&(name.to_string(), false)];
    let field = |key: &str| reference.0.get(key);
    assert_eq!(field("eflag").unwrap()[0], 1.0, "{name}: MIPS converges");
    assert!(solution.converged(), "{name}: {:?}", solution.status);
    assert_eq!(
        solution.iterations,
        expected_iterations(name),
        "{name}: iterations"
    );
    close_within(
        &[solution.f],
        field("f").unwrap(),
        OBJECTIVE,
        &format!("{name}: f"),
    );
    let mut worst = close_within(
        &solution.x,
        field("x").unwrap(),
        AGREEMENT,
        &format!("{name}: x"),
    );
    let m = &solution.multipliers;
    for (key, ours) in [
        ("mu_l", &m.mu_l),
        ("mu_u", &m.mu_u),
        ("lower", &m.lower),
        ("upper", &m.upper),
    ] {
        let theirs = field(key).unwrap();
        worst = worst.max(close_within(
            ours,
            theirs,
            AGREEMENT,
            &format!("{name}: {key}"),
        ));
    }
    let residuals = kkt_residuals(problem, constraints, &solution.x, m).unwrap();
    assert!(residuals.max() < RESIDUAL, "{name}: {residuals:?}");
    eprintln!(
        "{name}: {} iterations (MIPS {}), largest difference from MIPS {worst:.1e}, residual {:.1e}",
        solution.iterations,
        field("iterations").unwrap()[0],
        residuals.max()
    );
}

#[test]
fn unconstrained_problems_match_mips() {
    let options = MehrotraOptions::default();
    let none = Constraints::default();
    let s = mehrotra(&Banana, &[-1.9, 2.0], &none, &options).unwrap();
    check("banana", &Banana, &none, &s);
    known(&s.x, &[1.0, 1.0], 13, "banana x");

    let s = mehrotra(&Quadratic3, &[0.0; 3], &none, &options).unwrap();
    check("quadratic3", &Quadratic3, &none, &s);
    known(&s.x, &[3.0, 5.0, 7.0], 13, "quadratic3 x");
    known(&[s.f], &[-244.0], 13, "quadratic3 f");
}

#[test]
fn a_linearly_constrained_qp_matches_mips() {
    let qp = QuadraticProblem {
        h: Some(qp4_hessian()),
        c: vec![0.0; 4],
    };
    let s = mehrotra(
        &qp,
        &[1.0, 0.0, 0.0, 1.0],
        &qp4_constraints(),
        &MehrotraOptions::default(),
    )
    .unwrap();
    check("qp4", &qp, &qp4_constraints(), &s);
    let third = |v: &[f64]| v.iter().map(|v| v / 3.0).collect::<Vec<_>>();
    known(&s.x, &third(&[0.0, 2.8, 0.2, 0.0]), 6, "qp4 x");
    known(&[s.f], &[3.29 / 3.0], 6, "qp4 f");
    known(&s.multipliers.mu_l, &third(&[6.58, 0.0]), 6, "qp4 mu_l");
    known(
        &s.multipliers.lower,
        &[2.24, 0.0, 0.0, 1.7667],
        4,
        "qp4 lower",
    );
}

#[test]
fn quadratic_and_linear_programs_match_mips() {
    let mut opt = MehrotraOptions::default();
    opt.comptol = 1e-8;

    let lp = QuadraticProblem {
        h: None,
        c: vec![-5.0, -4.0, -6.0],
    };
    let constraints = Constraints {
        a: Some(dense(&[
            &[1.0, -1.0, 1.0],
            &[-3.0, -2.0, -4.0],
            &[3.0, 2.0, 0.0],
        ])),
        l: Some(vec![f64::NEG_INFINITY, -42.0, f64::NEG_INFINITY]),
        u: Some(vec![20.0, f64::INFINITY, 30.0]),
        xmin: Some(vec![0.0; 3]),
        xmax: None,
    };
    let s = qps_mehrotra(&lp, &constraints, None, &opt).unwrap();
    check("lp3", &lp, &constraints, &s);
    known(&s.x, &[0.0, 15.0, 3.0], 6, "lp3 x");
    known(&s.multipliers.mu_l, &[0.0, 1.5, 0.0], 9, "lp3 mu_l");
    known(&s.multipliers.mu_u, &[0.0, 0.0, 0.5], 9, "lp3 mu_u");

    let qp = QuadraticProblem {
        h: Some(Quadratic3.objective_hessian(&[])),
        c: vec![2.0, -35.0, -47.0],
    };
    let s = qps_mehrotra(&qp, &Constraints::default(), Some(&[0.0; 3]), &opt).unwrap();
    check("qps_quadratic3", &qp, &Constraints::default(), &s);
    known(&s.x, &[3.0, 5.0, 7.0], 8, "qps_quadratic3 x");
    known(&[s.f], &[-249.0], 13, "qps_quadratic3 f");

    let qp = QuadraticProblem {
        h: Some(dense(&[&[1.0, -1.0], &[-1.0, 2.0]])),
        c: vec![-2.0, -6.0],
    };
    let constraints = Constraints {
        a: Some(dense(&[&[1.0, 1.0], &[-1.0, 2.0], &[2.0, 1.0]])),
        l: None,
        u: Some(vec![2.0, 2.0, 3.0]),
        xmin: Some(vec![0.0; 2]),
        xmax: None,
    };
    let s = qps_mehrotra(&qp, &constraints, None, &opt).unwrap();
    check("qp2", &qp, &constraints, &s);
    known(&s.x, &[2.0 / 3.0, 4.0 / 3.0], 7, "qp2 x");
    known(&[s.f], &[-74.0 / 9.0], 6, "qp2 f");
    known(
        &s.multipliers.mu_u,
        &[28.0 / 9.0, 4.0 / 9.0, 0.0],
        4,
        "qp2 mu_u",
    );

    let qp = QuadraticProblem {
        h: Some(qp4_hessian()),
        c: vec![0.0; 4],
    };
    let s = qps_mehrotra(&qp, &qp4_constraints(), Some(&[1.0, 0.0, 0.0, 1.0]), &opt).unwrap();
    check("qps_qp4", &qp, &qp4_constraints(), &s);
    known(&[s.f], &[3.29 / 3.0], 6, "qps_qp4 f");
}

#[test]
fn an_infeasible_lp_does_not_converge() {
    // x₁ + x₂ ≤ −1 with x ≥ 0, which MIPS fails on numerically
    let lp = QuadraticProblem {
        h: None,
        c: vec![1.0, 1.0],
    };
    let constraints = Constraints {
        a: Some(dense(&[&[1.0, 1.0]])),
        l: None,
        u: Some(vec![-1.0]),
        xmin: Some(vec![0.0; 2]),
        xmax: None,
    };
    let s = qps_mehrotra(&lp, &constraints, None, &MehrotraOptions::default()).unwrap();
    assert!(
        !s.converged(),
        "{:?} after {} iterations",
        s.status,
        s.iterations
    );
    eprintln!(
        "infeasible_lp: {:?} after {} iterations",
        s.status, s.iterations
    );
}

#[test]
fn fixed_variables_and_two_sided_bounds_are_handled() {
    // min (x₁ − 3)² + (x₂ + 1)² + (x₃ − 0.5)² with x₁ fixed at 1, −0.5 ≤ x₂ ≤ 2 and
    //   0 ≤ x₃ ≤ 0.25: the optimum is (1, −0.5, 0.25). x₁'s fixed value holds it below its
    //   unconstrained optimum (∂f/∂x₁ = −4, so its upper multiplier is 4), x₂'s lower bound
    //   binds (∂f/∂x₂ = 1, its lower multiplier 1), and x₃'s upper one (∂f/∂x₃ = −0.5, its
    //   upper multiplier 0.5)
    let qp = QuadraticProblem {
        h: Some(dense(&[
            &[2.0, 0.0, 0.0],
            &[0.0, 2.0, 0.0],
            &[0.0, 0.0, 2.0],
        ])),
        c: vec![-6.0, 2.0, -1.0],
    };
    let constraints = Constraints {
        xmin: Some(vec![1.0, -0.5, 0.0]),
        xmax: Some(vec![1.0, 2.0, 0.25]),
        ..Constraints::default()
    };
    let s = qps_mehrotra(
        &qp,
        &constraints,
        Some(&[5.0, 5.0, 5.0]),
        &MehrotraOptions::default(),
    )
    .unwrap();
    assert!(s.converged(), "{:?}", s.status);
    known(&s.x, &[1.0, -0.5, 0.25], 6, "x");
    known(&s.multipliers.upper, &[4.0, 0.0, 0.5], 5, "upper");
    known(&s.multipliers.lower, &[0.0, 1.0, 0.0], 5, "lower");
    let residuals = kkt_residuals(&qp, &constraints, &s.x, &s.multipliers).unwrap();
    assert!(residuals.max() < RESIDUAL, "{residuals:?}");
}

#[test]
fn nonlinear_constraints_are_not_supported_yet() {
    let bounds = Constraints {
        xmin: Some(vec![0.0; 2]),
        ..Constraints::default()
    };
    assert_eq!(
        mehrotra(
            &Nonlinear2,
            &[1.1, 0.0],
            &bounds,
            &MehrotraOptions::default()
        )
        .unwrap_err(),
        InteriorPointError::Unsupported {
            what: "nonlinear constraints"
        }
    );
}
