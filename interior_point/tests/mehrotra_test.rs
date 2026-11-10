// Tests the predictor-corrector method (interior_point::mehrotra) on the test problems of
//   MATPOWER's MIPS tests (tests/common), with and without nonlinear constraints: it must converge
//   wherever MATPOWER's MIPS does, with the objective within 1e-6 (relative) of MIPS's, the
//   solution and multipliers close to MIPS's and to each problem's known answer, and small KKT
//   residuals; and fail where MIPS fails. On a nonconvex problem where the two may stop at
//   different local optima, it must reach a verified KKT point at a known one. Its iteration
//   counts are printed beside MIPS's.

mod common;

use common::{
    dense, known, qp4_constraints, qp4_hessian, references, Banana, Hs71, Nonlinear2, Nonlinear3,
    Quadratic3,
};
use interior_point::mehrotra::{mehrotra, qps_mehrotra, MehrotraOptions};
use interior_point::problem::{
    Constraints, NonlinearConstraints, NonlinearMultipliers, Problem, QuadraticProblem,
};
use interior_point::solution::kkt_residuals;
use interior_point::solution::{Solution, Status};
use sparse_mat::compressed_matrix::CsrMatrix;

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
        // With nonlinear constraints the corrector's second-order term is damped (phase O3's
        //   PR 4), which changed nonlinear2 from 12 and hs71 from 11
        "nonlinear2" => 8,
        "nonlinear3" => 6,
        "hs71" => 9,
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
fn check<P: Problem>(name: &str, problem: &P, constraints: &Constraints, solution: &Solution) {
    let runs = references();
    let reference = &runs[&(name.to_string(), false)];
    let field = |key: &str| reference.0.get(key);
    assert_eq!(field("eflag").unwrap()[0], 1.0, "{name}: MIPS converges");
    assert!(solution.converged(), "{name}: {:?}", solution.status);
    // One analysis serves the whole solve where the KKT pattern holds. nonlinear3 builds its
    //   Jacobian with from_dense, which stores no zeros, and starts at x₃ = 0, so an entry
    //   appears once x₃ moves (nonlinear2's x₂ = 0 is a bound, which this method moves off
    //   first, so its pattern holds)
    let expected = if name == "nonlinear3" { 2 } else { 1 };
    assert_eq!(solution.analyses, expected, "{name}: analyses");
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
        ("eqnonlin", &m.eq_nonlin),
        ("ineqnonlin", &m.ineq_nonlin),
        ("mu_l", &m.mu_l),
        ("mu_u", &m.mu_u),
        ("lower", &m.lower),
        ("upper", &m.upper),
    ] {
        match field(key) {
            Some(theirs) => {
                let what = format!("{name}: {key}");
                worst = worst.max(close_within(ours, theirs, AGREEMENT, &what));
            }
            // MIPS leaves out the nonlinear multipliers when there are no such constraints
            None => assert!(ours.is_empty(), "{name}: {key} should be empty"),
        }
    }
    let residuals = kkt_residuals(problem, constraints, &solution.x, m).unwrap();
    assert!(residuals.max() < RESIDUAL, "{name}: {residuals:?}");
    eprintln!(
        "{name}: {} iterations (MIPS {}), largest difference from MIPS {worst:.1e}, residual {:.1e}",
        solution.iterations,
        field("iterations").unwrap()[0],
        residuals.max()
    );
    assert_eq!(
        solution.iterations,
        expected_iterations(name),
        "{name}: iterations"
    );
}

#[test]
fn unconstrained_problems_match_mips() {
    let options = MehrotraOptions::default();
    let none = Constraints::default();
    let s = mehrotra(&Banana, &none, &[-1.9, 2.0], &options).unwrap();
    check("banana", &Banana, &none, &s);
    known(&s.x, &[1.0, 1.0], 13, "banana x");

    let s = mehrotra(&Quadratic3, &none, &[0.0; 3], &options).unwrap();
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
        &qp4_constraints(),
        &[1.0, 0.0, 0.0, 1.0],
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

/// hs71's bounds, 1 ≤ x ≤ 5
fn hs71_bounds() -> Constraints {
    Constraints {
        xmin: Some(vec![1.0; 4]),
        xmax: Some(vec![5.0; 4]),
        ..Constraints::default()
    }
}

/// hs71's known KKT points:
///
/// - **the local point:** where MIPS without step control stops, from (4.5, 4.5, 4.5, 4.5);
/// - **the optimum:** the one Hock and Schittkowski give;
/// - **a constrained local maximum:** x ≈ (5, 1.388, 1.472, 3.303), with x₁ at its upper bound
///   and the product constraint slack. The Lagrangian's Hessian, reduced to the tangent space
///   of the active constraints (the sphere and x₁'s bound), has eigenvalues −18.2 and −11.9.
///   The method converges to it from (2, 2, 2, 2) since its corrector was damped (phase O3's
///   PR 4): with no inertia correction, an interior-point method can't tell a maximum from a
///   minimum (left to phase O4's LDLᵀ). MIPS reaches its iteration limit from there.
const HS71_OPTIMA: [(&str, f64); 3] = [
    ("the local point", 35.278_137_612_159_04),
    ("the optimum", 17.0140173),
    ("a constrained local maximum", 131.261_236_977_572_6),
];

/// Which of hs71's known KKT points `f` is, within the objective bar
fn hs71_optimum(f: f64) -> Option<&'static str> {
    HS71_OPTIMA
        .iter()
        .find(|(_, optimum)| (f - optimum).abs() <= OBJECTIVE * optimum.abs())
        .map(|(name, _)| *name)
}

#[test]
fn nonlinear_programs_match_mips() {
    let options = MehrotraOptions::default();
    let bounds = Constraints {
        xmin: Some(vec![0.0; 2]),
        ..Constraints::default()
    };
    let s = mehrotra(&Nonlinear2, &bounds, &[1.1, 0.0], &options).unwrap();
    check("nonlinear2", &Nonlinear2, &bounds, &s);
    known(&s.x, &[1.0, 1.0], 6, "nonlinear2 x");
    known(&s.multipliers.ineq_nonlin, &[0.0, 0.5], 6, "nonlinear2 mu");

    let none = Constraints::default();
    let s = mehrotra(&Nonlinear3, &none, &[1.0, 1.0, 0.0], &options).unwrap();
    check("nonlinear3", &Nonlinear3, &none, &s);
    known(
        &s.x,
        &[1.58113883, 2.23606798, 1.58113883],
        6,
        "nonlinear3 x",
    );
    known(&[s.f], &[-5.0 * 2f64.sqrt()], 6, "nonlinear3 f");
    known(
        &s.multipliers.ineq_nonlin,
        &[0.0, 2f64.sqrt() / 2.0],
        7,
        "nonlinear3 mu",
    );

    let s = mehrotra(&Hs71, &hs71_bounds(), &[1.0, 5.0, 5.0, 1.0], &options).unwrap();
    check("hs71", &Hs71, &hs71_bounds(), &s);
    known(&s.x, &[1.0, 4.7429994, 3.8211503, 1.3794082], 6, "hs71 x");
    known(&[s.f], &[17.0140173], 6, "hs71 f");
    known(&s.multipliers.eq_nonlin, &[0.1614686], 5, "hs71 lambda");
    known(&s.multipliers.ineq_nonlin, &[0.55229366], 5, "hs71 mu");
    known(
        &s.multipliers.lower,
        &[1.08787121024, 0.0, 0.0, 0.0],
        5,
        "hs71 lower",
    );
}

/// What this method does on hs71 from each other start (2026-11-07), pinned as the iteration
/// counts are: the status, the iterations, and which known optimum it reached
fn expected_hs71(label: &str) -> (Status, usize, Option<&'static str>) {
    match label {
        // MIPS without step control stops at the local point (11 iterations)
        "4.5" => (Status::Converged, 11, Some("the optimum")),
        // MIPS reaches its iteration limit (150). With the damped corrector this method
        //   converges to a constrained local maximum (see HS71_OPTIMA); with the full one, it
        //   reached the optimum in 13
        "2" => (Status::Converged, 8, Some("a constrained local maximum")),
        _ => unreachable!(),
    }
}

#[test]
fn hs71_from_other_starts_reaches_a_known_kkt_point() {
    let options = MehrotraOptions::default();
    for (x0, label) in [([4.5; 4], "4.5"), ([2.0; 4], "2")] {
        let s = mehrotra(&Hs71, &hs71_bounds(), &x0, &options).unwrap();
        let runs = references();
        let mips = &runs[&(format!("hs71_from_{label}"), false)];
        let optimum = if s.converged() {
            hs71_optimum(s.f)
        } else {
            None
        };
        eprintln!(
            "hs71 from {label}: {:?} after {} iterations, f = {}, {} (MIPS: eflag {}, {} \
             iterations, f = {})",
            s.status,
            s.iterations,
            s.f,
            optimum.unwrap_or("not a known optimum"),
            mips.0["eflag"][0],
            mips.0["iterations"][0],
            mips.0["f"][0],
        );
        if s.converged() {
            let residuals = kkt_residuals(&Hs71, &hs71_bounds(), &s.x, &s.multipliers).unwrap();
            assert!(
                residuals.max() < RESIDUAL,
                "hs71 from {label}: {residuals:?}"
            );
            assert!(
                optimum.is_some(),
                "hs71 from {label}: f = {} is no known KKT point",
                s.f
            );
        }
        // Where MIPS converges, so must this method
        if mips.0["eflag"][0] == 1.0 {
            assert!(s.converged(), "hs71 from {label}: {:?}", s.status);
        }
        assert_eq!(
            (s.status, s.iterations, optimum),
            expected_hs71(label),
            "hs71 from {label}"
        );
    }
}

#[test]
fn the_cost_scale_comes_off_the_multipliers() {
    // The AC OPF solves with cost_mult 1e-4: the same point, objective and multipliers
    let x0 = [1.0, 5.0, 5.0, 1.0];
    let plain = mehrotra(&Hs71, &hs71_bounds(), &x0, &MehrotraOptions::default()).unwrap();
    let mut options = MehrotraOptions::default();
    options.cost_mult = 1e-2;
    let scaled = mehrotra(&Hs71, &hs71_bounds(), &x0, &options).unwrap();
    assert!(scaled.converged(), "{:?}", scaled.status);
    let tolerance = 1e-6;
    close_within(&scaled.x, &plain.x, tolerance, "x");
    close_within(&[scaled.f], &[plain.f], tolerance, "f");
    let (a, b) = (&scaled.multipliers, &plain.multipliers);
    close_within(&a.eq_nonlin, &b.eq_nonlin, tolerance, "eq_nonlin");
    close_within(&a.ineq_nonlin, &b.ineq_nonlin, tolerance, "ineq_nonlin");
    close_within(&a.lower, &b.lower, tolerance, "lower");
    close_within(&a.upper, &b.upper, tolerance, "upper");
}

/// min x₂ − x₁ subject to x₁² − x₂ = 0 and 0 ≤ x₁ ≤ 2: on the parabola, x₁² − x₁ is least at
/// x₁ = 1/2, so the optimum is (1/2, 1/4) with f = −1/4 and λ = 1 (∇f + λ ∇g = (−1 + 2λx₁,
/// 1 − λ) = 0). λ starts at 0, so the Lagrangian's Hessian λ ∇²g is zero at first, though the
/// constraint's Jacobian still couples x to λ
struct Parabola;

impl Problem for Parabola {
    fn dimension(&self) -> usize {
        2
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        (x[1] - x[0], vec![-1.0, 1.0])
    }
    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        CsrMatrix::from_triplets(2, 2, vec![]).unwrap()
    }
    fn has_nonlinear_constraints(&self) -> bool {
        true
    }
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        NonlinearConstraints {
            h: vec![],
            g: vec![x[0] * x[0] - x[1]],
            dh: CsrMatrix::from_triplets(0, 2, vec![]).unwrap(),
            dg: dense(&[&[2.0 * x[0], -1.0]]),
        }
    }
    fn lagrangian_hessian(&self, _: &[f64], m: NonlinearMultipliers<'_>, _: f64) -> CsrMatrix<f64> {
        CsrMatrix::from_triplets(2, 2, vec![(0, 0, 2.0 * m.eq[0])]).unwrap()
    }
}

#[test]
fn a_nonlinear_equality_couples_the_steps_from_the_start() {
    let bounds = Constraints {
        xmin: Some(vec![0.0, f64::NEG_INFINITY]),
        xmax: Some(vec![2.0, f64::INFINITY]),
        ..Constraints::default()
    };
    let s = mehrotra(&Parabola, &bounds, &[1.5, 0.0], &MehrotraOptions::default()).unwrap();
    eprintln!("parabola: {:?} after {} iterations", s.status, s.iterations);
    assert!(s.converged(), "{:?}", s.status);
    // Pinned (2026-11-07): with separate step lengths while the Hessian is zero, it takes 5
    assert_eq!(s.iterations, 6, "iterations");
    known(&s.x, &[0.5, 0.25], 6, "x");
    known(&[s.f], &[-0.25], 6, "f");
    known(&s.multipliers.eq_nonlin, &[1.0], 6, "lambda");
    let residuals = kkt_residuals(&Parabola, &bounds, &s.x, &s.multipliers).unwrap();
    assert!(residuals.max() < RESIDUAL, "{residuals:?}");
}
