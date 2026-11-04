// Tests MIPS against MATPOWER's MIPS on the test problems of MATPOWER's own MIPS tests
//   (mips/lib/t/t_mips.m and t_qps_mips.m), restated here: the same exit status and iteration
//   count, and the solution, objective and multipliers close to MATPOWER's, whose results
//   scripts/generate-mips-reference.sh records in tests/data/mips_reference.txt; and each
//   problem's known answer to the precision t_mips checks it.
//
// The problems, starting points and expected answers are MIPS's, from t_mips.m and
//   t_qps_mips.m, restated in Rust:
//   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
//   by Ray Zimmerman, PSERC Cornell
//   Covered by the 3-clause BSD License; see LICENSE-MIPS at the repository's root.
// Two of them are classics: Rosenbrock's function (H. H. Rosenbrock, "An automatic method for
//   finding the greatest or least value of a function," The Computer Journal 3(3):175-184,
//   1960, doi:10.1093/comjnl/3.3.175), through MATLAB's bandem.m, and problem 71 of W. Hock
//   and K. Schittkowski, Test Examples for Nonlinear Programming Codes, Lecture Notes in
//   Economics and Mathematical Systems 187, Springer, 1981, doi:10.1007/978-3-642-48320-2.
//   t_mips.m cites the sources of the others.

use std::collections::HashMap;

use interior_point::mips::{mips, qps_mips, MipsOptions, MipsSolution, MipsStatus};
use interior_point::problem::{
    Constraints, NonlinearConstraints, NonlinearMultipliers, Problem, QuadraticProblem,
};
use sparse_mat::compressed_matrix::CsrMatrix;

/// How far our results may be from MATPOWER's, relative to the larger magnitude (at least 1).
/// Both solve the same small systems, but by different LU factorizations (UMFPACK's in
/// MATPOWER), so they agree to rounding, compounded over the iterations: the largest
/// difference measured on a converged run is 6.5e-14 (hs71's lower-bound multiplier).
const AGREEMENT: f64 = 1e-12;

/// The same for the multipliers of a run that fails numerically, which have diverged (to about
/// 3e24 on the infeasible LP) and amplify rounding: 3.3e-7 measured there.
const DIVERGED: f64 = 1e-6;

fn dense(rows: &[&[f64]]) -> CsrMatrix<f64> {
    CsrMatrix::from_dense(rows).unwrap()
}

/// A recorded MATPOWER run: each field's values, by name.
struct Reference(HashMap<String, Vec<f64>>);

fn references() -> HashMap<(String, bool), Reference> {
    let text = include_str!("data/mips_reference.txt");
    let mut runs = HashMap::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let mut parts = line.splitn(3, ',');
        let name = parts.next().unwrap().to_string();
        let step_control = parts.next().unwrap() == "1";
        let mut fields = HashMap::new();
        for field in parts.next().unwrap().split('|') {
            let (key, values) = field.split_once('=').unwrap();
            let values = values
                .split_whitespace()
                .map(|v| v.parse::<f64>().unwrap())
                .collect();
            fields.insert(key.to_string(), values);
        }
        runs.insert((name, step_control), Reference(fields));
    }
    runs
}

fn close(ours: &[f64], theirs: &[f64], what: &str) {
    close_within(ours, theirs, AGREEMENT, what);
}

fn close_within(ours: &[f64], theirs: &[f64], tolerance: f64, what: &str) {
    assert_eq!(ours.len(), theirs.len(), "{what}: lengths");
    for (k, (a, b)) in ours.iter().zip(theirs).enumerate() {
        let scale = a.abs().max(b.abs()).max(1.0);
        assert!(
            (a - b).abs() <= tolerance * scale,
            "{what}[{k}]: ours {a:e}, MATPOWER's {b:e}"
        );
    }
}

/// Checks a solution against MATPOWER's run of the same problem.
fn check(name: &str, step_control: bool, solution: &MipsSolution) {
    let runs = references();
    let reference = &runs[&(name.to_string(), step_control)];
    let field = |key: &str| reference.0.get(key);
    let eflag = field("eflag").unwrap()[0];
    let status = match eflag as i64 {
        1 => MipsStatus::Converged,
        0 => MipsStatus::MaxIterations,
        _ => MipsStatus::NumericalFailure,
    };
    let what = format!("{name} (step control {step_control})");
    assert_eq!(solution.status, status, "{what}: status");
    assert_eq!(
        solution.iterations,
        field("iterations").unwrap()[0] as usize,
        "{what}: iterations"
    );
    close(&[solution.f], field("f").unwrap(), &format!("{what}: f"));
    close(&solution.x, field("x").unwrap(), &format!("{what}: x"));
    let m = &solution.multipliers;
    for (key, ours) in [
        ("eqnonlin", &m.eq_nonlin),
        ("ineqnonlin", &m.ineq_nonlin),
        ("mu_l", &m.mu_l),
        ("mu_u", &m.mu_u),
        ("lower", &m.lower),
        ("upper", &m.upper),
    ] {
        let tolerance = if status == MipsStatus::NumericalFailure {
            DIVERGED
        } else {
            AGREEMENT
        };
        match field(key) {
            Some(theirs) => close_within(ours, theirs, tolerance, &format!("{what}: {key}")),
            // MIPS leaves out the nonlinear multipliers when there are no such constraints
            None => assert!(ours.is_empty(), "{what}: {key} should be empty"),
        }
    }
}

/// t_mips's check: within 10^-digits of the known answer.
fn known(ours: &[f64], answer: &[f64], digits: i32, what: &str) {
    for (k, (a, b)) in ours.iter().zip(answer).enumerate() {
        assert!(
            (a - b).abs() < 10f64.powi(-digits),
            "{what}[{k}]: {a} against the known {b}"
        );
    }
}

fn options(step_control: bool) -> MipsOptions {
    let mut options = MipsOptions::default();
    options.step_control = step_control;
    options
}

// The problems of t_mips.m

/// Rosenbrock's banana function, unconstrained
struct Banana;

impl Problem for Banana {
    fn dimension(&self) -> usize {
        2
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        let a = 100.0;
        let f = a * (x[1] - x[0].powi(2)).powi(2) + (1.0 - x[0]).powi(2);
        let df = vec![
            4.0 * a * (x[0].powi(3) - x[0] * x[1]) + 2.0 * x[0] - 2.0,
            2.0 * a * (x[1] - x[0].powi(2)),
        ];
        (f, df)
    }
    fn objective_hessian(&self, x: &[f64]) -> CsrMatrix<f64> {
        let a = 100.0;
        dense(&[
            &[
                4.0 * a * (3.0 * x[0].powi(2) - x[1] + 1.0 / (2.0 * a)),
                -4.0 * a * x[0],
            ],
            &[-4.0 * a * x[0], 4.0 * a * 0.5],
        ])
    }
}

/// An unconstrained 3-d quadratic, with a constant term
struct Quadratic3;

impl Problem for Quadratic3 {
    fn dimension(&self) -> usize {
        3
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        let qp = QuadraticProblem {
            h: Some(Quadratic3.objective_hessian(x)),
            c: vec![2.0, -35.0, -47.0],
        };
        let (f, df) = qp.objective(x);
        (f + 5.0, df)
    }
    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        dense(&[&[5.0, -2.0, -1.0], &[-2.0, 4.0, 3.0], &[-1.0, 3.0, 5.0]])
    }
}

fn qp4_hessian() -> CsrMatrix<f64> {
    dense(&[
        &[1003.1, 4.3, 6.3, 5.9],
        &[4.3, 2.2, 2.1, 3.9],
        &[6.3, 2.1, 3.5, 4.8],
        &[5.9, 3.9, 4.8, 10.0],
    ])
}

fn qp4_constraints() -> Constraints {
    Constraints {
        a: Some(dense(&[&[1.0, 1.0, 1.0, 1.0], &[0.17, 0.11, 0.10, 0.18]])),
        l: Some(vec![1.0, 0.10]),
        u: Some(vec![1.0, f64::INFINITY]),
        xmin: Some(vec![0.0; 4]),
        xmax: None,
    }
}

/// A linear objective with two quadratic inequalities (a 2-d nonlinear program)
struct Nonlinear2;

impl Problem for Nonlinear2 {
    fn dimension(&self) -> usize {
        2
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        (-x[0] - x[1], vec![-1.0, -1.0])
    }
    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        CsrMatrix::from_triplets(2, 2, vec![]).unwrap()
    }
    fn has_nonlinear_constraints(&self) -> bool {
        true
    }
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        let s = x[0].powi(2) + x[1].powi(2);
        NonlinearConstraints {
            h: vec![-s + 1.0, s - 2.0],
            g: vec![],
            dh: dense(&[&[-2.0 * x[0], -2.0 * x[1]], &[2.0 * x[0], 2.0 * x[1]]]),
            dg: CsrMatrix::from_triplets(0, 2, vec![]).unwrap(),
        }
    }
    fn lagrangian_hessian(&self, _: &[f64], m: NonlinearMultipliers<'_>, _: f64) -> CsrMatrix<f64> {
        let d = 2.0 * (-m.ineq[0] + m.ineq[1]);
        dense(&[&[d, 0.0], &[0.0, d]])
    }
}

/// A bilinear objective with two quadratic inequalities (a 3-d nonlinear program)
struct Nonlinear3;

impl Problem for Nonlinear3 {
    fn dimension(&self) -> usize {
        3
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        (
            -x[0] * x[1] - x[1] * x[2],
            vec![-x[1], -(x[0] + x[2]), -x[1]],
        )
    }
    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        dense(&[&[0.0, -1.0, 0.0], &[-1.0, 0.0, -1.0], &[0.0, -1.0, 0.0]])
    }
    fn has_nonlinear_constraints(&self) -> bool {
        true
    }
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        let sq: Vec<f64> = x.iter().map(|v| v * v).collect();
        NonlinearConstraints {
            h: vec![sq[0] - sq[1] + sq[2] - 2.0, sq[0] + sq[1] + sq[2] - 10.0],
            g: vec![],
            dh: dense(&[
                &[2.0 * x[0], -2.0 * x[1], 2.0 * x[2]],
                &[2.0 * x[0], 2.0 * x[1], 2.0 * x[2]],
            ]),
            dg: CsrMatrix::from_triplets(0, 3, vec![]).unwrap(),
        }
    }
    fn lagrangian_hessian(
        &self,
        x: &[f64],
        m: NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        let mu = m.ineq;
        let objective = self.objective_hessian(x).scale(cost_mult);
        let diagonal = dense(&[
            &[2.0 * (mu[0] + mu[1]), 0.0, 0.0],
            &[0.0, 2.0 * (-mu[0] + mu[1]), 0.0],
            &[0.0, 0.0, 2.0 * (mu[0] + mu[1])],
        ]);
        objective.try_add(&diagonal).unwrap()
    }
}

/// Hock and Schittkowski's problem 71: one nonlinear equality, one nonlinear inequality, and
/// bounds
struct Hs71;

impl Hs71 {
    fn hessian_of_f(x: &[f64]) -> CsrMatrix<f64> {
        let s = 2.0 * x[0] + x[1] + x[2];
        dense(&[
            &[2.0 * x[3], x[3], x[3], s],
            &[x[3], 0.0, 0.0, x[0]],
            &[x[3], 0.0, 0.0, x[0]],
            &[s, x[0], x[0], 0.0],
        ])
    }
}

impl Problem for Hs71 {
    fn dimension(&self) -> usize {
        4
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        let sum = x[0] + x[1] + x[2];
        let f = x[0] * x[3] * sum + x[2];
        let df = vec![
            x[0] * x[3] + x[3] * sum,
            x[0] * x[3],
            x[0] * x[3] + 1.0,
            x[0] * sum,
        ];
        (f, df)
    }
    fn objective_hessian(&self, x: &[f64]) -> CsrMatrix<f64> {
        Hs71::hessian_of_f(x)
    }
    fn has_nonlinear_constraints(&self) -> bool {
        true
    }
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        let product: f64 = x.iter().product();
        NonlinearConstraints {
            h: vec![-product + 25.0],
            g: vec![x.iter().map(|v| v * v).sum::<f64>() - 40.0],
            dh: dense(&[&x.iter().map(|v| -product / v).collect::<Vec<_>>()]),
            dg: dense(&[&x.iter().map(|v| 2.0 * v).collect::<Vec<_>>()]),
        }
    }
    fn lagrangian_hessian(
        &self,
        x: &[f64],
        m: NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        let (lambda, mu) = (m.eq[0], m.ineq[0]);
        let identity = dense(&[
            &[2.0, 0.0, 0.0, 0.0],
            &[0.0, 2.0, 0.0, 0.0],
            &[0.0, 0.0, 2.0, 0.0],
            &[0.0, 0.0, 0.0, 2.0],
        ]);
        let products = dense(&[
            &[0.0, x[2] * x[3], x[1] * x[3], x[1] * x[2]],
            &[x[2] * x[3], 0.0, x[0] * x[3], x[0] * x[2]],
            &[x[1] * x[3], x[0] * x[3], 0.0, x[0] * x[1]],
            &[x[1] * x[2], x[0] * x[2], x[0] * x[1], 0.0],
        ]);
        Hs71::hessian_of_f(x)
            .scale(cost_mult)
            .try_add(&identity.scale(lambda))
            .unwrap()
            .try_sub(&products.scale(mu))
            .unwrap()
    }
}

#[test]
fn unconstrained_problems_match_matpower() {
    for sc in [false, true] {
        let s = mips(&Banana, &[-1.9, 2.0], &Constraints::default(), &options(sc)).unwrap();
        check("banana", sc, &s);
        known(&s.x, &[1.0, 1.0], 13, "banana x");
        known(&[s.f], &[0.0], 13, "banana f");

        let s = mips(
            &Quadratic3,
            &[0.0; 3],
            &Constraints::default(),
            &options(sc),
        )
        .unwrap();
        check("quadratic3", sc, &s);
        known(&s.x, &[3.0, 5.0, 7.0], 13, "quadratic3 x");
        known(&[s.f], &[-244.0], 13, "quadratic3 f");
    }
}

#[test]
fn a_linearly_constrained_qp_matches_matpower() {
    for sc in [false, true] {
        let qp = QuadraticProblem {
            h: Some(qp4_hessian()),
            c: vec![0.0; 4],
        };
        let s = mips(&qp, &[1.0, 0.0, 0.0, 1.0], &qp4_constraints(), &options(sc)).unwrap();
        check("qp4", sc, &s);
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
}

#[test]
fn nonlinear_programs_match_matpower() {
    for sc in [false, true] {
        let bounds = Constraints {
            xmin: Some(vec![0.0; 2]),
            ..Constraints::default()
        };
        let s = mips(&Nonlinear2, &[1.1, 0.0], &bounds, &options(sc)).unwrap();
        check("nonlinear2", sc, &s);
        known(&s.x, &[1.0, 1.0], 6, "nonlinear2 x");
        known(&s.multipliers.ineq_nonlin, &[0.0, 0.5], 6, "nonlinear2 mu");

        let s = mips(
            &Nonlinear3,
            &[1.0, 1.0, 0.0],
            &Constraints::default(),
            &options(sc),
        )
        .unwrap();
        check("nonlinear3", sc, &s);
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

        let bounds = Constraints {
            xmin: Some(vec![1.0; 4]),
            xmax: Some(vec![5.0; 4]),
            ..Constraints::default()
        };
        let s = mips(&Hs71, &[1.0, 5.0, 5.0, 1.0], &bounds, &options(sc)).unwrap();
        check("hs71", sc, &s);
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
}

#[test]
fn step_control_and_the_iteration_limit_match_matpower() {
    let bounds = Constraints {
        xmin: Some(vec![1.0; 4]),
        xmax: Some(vec![5.0; 4]),
        ..Constraints::default()
    };
    // From (4.5, 4.5, 4.5, 4.5), step control changes the run: without it MIPS stops at a
    //   local point (f ≈ 35.28) in 11 iterations, with it reaches the optimum in 12
    let plain = mips(&Hs71, &[4.5; 4], &bounds, &options(false)).unwrap();
    let controlled = mips(&Hs71, &[4.5; 4], &bounds, &options(true)).unwrap();
    check("hs71_from_4.5", false, &plain);
    check("hs71_from_4.5", true, &controlled);
    assert_ne!(plain.iterations, controlled.iterations);
    known(
        &[controlled.f],
        &[17.0140173],
        6,
        "hs71 from 4.5, step control, f",
    );
    // From (2, 2, 2, 2), MIPS reaches its iteration limit
    for sc in [false, true] {
        let s = mips(&Hs71, &[2.0; 4], &bounds, &options(sc)).unwrap();
        check("hs71_from_2", sc, &s);
        assert_eq!(s.status, MipsStatus::MaxIterations);
        assert_eq!(s.history.len(), 151);
    }
}

#[test]
fn quadratic_and_linear_programs_match_matpower() {
    for sc in [false, true] {
        let mut opt = options(sc);
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
        let s = qps_mips(&lp, &constraints, None, &opt).unwrap();
        check("lp3", sc, &s);
        known(&s.x, &[0.0, 15.0, 3.0], 6, "lp3 x");
        known(&s.multipliers.mu_l, &[0.0, 1.5, 0.0], 9, "lp3 mu_l");
        known(&s.multipliers.mu_u, &[0.0, 0.0, 0.5], 9, "lp3 mu_u");

        let qp = QuadraticProblem {
            h: Some(Quadratic3.objective_hessian(&[])),
            c: vec![2.0, -35.0, -47.0],
        };
        let s = qps_mips(&qp, &Constraints::default(), Some(&[0.0; 3]), &opt).unwrap();
        check("qps_quadratic3", sc, &s);
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
        let s = qps_mips(&qp, &constraints, None, &opt).unwrap();
        check("qp2", sc, &s);
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
        let s = qps_mips(&qp, &qp4_constraints(), Some(&[1.0, 0.0, 0.0, 1.0]), &opt).unwrap();
        check("qps_qp4", sc, &s);
        known(&[s.f], &[3.29 / 3.0], 6, "qps_qp4 f");

        // Infeasible: x₁ + x₂ ≤ −1 with x ≥ 0. MIPS fails numerically, as MATPOWER's does
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
        let s = qps_mips(&lp, &constraints, None, &opt).unwrap();
        check("infeasible_lp", sc, &s);
        assert!(!s.converged());
    }
}
