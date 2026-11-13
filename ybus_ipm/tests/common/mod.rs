//! Shared by the solver tests: the test problems of MATPOWER's MIPS tests (mips/lib/t/t_mips.m
//!   and t_qps_mips.m), restated here, MATPOWER's recorded results on them
//!   (`data/mips_reference.txt`, from scripts/generate-mips-reference.sh), and t_mips's check
//!   against a problem's known answer.
//!
//! The problems, starting points and expected answers are MIPS's, from t_mips.m and
//!   t_qps_mips.m, restated in Rust:
//!   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
//!   by Ray Zimmerman, PSERC Cornell
//!   Covered by the 3-clause BSD License; see LICENSE-MIPS at the repository's root.
//! Two of them are classics: Rosenbrock's function (H. H. Rosenbrock, "An automatic method for
//!   finding the greatest or least value of a function," The Computer Journal 3(3):175-184,
//!   1960, doi:10.1093/comjnl/3.3.175), through MATLAB's bandem.m, and problem 71 of W. Hock
//!   and K. Schittkowski, Test Examples for Nonlinear Programming Codes, Lecture Notes in
//!   Economics and Mathematical Systems 187, Springer, 1981, doi:10.1007/978-3-642-48320-2.
//!   t_mips.m cites the sources of the others.

// Each test file compiles this module separately and uses only part of it
#![allow(dead_code)]

use std::collections::HashMap;

use ybus_ipm::problem::{
    Constraints, NonlinearConstraints, NonlinearMultipliers, Problem, QuadraticProblem,
};
use ybus_sparse::compressed_matrix::CsrMatrix;

pub fn dense(rows: &[&[f64]]) -> CsrMatrix<f64> {
    CsrMatrix::from_dense(rows).unwrap()
}

/// A recorded MATPOWER run: each field's values, by name.
pub struct Reference(pub HashMap<String, Vec<f64>>);

pub fn references() -> HashMap<(String, bool), Reference> {
    let text = include_str!("../data/mips_reference.txt");
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

/// t_mips's check: within 10^-digits of the known answer.
pub fn known(ours: &[f64], answer: &[f64], digits: i32, what: &str) {
    for (k, (a, b)) in ours.iter().zip(answer).enumerate() {
        assert!(
            (a - b).abs() < 10f64.powi(-digits),
            "{what}[{k}]: {a} against the known {b}"
        );
    }
}

// The problems of t_mips.m

/// Rosenbrock's banana function, unconstrained
pub struct Banana;

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
pub struct Quadratic3;

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

pub fn qp4_hessian() -> CsrMatrix<f64> {
    dense(&[
        &[1003.1, 4.3, 6.3, 5.9],
        &[4.3, 2.2, 2.1, 3.9],
        &[6.3, 2.1, 3.5, 4.8],
        &[5.9, 3.9, 4.8, 10.0],
    ])
}

pub fn qp4_constraints() -> Constraints {
    Constraints {
        a: Some(dense(&[&[1.0, 1.0, 1.0, 1.0], &[0.17, 0.11, 0.10, 0.18]])),
        l: Some(vec![1.0, 0.10]),
        u: Some(vec![1.0, f64::INFINITY]),
        xmin: Some(vec![0.0; 4]),
        xmax: None,
    }
}

/// A linear objective with two quadratic inequalities (a 2-d nonlinear program)
pub struct Nonlinear2;

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
pub struct Nonlinear3;

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
pub struct Hs71;

impl Hs71 {
    pub fn hessian_of_f(x: &[f64]) -> CsrMatrix<f64> {
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
