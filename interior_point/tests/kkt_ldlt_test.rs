// Tests the solvers with the KKT matrix factored by the pivoted LDLᵀ (KktFactorization::Ldlt), on
//   MATPOWER's MIPS test problems (tests/common): each solver must end as it does with LU, with
//   the objective within 1e-6 (relative) and small KKT residuals, the bar of phase O3. The
//   iterations are printed beside LU's, not pinned to them. Then the inertia correction of the
//   predictor-corrector method, on a nonconvex quadratic program and on hs71.

mod common;

use common::{qp4_constraints, qp4_hessian, Banana, Hs71, Nonlinear2, Nonlinear3, Quadratic3};
use interior_point::kkt::KktFactorization;
use interior_point::mehrotra::{mehrotra, qps_mehrotra, MehrotraIteration, MehrotraOptions};
use interior_point::mips::{mips, MipsOptions};
use interior_point::problem::{Constraints, Problem, QuadraticProblem};
use interior_point::solution::{kkt_residuals, History, Multipliers, Solution, Status};
use sparse_mat::compressed_matrix::CsrMatrix;

/// The bar for the objective, relative to the larger magnitude (at least 1)
const OBJECTIVE: f64 = 1e-6;

/// The largest KKT residual accepted, unscaled, as in the predictor-corrector method's tests
const RESIDUAL: f64 = 1e-5;

/// One solve's outcome
struct Outcome {
    status: Status,
    iterations: usize,
    f: f64,
    x: Vec<f64>,
    multipliers: Multipliers,
    refinement_steps: usize,
}

fn mips_options(kkt: KktFactorization) -> MipsOptions {
    let mut options = MipsOptions::default();
    options.kkt_factorization = kkt;
    options
}

fn mehrotra_options(kkt: KktFactorization) -> MehrotraOptions {
    let mut options = MehrotraOptions::default();
    options.kkt_factorization = kkt;
    options
}

/// Both solvers on `problem` from `x0`, with LU and with LDLᵀ: LDLᵀ must end as LU does, and,
/// where they converge, agree on the objective with small KKT residuals
fn compare<P: Problem>(name: &str, problem: &P, x0: &[f64], constraints: &Constraints) {
    let run = |solver: &str, kkt: KktFactorization| -> Outcome {
        if solver == "MIPS" {
            let s = mips(problem, constraints, x0, &mips_options(kkt)).unwrap();
            Outcome {
                status: s.status,
                iterations: s.iterations,
                f: s.f,
                x: s.x,
                multipliers: s.multipliers,
                refinement_steps: s.refinement_steps,
            }
        } else {
            let s = mehrotra(problem, constraints, x0, &mehrotra_options(kkt)).unwrap();
            Outcome {
                status: s.status,
                iterations: s.iterations,
                f: s.f,
                x: s.x,
                multipliers: s.multipliers,
                refinement_steps: s.refinement_steps,
            }
        }
    };
    for solver in ["MIPS", "predictor-corrector"] {
        let lu = run(solver, KktFactorization::Lu);
        let ldlt = run(solver, KktFactorization::Ldlt);
        assert_eq!(ldlt.status, lu.status, "{name}, {solver}: status");
        assert_eq!(lu.refinement_steps, 0, "{name}, {solver}: LU isn't refined");
        let objective = (ldlt.f - lu.f).abs() / lu.f.abs().max(ldlt.f.abs()).max(1.0);
        let residual = kkt_residuals(problem, constraints, &ldlt.x, &ldlt.multipliers)
            .unwrap()
            .max();
        eprintln!(
            "{name}, {solver}: {:?} in {} iterations (LU {}), objective {objective:.1e} from LU's, \
             residual {residual:.1e}, {} refinement steps",
            ldlt.status, ldlt.iterations, lu.iterations, ldlt.refinement_steps
        );
        if lu.status == Status::Converged {
            assert!(
                objective <= OBJECTIVE,
                "{name}, {solver}: objective {objective:e}"
            );
            assert!(
                residual <= RESIDUAL,
                "{name}, {solver}: residual {residual:e}"
            );
        }
    }
}

/// A predictor-corrector solve's iteration records
fn mehrotra_history(s: &Solution) -> &[MehrotraIteration] {
    match &s.history {
        History::Mehrotra(records) => records,
        History::Mips(_) => panic!("a MIPS history from the predictor-corrector method"),
    }
}

#[test]
fn both_solvers_end_as_they_do_with_lu() {
    let none = Constraints::default();
    compare("banana", &Banana, &[-1.9, 2.0], &none);
    compare("quadratic3", &Quadratic3, &[0.0; 3], &none);
    let qp4 = QuadraticProblem {
        h: Some(qp4_hessian()),
        c: vec![0.0; 4],
    };
    compare("qp4", &qp4, &[1.0, 0.0, 0.0, 1.0], &qp4_constraints());
    let nonnegative = |n| Constraints {
        xmin: Some(vec![0.0; n]),
        ..Constraints::default()
    };
    compare("nonlinear2", &Nonlinear2, &[1.1, 0.0], &nonnegative(2));
    compare("nonlinear3", &Nonlinear3, &[1.0, 1.0, 0.0], &none);
    let hs71_bounds = Constraints {
        xmin: Some(vec![1.0; 4]),
        xmax: Some(vec![5.0; 4]),
        ..Constraints::default()
    };
    compare("hs71", &Hs71, &[1.0, 5.0, 5.0, 1.0], &hs71_bounds);
    compare("hs71 from 4.5", &Hs71, &[4.5; 4], &hs71_bounds);
}

/// `min ½ (x₂² − x₁²)` over the box `−1 ≤ x ≤ 1`: a saddle at the origin, minima at
/// `x₁ = ±1, x₂ = 0` with `f = −½`
fn saddle() -> (QuadraticProblem, Constraints) {
    let h = CsrMatrix::from_triplets(2, 2, vec![(0, 0, -1.0), (1, 1, 1.0)]).unwrap();
    (
        QuadraticProblem {
            h: Some(h),
            c: vec![0.0; 2],
        },
        Constraints {
            xmin: Some(vec![-1.0; 2]),
            xmax: Some(vec![1.0; 2]),
            ..Constraints::default()
        },
    )
}

#[test]
fn the_inertia_correction_leaves_a_saddle_for_a_minimum() {
    let (qp, constraints) = saddle();
    let x0 = [0.1, 0.1];
    let lu = qps_mehrotra(
        &qp,
        &constraints,
        Some(&x0),
        &mehrotra_options(KktFactorization::Lu),
    )
    .unwrap();
    let mut corrected = mehrotra_options(KktFactorization::Ldlt);
    corrected.inertia_correction = true;
    let ldlt = qps_mehrotra(&qp, &constraints, Some(&x0), &corrected).unwrap();
    eprintln!(
        "saddle: LU {:?} at {:?}, f = {}, in {}; LDLᵀ {:?} at {:?}, f = {}, in {}, {} corrections",
        lu.status,
        lu.x,
        lu.f,
        lu.iterations,
        ldlt.status,
        ldlt.x,
        ldlt.f,
        ldlt.iterations,
        ldlt.inertia_corrections
    );
    assert_eq!(lu.inertia_corrections, 0);
    assert!(ldlt.inertia_corrections > 0);
    assert!(mehrotra_history(&ldlt).iter().any(|h| h.delta_w > 0.0));
    assert_eq!(ldlt.status, Status::Converged);
    assert!((ldlt.f + 0.5).abs() < 1e-6, "f = {}", ldlt.f);
    assert!((ldlt.x[0].abs() - 1.0).abs() < 1e-6 && ldlt.x[1].abs() < 1e-6);
}

#[test]
fn without_the_correction_ldlt_stops_at_the_saddle_as_lu_does() {
    let (qp, constraints) = saddle();
    // The default
    let options = mehrotra_options(KktFactorization::Ldlt);
    assert!(!options.inertia_correction);
    let s = qps_mehrotra(&qp, &constraints, Some(&[0.1, 0.1]), &options).unwrap();
    eprintln!(
        "saddle, LDLᵀ uncorrected: {:?} at {:?}, f = {}, in {}",
        s.status, s.x, s.f, s.iterations
    );
    assert_eq!(s.status, Status::Converged);
    assert_eq!(s.inertia_corrections, 0);
    assert!(mehrotra_history(&s).iter().all(|h| h.delta_w == 0.0));
    assert!(s.f.abs() < 1e-6 && s.x[0].abs() < 1e-6, "f = {}", s.f);
}

#[test]
fn hs71_from_2_is_refused_as_a_maximum() {
    // With LU the method converges here to a constrained local maximum (f = 131.26, see
    //   mehrotra_test). The inertia correction refuses it, but with no line search the method
    //   then fails numerically (2026-11-07), as MIPS fails here too (150 iterations)
    let bounds = Constraints {
        xmin: Some(vec![1.0; 4]),
        xmax: Some(vec![5.0; 4]),
        ..Constraints::default()
    };
    let mut corrected = mehrotra_options(KktFactorization::Ldlt);
    corrected.inertia_correction = true;
    let s = mehrotra(&Hs71, &bounds, &[2.0; 4], &corrected).unwrap();
    eprintln!(
        "hs71 from 2: {:?} after {} iterations, f = {}, {} corrections",
        s.status, s.iterations, s.f, s.inertia_corrections
    );
    assert_eq!(s.status, Status::NumericalFailure);
    assert!(s.inertia_corrections > 0);
}
