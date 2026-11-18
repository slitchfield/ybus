// Warm starts (ybus_ipm::start), on both solvers: the inverse mapping of the multipliers, tested
//   by iteration 0's measures on a QP with every kind of row; convergence from a problem's own
//   solution and from a neighbour's, within O3's bar; the cold path, bit for bit; and errors.

mod common;

use common::{qp4_constraints, qp4_hessian, Banana, Hs71, Nonlinear2, Nonlinear3, Quadratic3};
use ybus_ipm::interior_point_error::InteriorPointError;
use ybus_ipm::mehrotra::{mehrotra, mehrotra_from, qps_mehrotra_from, MehrotraOptions};
use ybus_ipm::mips::{mips, mips_from, qps_mips, qps_mips_from, MipsOptions};
use ybus_ipm::problem::{Constraints, Problem, QuadraticProblem};
use ybus_ipm::solution::{History, Solution, Status};
use ybus_ipm::start::{Start, WarmStart};
use ybus_sparse::compressed_matrix::CsrMatrix;

/// O3's bar for the objective, relative to the larger magnitude (at least 1)
const OBJECTIVE: f64 = 1e-6;

const INF: f64 = f64::INFINITY;

fn matrix(rows: usize, cols: usize, entries: &[(usize, usize, f64)]) -> CsrMatrix<f64> {
    CsrMatrix::from_triplets(rows, cols, entries.to_vec()).unwrap()
}

/// A QP with every kind of row: bounds lower-only (x0), upper-only (x1), boxed (x2) and fixed
///   (x3); rows of A an equality, upper-only, lower-only and boxed; and several of them active
fn every_row_kind(c: [f64; 4]) -> (QuadraticProblem, Constraints) {
    let qp = QuadraticProblem {
        h: Some(matrix(
            4,
            4,
            &[
                (0, 0, 1.0),
                (1, 1, 2.0),
                (2, 2, 1.0),
                (3, 3, 1.0),
                (0, 1, 0.5),
                (1, 0, 0.5),
            ],
        )),
        c: c.to_vec(),
    };
    let constraints = Constraints {
        a: Some(matrix(
            4,
            4,
            &[
                (0, 0, 1.0),
                (0, 1, 1.0),
                (0, 2, 1.0),
                (1, 0, 1.0),
                (1, 2, -1.0),
                (2, 1, 1.0),
                (2, 3, 1.0),
                (3, 0, 1.0),
                (3, 3, 1.0),
            ],
        )),
        l: Some(vec![2.0, -INF, -4.0, -2.0]),
        u: Some(vec![2.0, 1.5, INF, 2.0]),
        xmin: Some(vec![0.0, -INF, -1.0, 0.5]),
        xmax: Some(vec![INF, 5.0, 1.0, 0.5]),
    };
    (qp, constraints)
}

/// The predictor-corrector method's options with `warm_shift` and `max_it`
fn mehrotra_options(shift: f64, max_it: usize) -> MehrotraOptions {
    let mut options = MehrotraOptions::default();
    options.warm_shift = shift;
    options.max_it = max_it;
    options
}

fn mips_options(shift: f64, max_it: usize) -> MipsOptions {
    let mut options = MipsOptions::default();
    options.warm_shift = shift;
    options.max_it = max_it;
    options
}

/// Iteration `k`'s gradient and feasibility measures, either solver
fn measures(s: &Solution, k: usize) -> (f64, f64) {
    match &s.history {
        History::Mips(records) => (records[k].gradcond, records[k].feascond),
        History::Mehrotra(records) => (records[k].gradcond, records[k].feascond),
    }
}

fn last(s: &Solution) -> usize {
    match &s.history {
        History::Mips(records) => records.len() - 1,
        History::Mehrotra(records) => records.len() - 1,
    }
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= OBJECTIVE * a.abs().max(b.abs()).max(1.0)
}

#[test]
fn the_multipliers_map_back_onto_every_kind_of_row() {
    // Warm started from a solution with a negligible shift, and no iteration taken, the start's
    //   gradient measure is the solution's: any multiplier on the wrong row, or of the wrong sign,
    //   would make it large
    let (qp, constraints) = every_row_kind([-5.0, 3.0, -4.0, 0.0]);
    for solver in ["MIPS", "predictor-corrector"] {
        let solve = |start: Start<'_>, shift: f64, max_it: usize| {
            match solver {
                "MIPS" => qps_mips_from(&qp, &constraints, start, &mips_options(shift, max_it)),
                _ => qps_mehrotra_from(&qp, &constraints, start, &mehrotra_options(shift, max_it)),
            }
            .unwrap()
        };
        let cold = solve(Start::Cold(&[0.0; 4]), 1e-4, 150);
        assert_eq!(cold.status, Status::Converged, "{solver}");
        let warm = WarmStart::from_solution(&cold);
        let start = solve(Start::Warm(&warm), 1e-12, 0);
        let (gradcond, feascond) = measures(&start, 0);
        let (final_gradcond, final_feascond) = measures(&cold, last(&cold));
        eprintln!(
            "{solver}: the start's gradcond {gradcond:.1e} (solution's {final_gradcond:.1e}), \
             feascond {feascond:.1e} (solution's {final_feascond:.1e})"
        );
        assert!(
            gradcond <= (10.0 * final_gradcond).max(1e-6),
            "{solver}: gradcond {gradcond:e}"
        );
        assert!(
            feascond <= (10.0 * final_feascond).max(1e-6),
            "{solver}: feascond {feascond:e}"
        );
    }
}

/// A solver on `problem` from `start`
fn run<P: Problem>(
    solver: &str,
    problem: &P,
    constraints: &Constraints,
    start: Start<'_>,
) -> Solution {
    match solver {
        "MIPS" => mips_from(problem, constraints, start, &MipsOptions::default()),
        _ => mehrotra_from(problem, constraints, start, &MehrotraOptions::default()),
    }
    .unwrap()
}

fn from_own_solution<P: Problem>(name: &str, problem: &P, x0: &[f64], constraints: &Constraints) {
    for solver in ["MIPS", "predictor-corrector"] {
        let cold = run(solver, problem, constraints, Start::Cold(x0));
        if cold.status != Status::Converged {
            continue;
        }
        let warm = run(
            solver,
            problem,
            constraints,
            Start::Warm(&WarmStart::from_solution(&cold)),
        );
        eprintln!(
            "{name}, {solver}: cold {} iterations, warm {}",
            cold.iterations, warm.iterations
        );
        assert_eq!(warm.status, Status::Converged, "{name}, {solver}");
        assert!(
            close(warm.f, cold.f),
            "{name}, {solver}: {} vs {}",
            warm.f,
            cold.f
        );
        assert!(
            warm.iterations < cold.iterations,
            "{name}, {solver}: no fewer iterations"
        );
    }
}

#[test]
fn a_warm_start_from_its_own_solution_converges_sooner() {
    let none = Constraints::default();
    from_own_solution("banana", &Banana, &[-1.9, 2.0], &none);
    from_own_solution("quadratic3", &Quadratic3, &[0.0; 3], &none);
    let nonnegative = Constraints {
        xmin: Some(vec![0.0; 2]),
        ..Constraints::default()
    };
    from_own_solution("nonlinear2", &Nonlinear2, &[1.1, 0.0], &nonnegative);
    from_own_solution("nonlinear3", &Nonlinear3, &[1.0, 1.0, 0.0], &none);
    let hs71 = Constraints {
        xmin: Some(vec![1.0; 4]),
        xmax: Some(vec![5.0; 4]),
        ..Constraints::default()
    };
    from_own_solution("hs71", &Hs71, &[1.0, 5.0, 5.0, 1.0], &hs71);
    let qp4 = QuadraticProblem {
        h: Some(qp4_hessian()),
        c: vec![0.0; 4],
    };
    from_own_solution("qp4", &qp4, &[1.0, 0.0, 0.0, 1.0], &qp4_constraints());
    let (qp, constraints) = every_row_kind([-5.0, 3.0, -4.0, 0.0]);
    from_own_solution("every row kind", &qp, &[0.0; 4], &constraints);
}

#[test]
fn a_warm_start_from_a_neighbour_converges_to_the_neighbours_answer() {
    let (qp, constraints) = every_row_kind([-5.0, 3.0, -4.0, 0.0]);
    let (moved, _) = every_row_kind([-5.05, 3.03, -4.04, 0.0]);
    for solver in ["MIPS", "predictor-corrector"] {
        let first = run(solver, &qp, &constraints, Start::Cold(&[0.0; 4]));
        let cold = run(solver, &moved, &constraints, Start::Cold(&[0.0; 4]));
        let warm = run(
            solver,
            &moved,
            &constraints,
            Start::Warm(&WarmStart::from_solution(&first)),
        );
        eprintln!(
            "{solver}: the neighbour cold {} iterations, warm {}",
            cold.iterations, warm.iterations
        );
        assert_eq!(warm.status, Status::Converged, "{solver}");
        assert!(close(warm.f, cold.f), "{solver}: {} vs {}", warm.f, cold.f);
    }
}

#[test]
fn a_cold_start_is_the_original_path_bit_for_bit() {
    let (qp, constraints) = every_row_kind([-5.0, 3.0, -4.0, 0.0]);
    let bits = |s: &Solution| {
        format!(
            "{:?} {:?} {} {}",
            s.x,
            s.multipliers,
            s.f.to_bits(),
            s.iterations
        )
    };
    let x0 = [0.0; 4];
    let a = mips(&qp, &constraints, &x0, &MipsOptions::default()).unwrap();
    let b = mips_from(&qp, &constraints, Start::Cold(&x0), &MipsOptions::default()).unwrap();
    let c = qps_mips(&qp, &constraints, Some(&x0), &MipsOptions::default()).unwrap();
    assert_eq!(bits(&a), bits(&b));
    assert_eq!(bits(&a), bits(&c));
    let a = mehrotra(&qp, &constraints, &x0, &MehrotraOptions::default()).unwrap();
    let b = mehrotra_from(
        &qp,
        &constraints,
        Start::Cold(&x0),
        &MehrotraOptions::default(),
    )
    .unwrap();
    assert_eq!(bits(&a), bits(&b));
}

#[test]
fn a_warm_start_that_doesnt_fit_is_an_error() {
    let (qp, constraints) = every_row_kind([-5.0, 3.0, -4.0, 0.0]);
    let solution = mips(&qp, &constraints, &[0.0; 4], &MipsOptions::default()).unwrap();
    let mut short = WarmStart::from_solution(&solution);
    short.multipliers.mu_l.pop();
    for result in [
        mips_from(
            &qp,
            &constraints,
            Start::Warm(&short),
            &MipsOptions::default(),
        ),
        mehrotra_from(
            &qp,
            &constraints,
            Start::Warm(&short),
            &MehrotraOptions::default(),
        ),
    ] {
        assert!(matches!(result, Err(InteriorPointError::Dimension { .. })));
    }
    let warm = WarmStart::from_solution(&solution);
    assert!(matches!(
        mips_from(
            &qp,
            &constraints,
            Start::Warm(&warm),
            &mips_options(-1.0, 150)
        ),
        Err(InteriorPointError::InvalidOption { name: "warm_shift" })
    ));
    assert!(matches!(
        mehrotra_from(
            &qp,
            &constraints,
            Start::Warm(&warm),
            &mehrotra_options(f64::NAN, 150)
        ),
        Err(InteriorPointError::InvalidOption { name: "warm_shift" })
    ));
}
