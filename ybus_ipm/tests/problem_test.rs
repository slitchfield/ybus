// Tests for the problem interface: QuadraticProblem's objective, gradient and Hessian.

use ybus_ipm::problem::{Problem, QuadraticProblem};
use ybus_sparse::compressed_matrix::CsrMatrix;

#[test]
fn a_quadratic_problem_evaluates_its_objective() {
    // f(x) = ½ xᵀ H x + cᵀ x with H = [[2, 1], [1, 4]], c = (1, −1), at x = (1, 2):
    //   ½ (2 + 4 + 16) + (1 − 2) = 10, ∇f = H x + c = (5, 8)
    let h = CsrMatrix::from_triplets(
        2,
        2,
        vec![(0, 0, 2.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 4.0)],
    )
    .unwrap();
    let qp = QuadraticProblem {
        h: Some(h.clone()),
        c: vec![1.0, -1.0],
    };
    assert_eq!(qp.dimension(), 2);
    assert_eq!(qp.objective(&[1.0, 2.0]), (10.0, vec![5.0, 8.0]));
    assert_eq!(qp.objective_hessian(&[1.0, 2.0]), h);
    assert!(!qp.has_nonlinear_constraints());

    // A linear program: no H
    let lp = QuadraticProblem {
        h: None,
        c: vec![3.0, -2.0],
    };
    assert_eq!(lp.objective(&[1.0, 2.0]), (-1.0, vec![3.0, -2.0]));
    assert_eq!(lp.objective_hessian(&[0.0, 0.0]).num_stored(), 0);
}
