// Tests that both solvers keep the KKT matrix's analysis while its pattern holds and analyze it
//   afresh when it changes: a problem whose Hessian's pattern changes part way through the solve
//   still solves correctly, with a second analysis.

use interior_point::mehrotra::{mehrotra, MehrotraOptions};
use interior_point::mips::{mips, MipsOptions};
use interior_point::problem::{Constraints, Problem};
use sparse_mat::compressed_matrix::CsrMatrix;

/// min (x₁ − 3)² + (x₂ − 3)² subject to x ≤ 2, whose optimum is (2, 2) with upper multipliers 2.
/// Its Hessian is 2I, but once x₁ passes 1.5 it also stores its off-diagonal zeros, so the KKT
/// matrix's pattern changes there while every value stays the same
struct ChangingPattern;

impl Problem for ChangingPattern {
    fn dimension(&self) -> usize {
        2
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        let d = [x[0] - 3.0, x[1] - 3.0];
        (d[0] * d[0] + d[1] * d[1], vec![2.0 * d[0], 2.0 * d[1]])
    }
    fn objective_hessian(&self, x: &[f64]) -> CsrMatrix<f64> {
        let mut entries = vec![(0, 0, 2.0), (1, 1, 2.0)];
        if x[0] > 1.5 {
            entries.extend([(0, 1, 0.0), (1, 0, 0.0)]);
        }
        CsrMatrix::from_triplets(2, 2, entries).unwrap()
    }
}

fn bounds() -> Constraints {
    Constraints {
        xmax: Some(vec![2.0; 2]),
        ..Constraints::default()
    }
}

fn assert_solved(what: &str, x: &[f64], upper: &[f64]) {
    for (k, (&x, &upper)) in x.iter().zip(upper).enumerate() {
        assert!((x - 2.0).abs() < 1e-6, "{what}: x[{k}] = {x}");
        assert!((upper - 2.0).abs() < 1e-5, "{what}: upper[{k}] = {upper}");
    }
}

#[test]
fn mips_reanalyzes_when_the_pattern_changes() {
    let s = mips(
        &ChangingPattern,
        &[0.0, 0.0],
        &bounds(),
        &MipsOptions::default(),
    )
    .unwrap();
    assert!(s.converged(), "{:?}", s.status);
    assert_solved("MIPS", &s.x, &s.multipliers.upper);
    assert_eq!(
        s.analyses, 2,
        "one analysis before x₁ passes 1.5, one after"
    );
}

#[test]
fn mehrotra_reanalyzes_when_the_pattern_changes() {
    let s = mehrotra(
        &ChangingPattern,
        &[0.0, 0.0],
        &bounds(),
        &MehrotraOptions::default(),
    )
    .unwrap();
    assert!(s.converged(), "{:?}", s.status);
    assert_solved("Mehrotra", &s.x, &s.multipliers.upper);
    assert_eq!(
        s.analyses, 2,
        "one analysis before x₁ passes 1.5, one after"
    );
}
