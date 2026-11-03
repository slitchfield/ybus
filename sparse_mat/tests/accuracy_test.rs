// Tests for Phase 5's accuracy tools: Lu::solve_transpose, iterative refinement (solve_refined)
//   and the 1-norm condition estimate, against dense computations.

use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix, LuOptions, Ordering};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::Complex;

/// Deterministic pseudo-random numbers (a linear congruential generator).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// A pseudo-random SPD matrix with a dominant diagonal.
fn random_spd(n: usize, fill_percent: u64, seed: u64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    let mut row_sums = vec![0.0; n];
    for i in 0..n {
        for j in 0..i {
            if rng.next() % 100 < fill_percent {
                let value = -((rng.next() % 1000) as f64 + 1.0) / 997.0;
                triplets.extend([(i, j, value), (j, i, value)]);
                row_sums[i] -= value;
                row_sums[j] -= value;
            }
        }
    }
    for (i, sum) in row_sums.iter().enumerate() {
        triplets.push((i, i, sum + 1.0 + i as f64 / n as f64));
    }
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// A pseudo-random unsymmetric matrix whose diagonal is `diagonal_scale` times its row's
/// largest entry, so a small scale forces pivoting off the diagonal.
fn random_square(n: usize, fill_percent: u64, diagonal_scale: f64, seed: u64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for i in 0..n {
        let mut largest: f64 = 1.0;
        for j in 0..n {
            if i != j && rng.next() % 100 < fill_percent {
                let value = ((rng.next() % 2001) as f64 - 1000.0) / 313.0;
                largest = largest.max(value.abs());
                triplets.push((i, j, value));
            }
        }
        let sign = if rng.next().is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        triplets.push((
            i,
            i,
            sign * diagonal_scale * largest * (1.0 + i as f64 / n as f64),
        ));
    }
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// The n × n Hilbert matrix, 1 / (i + j + 1): dense, SPD and badly conditioned.
fn hilbert(n: usize) -> CsrMatrix<f64> {
    let triplets: Vec<(usize, usize, f64)> = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j, 1.0 / (i + j + 1) as f64)))
        .collect();
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// The exact ‖A⁻¹‖₁, the largest column sum of the inverse, by solving for each column.
fn exact_inverse_norm_1(n: usize, solve: impl Fn(&[f64]) -> Vec<f64>) -> f64 {
    (0..n)
        .map(|j| {
            let mut unit = vec![0.0; n];
            unit[j] = 1.0;
            solve(&unit).iter().map(|v| v.abs()).sum::<f64>()
        })
        .fold(0.0, f64::max)
}

#[test]
fn solve_transpose_solves_with_the_transpose() {
    for (n, fill, scale, seed) in [(30, 10, 1.5, 1), (80, 4, 0.05, 2), (120, 2, 0.5, 3)] {
        let a = random_square(n, fill, scale, seed).to_csc();
        let lu = a.analyze_lu(Ordering::Amd).unwrap();
        let lu = lu.factor(&a, LuOptions::default()).unwrap();
        let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.7).sin()).collect();
        let x = lu.solve_transpose(&b).unwrap();
        // Aᵀ x = b, checked by the product with the transpose
        let at_x = a.clone().transpose().mul_vec(&x).unwrap();
        let scale = b.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        for (got, want) in at_x.iter().zip(&b) {
            assert!((got - want).abs() <= 1e-10 * scale, "{got} vs {want}");
        }
    }
    // Complex: the plain transpose, not the conjugate
    let a: CscMatrix<Complex<f64>> = random_square(25, 12, 0.3, 4)
        .to_csc()
        .map(|v| Complex::new(v, 0.4 * v - 0.1));
    let lu = a.lu().unwrap();
    let b: Vec<Complex<f64>> = (0..25).map(|i| Complex::new(i as f64, 1.0)).collect();
    let x = lu.solve_transpose(&b).unwrap();
    let at_x = a.clone().transpose().mul_vec(&x).unwrap();
    for (got, want) in at_x.iter().zip(&b) {
        assert!((got - want).norm() <= 1e-10 * 25.0);
    }
    assert_eq!(
        lu.solve_transpose(&b[..24]),
        Err(SparseError::ShapeMismatch {
            left: (25, 25),
            right: (24, 1)
        })
    );
}

#[test]
fn the_condition_estimate_is_a_lower_bound_and_usually_exact() {
    let (mut tried, mut exact) = (0, 0);
    for seed in 0..20 {
        for (n, fill) in [(10, 30), (40, 8), (100, 3)] {
            // Cholesky, on an SPD matrix
            let a = random_spd(n, fill, seed);
            let factor = a.cholesky().unwrap();
            let truth = exact_inverse_norm_1(n, |rhs| factor.solve(rhs).unwrap());
            let estimate = factor.inverse_norm_1_estimate();
            assert!(estimate <= truth * (1.0 + 1e-12), "{estimate} > {truth}");
            assert!(estimate >= truth / 3.0, "{estimate} far below {truth}");
            tried += 1;
            exact += usize::from(estimate >= truth * (1.0 - 1e-12));

            // LU, on an unsymmetric matrix with pivoting
            let a = random_square(n, fill, 0.2, seed + 100).to_csc();
            let lu = a.lu().unwrap();
            let truth = exact_inverse_norm_1(n, |rhs| lu.solve(rhs).unwrap());
            let estimate = lu.inverse_norm_1_estimate();
            assert!(estimate <= truth * (1.0 + 1e-12), "{estimate} > {truth}");
            assert!(estimate >= truth / 3.0, "{estimate} far below {truth}");
            tried += 1;
            exact += usize::from(estimate >= truth * (1.0 - 1e-12));
        }
    }
    eprintln!("the estimate was exact on {exact} of {tried} matrices");
    assert!(exact * 10 >= tried * 8, "exact on only {exact} of {tried}");

    // A diagonal matrix's is exact, and so is the condition number
    let diagonal =
        CsrMatrix::from_triplets(3, 3, vec![(0, 0, 2.0), (1, 1, -0.5), (2, 2, 4.0)]).unwrap();
    let factor = diagonal.to_csc().lu().unwrap();
    assert_eq!(factor.inverse_norm_1_estimate(), 2.0);
    assert_eq!(factor.condition_estimate_1(&diagonal.to_csc()), 8.0);
    // The Hilbert matrix: κ₁(H₆) is about 2.9e7
    let h = hilbert(6);
    let factor = h.cholesky().unwrap();
    let kappa = factor.condition_estimate_1(&h);
    let truth = h.norm_1() * exact_inverse_norm_1(6, |rhs| factor.solve(rhs).unwrap());
    assert!(kappa <= truth * (1.0 + 1e-8) && kappa >= truth / 3.0);
    assert!(kappa > 1e7 && kappa < 1e8, "{kappa}");
}

#[test]
fn refinement_brings_the_backward_error_to_roundoff() {
    let unit_roundoff = f64::EPSILON / 2.0;
    // An unsymmetric matrix whose small diagonal the threshold pivoting keeps, and a badly
    //   conditioned dense one: refinement can only improve the backward error, and here
    //   brings it to the unit roundoff's order
    let a = random_square(150, 3, 0.11, 9).to_csc();
    let lu = a.lu().unwrap();
    // No zero entries: where bᵢ = 0 on a row holding only its diagonal, the exact xᵢ is 0, and
    //   any rounding left in xᵢ makes that row's ratio exactly 1
    let b: Vec<f64> = (0..150).map(|i| 1.0 + ((i * 37) % 11) as f64).collect();
    let first = lu.solve_refined(&a, &b, 0).unwrap();
    let refined = lu.solve_refined(&a, &b, 10).unwrap();
    assert_eq!(first.steps, 0);
    assert_eq!(first.x, lu.solve(&b).unwrap());
    assert!(refined.backward_error <= first.backward_error);
    assert!(
        refined.backward_error <= 4.0 * unit_roundoff,
        "{:e}",
        refined.backward_error
    );
    eprintln!(
        "LU: componentwise backward error {:e} → {:e} in {} steps",
        first.backward_error, refined.backward_error, refined.steps
    );

    let h = hilbert(10);
    let factor = h.cholesky().unwrap();
    let b: Vec<f64> = (0..10).map(|i| 1.0 + i as f64).collect();
    let first = factor.solve_refined(&h, &b, 0).unwrap();
    let refined = factor.solve_refined(&h, &b, 10).unwrap();
    assert!(refined.backward_error <= first.backward_error);
    assert!(
        refined.backward_error <= 4.0 * unit_roundoff,
        "{:e}",
        refined.backward_error
    );
    eprintln!(
        "Hilbert(10), κ₁ ≈ {:.1e}: componentwise backward error {:e} → {:e} in {} steps",
        factor.condition_estimate_1(&h),
        first.backward_error,
        refined.backward_error,
        refined.steps
    );

    // Shapes are checked
    assert!(matches!(
        factor.solve_refined(&h, &b[..9], 3),
        Err(SparseError::ShapeMismatch { .. })
    ));
    assert!(matches!(
        lu.solve_refined(&h.to_csc(), &b, 3),
        Err(SparseError::ShapeMismatch { .. })
    ));
}

#[test]
fn an_exact_solve_needs_no_refinement() {
    // LU solves a diagonal system exactly (one division per entry), so the backward error
    //   starts at zero. (Cholesky wouldn't: its factor holds square roots.)
    let a = CscMatrix::from_triplets(3, 3, vec![(0, 0, 2.0), (1, 1, 4.0), (2, 2, 8.0)]).unwrap();
    let refined = a
        .lu()
        .unwrap()
        .solve_refined(&a, &[2.0, 4.0, 8.0], 5)
        .unwrap();
    assert_eq!(
        (refined.x, refined.steps, refined.backward_error),
        (vec![1.0; 3], 0, 0.0)
    );
}
