// Tests for the static-pivot LDLᵀ factorization (SymbolicCholesky::factor_ldlt, CsrMatrix::ldlt):
//   against the Cholesky on positive definite matrices, against a known D on matrices built as
//   L₀ D₀ L₀ᵀ, on quasi-definite KKT-shaped matrices under every ordering, its inertia, a
//   complex Hermitian matrix, and its errors; then the tools for solving with it on regularized
//   matrices: refinement against the unregularized matrix (componentwise and normwise), pivot
//   replacement (factor_ldlt_replacing) and symmetric scaling.

use sparse_mat::compressed_matrix::{CsrMatrix, Inertia, Ordering, PivotReplacement};
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

    /// A value in [-1, 1]
    fn unit(&mut self) -> f64 {
        (self.next() % 2001) as f64 / 1000.0 - 1.0
    }
}

/// A pseudo-random SPD matrix with a dominant diagonal, as accuracy_test.rs builds them.
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

/// The quasi-definite KKT-shaped matrix `[H, Bᵀ; B, −δ I]`: `H` random SPD (n × n), `B` random
/// (m × n), and `δ > 0`
fn quasi_definite(n: usize, m: usize, delta: f64, seed: u64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed ^ 0x9e37);
    let mut triplets: Vec<(usize, usize, f64)> = random_spd(n, 15, seed).iter().collect();
    for r in 0..m {
        // Each constraint row touches a few variables
        for _ in 0..3 {
            let c = (rng.next() as usize) % n;
            let v = rng.unit() + 0.1;
            triplets.extend([(n + r, c, v), (c, n + r, v)]);
        }
        triplets.push((n + r, n + r, -delta));
    }
    CsrMatrix::from_triplets(n + m, n + m, triplets).unwrap()
}

/// `‖A x − b‖∞ / (‖A‖∞ ‖x‖∞ + ‖b‖∞)`, a scaled residual
fn residual(a: &CsrMatrix<f64>, x: &[f64], b: &[f64]) -> f64 {
    let ax = a.mul_vec(x).unwrap();
    let r = ax
        .iter()
        .zip(b)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
    let norm_a = (0..a.shape().0)
        .map(|i| a.outer(i).unwrap().1.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0f64, f64::max);
    let norm = |v: &[f64]| v.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    r / (norm_a * norm(x) + norm(b))
}

/// The scaled residual allowed for a quasi-definite matrix regularized by `δ`. Static pivots
/// lose accuracy as the regularization shrinks: measured over 20 seeds and three orderings of
/// `quasi_definite(50, 20, δ, ·)`, the worst residual is 3.1e-10 at δ = 1e-8, 3.1e-12 at 1e-6,
/// 2.9e-14 at 1e-4, and near machine precision from 1e-2 up (2026-11-07), growing as 1/δ
fn residual_bar(delta: f64) -> f64 {
    1e-14 * (1e-3 / delta).max(1.0)
}

fn right_hand_side(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed);
    (0..n).map(|_| rng.unit()).collect()
}

#[test]
fn on_positive_definite_matrices_it_is_the_cholesky_scaled() {
    for seed in 1..6 {
        let a = random_spd(40, 10, seed);
        let symbolic = a.analyze_cholesky(Ordering::Amd).unwrap();
        let cholesky = symbolic.factor(&a).unwrap();
        let ldlt = symbolic.factor_ldlt(&a).unwrap();
        assert_eq!(ldlt.permutation(), cholesky.permutation());
        assert_eq!(
            ldlt.inertia(),
            Inertia {
                positive: 40,
                negative: 0
            }
        );
        // The Cholesky's L is L sqrt(D), column by column
        let (lc, ll) = (cholesky.l(), ldlt.l());
        assert_eq!(lc.pointers(), ll.pointers());
        assert_eq!(lc.indices(), ll.indices());
        for i in 0..40 {
            let (cols, values) = ll.outer(i).unwrap();
            let (_, chol) = lc.outer(i).unwrap();
            for ((&j, &l), &c) in cols.iter().zip(values).zip(chol) {
                let scaled = l * ldlt.d()[j].sqrt();
                assert!(
                    (scaled - c).abs() <= 1e-12 * c.abs().max(1.0),
                    "seed {seed}: L({i}, {j}) √d = {scaled}, Cholesky {c}"
                );
            }
        }
        let b = right_hand_side(40, seed);
        let (xl, xc) = (ldlt.solve(&b).unwrap(), cholesky.solve(&b).unwrap());
        for (l, c) in xl.iter().zip(&xc) {
            assert!((l - c).abs() <= 1e-12 * c.abs().max(1.0), "seed {seed}");
        }
    }
}

#[test]
fn a_known_factorization_is_recovered_with_its_inertia() {
    // A = L₀ D₀ L₀ᵀ with a sparse unit lower L₀ and D₀ of mixed signs: under the natural
    //   ordering, the static-pivot LDLᵀ is unique, so it must give back D₀ and L₀
    for seed in 1..6 {
        let n = 30;
        let mut rng = Lcg(seed);
        let mut l0 = vec![vec![0.0; n]; n];
        let mut d0 = vec![0.0; n];
        for i in 0..n {
            l0[i][i] = 1.0;
            for row in l0[i].iter_mut().take(i) {
                if rng.next() % 100 < 15 {
                    *row = rng.unit();
                }
            }
            let magnitude = 0.5 + (rng.next() % 1500) as f64 / 1000.0;
            d0[i] = if rng.next().is_multiple_of(3) {
                -magnitude
            } else {
                magnitude
            };
        }
        // The lower triangle, mirrored, so the matrix is exactly symmetric
        let mut triplets = vec![];
        for i in 0..n {
            for j in 0..=i {
                let v: f64 = (0..=j).map(|k| l0[i][k] * d0[k] * l0[j][k]).sum();
                if v != 0.0 {
                    triplets.push((i, j, v));
                    if i != j {
                        triplets.push((j, i, v));
                    }
                }
            }
        }
        let a = CsrMatrix::from_triplets(n, n, triplets).unwrap();
        let ldlt = a
            .analyze_cholesky(Ordering::Natural)
            .unwrap()
            .factor_ldlt(&a)
            .unwrap();
        for (i, (&d, &expected)) in ldlt.d().iter().zip(&d0).enumerate() {
            assert!(
                (d - expected).abs() <= 1e-10,
                "seed {seed}: d[{i}] = {d}, {expected}"
            );
        }
        let negative = d0.iter().filter(|&&d| d < 0.0).count();
        assert_eq!(
            ldlt.inertia(),
            Inertia {
                positive: n - negative,
                negative
            },
            "seed {seed}"
        );
        let b = right_hand_side(n, seed);
        let x = ldlt.solve(&b).unwrap();
        assert!(residual(&a, &x, &b) < 1e-13, "seed {seed}");
    }
}

#[test]
fn quasi_definite_matrices_factor_under_every_ordering() {
    // Vanderbei: [H, Bᵀ; B, −G] with H and G positive definite factors under any symmetric
    //   ordering, with D positive on H's rows and negative on G's
    for seed in 1..6 {
        let (n, m) = (50, 20);
        for delta in [1e-8, 1e-4, 1.0] {
            let a = quasi_definite(n, m, delta, seed);
            for ordering in [Ordering::Natural, Ordering::Rcm, Ordering::Amd] {
                let ldlt = a
                    .analyze_cholesky(ordering)
                    .unwrap()
                    .factor_ldlt(&a)
                    .unwrap_or_else(|e| panic!("seed {seed}, δ {delta}, {ordering:?}: {e}"));
                assert_eq!(
                    ldlt.inertia(),
                    Inertia {
                        positive: n,
                        negative: m
                    },
                    "seed {seed}, δ {delta}, {ordering:?}"
                );
                let b = right_hand_side(n + m, seed);
                let x = ldlt.solve(&b).unwrap();
                let r = residual(&a, &x, &b);
                assert!(
                    r < residual_bar(delta),
                    "seed {seed}, δ {delta}, {ordering:?}: residual {r:e}"
                );
            }
        }
    }
}

#[test]
fn a_matrix_that_needs_pivoting_fails_until_regularized() {
    let swap = CsrMatrix::from_triplets(2, 2, vec![(0, 1, 1.0), (1, 0, 1.0)]).unwrap();
    let symbolic = swap.analyze_cholesky(Ordering::Natural).unwrap();
    assert_eq!(
        symbolic.factor_ldlt(&swap).unwrap_err(),
        SparseError::SingularMatrix { pivot: 0 }
    );
    let delta = 1e-6;
    let regularized = CsrMatrix::from_triplets(
        2,
        2,
        vec![(0, 0, delta), (0, 1, 1.0), (1, 0, 1.0), (1, 1, -delta)],
    )
    .unwrap();
    let ldlt = regularized.ldlt().unwrap();
    assert_eq!(
        ldlt.inertia(),
        Inertia {
            positive: 1,
            negative: 1
        }
    );
    let x = ldlt.solve(&[1.0, 2.0]).unwrap();
    // With the tiny pivot first, L's entry is 1/δ, and the residual grows accordingly
    assert!(residual(&regularized, &x, &[1.0, 2.0]) < residual_bar(delta));
}

#[test]
fn errors_match_the_choleskys() {
    let a = random_spd(10, 20, 3);
    let symbolic = a.analyze_cholesky(Ordering::Amd).unwrap();
    let other = random_spd(10, 40, 4);
    assert_eq!(
        symbolic.factor_ldlt(&other).unwrap_err(),
        SparseError::PatternMismatch
    );
    let unsymmetric = CsrMatrix::from_triplets(
        2,
        2,
        vec![(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 1.0)],
    )
    .unwrap();
    assert_eq!(unsymmetric.ldlt().unwrap_err(), SparseError::NotSymmetric);
    let ldlt = a.ldlt().unwrap();
    assert!(matches!(
        ldlt.solve(&[1.0; 9]),
        Err(SparseError::ShapeMismatch { .. })
    ));
}

#[test]
fn a_hermitian_indefinite_matrix_factors() {
    // [[2, 1 + i], [1 − i, −3]]: Hermitian, with one eigenvalue of each sign
    let c = |re: f64, im: f64| Complex::new(re, im);
    let a = CsrMatrix::from_triplets(
        2,
        2,
        vec![
            (0, 0, c(2.0, 0.0)),
            (0, 1, c(1.0, 1.0)),
            (1, 0, c(1.0, -1.0)),
            (1, 1, c(-3.0, 0.0)),
        ],
    )
    .unwrap();
    let ldlt = a
        .analyze_cholesky(Ordering::Natural)
        .unwrap()
        .factor_ldlt(&a)
        .unwrap();
    assert_eq!(
        ldlt.inertia(),
        Inertia {
            positive: 1,
            negative: 1
        }
    );
    // d₂ = −3 − |1 + i|² / 2 = −4
    assert!((ldlt.d()[1] - c(-4.0, 0.0)).norm() < 1e-14);
    let b = [c(1.0, 2.0), c(-1.0, 0.5)];
    let x = ldlt.solve(&b).unwrap();
    let ax = a.mul_vec(&x).unwrap();
    for (ax, b) in ax.iter().zip(&b) {
        assert!((ax - b).norm() < 1e-14);
    }
}

#[test]
fn solve_many_is_solve_bit_for_bit() {
    let a = quasi_definite(30, 10, 1e-6, 7);
    let ldlt = a.ldlt().unwrap();
    let rhs: Vec<Vec<f64>> = (0..5).map(|s| right_hand_side(40, s)).collect();
    let many = ldlt.solve_many(&rhs).unwrap();
    for (b, x) in rhs.iter().zip(&many) {
        assert_eq!(&ldlt.solve(b).unwrap(), x);
        let mut in_place = b.clone();
        ldlt.solve_in_place(&mut in_place).unwrap();
        assert_eq!(&in_place, x);
    }
    #[cfg(feature = "rayon")]
    assert_eq!(ldlt.par_solve_many(&rhs).unwrap(), many);
}

#[test]
fn a_larger_hermitian_quasi_definite_matrix_factors() {
    // [H, Bᴴ; B, −δI] with complex off-diagonal entries, so the updates multiply complex values
    //   (the 2 × 2 case above has none)
    let (n, m, delta) = (24, 8, 1e-2);
    let mut rng = Lcg(11);
    let mut triplets: Vec<(usize, usize, Complex<f64>)> = vec![];
    let mut row_sums = vec![0.0; n];
    for i in 0..n {
        for j in 0..i {
            if rng.next() % 100 < 25 {
                let v = Complex::new(rng.unit(), rng.unit());
                triplets.extend([(i, j, v), (j, i, v.conj())]);
                row_sums[i] += v.norm();
                row_sums[j] += v.norm();
            }
        }
    }
    for (i, sum) in row_sums.iter().enumerate() {
        triplets.push((i, i, Complex::new(sum + 1.0, 0.0)));
    }
    for r in 0..m {
        for _ in 0..3 {
            let c = (rng.next() as usize) % n;
            let v = Complex::new(rng.unit() + 0.1, rng.unit());
            triplets.extend([(n + r, c, v), (c, n + r, v.conj())]);
        }
        triplets.push((n + r, n + r, Complex::new(-delta, 0.0)));
    }
    let a = CsrMatrix::from_triplets(n + m, n + m, triplets).unwrap();
    for ordering in [Ordering::Natural, Ordering::Amd] {
        let ldlt = a
            .analyze_cholesky(ordering)
            .unwrap()
            .factor_ldlt(&a)
            .unwrap();
        assert_eq!(
            ldlt.inertia(),
            Inertia {
                positive: n,
                negative: m
            },
            "{ordering:?}"
        );
        let b: Vec<Complex<f64>> = (0..n + m)
            .map(|_| Complex::new(rng.unit(), rng.unit()))
            .collect();
        let x = ldlt.solve(&b).unwrap();
        let ax = a.mul_vec(&x).unwrap();
        let worst = ax
            .iter()
            .zip(&b)
            .fold(0.0f64, |w, (ax, b)| w.max((ax - b).norm()));
        assert!(worst < 1e-12, "{ordering:?}: residual {worst:e}");
    }
}

#[test]
fn refinement_against_the_unregularized_matrix_recovers_its_solution() {
    // A KKT-shaped matrix with a zero lower-right block, which static pivots can't factor under
    //   every ordering, factored regularized and refined against itself
    let unit_roundoff = f64::EPSILON / 2.0;
    let (n, m) = (50, 20);
    let exact = quasi_definite(n, m, 0.0, 41);
    let b = right_hand_side(n + m, 42);
    let reference = exact.to_csc().lu().unwrap().solve(&b).unwrap();
    let distance = |x: &[f64]| {
        x.iter()
            .zip(&reference)
            .fold(0.0f64, |d, (x, r)| d.max((x - r).abs()))
            / reference.iter().fold(0.0f64, |d, r| d.max(r.abs()))
    };
    for delta in [1e-8, 1e-6, 1e-4] {
        let regularized = quasi_definite(n, m, delta, 41);
        let ldlt = regularized.ldlt().unwrap();
        let first = ldlt.solve_refined(&exact, &b, 0).unwrap();
        let refined = ldlt.solve_refined(&exact, &b, 10).unwrap();
        eprintln!(
            "δ = {delta:e}: backward error {:e} → {:e} in {} steps; distance from the exact \
             solution {:e} → {:e}",
            first.backward_error,
            refined.backward_error,
            refined.steps,
            distance(&first.x),
            distance(&refined.x)
        );
        assert_eq!(first.x, ldlt.solve(&b).unwrap());
        assert!(refined.steps >= 1);
        assert!(
            refined.backward_error <= 4.0 * unit_roundoff,
            "{:e}",
            refined.backward_error
        );
        assert!(distance(&refined.x) < 1e-12, "{:e}", distance(&refined.x));
        // Refining against the regularized matrix instead converges to its own solution,
        //   which differs from the exact one by about δ
        let own = ldlt.solve_refined(&regularized, &b, 10).unwrap();
        assert!(distance(&own.x) > 0.1 * delta, "{:e}", distance(&own.x));
    }
    // Shapes are checked
    let ldlt = quasi_definite(n, m, 1e-6, 41).ldlt().unwrap();
    assert!(matches!(
        ldlt.solve_refined(&exact, &b[..n], 3),
        Err(SparseError::ShapeMismatch { .. })
    ));
}

/// The relative residual `‖b − A x‖∞ / ‖b‖∞`
fn relative_residual(a: &CsrMatrix<f64>, x: &[f64], b: &[f64]) -> f64 {
    let ax = a.mul_vec(x).unwrap();
    let norm = |v: &mut dyn Iterator<Item = f64>| v.fold(0.0f64, |m, v| m.max(v.abs()));
    norm(&mut ax.iter().zip(b).map(|(a, b)| a - b)) / norm(&mut b.iter().copied())
}

#[test]
fn normwise_refinement_reaches_its_tolerance_against_the_unregularized_matrix() {
    let (n, m) = (50, 20);
    let exact = quasi_definite(n, m, 0.0, 41);
    let b = right_hand_side(n + m, 42);
    for delta in [1e-8, 1e-6, 1e-4] {
        let ldlt = quasi_definite(n, m, delta, 41).ldlt().unwrap();
        let first = ldlt.solve_refined_normwise(&exact, &b, 0, 0.0).unwrap();
        assert_eq!(first.x, ldlt.solve(&b).unwrap());
        assert_eq!(
            first.relative_residual,
            relative_residual(&exact, &first.x, &b)
        );
        let refined = ldlt.solve_refined_normwise(&exact, &b, 20, 1e-14).unwrap();
        eprintln!(
            "δ = {delta:e}: relative residual {:e} → {:e} in {} steps",
            first.relative_residual, refined.relative_residual, refined.steps
        );
        assert!(refined.steps >= 1);
        assert!(
            refined.relative_residual <= 1e-14,
            "{:e}",
            refined.relative_residual
        );
        assert_eq!(
            refined.relative_residual,
            relative_residual(&exact, &refined.x, &b)
        );
        // A tolerance already met takes no step
        let loose = ldlt
            .solve_refined_normwise(&exact, &b, 20, first.relative_residual)
            .unwrap();
        assert_eq!((loose.steps, loose.x), (0, first.x));
    }
    // Refinement against a matrix too far from the factored one diverges, and the best iterate,
    //   the first, is returned: factored [0.25] against [1], each step's error is 3 times the
    //   last (x = 4, then −8)
    let far = CsrMatrix::from_triplets(1, 1, vec![(0, 0, 0.25)])
        .unwrap()
        .ldlt()
        .unwrap();
    let one = CsrMatrix::from_triplets(1, 1, vec![(0, 0, 1.0)]).unwrap();
    let diverged = far.solve_refined_normwise(&one, &[1.0], 5, 0.0).unwrap();
    assert_eq!(
        (diverged.x, diverged.steps, diverged.relative_residual),
        (vec![4.0], 0, 3.0)
    );
    let ldlt = quasi_definite(n, m, 1e-6, 41).ldlt().unwrap();
    assert!(matches!(
        ldlt.solve_refined_normwise(&exact, &b[..n], 3, 1e-14),
        Err(SparseError::ShapeMismatch { .. })
    ));
}

#[test]
fn replaced_pivots_take_the_expected_sign_and_refinement_recovers_the_solution() {
    // [[0, 1], [1, 0]] under the natural ordering: the first pivot is exactly zero, and is
    //   replaced by +1e-6 (row 0 expected positive); the second is then 0 − 1/1e-6 = −1e6
    let swap: CsrMatrix<f64> =
        CsrMatrix::from_triplets(2, 2, vec![(0, 1, 1.0), (1, 0, 1.0)]).unwrap();
    let symbolic = swap.analyze_cholesky(Ordering::Natural).unwrap();
    let replacement = PivotReplacement {
        positive: &[true, false],
        threshold: 1e-13,
        value: 1e-6,
    };
    let ldlt = symbolic.factor_ldlt_replacing(&swap, &replacement).unwrap();
    assert_eq!(ldlt.replaced(), 1);
    assert_eq!(ldlt.d(), &[1e-6, -1e6]);
    assert_eq!(
        ldlt.inertia(),
        Inertia {
            positive: 1,
            negative: 1
        }
    );
    let refined = ldlt
        .solve_refined_normwise(&swap, &[1.0, 2.0], 10, 1e-15)
        .unwrap();
    assert!(
        refined.relative_residual <= 1e-15,
        "{:e}",
        refined.relative_residual
    );
    assert!((refined.x[0] - 2.0).abs() < 1e-14 && (refined.x[1] - 1.0).abs() < 1e-14);

    // A row expected negative takes −value
    let flipped = PivotReplacement {
        positive: &[false, true],
        ..replacement
    };
    let ldlt = symbolic.factor_ldlt_replacing(&swap, &flipped).unwrap();
    assert_eq!(ldlt.d(), &[-1e-6, 1e6]);

    // Without replacement, factor_ldlt is unchanged, and counts none
    let plain = quasi_definite(30, 10, 1e-6, 5);
    let symbolic = plain.analyze_cholesky(Ordering::Amd).unwrap();
    let positive: Vec<bool> = (0..40).map(|i| i < 30).collect();
    let none_small = PivotReplacement {
        positive: &positive,
        threshold: 1e-13,
        value: 1e-6,
    };
    let (a, b) = (
        symbolic.factor_ldlt(&plain).unwrap(),
        symbolic.factor_ldlt_replacing(&plain, &none_small).unwrap(),
    );
    assert_eq!((a.replaced(), b.replaced()), (0, 0));
    assert_eq!(a.d(), b.d());
    assert_eq!(a.l(), b.l());

    // A large pivot of the unexpected sign is left alone: here every row is expected positive,
    //   and the inertia still shows the negative block
    let all_positive = vec![true; 40];
    let unexpected = PivotReplacement {
        positive: &all_positive,
        ..none_small
    };
    let ldlt = symbolic.factor_ldlt_replacing(&plain, &unexpected).unwrap();
    assert_eq!(ldlt.replaced(), 0);
    assert_eq!(
        ldlt.inertia(),
        Inertia {
            positive: 30,
            negative: 10
        }
    );

    // Errors: one sign per row, and a NaN still fails
    assert!(matches!(
        symbolic.factor_ldlt_replacing(
            &plain,
            &PivotReplacement {
                positive: &[true],
                ..none_small
            }
        ),
        Err(SparseError::ShapeMismatch { .. })
    ));
    let nan = CsrMatrix::from_triplets(1, 1, vec![(0, 0, f64::NAN)]).unwrap();
    let nan_symbolic = nan.analyze_cholesky(Ordering::Natural);
    if let Ok(nan_symbolic) = nan_symbolic {
        assert!(matches!(
            nan_symbolic.factor_ldlt_replacing(
                &nan,
                &PivotReplacement {
                    positive: &[true],
                    ..none_small
                }
            ),
            Err(SparseError::SingularMatrix { pivot: 0 } | SparseError::NotSymmetric)
        ));
    }
}

#[test]
fn symmetric_scaling_equilibrates_and_keeps_the_solution() {
    // A quasi-definite matrix with rows and columns scaled by powers of ten from 1e-6 to 1e6
    let (n, m) = (40, 15);
    let base = quasi_definite(n, m, 1e-6, 17);
    let mut rng = Lcg(23);
    let skew: Vec<f64> = (0..n + m)
        .map(|_| 10f64.powi((rng.next() % 13) as i32 - 6))
        .collect();
    let skewed = CsrMatrix::from_triplets(
        n + m,
        n + m,
        base.iter().map(|(r, c, v)| (r, c, v * (skew[r] * skew[c]))),
    )
    .unwrap();
    let largest_per_row = |a: &CsrMatrix<f64>| -> Vec<f64> {
        (0..a.shape().0)
            .map(|r| {
                a.outer(r)
                    .unwrap()
                    .1
                    .iter()
                    .fold(0.0f64, |m, v| m.max(v.abs()))
            })
            .collect()
    };
    let spread = |v: &[f64]| {
        v.iter().fold(0.0f64, |m, v| m.max(*v)) / v.iter().fold(f64::MAX, |m, v| m.min(*v))
    };
    let s = skewed.symmetric_scaling(10).unwrap();
    assert!(s.iter().all(|s| *s > 0.0));
    let scaled = CsrMatrix::from_triplets(
        n + m,
        n + m,
        skewed.iter().map(|(r, c, v)| (r, c, v * (s[r] * s[c]))),
    )
    .unwrap();
    assert!(scaled.is_hermitian());
    let before = spread(&largest_per_row(&skewed));
    let after = largest_per_row(&scaled);
    eprintln!(
        "rows' largest magnitudes: spread {before:e} before, {:e} after, in [{:e}, {:e}]",
        spread(&after),
        after.iter().fold(f64::MAX, |m, v| m.min(*v)),
        after.iter().fold(0.0f64, |m, v| m.max(*v))
    );
    assert!(before > 1e6);
    // Each sweep about halves the logarithm's distance from 1: after 10, measured, the rows'
    //   largest magnitudes are in [0.977, 1] (2026-11-07)
    assert!(
        after.iter().all(|v| (0.95..=1.0 + 1e-12).contains(v)),
        "{after:?}"
    );

    // Solving the scaled system and unscaling gives the original's solution
    let b = right_hand_side(n + m, 29);
    let sb: Vec<f64> = b.iter().zip(&s).map(|(b, s)| b * s).collect();
    let y = scaled.ldlt().unwrap().solve(&sb).unwrap();
    let x: Vec<f64> = y.iter().zip(&s).map(|(y, s)| y * s).collect();
    let reference = skewed.to_csc().lu().unwrap().solve(&b).unwrap();
    let worst = x.iter().zip(&reference).fold(0.0f64, |w, (x, r)| {
        w.max((x - r).abs() / r.abs().max(1e-300))
    });
    assert!(worst < 1e-6, "{worst:e}");

    // Zero sweeps is the identity; errors as for the factorization
    assert_eq!(skewed.symmetric_scaling(0).unwrap(), vec![1.0; n + m]);
    let rectangular = CsrMatrix::from_triplets(2, 3, vec![(0, 0, 1.0)]).unwrap();
    assert!(matches!(
        rectangular.symmetric_scaling(3),
        Err(SparseError::NotSquare { .. })
    ));
    let asymmetric = CsrMatrix::from_triplets(2, 2, vec![(0, 1, 1.0), (1, 0, 2.0)]).unwrap();
    assert_eq!(
        asymmetric.symmetric_scaling(3).unwrap_err(),
        SparseError::NotSymmetric
    );
}
