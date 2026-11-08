// Tests for the static-pivot LDLᵀ factorization (SymbolicCholesky::factor_ldlt, CsrMatrix::ldlt):
//   against the Cholesky on positive definite matrices, against a known D on matrices built as
//   L₀ D₀ L₀ᵀ, on quasi-definite KKT-shaped matrices under every ordering, its inertia, a
//   complex Hermitian matrix, and its errors.

use sparse_mat::compressed_matrix::{CsrMatrix, Inertia, Ordering};
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
