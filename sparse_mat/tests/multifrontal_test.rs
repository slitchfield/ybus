// Tests for the multifrontal pivoted LDLᵀ (SymbolicCholesky::factor_indefinite,
//   CsrMatrix::ldlt_indefinite): against the dense kernel and the LU on random symmetric
//   indefinite matrices, on KKT-shaped matrices whose zero diagonals static pivots fail on, with
//   delayed pivots, against the Cholesky on positive definite matrices, complex Hermitian, its
//   refinement, and its errors.

use sparse_mat::compressed_matrix::{
    Amalgamation, CsrMatrix, DenseLdlt, IndefiniteOptions, Inertia, Ordering, SymbolicIndefinite,
};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::Complex;

type C = Complex<f64>;

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

/// A random sparse symmetric matrix with entries in [-1, 1], about `fill_percent`% of the
/// off-diagonal positions stored, and every diagonal stored (zero in about a third of them, each
/// such row also joined to the next row, so that few draws are structurally singular)
fn random_symmetric(n: usize, fill_percent: u64, seed: u64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for i in 0..n {
        for j in 0..i {
            if rng.next() % 100 < fill_percent {
                let v = rng.unit();
                triplets.extend([(i, j, v), (j, i, v)]);
            }
        }
        let zero = rng.next().is_multiple_of(3) && n > 1;
        let d = if zero { 0.0 } else { rng.unit() * 3.0 };
        triplets.push((i, i, d));
        if zero {
            let next = (i + 1) % n;
            let v = 0.5 + rng.unit().abs();
            triplets.extend([(i, next, v), (next, i, v)]);
        }
    }
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// `[W Bᵀ; B 0]`: `W` (n × n, n ≥ 2m) diagonal, zero at the even rows below 2m and in [1, 1001]
/// elsewhere, and `B` (m × n) sparse, row r touching variable 2r and three others, every
/// constraint row's zero diagonal stored. Each zero-diagonal variable has its own row of `B`, so
/// `W` is positive definite on `B`'s null space (for generic values) and the inertia is (n, m).
fn kkt(n: usize, m: usize, seed: u64) -> CsrMatrix<f64> {
    assert!(n >= 2 * m);
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for i in 0..n {
        let d = if i.is_multiple_of(2) && i < 2 * m {
            0.0
        } else {
            1.0 + rng.unit().abs() * 1e3
        };
        triplets.push((i, i, d));
    }
    for r in 0..m {
        // Each zero-diagonal variable is touched by some row, and each row touches a few more
        let mut cols = vec![2 * r];
        for _ in 0..3 {
            cols.push((rng.next() as usize) % n);
        }
        for c in cols {
            let v = rng.unit() + 2.0;
            triplets.extend([(n + r, c, v), (c, n + r, v)]);
        }
        triplets.push((n + r, n + r, 0.0));
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

fn right_hand_side(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed);
    (0..n).map(|_| rng.unit()).collect()
}

/// The scaled residual allowed: a backward-stable solve stays within a modest multiple of n·ε,
/// L's entries being up to 1/u = 100. The worst measured over these tests' matrices was 7.1e-15,
/// on the random ones (2026-11-07)
const RESIDUAL_BAR: f64 = 1e-13;

#[test]
fn random_symmetric_indefinite_matrices_match_the_dense_factorization() {
    let mut worst = 0.0f64;
    let (mut compared, mut skipped) = (0, 0);
    for seed in 0..30 {
        for (n, fill) in [(1, 0), (6, 30), (25, 10), (60, 5)] {
            let a = random_symmetric(n, fill, seed * 13 + n as u64);
            let dense = a.to_dense();
            // A draw singular to working precision has no well-defined inertia (the sign of a
            //   pivot of 1e-17 is rounding's): only well-conditioned ones are compared
            let Ok(reference) = DenseLdlt::factor(&dense, 0.01) else {
                skipped += 1;
                continue;
            };
            if inverse_norm_1(&reference, n) > 1e8 {
                skipped += 1;
                continue;
            }
            compared += 1;
            for ordering in [Ordering::Natural, Ordering::Amd, Ordering::Rcm] {
                let f = a
                    .analyze_cholesky(ordering)
                    .unwrap()
                    .factor_indefinite(&a, IndefiniteOptions::default())
                    .unwrap();
                assert_eq!(f.inertia(), reference.inertia(), "n {n}, seed {seed}");
                let b = right_hand_side(n, seed);
                let x = f.solve(&b).unwrap();
                worst = worst.max(residual(&a, &x, &b));
            }
        }
    }
    eprintln!("worst scaled residual {worst:e}; {compared} draws compared, {skipped} too badly conditioned");
    assert!(worst < RESIDUAL_BAR, "{worst:e}");
    assert!(
        skipped * 10 < compared,
        "{skipped} skipped of {}",
        compared + skipped
    );
}

/// `‖A⁻¹‖₁` from the dense factorization, column by column
fn inverse_norm_1(f: &DenseLdlt<f64>, n: usize) -> f64 {
    (0..n)
        .map(|j| {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            f.solve(&e).unwrap().iter().map(|v| v.abs()).sum::<f64>()
        })
        .fold(0.0, f64::max)
}

#[test]
fn kkt_matrices_with_zero_diagonals_factor_under_every_ordering() {
    let mut worst = 0.0f64;
    for seed in 0..10 {
        let (n, m) = (60, 20);
        let a = kkt(n, m, seed);
        // Static pivots meet a zero (the natural order starts on one)
        assert!(a
            .analyze_cholesky(Ordering::Natural)
            .unwrap()
            .factor_ldlt(&a)
            .is_err());
        for ordering in [Ordering::Natural, Ordering::Amd, Ordering::Rcm] {
            let f = a
                .analyze_cholesky(ordering)
                .unwrap()
                .factor_indefinite(&a, IndefiniteOptions::default())
                .unwrap();
            assert_eq!(
                f.inertia(),
                Inertia {
                    positive: n,
                    negative: m
                },
                "seed {seed}, {ordering:?}"
            );
            assert!(f.two_by_two() > 0);
            let b = right_hand_side(n + m, seed);
            let x = f.solve(&b).unwrap();
            worst = worst.max(residual(&a, &x, &b));
        }
    }
    eprintln!("KKT: worst scaled residual {worst:e}");
    assert!(worst < RESIDUAL_BAR, "{worst:e}");
}

#[test]
fn small_leaf_pivots_are_delayed_to_the_root() {
    // A star: k leaves, each with a diagonal of 1e-6, joined only to the centre. Each leaf is its
    //   own front, whose one fully summed pivot fails the test against its entry of 1, so it's
    //   delayed; the root pivots them in 2×2 blocks with the centre or with each other. That
    //   needs the fundamental supernodes (no amalgamation), which keep each leaf apart
    let k = 8;
    let mut triplets = vec![(k, k, 1.0)];
    for leaf in 0..k {
        triplets.extend([(leaf, leaf, 1e-6), (leaf, k, 1.0), (k, leaf, 1.0)]);
    }
    let a = CsrMatrix::from_triplets(k + 1, k + 1, triplets).unwrap();
    let unmerged = SymbolicIndefinite::with_amalgamation(
        a.analyze_cholesky(Ordering::Natural).unwrap(),
        Amalgamation::NONE,
    );
    let f = unmerged.factor(&a, IndefiniteOptions::default()).unwrap();
    assert_eq!(f.delayed(), k);
    let reference = DenseLdlt::factor(&a.to_dense(), 0.01).unwrap();
    assert_eq!(f.inertia(), reference.inertia());
    let b = right_hand_side(k + 1, 1);
    let x = f.solve(&b).unwrap();
    assert!(
        residual(&a, &x, &b) < RESIDUAL_BAR,
        "{:e}",
        residual(&a, &x, &b)
    );
    // The delays enlarge the root: more of L is stored than the Cholesky pattern's 8 entries
    assert!(f.l_entries() > k, "{}", f.l_entries());

    // A threshold small enough lets the leaves pivot where they are
    let loose = unmerged
        .factor(&a, IndefiniteOptions { threshold: 1e-7 })
        .unwrap();
    assert_eq!((loose.delayed(), loose.l_entries()), (0, k));

    // With the default amalgamation the leaves and the centre are one front (9 columns, under
    //   16, with fewer than 80% zeros), where the pivots are found without any delay
    let merged = a
        .analyze_cholesky(Ordering::Natural)
        .unwrap()
        .factor_indefinite(&a, IndefiniteOptions::default())
        .unwrap();
    assert_eq!(merged.delayed(), 0);
    assert_eq!(merged.inertia(), reference.inertia());
    let x = merged.solve(&b).unwrap();
    assert!(residual(&a, &x, &b) < RESIDUAL_BAR);
}

#[test]
fn positive_definite_matrices_factor_as_the_cholesky_does() {
    // A diagonally dominant SPD matrix: every 1×1 pivot passes, nothing is delayed, and without
    //   amalgamation L has the Cholesky pattern
    let mut rng = Lcg(77);
    let n = 80;
    let mut triplets = vec![];
    let mut sums = vec![0.0; n];
    for i in 0..n {
        for j in 0..i {
            if rng.next() % 100 < 6 {
                let v = rng.unit();
                triplets.extend([(i, j, v), (j, i, v)]);
                sums[i] += v.abs();
                sums[j] += v.abs();
            }
        }
    }
    for (i, s) in sums.iter().enumerate() {
        triplets.push((i, i, s + 1.0));
    }
    let a = CsrMatrix::from_triplets(n, n, triplets).unwrap();
    let symbolic = a.analyze_cholesky(Ordering::Amd).unwrap();
    let f = SymbolicIndefinite::with_amalgamation(symbolic.clone(), Amalgamation::NONE)
        .factor(&a, IndefiniteOptions::default())
        .unwrap();
    let cholesky = symbolic.factor(&a).unwrap();
    assert_eq!((f.delayed(), f.two_by_two()), (0, 0));
    assert_eq!(f.l_entries(), cholesky.l().values().len() - n);
    // Amalgamation stores explicit zeros too, and still solves alike
    let merged = symbolic
        .factor_indefinite(&a, IndefiniteOptions::default())
        .unwrap();
    assert!(merged.l_entries() > f.l_entries());
    let b = right_hand_side(n, 5);
    let (x, y) = (merged.solve(&b).unwrap(), cholesky.solve(&b).unwrap());
    let worst = x
        .iter()
        .zip(&y)
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    assert!(worst < 1e-13, "{worst:e}");
    assert_eq!(
        f.inertia(),
        Inertia {
            positive: n,
            negative: 0
        }
    );
    let b = right_hand_side(n, 5);
    let (x, y) = (f.solve(&b).unwrap(), cholesky.solve(&b).unwrap());
    let worst = x
        .iter()
        .zip(&y)
        .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()));
    assert!(worst < 1e-13, "{worst:e}");
}

#[test]
fn complex_hermitian_kkt_matrices_factor() {
    // The KKT shape with complex couplings: Hermitian, inertia (n, m)
    let (n, m) = (30, 10);
    let real = kkt(n, m, 9);
    let mut rng = Lcg(11);
    let triplets: Vec<(usize, usize, C)> = real
        .iter()
        .filter(|(i, j, _)| i >= j)
        .flat_map(|(i, j, v)| {
            let z = if i == j {
                C::new(v, 0.0)
            } else {
                C::new(v, rng.unit())
            };
            if i == j {
                vec![(i, j, z)]
            } else {
                vec![(i, j, z), (j, i, z.conj())]
            }
        })
        .collect();
    let a = CsrMatrix::from_triplets(n + m, n + m, triplets).unwrap();
    let f = a.ldlt_indefinite().unwrap();
    assert_eq!(
        f.inertia(),
        Inertia {
            positive: n,
            negative: m
        }
    );
    let b: Vec<C> = (0..n + m).map(|_| C::new(rng.unit(), rng.unit())).collect();
    let x = f.solve(&b).unwrap();
    let ax = a.mul_vec(&x).unwrap();
    let worst = ax
        .iter()
        .zip(&b)
        .fold(0.0f64, |w, (a, b)| w.max((a - b).norm()));
    let scale = x.iter().fold(0.0f64, |w, x| w.max(x.norm()));
    eprintln!("complex: worst residual {worst:e}, ‖x‖∞ {scale:e}");
    assert!(worst < 1e-12 * (1.0 + scale), "{worst:e}");
}

#[test]
fn one_analysis_serves_matrices_with_the_same_pattern() {
    let a = kkt(40, 15, 3);
    let symbolic = a.analyze_cholesky(Ordering::Amd).unwrap();
    // The same pattern with other values
    let mut rng = Lcg(4);
    let mut other = a.clone();
    for v in other.values_mut() {
        if *v != 0.0 {
            *v *= 1.0 + 0.5 * rng.unit().abs();
        }
    }
    // Keep it exactly symmetric: rebuild from the lower triangle
    let other = CsrMatrix::from_triplets(
        55,
        55,
        other
            .iter()
            .filter(|(i, j, _)| i >= j)
            .flat_map(|(i, j, v)| {
                if i == j {
                    vec![(i, j, v)]
                } else {
                    vec![(i, j, v), (j, i, v)]
                }
            }),
    )
    .unwrap();
    for matrix in [&a, &other] {
        let f = symbolic
            .factor_indefinite(matrix, IndefiniteOptions::default())
            .unwrap();
        let b = right_hand_side(55, 8);
        assert!(residual(matrix, &f.solve(&b).unwrap(), &b) < RESIDUAL_BAR);
    }
}

#[test]
fn refinement_and_many_right_hand_sides() {
    let a = kkt(50, 20, 21);
    let f = a.ldlt_indefinite().unwrap();
    let b = right_hand_side(70, 2);
    let refined = f.solve_refined_normwise(&a, &b, 5, 1e-15).unwrap();
    assert!(refined.relative_residual <= 1e-15 || refined.steps <= 5);
    assert!(
        refined.relative_residual
            <= f.solve_refined_normwise(&a, &b, 0, 0.0)
                .unwrap()
                .relative_residual
    );
    let componentwise = f.solve_refined(&a, &b, 5).unwrap();
    assert!(
        componentwise.backward_error < 1e-15,
        "{:e}",
        componentwise.backward_error
    );
    let rhs: Vec<Vec<f64>> = (0..4).map(|s| right_hand_side(70, s)).collect();
    let many = f.solve_many(&rhs).unwrap();
    for (b, x) in rhs.iter().zip(&many) {
        assert_eq!(x, &f.solve(b).unwrap());
    }
    #[cfg(feature = "rayon")]
    assert_eq!(f.par_solve_many(&rhs).unwrap(), many);
    assert!(matches!(
        f.solve(&b[..69]),
        Err(SparseError::ShapeMismatch { .. })
    ));
    assert!(matches!(
        f.solve_many(&[b.clone(), b[..3].to_vec()]),
        Err(SparseError::ShapeMismatch { .. })
    ));
}

#[test]
fn errors() {
    let a = kkt(10, 4, 1);
    let symbolic = a.analyze_cholesky(Ordering::Amd).unwrap();
    // Another pattern
    let other = kkt(10, 4, 2);
    assert_eq!(
        symbolic
            .factor_indefinite(&other, IndefiniteOptions::default())
            .unwrap_err(),
        SparseError::PatternMismatch
    );
    // Not symmetric
    let mut asymmetric = a.clone();
    asymmetric.values_mut()[1] += 1.0;
    assert_eq!(
        symbolic
            .factor_indefinite(&asymmetric, IndefiniteOptions::default())
            .unwrap_err(),
        SparseError::NotSymmetric
    );
    // Singular: two identical rows, or a zero matrix
    let singular = CsrMatrix::from_dense(&[
        vec![1.0, 1.0, 0.0],
        vec![1.0, 1.0, 0.0],
        vec![0.0, 0.0, 2.0],
    ])
    .unwrap();
    assert!(matches!(
        singular.ldlt_indefinite(),
        Err(SparseError::SingularMatrix { .. })
    ));
    let zero = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 0.0), (1, 1, 0.0)]).unwrap();
    assert!(matches!(
        zero.ldlt_indefinite(),
        Err(SparseError::SingularMatrix { .. })
    ));
}

#[test]
#[should_panic(expected = "the threshold must be in (0, 0.5]")]
fn a_threshold_above_one_half_is_refused() {
    let a = kkt(4, 2, 1);
    let _ = a
        .analyze_cholesky(Ordering::Amd)
        .unwrap()
        .factor_indefinite(&a, IndefiniteOptions { threshold: 0.75 });
}

/// Everything a factorization reports, to compare two bit for bit: the solve of a fixed
/// right-hand side, the inertia and the counts
fn fingerprint(
    f: &sparse_mat::compressed_matrix::IndefiniteLdlt<f64>,
    n: usize,
) -> (Vec<u64>, Inertia, usize, usize, usize) {
    let b = right_hand_side(n, 99);
    let x = f.solve(&b).unwrap();
    (
        x.iter().map(|v| v.to_bits()).collect(),
        f.inertia(),
        f.delayed(),
        f.two_by_two(),
        f.l_entries(),
    )
}

#[test]
fn one_kept_analysis_factors_as_fresh_ones_do() {
    // Matrices with the same pattern and other values, factored in turn by one kept analysis
    //   and each by a fresh one, and by SymbolicCholesky::factor_indefinite: the same bit for
    //   bit, so nothing is left behind in the analysis or reused buffers
    let base = kkt(60, 25, 13);
    let kept = base.analyze_indefinite(Ordering::Amd).unwrap();
    let mut rng = Lcg(14);
    for round in 0..6 {
        let lower: Vec<(usize, usize, f64)> = base
            .iter()
            .filter(|(i, j, _)| i >= j)
            .map(|(i, j, v)| {
                let scaled = if v == 0.0 {
                    0.0
                } else {
                    v * (0.5 + rng.unit().abs())
                };
                (i, j, scaled)
            })
            .collect();
        let a = CsrMatrix::from_triplets(
            85,
            85,
            lower.iter().flat_map(|&(i, j, v)| {
                if i == j {
                    vec![(i, j, v)]
                } else {
                    vec![(i, j, v), (j, i, v)]
                }
            }),
        )
        .unwrap();
        let fresh = a.analyze_indefinite(Ordering::Amd).unwrap();
        let options = IndefiniteOptions {
            threshold: [0.01, 1e-4, 0.5][round % 3],
        };
        let by_kept = fingerprint(&kept.factor(&a, options).unwrap(), 85);
        assert_eq!(
            by_kept,
            fingerprint(&fresh.factor(&a, options).unwrap(), 85),
            "round {round}"
        );
        assert_eq!(
            by_kept,
            fingerprint(
                &fresh.cholesky().factor_indefinite(&a, options).unwrap(),
                85
            ),
            "round {round}"
        );
    }
}

#[test]
fn the_symmetry_check_refuses_what_is_hermitian_refuses() {
    // Each matrix is analyzed on its own pattern, so only the values decide
    let real =
        |triplets: Vec<(usize, usize, f64)>| CsrMatrix::from_triplets(3, 3, triplets).unwrap();
    let diagonal = vec![(0, 0, 2.0), (1, 1, 3.0), (2, 2, 4.0)];
    let with = |extra: &[(usize, usize, f64)]| {
        let mut t = diagonal.clone();
        t.extend_from_slice(extra);
        real(t)
    };
    let cases = [
        ("symmetric", with(&[(0, 1, 1.0), (1, 0, 1.0)])),
        ("asymmetric values", with(&[(0, 1, 1.0), (1, 0, 1.5)])),
        // (0, 1) stored, (1, 0) not: symmetric while the stored value is zero
        ("a stored zero without its mirror", with(&[(0, 1, 0.0)])),
        ("a nonzero without its mirror", with(&[(0, 1, 1.0)])),
        (
            "a nonzero against a stored zero",
            with(&[(0, 1, 1.0), (1, 0, 0.0)]),
        ),
        (
            "negative zero against zero",
            with(&[(0, 1, -0.0), (1, 0, 0.0)]),
        ),
        ("a NaN pair", with(&[(0, 1, f64::NAN), (1, 0, f64::NAN)])),
    ];
    for (name, a) in &cases {
        let refused = matches!(
            a.analyze_indefinite(Ordering::Natural)
                .unwrap()
                .factor(a, IndefiniteOptions::default()),
            Err(SparseError::NotSymmetric)
        );
        assert_eq!(refused, !a.is_hermitian(), "{name}");
    }
    let complex =
        |triplets: Vec<(usize, usize, C)>| CsrMatrix::from_triplets(2, 2, triplets).unwrap();
    let complex_cases = [
        (
            "Hermitian",
            complex(vec![
                (0, 0, C::new(2.0, 0.0)),
                (1, 1, C::new(3.0, 0.0)),
                (0, 1, C::new(1.0, 1.0)),
                (1, 0, C::new(1.0, -1.0)),
            ]),
        ),
        (
            "symmetric but not Hermitian",
            complex(vec![
                (0, 0, C::new(2.0, 0.0)),
                (1, 1, C::new(3.0, 0.0)),
                (0, 1, C::new(1.0, 1.0)),
                (1, 0, C::new(1.0, 1.0)),
            ]),
        ),
        (
            "a complex diagonal",
            complex(vec![(0, 0, C::new(2.0, 0.5)), (1, 1, C::new(3.0, 0.0))]),
        ),
    ];
    for (name, a) in &complex_cases {
        let refused = matches!(
            a.analyze_indefinite(Ordering::Natural)
                .unwrap()
                .factor(a, IndefiniteOptions::default()),
            Err(SparseError::NotSymmetric)
        );
        assert_eq!(refused, !a.is_hermitian(), "{name}");
    }
    // Each kind occurs: some refused, some accepted
    assert!(
        cases.iter().any(|(_, a)| a.is_hermitian()) && cases.iter().any(|(_, a)| !a.is_hermitian())
    );
}

#[test]
fn amalgamation_postorders_the_tree_and_merges_without_changing_the_answers() {
    for seed in 0..8 {
        let a = kkt(80, 30, 40 + seed);
        for ordering in [Ordering::Amd, Ordering::Rcm, Ordering::Natural] {
            let cholesky = a.analyze_cholesky(ordering).unwrap();
            let merged = SymbolicIndefinite::new(cholesky.clone());
            // The same tree, renumbered in postorder: the same pattern's size, and a postorder
            //   that is the identity
            let postordered = merged.cholesky();
            assert!(postordered
                .postorder()
                .iter()
                .enumerate()
                .all(|(k, &node)| k == node));
            assert_eq!(
                postordered.column_counts().iter().sum::<usize>(),
                cholesky.column_counts().iter().sum::<usize>()
            );
            let unmerged = SymbolicIndefinite::with_amalgamation(cholesky, Amalgamation::NONE);
            for threshold in [0.01, 1e-4] {
                let options = IndefiniteOptions { threshold };
                let (f, g) = (
                    merged.factor(&a, options).unwrap(),
                    unmerged.factor(&a, options).unwrap(),
                );
                assert!(f.l_entries() >= g.l_entries());
                assert_eq!(f.inertia(), g.inertia());
                let b = right_hand_side(110, seed);
                for factor in [&f, &g] {
                    let x = factor.solve(&b).unwrap();
                    assert!(
                        residual(&a, &x, &b) < RESIDUAL_BAR,
                        "{:e}",
                        residual(&a, &x, &b)
                    );
                }
            }
        }
    }
}

#[test]
fn merging_everything_merges_only_adjacent_fronts() {
    // A rule that merges whatever it may: each child joins its parent only when its columns come
    //   just before the parent's, so every tree ends as one front, dense, and the factorization
    //   is still exact
    let everything = Amalgamation {
        columns: [usize::MAX; 3],
        zeros: [2.0; 3],
    };
    for seed in 0..4 {
        let a = kkt(40, 15, 70 + seed);
        let cholesky = a.analyze_cholesky(Ordering::Amd).unwrap();
        let f = SymbolicIndefinite::with_amalgamation(cholesky, everything)
            .factor(&a, IndefiniteOptions::default())
            .unwrap();
        // One front per tree of the forest: every pivot found within it, nothing delayed
        assert_eq!(f.delayed(), 0);
        assert_eq!(
            f.inertia(),
            Inertia {
                positive: 40,
                negative: 15
            }
        );
        let b = right_hand_side(55, seed);
        let x = f.solve(&b).unwrap();
        assert!(
            residual(&a, &x, &b) < RESIDUAL_BAR,
            "{:e}",
            residual(&a, &x, &b)
        );
    }
}

#[test]
fn amalgamation_keeps_the_fill_bounded() {
    // The default rule allows explicit zeros only in small fronts, or a few percent in large
    //   ones, so the stored entries stay within a small multiple of the fundamental supernodes'.
    //   Measured on these matrices: up to 2.4× (2026-11-08). (Counting a merged front's zeros as
    //   true entries once let merges cascade into one dense front: 480× on an OPF KKT matrix.)
    let mut worst = 0.0f64;
    for seed in 0..4 {
        let a = kkt(400, 150, 90 + seed);
        let cholesky = a.analyze_cholesky(Ordering::Amd).unwrap();
        let options = IndefiniteOptions::default();
        let merged = SymbolicIndefinite::new(cholesky.clone())
            .factor(&a, options)
            .unwrap();
        let unmerged = SymbolicIndefinite::with_amalgamation(cholesky, Amalgamation::NONE)
            .factor(&a, options)
            .unwrap();
        worst = worst.max(merged.l_entries() as f64 / unmerged.l_entries() as f64);
    }
    eprintln!("stored entries with amalgamation: up to {worst:.2}× without");
    assert!(worst < 4.0, "{worst:.2}");
}
