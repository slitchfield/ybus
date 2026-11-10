// Tests for the dense threshold Bunch–Kaufman factorization (DenseLdlt): on matrices built with a
//   known inertia, real and complex Hermitian, on KKT-shaped matrices with zero diagonals that
//   static pivots fail on, its bound on L, the threshold's effect, and its errors.

// Dense matrices' entries are reached by (i, j) here, which reads more plainly than iterators
#![allow(clippy::needless_range_loop)]

use sparse_mat::compressed_matrix::{CsrMatrix, DBlock, DenseLdlt, Inertia, Ordering};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::Complex;
use sparse_mat::Factorization;

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

fn to_complex(rows: &[Vec<f64>]) -> Vec<Vec<C>> {
    rows.iter()
        .map(|row| row.iter().map(|&v| C::new(v, 0.0)).collect())
        .collect()
}

fn matmul(a: &[Vec<C>], b: &[Vec<C>]) -> Vec<Vec<C>> {
    let (n, k, m) = (a.len(), b.len(), b[0].len());
    (0..n)
        .map(|i| {
            (0..m)
                .map(|j| (0..k).map(|p| a[i][p] * b[p][j]).sum())
                .collect()
        })
        .collect()
}

fn adjoint(a: &[Vec<C>]) -> Vec<Vec<C>> {
    (0..a[0].len())
        .map(|j| (0..a.len()).map(|i| a[i][j].conj()).collect())
        .collect()
}

/// A Hermitian matrix `Q Λ Qᴴ` with the given eigenvalues, `Q` a product of three random
/// Householder reflectors (complex ones if `complex`), made exactly Hermitian
fn with_eigenvalues(eigenvalues: &[f64], complex: bool, seed: u64) -> Vec<Vec<C>> {
    let n = eigenvalues.len();
    let mut rng = Lcg(seed);
    let mut q: Vec<Vec<C>> = (0..n)
        .map(|i| (0..n).map(|j| C::new((i == j) as u8 as f64, 0.0)).collect())
        .collect();
    for _ in 0..3 {
        let v: Vec<C> = (0..n)
            .map(|_| C::new(rng.unit(), if complex { rng.unit() } else { 0.0 }))
            .collect();
        let norm2: f64 = v.iter().map(|v| v.norm_sqr()).sum();
        let h: Vec<Vec<C>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        C::new((i == j) as u8 as f64, 0.0) - v[i] * v[j].conj() * (2.0 / norm2)
                    })
                    .collect()
            })
            .collect();
        q = matmul(&q, &h);
    }
    let lambda: Vec<Vec<C>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| C::new(if i == j { eigenvalues[i] } else { 0.0 }, 0.0))
                .collect()
        })
        .collect();
    let mut a = matmul(&matmul(&q, &lambda), &adjoint(&q));
    for i in 0..n {
        a[i][i] = C::new(a[i][i].re, 0.0);
        for j in 0..i {
            a[j][i] = a[i][j].conj();
        }
    }
    a
}

/// `D`, dense, from its blocks
fn d_dense(blocks: &[DBlock<C>], n: usize) -> Vec<Vec<C>> {
    let mut d = vec![vec![C::new(0.0, 0.0); n]; n];
    let mut at = 0;
    for block in blocks {
        match *block {
            DBlock::One(v) => d[at][at] = v,
            DBlock::Two { a, b, c } => {
                d[at][at] = a;
                d[at + 1][at] = b;
                d[at][at + 1] = b.conj();
                d[at + 1][at + 1] = c;
            }
        }
        at += block.size();
    }
    d
}

fn max_abs(a: &[Vec<C>]) -> f64 {
    a.iter().flatten().fold(0.0, |m, v| m.max(v.norm()))
}

/// What every factorization must satisfy: `P A Pᵀ = L D Lᴴ` to rounding, `L` unit lower
/// triangular with every entry at most `1/u`, zero where a 2×2 block sits, and a
/// backward-stable solve. Returns the reconstruction's and the solve's scaled errors.
fn check(a: &[Vec<C>], f: &DenseLdlt<C>, u: f64) -> (f64, f64) {
    let n = a.len();
    let l = f.l();
    let perm = f.permutation();
    let mut sorted = perm.to_vec();
    sorted.sort_unstable();
    assert_eq!(sorted, (0..n).collect::<Vec<_>>(), "a permutation");
    let mut at = 0;
    for block in f.d_blocks() {
        if block.size() == 2 {
            assert_eq!(l[at + 1][at], C::new(0.0, 0.0), "L inside a 2×2 block");
        }
        at += block.size();
    }
    assert_eq!(at, n, "the blocks cover the matrix");
    for i in 0..n {
        assert_eq!(l[i][i], C::new(1.0, 0.0));
        for j in i + 1..n {
            assert_eq!(l[i][j], C::new(0.0, 0.0));
        }
        for j in 0..i {
            assert!(
                l[i][j].norm() <= (1.0 + 1e-12) / u,
                "|L({i}, {j})| = {} above 1/u = {}",
                l[i][j].norm(),
                1.0 / u
            );
        }
    }
    let ldl = matmul(&matmul(&l, &d_dense(f.d_blocks(), n)), &adjoint(&l));
    let pap: Vec<Vec<C>> = (0..n)
        .map(|i| (0..n).map(|j| a[perm[i]][perm[j]]).collect())
        .collect();
    let difference: Vec<Vec<C>> = (0..n)
        .map(|i| (0..n).map(|j| pap[i][j] - ldl[i][j]).collect())
        .collect();
    let reconstruction = max_abs(&difference) / max_abs(a);

    let mut rng = Lcg(n as u64 * 31 + 7);
    let b: Vec<C> = (0..n).map(|_| C::new(rng.unit(), rng.unit())).collect();
    let x = f.solve(&b).unwrap();
    let ax: Vec<C> = (0..n)
        .map(|i| (0..n).map(|j| a[i][j] * x[j]).sum())
        .collect();
    let norm_a = (0..n)
        .map(|i| a[i].iter().map(|v| v.norm()).sum::<f64>())
        .fold(0.0, f64::max);
    let inf = |v: &[C]| v.iter().fold(0.0f64, |m, v| m.max(v.norm()));
    let r: Vec<C> = ax.iter().zip(&b).map(|(a, b)| a - b).collect();
    let solve = inf(&r) / (norm_a * inf(&x) + inf(&b));
    (reconstruction, solve)
}

fn inertia_of(eigenvalues: &[f64]) -> Inertia {
    Inertia {
        positive: eigenvalues.iter().filter(|e| **e > 0.0).count(),
        negative: eigenvalues.iter().filter(|e| **e < 0.0).count(),
    }
}

/// Eigenvalues of mixed signs and magnitudes from 1e-3 to 1e3
fn mixed_eigenvalues(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed);
    (0..n)
        .map(|_| {
            let sign = if rng.next().is_multiple_of(2) {
                1.0
            } else {
                -1.0
            };
            sign * 10f64.powf(rng.unit() * 3.0)
        })
        .collect()
}

/// The worst scaled errors allowed, about ten times the worst measured over these tests'
/// matrices (n up to 40, real and complex, thresholds 0.01 to 0.5): 5.3e-15 for the
/// reconstruction and 5.7e-16 for the solve (2026-11-07)
const RECONSTRUCTION_BAR: f64 = 5e-14;
const SOLVE_BAR: f64 = 5e-15;

#[test]
fn a_zero_diagonal_takes_a_2x2_pivot() {
    let swap = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
    let f = DenseLdlt::factor(&swap, DenseLdlt::<f64>::DEFAULT_THRESHOLD).unwrap();
    assert_eq!(
        f.d_blocks(),
        &[DBlock::Two {
            a: 0.0,
            b: 1.0,
            c: 0.0
        }]
    );
    assert_eq!(
        f.inertia(),
        Inertia {
            positive: 1,
            negative: 1
        }
    );
    assert_eq!(f.solve(&[1.0, 2.0]).unwrap(), vec![2.0, 1.0]);
    // The static pivots of Ldlt fail on it
    let csr = CsrMatrix::from_dense(&swap).unwrap();
    assert!(csr.ldlt().is_err());
}

#[test]
fn matrices_of_known_inertia_factor_stably_at_every_threshold() {
    let mut worst = (0.0f64, 0.0f64);
    for seed in 0..20 {
        for n in [1, 2, 5, 17, 40] {
            let eigenvalues = mixed_eigenvalues(n, seed * 101 + n as u64);
            let a = with_eigenvalues(&eigenvalues, false, seed);
            let real: Vec<Vec<f64>> = a.iter().map(|r| r.iter().map(|v| v.re).collect()).collect();
            for u in [0.01, 0.1, 0.5] {
                let f = DenseLdlt::factor(&real, u).unwrap();
                assert_eq!(
                    f.inertia(),
                    inertia_of(&eigenvalues),
                    "n {n}, seed {seed}, u {u}"
                );
                // Check through the complex helpers, with the real factors converted
                let fc = DenseLdlt::factor(&a, u).unwrap();
                assert_eq!(fc.permutation(), f.permutation());
                let (reconstruction, solve) = check(&a, &fc, u);
                worst = (worst.0.max(reconstruction), worst.1.max(solve));
            }
        }
    }
    eprintln!(
        "worst scaled reconstruction error {:e}, solve residual {:e}",
        worst.0, worst.1
    );
    assert!(
        worst.0 < RECONSTRUCTION_BAR && worst.1 < SOLVE_BAR,
        "{worst:?}"
    );
}

#[test]
fn complex_hermitian_matrices_factor_with_their_inertia() {
    let mut worst = (0.0f64, 0.0f64);
    let mut two_by_two = 0;
    for seed in 0..20 {
        for n in [2, 9, 30] {
            let eigenvalues = mixed_eigenvalues(n, seed * 7 + 3);
            let a = with_eigenvalues(&eigenvalues, true, seed + 1000);
            for u in [0.01, 0.5] {
                let f = DenseLdlt::factor(&a, u).unwrap();
                assert_eq!(f.inertia(), inertia_of(&eigenvalues));
                for block in f.d_blocks() {
                    match *block {
                        DBlock::One(d) => assert_eq!(d.im, 0.0),
                        DBlock::Two { a, c, .. } => {
                            assert_eq!((a.im, c.im), (0.0, 0.0));
                            two_by_two += 1;
                        }
                    }
                }
                let (reconstruction, solve) = check(&a, &f, u);
                worst = (worst.0.max(reconstruction), worst.1.max(solve));
            }
        }
    }
    eprintln!(
        "complex: worst reconstruction {:e}, solve {:e}, {two_by_two} 2×2 blocks",
        worst.0, worst.1
    );
    assert!(two_by_two > 0, "the complex 2×2 path is exercised");
    assert!(
        worst.0 < RECONSTRUCTION_BAR && worst.1 < SOLVE_BAR,
        "{worst:?}"
    );
}

#[test]
fn kkt_matrices_with_zero_diagonals_factor_where_static_pivots_fail() {
    // [W Bᵀ; B 0] with W diagonal, half its entries zero, and B dense: W is positive definite
    //   on B's null space, so the inertia is (n, m)
    for seed in 0..10 {
        let (n, m) = (24, 12);
        let mut rng = Lcg(seed + 500);
        let mut a = vec![vec![0.0; n + m]; n + m];
        for i in 0..n {
            a[i][i] = if i % 2 == 0 {
                0.0
            } else {
                1.0 + rng.unit().abs() * 1e3
            };
        }
        for r in 0..m {
            for c in 0..n {
                let v = rng.unit();
                a[n + r][c] = v;
                a[c][n + r] = v;
            }
        }
        let f = DenseLdlt::factor(&a, DenseLdlt::<f64>::DEFAULT_THRESHOLD).unwrap();
        assert_eq!(
            f.inertia(),
            Inertia {
                positive: n,
                negative: m
            }
        );
        let (reconstruction, solve) = check(
            &to_complex(&a),
            &DenseLdlt::factor(&to_complex(&a), 0.01).unwrap(),
            0.01,
        );
        assert!(
            reconstruction < RECONSTRUCTION_BAR && solve < SOLVE_BAR,
            "{reconstruction:e} {solve:e}"
        );
        // Static pivots in the natural order meet the zero diagonal first
        let csr = CsrMatrix::from_triplets(
            n + m,
            n + m,
            (0..n + m)
                .flat_map(|i| (0..n + m).map(move |j| (i, j)))
                .map(|(i, j)| (i, j, a[i][j])),
        )
        .unwrap();
        let symbolic = csr.analyze_cholesky(Ordering::Natural).unwrap();
        assert_eq!(
            symbolic.factor_ldlt(&csr).unwrap_err(),
            SparseError::SingularMatrix { pivot: 0 }
        );
    }
}

#[test]
fn the_threshold_trades_the_diagonal_order_for_smaller_entries() {
    // A larger threshold takes more 2×2 pivots
    let (mut small, mut large) = (0, 0);
    for seed in 0..20 {
        let eigenvalues = mixed_eigenvalues(30, seed);
        let a: Vec<Vec<f64>> = with_eigenvalues(&eigenvalues, false, seed)
            .iter()
            .map(|r| r.iter().map(|v| v.re).collect())
            .collect();
        let count = |u| {
            DenseLdlt::factor(&a, u)
                .unwrap()
                .d_blocks()
                .iter()
                .filter(|b| b.size() == 2)
                .count()
        };
        small += count(0.01);
        large += count(0.5);
    }
    eprintln!("2×2 blocks: {small} at u = 0.01, {large} at u = 0.5");
    assert!(large > small);

    // A diagonally dominant matrix keeps its order and its diagonal pivots, at any threshold
    let mut rng = Lcg(3);
    let n = 20;
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..i {
            let v = rng.unit();
            a[i][j] = v;
            a[j][i] = v;
        }
        a[i][i] = if i.is_multiple_of(3) { -1.0 } else { 1.0 } * (n as f64 + 1.0);
    }
    for u in [0.01, 0.5] {
        let f = DenseLdlt::factor(&a, u).unwrap();
        assert_eq!(f.permutation(), (0..n).collect::<Vec<_>>());
        assert!(f.d_blocks().iter().all(|b| b.size() == 1));
        assert_eq!(
            f.inertia().negative,
            (0..n).filter(|i| i.is_multiple_of(3)).count()
        );
    }
}

#[test]
fn errors() {
    // A singular matrix: no pivot exists once the zero row is reached
    let singular = vec![
        vec![1.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0],
        vec![0.0, 0.0, 2.0],
    ];
    assert_eq!(
        DenseLdlt::factor(&singular, 0.01).unwrap_err(),
        SparseError::SingularMatrix { pivot: 1 }
    );
    let rank_one = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
    assert!(matches!(
        DenseLdlt::factor(&rank_one, 0.5),
        Err(SparseError::SingularMatrix { .. })
    ));
    assert_eq!(
        DenseLdlt::factor(&[vec![1.0, 2.0], vec![3.0, 1.0]], 0.01).unwrap_err(),
        SparseError::NotSymmetric
    );
    let complex_diagonal = vec![vec![C::new(1.0, 1.0)]];
    assert_eq!(
        DenseLdlt::factor(&complex_diagonal, 0.01).unwrap_err(),
        SparseError::NotSymmetric
    );
    assert!(matches!(
        DenseLdlt::factor(&[vec![1.0, 2.0]], 0.01),
        Err(SparseError::NotSquare { .. })
    ));
    assert!(matches!(
        DenseLdlt::factor(&[vec![1.0, 2.0], vec![3.0]], 0.01),
        Err(SparseError::RaggedRows { .. })
    ));
    let f = DenseLdlt::factor(&[vec![2.0]], 0.01).unwrap();
    assert!(matches!(
        f.solve(&[1.0, 2.0]),
        Err(SparseError::ShapeMismatch { .. })
    ));
    let empty = DenseLdlt::<f64>::factor(&[], 0.01).unwrap();
    assert_eq!(empty.solve(&[]).unwrap(), Vec::<f64>::new());
}

#[test]
#[should_panic(expected = "the threshold must be in (0, 0.5]")]
fn a_threshold_above_one_half_is_refused() {
    let _ = DenseLdlt::factor(&[vec![1.0]], 0.6);
}
