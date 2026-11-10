// Tests for Phase 5's sparse right-hand sides: `SparseVector::new`, and the sparse forward and
//   full solves of `Cholesky` and `Lu`, against their dense solves and a brute-force reach.

use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix, Ordering};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::sparse_matrix::SparseVector;
use sparse_mat::Complex;
use sparse_mat::Factorization;

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

/// A pseudo-random SPD matrix: about `fill_percent` of the off-diagonal positions in a
/// symmetric pattern, with a dominant diagonal.
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

/// Right-hand sides with one, two and a few entries, spread over `0..n`.
fn sparse_rhs(n: usize) -> Vec<SparseVector<f64>> {
    vec![
        SparseVector::new(n, vec![0], vec![1.0]).unwrap(),
        SparseVector::new(n, vec![n / 2], vec![-2.5]).unwrap(),
        SparseVector::new(n, vec![n / 3, n - 1], vec![1.0, -1.0]).unwrap(),
        SparseVector::new(
            n,
            vec![1, n / 4, n / 2 + 1, n - 2],
            vec![0.5, 3.0, -1.0, 2.0],
        )
        .unwrap(),
    ]
}

/// The structural reach by brute force: every position `L y = c` can make nonzero, where `L`
/// is lower triangular and `starts` are `c`'s positions, by propagating down the rows.
fn brute_force_reach(l_lower: &CsrMatrix<f64>, starts: &[usize]) -> Vec<usize> {
    let n = l_lower.shape().0;
    let mut reached = vec![false; n];
    for &start in starts {
        reached[start] = true;
    }
    for i in 0..n {
        let (cols, _) = l_lower.outer(i).unwrap();
        if cols.iter().any(|&j| j < i && reached[j]) {
            reached[i] = true;
        }
    }
    (0..n).filter(|&i| reached[i]).collect()
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |m, v| m.max(v.abs()))
}

#[test]
fn sparse_vector_new_validates() {
    let v = SparseVector::new(5, vec![1, 3], vec![2.0, 0.0]).unwrap();
    assert_eq!((v.len(), v.nnz()), (5, 2));
    assert_eq!(v.to_dense(), vec![0.0, 2.0, 0.0, 0.0, 0.0]);
    assert_eq!(
        SparseVector::new(5, vec![1, 3], vec![2.0]),
        Err(SparseError::LengthMismatch {
            expected: 2,
            found: 1
        })
    );
    assert_eq!(
        SparseVector::new(5, vec![1, 5], vec![2.0, 1.0]),
        Err(SparseError::OutOfRange { index: 5, len: 5 })
    );
    assert!(matches!(
        SparseVector::new(5, vec![3, 1], vec![2.0, 1.0]),
        Err(SparseError::InvalidCompressed { .. })
    ));
    assert!(matches!(
        SparseVector::new(5, vec![2, 2], vec![2.0, 1.0]),
        Err(SparseError::InvalidCompressed { .. })
    ));
    assert!(SparseVector::<f64>::new(0, vec![], vec![])
        .unwrap()
        .is_empty());
}

#[test]
fn cholesky_sparse_solves_match_the_dense_solve() {
    for (n, fill, seed) in [(40, 5, 1), (80, 3, 2), (120, 2, 3)] {
        let a = random_spd(n, fill, seed);
        let factor = a
            .analyze_cholesky(Ordering::Amd)
            .unwrap()
            .factor(&a)
            .unwrap();
        let perm = factor.permutation();
        for b in sparse_rhs(n) {
            let dense_b = b.to_dense();
            let expected = factor.solve(&dense_b).unwrap();
            let x = factor.solve_sparse(&b).unwrap();
            let scale = max_abs(&expected);
            for (got, want) in x.iter().zip(&expected) {
                assert!((got - want).abs() <= 1e-13 * scale, "{got} vs {want}");
            }

            // The forward solve: L y = P b, against the dense triangular solve
            let permuted: Vec<f64> = perm.iter().map(|&row| dense_b[row]).collect();
            let dense_y = factor.l().solve_lower_triangular(&permuted).unwrap();
            let y = factor.solve_lower_sparse(&b).unwrap();
            let mut inverse = vec![0; n];
            for (new, &old) in perm.iter().enumerate() {
                inverse[old] = new;
            }
            let mut starts: Vec<usize> = b.indices().iter().map(|&row| inverse[row]).collect();
            starts.sort_unstable();
            assert_eq!(y.indices(), brute_force_reach(factor.l(), &starts));
            let y_scale = max_abs(&dense_y);
            for (i, want) in dense_y.iter().enumerate() {
                let got = y.get(i).unwrap_or(0.0);
                assert!(
                    (got - want).abs() <= 1e-14 * y_scale,
                    "{i}: {got} vs {want}"
                );
            }
            // ‖y‖² is bᵀ A⁻¹ b
            let quadratic: f64 = y.values().iter().map(|v| v * v).sum();
            let direct: f64 = dense_b.iter().zip(&expected).map(|(b, x)| b * x).sum();
            assert!((quadratic - direct).abs() <= 1e-12 * direct.abs());
        }
    }
}

#[test]
fn a_sparse_forward_solve_can_touch_little_of_the_factor() {
    // A path graph's Laplacian, plus 1 on the diagonal, ordered naturally: L is bidiagonal, so
    //   a right-hand side at the last row reaches only that row
    let n = 50;
    let mut triplets = vec![];
    for i in 0..n {
        let degree = if i == 0 || i == n - 1 { 1.0 } else { 2.0 };
        triplets.push((i, i, degree + 1.0));
        if i + 1 < n {
            triplets.extend([(i, i + 1, -1.0), (i + 1, i, -1.0)]);
        }
    }
    let a = CsrMatrix::from_triplets(n, n, triplets).unwrap();
    let factor = a
        .analyze_cholesky(Ordering::Natural)
        .unwrap()
        .factor(&a)
        .unwrap();
    let last = SparseVector::new(n, vec![n - 1], vec![1.0]).unwrap();
    assert_eq!(factor.solve_lower_sparse(&last).unwrap().indices(), [n - 1]);
    let first = SparseVector::new(n, vec![0], vec![1.0]).unwrap();
    assert_eq!(factor.solve_lower_sparse(&first).unwrap().nnz(), n);
}

#[test]
fn lu_sparse_solves_match_the_dense_solve() {
    for (n, fill, scale, seed) in [(40, 6, 1.5, 4), (80, 3, 0.05, 5), (120, 2, 0.5, 6)] {
        let a = random_square(n, fill, scale, seed).to_csc();
        let lu = a
            .analyze_lu(Ordering::Amd)
            .unwrap()
            .factor(&a, Default::default());
        let lu = lu.unwrap();
        for b in sparse_rhs(n) {
            let dense_b = b.to_dense();
            // Bit for bit, apart from the sign of a zero, which == ignores
            assert_eq!(lu.solve_sparse(&b).unwrap(), lu.solve(&dense_b).unwrap());

            let permuted: Vec<f64> = lu.permutation().iter().map(|&row| dense_b[row]).collect();
            let dense_y = lu.l().solve_lower_triangular(&permuted).unwrap();
            let y = lu.solve_lower_sparse(&b).unwrap();
            let mut starts: Vec<usize> = b
                .indices()
                .iter()
                .map(|&row| lu.permutation().iter().position(|&r| r == row).unwrap())
                .collect();
            starts.sort_unstable();
            assert_eq!(y.indices(), brute_force_reach(&lu.l().to_csr(), &starts));
            assert_eq!(y.to_dense(), dense_y);
        }
    }
}

#[test]
fn lu_sparse_solves_work_for_complex_matrices() {
    let n = 30;
    let a: CscMatrix<Complex<f64>> = random_square(n, 8, 0.3, 7)
        .to_csc()
        .map(|v| Complex::new(v, 0.25 * v));
    let lu = a.lu().unwrap();
    let b = SparseVector::new(
        n,
        vec![2, 17],
        vec![Complex::new(1.0, -1.0), Complex::new(0.0, 2.0)],
    )
    .unwrap();
    assert_eq!(
        lu.solve_sparse(&b).unwrap(),
        lu.solve(&b.to_dense()).unwrap()
    );
}

#[test]
fn sparse_solves_check_the_length() {
    let a = random_spd(10, 20, 8);
    let b = SparseVector::new(9, vec![0], vec![1.0]).unwrap();
    let mismatch = || {
        Some(SparseError::ShapeMismatch {
            left: (10, 10),
            right: (9, 1),
        })
    };
    let cholesky = a.cholesky().unwrap();
    assert_eq!(cholesky.solve_sparse(&b).err(), mismatch());
    assert_eq!(cholesky.solve_lower_sparse(&b).err(), mismatch());
    let lu = a.to_csc().lu().unwrap();
    assert_eq!(lu.solve_sparse(&b).err(), mismatch());
    assert_eq!(lu.solve_lower_sparse(&b).err(), mismatch());
}

/// Every value's bits, so that "bit for bit" also tells 0.0 from −0.0
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|v| v.to_bits()).collect()
}

/// More right-hand sides for the workspace tests: empty, and every position.
fn more_rhs(n: usize) -> Vec<SparseVector<f64>> {
    let mut rhs = sparse_rhs(n);
    rhs.push(SparseVector::new(n, vec![], vec![]).unwrap());
    rhs.push(
        SparseVector::new(
            n,
            (0..n).collect(),
            (0..n).map(|i| 1.0 + i as f64).collect(),
        )
        .unwrap(),
    );
    rhs
}

#[test]
fn solves_into_a_reused_workspace_match_solve_sparse_bit_for_bit() {
    use sparse_mat::compressed_matrix::SparseSolveWorkspace;
    // One workspace for every solve, through factors of different sizes and kinds, so stale
    //   marks or values from an earlier solve would show
    let mut workspace = SparseSolveWorkspace::new();
    for (n, fill, seed) in [(120, 2, 3), (40, 5, 1), (80, 3, 2)] {
        let spd = random_spd(n, fill, seed);
        let cholesky = spd
            .analyze_cholesky(Ordering::Amd)
            .unwrap()
            .factor(&spd)
            .unwrap();
        let square = random_square(n, fill, 0.5, seed).to_csc();
        let lu = square
            .analyze_lu(Ordering::Amd)
            .unwrap()
            .factor(&square, Default::default())
            .unwrap();
        for _ in 0..2 {
            for b in more_rhs(n) {
                let mut x = vec![f64::NAN; n];
                cholesky
                    .solve_sparse_into(&b, &mut workspace, &mut x)
                    .unwrap();
                assert_eq!(
                    bits(&x),
                    bits(&cholesky.solve_sparse(&b).unwrap()),
                    "Cholesky n={n}"
                );
                let mut x = vec![f64::NAN; n];
                lu.solve_sparse_into(&b, &mut workspace, &mut x).unwrap();
                assert_eq!(bits(&x), bits(&lu.solve_sparse(&b).unwrap()), "LU n={n}");
            }
        }
    }
}

#[test]
fn solves_into_check_both_lengths() {
    use sparse_mat::compressed_matrix::SparseSolveWorkspace;
    let a = random_spd(10, 20, 4);
    let cholesky = a
        .analyze_cholesky(Ordering::Amd)
        .unwrap()
        .factor(&a)
        .unwrap();
    let b = SparseVector::new(10, vec![3], vec![1.0]).unwrap();
    let mut workspace = SparseSolveWorkspace::new();
    let mut short = vec![0.0; 9];
    assert!(matches!(
        cholesky.solve_sparse_into(&b, &mut workspace, &mut short),
        Err(SparseError::ShapeMismatch { .. })
    ));
    let wrong = SparseVector::new(9, vec![3], vec![1.0]).unwrap();
    let mut x = vec![0.0; 10];
    assert!(matches!(
        cholesky.solve_sparse_into(&wrong, &mut workspace, &mut x),
        Err(SparseError::ShapeMismatch { .. })
    ));
}
