// Tests for COLAMD (Ordering::Colamd): valid column orderings, square and rectangular, with
//   empty and dense rows and columns; correct LU solves under it; AᵀA ordered about as well as
//   AMD orders it explicitly; less fill than AMD when the rows are out of order; and its
//   refusal for Cholesky.

use ybus_sparse::compressed_matrix::{CscMatrix, CsrMatrix, LuOptions, Ordering};
use ybus_sparse::sparse_error::SparseError;
use ybus_sparse::Complex;
use ybus_sparse::Factorization;

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

/// An `m × n` pattern with about `per_row` entries in each row, at random columns, with values.
fn random_pattern(m: usize, n: usize, per_row: usize, seed: u64) -> Vec<(usize, usize, f64)> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for row in 0..m {
        for _ in 0..per_row {
            let col = rng.next() as usize % n;
            triplets.push((row, col, ((rng.next() % 2001) as f64 - 1000.0) / 317.0));
        }
    }
    triplets
}

/// A square, nonsingular, structurally unsymmetric matrix: a dominant diagonal, plus random
/// entries placed without their transposes, more of them in some rows than others.
fn unsymmetric(n: usize, seed: u64) -> CscMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for row in 0..n {
        let count = 1 + (rng.next() % 6) as usize + if row % 7 == 0 { 8 } else { 0 };
        let mut sum = 1.0;
        for _ in 0..count {
            let col = rng.next() as usize % n;
            if col != row {
                let value = ((rng.next() % 2001) as f64 - 1000.0) / 500.0;
                sum += value.abs();
                triplets.push((row, col, value));
            }
        }
        triplets.push((row, row, sum));
    }
    CscMatrix::from_triplets(n, n, triplets).unwrap()
}

fn is_permutation(perm: &[usize], n: usize) -> bool {
    let mut seen = vec![false; n];
    perm.len() == n
        && perm
            .iter()
            .all(|&j| j < n && !std::mem::replace(&mut seen[j], true))
}

/// The LU's stored entries, L and U together, counting the diagonal once
fn fill(a: &CscMatrix<f64>, ordering: Ordering) -> usize {
    let lu = a
        .analyze_lu(ordering)
        .unwrap()
        .factor(a, LuOptions::default())
        .unwrap();
    lu.l().num_stored() + lu.u().num_stored() - a.shape().0
}

#[test]
fn colamd_gives_a_column_permutation_of_any_shape() {
    for (m, n, per_row, seed) in [
        (50, 50, 3, 1),
        (80, 30, 4, 2),
        (30, 80, 6, 3),
        (200, 200, 2, 4),
    ] {
        let a = CsrMatrix::from_triplets(m, n, random_pattern(m, n, per_row, seed)).unwrap();
        let perm = a.ordering(Ordering::Colamd).unwrap();
        assert!(is_permutation(&perm, n), "{m} × {n}");
        // The same in either orientation, and the same each time
        assert_eq!(a.to_csc().ordering(Ordering::Colamd).unwrap(), perm);
        assert_eq!(a.ordering(Ordering::Colamd).unwrap(), perm);
    }

    // Empty columns, a dense row and a dense column: the empty and dense columns go last
    let n = 300;
    let mut triplets: Vec<(usize, usize, f64)> = (0..n)
        .filter(|&i| i % 10 != 3)
        .map(|i| (i, i, 1.0))
        .collect();
    triplets.extend((0..n).map(|j| (5, j, 1.0))); // a dense row
    triplets.extend((0..n).map(|i| (i, 7, 1.0))); // a dense column
    triplets.extend((0..n - 1).map(|i| (i, i + 1, 2.0)));
    let a = CscMatrix::from_triplets(n, n, triplets).unwrap();
    let perm = a.ordering(Ordering::Colamd).unwrap();
    assert!(is_permutation(&perm, n));
    assert_eq!(*perm.last().unwrap(), 7, "the dense column is ordered last");

    // An empty matrix, and one with no entries
    assert_eq!(
        CscMatrix::<f64>::from_triplets(0, 0, vec![])
            .unwrap()
            .ordering(Ordering::Colamd)
            .unwrap(),
        Vec::<usize>::new()
    );
    let empty = CscMatrix::<f64>::from_triplets(4, 3, vec![]).unwrap();
    assert!(is_permutation(
        &empty.ordering(Ordering::Colamd).unwrap(),
        3
    ));
    // The other orderings still need a square matrix
    assert_eq!(
        empty.ordering(Ordering::Amd).unwrap_err(),
        SparseError::NotSquare { shape: (4, 3) }
    );
}

#[test]
fn lu_under_colamd_solves() {
    for (n, seed) in [(40, 5), (150, 6), (400, 7)] {
        let a = unsymmetric(n, seed);
        let lu = a
            .analyze_lu(Ordering::Colamd)
            .unwrap()
            .factor(&a, LuOptions::default())
            .unwrap();
        let truth: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
        let x = lu.solve(&a.mul_vec(&truth).unwrap()).unwrap();
        for (got, want) in x.iter().zip(&truth) {
            assert!((got - want).abs() <= 1e-10, "{got} vs {want}");
        }
        // Refactoring with new values reuses the analysis
        let mut lu = lu;
        lu.refactor(&a.map(|v| 1.5 * v)).unwrap();
    }

    // Complex
    let a: CscMatrix<Complex<f64>> = unsymmetric(60, 8).map(|v| Complex::new(v, 0.2 * v));
    let lu = a
        .analyze_lu(Ordering::Colamd)
        .unwrap()
        .factor(&a, LuOptions::default())
        .unwrap();
    let truth: Vec<Complex<f64>> = (0..60).map(|i| Complex::new(i as f64, -1.0)).collect();
    let x = lu.solve(&a.mul_vec(&truth).unwrap()).unwrap();
    assert!(x.iter().zip(&truth).all(|(g, w)| (g - w).norm() <= 1e-9));

    // A structurally singular matrix (an empty column) is reported
    let singular =
        CscMatrix::from_triplets(3, 3, vec![(0, 0, 1.0), (1, 0, 2.0), (2, 2, 1.0)]).unwrap();
    assert!(matches!(
        singular
            .analyze_lu(Ordering::Colamd)
            .unwrap()
            .factor(&singular, LuOptions::default()),
        Err(SparseError::SingularMatrix { .. })
    ));
}

/// The Cholesky factor's entries for `AᵀA` reordered by `perm`: what COLAMD minimizes
fn ata_fill(a: &CscMatrix<f64>, perm: &[usize]) -> usize {
    let at: CsrMatrix<f64> = a.clone().transpose();
    let ata = at.matmul(&a.to_csr()).unwrap();
    ata.permute_symmetric(perm)
        .unwrap()
        .analyze_cholesky(Ordering::Natural)
        .unwrap()
        .predicted_entries()
}

/// `a` with its rows in a pseudo-random order
fn shuffle_rows(a: &CscMatrix<f64>, seed: u64) -> CscMatrix<f64> {
    let n = a.shape().0;
    let mut rng = Lcg(seed);
    let mut perm: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        perm.swap(i, rng.next() as usize % (i + 1));
    }
    let triplets: Vec<(usize, usize, f64)> = a
        .iter()
        .map(|(row, col, value)| (perm[row], col, value))
        .collect();
    CscMatrix::from_triplets(n, a.shape().1, triplets).unwrap()
}

#[test]
fn colamd_orders_ata_as_well_as_amd_on_ata() {
    // COLAMD's aim: a sparse Cholesky factor of AᵀA, found without forming AᵀA. AMD on the
    //   AᵀA it forms explicitly is the yardstick
    let (mut colamd, mut amd) = (0, 0);
    for seed in 0..8 {
        let a = unsymmetric(300, 100 + seed);
        let at: CsrMatrix<f64> = a.clone().transpose();
        let ata = at.matmul(&a.to_csr()).unwrap();
        colamd += ata_fill(&a, &a.ordering(Ordering::Colamd).unwrap());
        amd += ata_fill(&a, &ata.ordering(Ordering::Amd).unwrap());
    }
    eprintln!("AᵀA Cholesky entries: COLAMD {colamd}, AMD on AᵀA {amd}");
    assert!(
        colamd as f64 <= 1.05 * amd as f64,
        "COLAMD {colamd} against {amd}"
    );
}

#[test]
fn colamd_beats_amd_when_rows_are_out_of_order() {
    // COLAMD depends on the rows' order only through ties, since AᵀA doesn't change; AMD on
    //   A + Aᵀ sees a scrambled pattern and a diagonal of mostly zeros. With a strong, zero-free
    //   diagonal in place, AMD (keeping the diagonal pivots) fills less: the usual reason to
    //   prefer it on structurally symmetric problems
    let (mut colamd, mut amd, mut colamd_in_order, mut amd_in_order) = (0, 0, 0, 0);
    for seed in 0..8 {
        let a = unsymmetric(300, 100 + seed);
        let b = shuffle_rows(&a, 900 + seed);
        colamd += fill(&b, Ordering::Colamd);
        amd += fill(&b, Ordering::Amd);
        colamd_in_order += fill(&a, Ordering::Colamd);
        amd_in_order += fill(&a, Ordering::Amd);
    }
    eprintln!(
        "LU entries, rows shuffled: COLAMD {colamd}, AMD {amd}; in order: COLAMD \
         {colamd_in_order}, AMD {amd_in_order}"
    );
    assert!(colamd < amd, "COLAMD {colamd} against AMD {amd}");
    // Close to the same fill either way round
    let ratio = colamd as f64 / colamd_in_order as f64;
    assert!((0.9..1.1).contains(&ratio), "{ratio}");
    assert!(amd_in_order < colamd_in_order);
}

#[test]
fn cholesky_refuses_colamd() {
    let a = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 2.0), (1, 1, 2.0)]).unwrap();
    assert!(matches!(
        a.analyze_cholesky(Ordering::Colamd),
        Err(SparseError::UnsupportedOrdering { .. })
    ));
    let wide = CscMatrix::<f64>::from_triplets(2, 3, vec![(0, 0, 1.0)]).unwrap();
    assert_eq!(
        wide.analyze_lu(Ordering::Colamd).unwrap_err(),
        SparseError::NotSquare { shape: (2, 3) }
    );
}
