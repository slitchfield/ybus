// Tests for Phase 3's factorization machinery: fill-reducing orderings, symbolic analysis and
//   the ordered factorizations.

use sparse_mat::compressed_matrix::{CsrMatrix, Ordering};
use sparse_mat::sparse_error::SparseError;

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

    /// A pseudo-random permutation of `0..n` (Fisher–Yates).
    fn permutation(&mut self, n: usize) -> Vec<usize> {
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            perm.swap(i, self.next() as usize % (i + 1));
        }
        perm
    }
}

/// The 5-point Laplacian of a `k x k` grid, in natural (row by row) order: SPD, with bandwidth
/// `k`.
fn grid_laplacian(k: usize) -> CsrMatrix<f64> {
    let mut triplets = vec![];
    for r in 0..k {
        for c in 0..k {
            let i = r * k + c;
            triplets.push((i, i, 4.0));
            if c + 1 < k {
                triplets.extend([(i, i + 1, -1.0), (i + 1, i, -1.0)]);
            }
            if r + 1 < k {
                triplets.extend([(i, i + k, -1.0), (i + k, i, -1.0)]);
            }
        }
    }
    CsrMatrix::from_triplets(k * k, k * k, triplets).unwrap()
}

/// An `n x n` banded SPD matrix with half-bandwidth `band`.
fn banded(n: usize, band: usize) -> CsrMatrix<f64> {
    let mut triplets = vec![];
    for i in 0..n {
        triplets.push((i, i, 2.0 * band as f64 + 1.0));
        for d in 1..=band.min(n - 1 - i) {
            triplets.extend([(i, i + d, -1.0), (i + d, i, -1.0)]);
        }
    }
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

fn is_permutation(perm: &[usize], n: usize) -> bool {
    let mut seen = vec![false; n];
    perm.len() == n
        && perm
            .iter()
            .all(|&i| i < n && !std::mem::replace(&mut seen[i], true))
}

// Bandwidth, profile and orderings

#[test]
fn bandwidth_and_profile() {
    assert_eq!(CsrMatrix::<f64>::identity(4).bandwidth(), 0);
    assert_eq!(CsrMatrix::<f64>::identity(4).profile(), 0);
    assert_eq!(CsrMatrix::<f64>::default().bandwidth(), 0);
    assert_eq!(banded(6, 2).bandwidth(), 2);
    // Rows 1.. reach two columns left, except row 1, which reaches one: 1 + 2 * 4
    assert_eq!(banded(6, 2).profile(), 9);
    assert_eq!(grid_laplacian(4).bandwidth(), 4);
    // A stored zero counts, and only the lower envelope enters the profile
    let a = CsrMatrix::from_triplets(3, 3, vec![(0, 0, 1.0), (2, 0, 0.0), (0, 2, 5.0)]).unwrap();
    assert_eq!(a.bandwidth(), 2);
    assert_eq!(a.profile(), 2);
    assert_eq!(a.to_csc().profile(), 2);
}

#[test]
fn natural_ordering_is_the_identity() {
    assert_eq!(
        grid_laplacian(3).ordering(Ordering::Natural).unwrap(),
        (0..9).collect::<Vec<_>>()
    );
}

#[test]
fn rcm_does_not_widen_banded_matrices() {
    for (n, band) in [(10, 1), (30, 2), (50, 5), (7, 6)] {
        let a = banded(n, band);
        let perm = a.ordering(Ordering::Rcm).unwrap();
        assert!(is_permutation(&perm, n));
        let reordered = a.permute_symmetric(&perm).unwrap();
        assert!(
            reordered.bandwidth() <= band,
            "n = {n}, band = {band}: {}",
            reordered.bandwidth()
        );
        assert!(reordered.profile() <= a.profile());
    }
}

#[test]
fn rcm_recovers_a_narrow_band_from_shuffled_grids() {
    let mut rng = Lcg(7);
    for k in [4, 10, 25] {
        let grid = grid_laplacian(k);
        let shuffled = grid.permute_symmetric(&rng.permutation(k * k)).unwrap();
        let perm = shuffled.ordering(Ordering::Rcm).unwrap();
        assert!(is_permutation(&perm, k * k));
        let reordered = shuffled.permute_symmetric(&perm).unwrap();
        // From a corner, the breadth-first levels are the grid's anti-diagonals, so RCM matches
        //   the natural order's bandwidth k, with a smaller envelope
        assert_eq!(
            reordered.bandwidth(),
            k,
            "k = {k} (shuffled: {})",
            shuffled.bandwidth()
        );
        assert!(reordered.profile() < grid.profile());
        assert!(reordered.profile() < shuffled.profile());
    }
}

#[test]
fn rcm_handles_components_and_isolated_nodes() {
    // Two paths, 0-2-4 and 1-3, plus the isolated node 5, interleaved
    let edges = [(0, 2), (2, 4), (1, 3)];
    let mut triplets: Vec<(usize, usize, f64)> = (0..6).map(|i| (i, i, 1.0)).collect();
    for (i, j) in edges {
        triplets.extend([(i, j, 1.0), (j, i, 1.0)]);
    }
    let a = CsrMatrix::from_triplets(6, 6, triplets).unwrap();
    let perm = a.ordering(Ordering::Rcm).unwrap();
    assert!(is_permutation(&perm, 6));
    // Each component is numbered contiguously, so every edge joins neighbours
    assert_eq!(a.permute_symmetric(&perm).unwrap().bandwidth(), 1);
}

#[test]
fn rcm_depends_only_on_the_pattern() {
    let mut rng = Lcg(3);
    let a = grid_laplacian(6)
        .permute_symmetric(&rng.permutation(36))
        .unwrap();
    let perm = a.ordering(Ordering::Rcm).unwrap();
    // The same in the other orientation, with other values, and with an unsymmetric pattern
    //   whose symmetrization is the same
    assert_eq!(a.to_csc().ordering(Ordering::Rcm).unwrap(), perm);
    assert_eq!(
        a.map(|v| v * 3.0 + 1.0).ordering(Ordering::Rcm).unwrap(),
        perm
    );
    let lower: Vec<(usize, usize, f64)> = a.iter().filter(|&(i, j, _)| j <= i).collect();
    let lower = CsrMatrix::from_triplets(36, 36, lower).unwrap();
    assert_eq!(lower.ordering(Ordering::Rcm).unwrap(), perm);
    // And deterministic
    assert_eq!(a.ordering(Ordering::Rcm).unwrap(), perm);
}

#[test]
fn ordering_edge_cases() {
    assert!(CsrMatrix::<f64>::default()
        .ordering(Ordering::Rcm)
        .unwrap()
        .is_empty());
    assert_eq!(
        CsrMatrix::<f64>::identity(1)
            .ordering(Ordering::Rcm)
            .unwrap(),
        [0]
    );
    let rectangular = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        rectangular.ordering(Ordering::Rcm),
        Err(SparseError::NotSquare { shape: (2, 3) })
    );
}
