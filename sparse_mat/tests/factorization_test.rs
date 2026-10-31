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

// Symbolic analysis for Cholesky

/// A pseudo-random symmetric positive definite matrix: about `fill_percent` of the
/// off-diagonal positions in a symmetric pattern, with a dominant diagonal.
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

/// The pattern of the Cholesky factor of `perm`'s reordering of `a`'s symmetrized pattern, by
/// dense symbolic elimination: `filled[i][j]` for `j <= i`.
fn brute_force_factor_pattern(a: &CsrMatrix<f64>, perm: &[usize]) -> Vec<Vec<bool>> {
    let n = perm.len();
    let mut inverse = vec![0; n];
    for (new, &old) in perm.iter().enumerate() {
        inverse[old] = new;
    }
    let mut filled = vec![vec![false; n]; n];
    for (i, j, _) in a.iter() {
        let (i, j) = (inverse[i], inverse[j]);
        filled[i.max(j)][i.min(j)] = true;
    }
    for (i, row) in filled.iter_mut().enumerate() {
        row[i] = true;
    }
    // Eliminating k connects every pair of later rows with an entry in column k
    for k in 0..n {
        let below: Vec<usize> = (k + 1..n).filter(|&i| filled[i][k]).collect();
        for (a, &i) in below.iter().enumerate() {
            for &j in &below[..=a] {
                filled[i][j] = true;
            }
        }
    }
    filled
}

#[test]
fn analysis_matches_brute_force_elimination() {
    for (n, fill, seed) in [(1, 0, 1), (12, 20, 2), (40, 8, 3), (60, 4, 4), (30, 50, 5)] {
        let a = random_spd(n, fill, seed);
        for ordering in [Ordering::Natural, Ordering::Rcm] {
            let symbolic = a.analyze_cholesky(ordering).unwrap();
            assert_eq!(symbolic.dim(), n);
            assert_eq!(symbolic.permutation(), a.ordering(ordering).unwrap());
            let filled = brute_force_factor_pattern(&a, symbolic.permutation());

            // The parent of j is the first row below j with an entry in column j
            let parent: Vec<Option<usize>> =
                (0..n).map(|j| (j + 1..n).find(|&i| filled[i][j])).collect();
            assert_eq!(symbolic.elimination_tree(), parent, "n = {n}, {ordering:?}");
            let counts: Vec<usize> = (0..n)
                .map(|j| (j..n).filter(|&i| filled[i][j]).count())
                .collect();
            assert_eq!(symbolic.column_counts(), counts, "n = {n}, {ordering:?}");
            assert_eq!(symbolic.predicted_entries(), counts.iter().sum::<usize>());

            // With values that never cancel exactly, the factor stores the whole pattern
            let reordered = a.permute_symmetric(symbolic.permutation()).unwrap();
            assert_eq!(
                reordered.cholesky().unwrap().num_stored(),
                symbolic.predicted_entries()
            );
        }
    }
}

#[test]
fn postorder_visits_descendants_first() {
    let a = grid_laplacian(7);
    for ordering in [Ordering::Natural, Ordering::Rcm] {
        let symbolic = a.analyze_cholesky(ordering).unwrap();
        let order = symbolic.postorder();
        assert!(is_permutation(order, 49));
        let mut position = vec![0; 49];
        for (k, &node) in order.iter().enumerate() {
            position[node] = k;
        }
        // Every node follows its children, and each subtree is contiguous: a node's subtree
        //   occupies the positions just before it, as many as the subtree has nodes
        let tree = symbolic.elimination_tree();
        let mut size = vec![1; 49];
        for node in 0..49 {
            if let Some(parent) = tree[node] {
                assert!(position[node] < position[parent]);
                assert!(parent > node);
                size[parent] += size[node];
            }
        }
        for node in 0..49 {
            let first = position[node] + 1 - size[node];
            for &member in &order[first..=position[node]] {
                let mut ancestor = Some(member);
                while ancestor.is_some_and(|a| a != node) {
                    ancestor = tree[ancestor.unwrap()];
                }
                assert_eq!(ancestor, Some(node), "{member} is not under {node}");
            }
        }
    }
}

#[test]
fn analysis_uses_the_symmetrized_pattern_and_counts_stored_zeros() {
    // Only the lower triangle stored, one entry of it a stored zero
    let a = CsrMatrix::from_triplets(
        3,
        3,
        vec![
            (0, 0, 4.0),
            (1, 1, 4.0),
            (2, 2, 4.0),
            (2, 0, 0.0),
            (2, 1, 1.0),
        ],
    )
    .unwrap();
    let symbolic = a.analyze_cholesky(Ordering::Natural).unwrap();
    assert_eq!(symbolic.elimination_tree(), [Some(2), Some(2), None]);
    assert_eq!(symbolic.column_counts(), [2, 2, 1]);
    // A forest: two components, so two roots
    let forest = CsrMatrix::<f64>::identity(3)
        .analyze_cholesky(Ordering::Rcm)
        .unwrap();
    assert_eq!(forest.elimination_tree(), [None, None, None]);
    assert_eq!(forest.predicted_entries(), 3);

    assert_eq!(
        CsrMatrix::<f64>::from_triplets(2, 3, Vec::new())
            .unwrap()
            .analyze_cholesky(Ordering::Rcm)
            .err(),
        Some(SparseError::NotSquare { shape: (2, 3) })
    );
    assert_eq!(
        CsrMatrix::<f64>::default()
            .analyze_cholesky(Ordering::Rcm)
            .unwrap()
            .dim(),
        0
    );
}
