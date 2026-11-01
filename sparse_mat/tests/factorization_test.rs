// Tests for Phase 3's factorization machinery: fill-reducing orderings, symbolic analysis and
//   the ordered factorizations.

mod common;

use sparse_mat::compressed_matrix::{CsrMatrix, Ordering};
use sparse_mat::sparse_error::SparseError;

use common::{reference_cholesky, reference_solve_spd};

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
                reference_cholesky(&reordered).unwrap().num_stored(),
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

// Ordered Cholesky

/// Whether two matrices hold the same values bit for bit, except that a zero may be stored in
/// one and absent from (or of the other sign in) the other.
fn same_bits(a: &CsrMatrix<f64>, b: &CsrMatrix<f64>) -> bool {
    a.to_dense()
        .iter()
        .flatten()
        .zip(b.to_dense().iter().flatten())
        .all(|(x, y)| x == y && (*x == 0.0 || x.to_bits() == y.to_bits()))
}

fn assert_close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!((a - e).abs() <= tolerance, "entry {i}: {a} vs {e}");
    }
}

#[test]
fn ordered_cholesky_reconstructs_the_reordered_matrix() {
    let mut rng = Lcg(11);
    let shuffled_grid = grid_laplacian(8)
        .permute_symmetric(&rng.permutation(64))
        .unwrap();
    for a in [random_spd(1, 0, 1), random_spd(25, 15, 2), shuffled_grid] {
        for ordering in [Ordering::Natural, Ordering::Rcm] {
            let symbolic = a.analyze_cholesky(ordering).unwrap();
            let factor = symbolic.factor(&a).unwrap();
            let l = factor.l();
            assert_eq!(factor.permutation(), symbolic.permutation());
            assert_eq!(l.num_stored(), symbolic.predicted_entries());
            assert!(l.iter().all(|(row, col, _)| col <= row));
            assert!(l.diagonal().values().iter().all(|&d| d > 0.0));
            let reordered = a.permute_symmetric(factor.permutation()).unwrap();
            let product = l * &l.clone().transpose().to_csr();
            for (got, want) in product.to_dense().iter().zip(reordered.to_dense()) {
                assert_close(got, &want, 1e-12);
            }
        }
    }
}

#[test]
fn natural_ordering_matches_the_unordered_cholesky_bit_for_bit() {
    let mut rng = Lcg(5);
    let shuffled_grid = grid_laplacian(9)
        .permute_symmetric(&rng.permutation(81))
        .unwrap();
    for a in [random_spd(40, 10, 6), random_spd(30, 40, 7), shuffled_grid] {
        let n = a.shape().0;
        let factor = a
            .analyze_cholesky(Ordering::Natural)
            .unwrap()
            .factor(&a)
            .unwrap();
        assert!(same_bits(factor.l(), &reference_cholesky(&a).unwrap()));
        let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.7).sin()).collect();
        let x = factor.solve(&b).unwrap();
        let reference = reference_solve_spd(&a, &b).unwrap();
        assert!(x
            .iter()
            .zip(&reference)
            .all(|(x, r)| x.to_bits() == r.to_bits()));
    }
}

#[test]
fn ordered_solves_match_the_true_solution() {
    let mut rng = Lcg(9);
    let a = grid_laplacian(12)
        .permute_symmetric(&rng.permutation(144))
        .unwrap();
    let x: Vec<f64> = (0..144).map(|i| (i as f64).cos()).collect();
    let b = a.mul_vec(&x).unwrap();
    assert_close(&a.solve_spd(&b).unwrap(), &x, 1e-12);
    assert_close(&a.cholesky().unwrap().solve(&b).unwrap(), &x, 1e-12);
    // RCM cuts the shuffled grid's fill (2,091 entries) below even that of the grid's own row
    //   by row numbering (1,739), whose factor fills its band
    let natural = a.analyze_cholesky(Ordering::Natural).unwrap();
    let rcm = a.analyze_cholesky(Ordering::Rcm).unwrap();
    let row_by_row = grid_laplacian(12)
        .analyze_cholesky(Ordering::Natural)
        .unwrap();
    assert!(rcm.predicted_entries() < row_by_row.predicted_entries());
    assert!(row_by_row.predicted_entries() < natural.predicted_entries());
    assert_eq!(
        a.cholesky().unwrap().solve(&b[..3]).err(),
        Some(SparseError::ShapeMismatch {
            left: (144, 144),
            right: (3, 1)
        })
    );
}

#[test]
fn one_analysis_factors_every_matrix_with_its_pattern() {
    let a = random_spd(20, 20, 3);
    let symbolic = a.analyze_cholesky(Ordering::Rcm).unwrap();
    // Same pattern, new values (here: a stored zero appears in place of a value)
    let mut b = a.clone();
    for (k, value) in b.values_mut().iter_mut().enumerate() {
        if k % 7 == 3 && *value < 0.0 {
            *value = 0.0;
        }
    }
    let mut fixed = b.clone();
    // Keep the matrix symmetric: zero the mirror images too
    for (row, col, value) in b.iter() {
        if value == 0.0 {
            let p = fixed.position(col, row).unwrap().unwrap();
            fixed.values_mut()[p] = 0.0;
        }
    }
    let factor = symbolic.factor(&fixed).unwrap();
    assert_eq!(factor.l().num_stored(), symbolic.predicted_entries());
    let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
    assert_close(
        &factor.solve(&fixed.mul_vec(&x).unwrap()).unwrap(),
        &x,
        1e-10,
    );

    // Another pattern is refused, even with the same shape and count
    let other = random_spd(20, 20, 4);
    assert_eq!(
        symbolic.factor(&other).err(),
        Some(SparseError::PatternMismatch)
    );
}

#[test]
fn cholesky_errors() {
    let not_symmetric = CsrMatrix::from_dense(&[[2.0, 1.0], [0.0, 2.0]]).unwrap();
    assert_eq!(
        not_symmetric.cholesky().err(),
        Some(SparseError::NotSymmetric)
    );
    // Symmetry is mathematical: a stored zero opposite a missing entry doesn't break it
    let stored_zero =
        CsrMatrix::from_triplets(2, 2, vec![(0, 0, 4.0), (0, 1, 0.0), (1, 1, 9.0)]).unwrap();
    let factor = stored_zero.cholesky().unwrap();
    assert_eq!(
        factor.l().num_stored(),
        3,
        "the stored zero's position is kept"
    );
    let rect = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        rect.cholesky().err(),
        Some(SparseError::NotSquare { shape: (2, 3) })
    );

    // Indefinite (eigenvalues 3 and -1): the failing pivot is named by its row of A, so the
    //   ordering decides which row that is
    let indefinite = CsrMatrix::from_dense(&[[1.0, 2.0], [2.0, 1.0]]).unwrap();
    for ordering in [Ordering::Natural, Ordering::Rcm] {
        let symbolic = indefinite.analyze_cholesky(ordering).unwrap();
        assert_eq!(
            symbolic.factor(&indefinite).err(),
            Some(SparseError::NotPositiveDefinite {
                pivot: symbolic.permutation()[1]
            })
        );
    }
    // Here only row 1 has a non-positive diagonal, wherever it's ordered
    let negative = CsrMatrix::from_diagonal(vec![1.0, -1.0, 2.0]);
    assert_eq!(
        negative.cholesky().err(),
        Some(SparseError::NotPositiveDefinite { pivot: 1 })
    );
    let zero = CsrMatrix::<f64>::from_triplets(2, 2, Vec::new()).unwrap();
    assert!(matches!(
        zero.cholesky().err(),
        Some(SparseError::NotPositiveDefinite { .. })
    ));
}

#[test]
fn cholesky_keeps_a_diagonal_matrix_diagonal() {
    let factor = CsrMatrix::from_diagonal(vec![4.0, 9.0, 16.0])
        .cholesky()
        .unwrap();
    assert_eq!(factor.l().num_stored(), 3);
    let mut diagonal = vec![0.0; 3];
    for (i, &row) in factor.permutation().iter().enumerate() {
        diagonal[row] = factor.l().diagonal().values()[i];
    }
    assert_eq!(diagonal, [2.0, 3.0, 4.0]);
}

// Approximate minimum degree

/// The SPD matrix of a graph's Laplacian plus the identity, from its edges.
fn graph_matrix(n: usize, edges: &[(usize, usize)]) -> CsrMatrix<f64> {
    let mut degree = vec![0.0; n];
    let mut triplets = vec![];
    for &(i, j) in edges {
        triplets.extend([(i, j, -1.0), (j, i, -1.0)]);
        degree[i] += 1.0;
        degree[j] += 1.0;
    }
    triplets.extend(degree.iter().enumerate().map(|(i, d)| (i, i, d + 1.0)));
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// Entries of a Cholesky factor that the matrix's lower triangle doesn't already hold.
fn fill_in(a: &CsrMatrix<f64>, ordering: Ordering) -> usize {
    let lower = a.iter().filter(|&(i, j, _)| j <= i).count();
    a.analyze_cholesky(ordering).unwrap().predicted_entries() - lower
}

#[test]
fn amd_returns_a_permutation_that_depends_only_on_the_pattern() {
    let mut rng = Lcg(21);
    for (n, fill, seed) in [(1, 0, 1), (2, 100, 2), (17, 30, 3), (80, 5, 4), (150, 2, 5)] {
        let a = random_spd(n, fill, seed);
        let perm = a.ordering(Ordering::Amd).unwrap();
        assert!(is_permutation(&perm, n), "n = {n}");
        assert_eq!(a.to_csc().ordering(Ordering::Amd).unwrap(), perm);
        assert_eq!(a.map(|v| v * -2.0).ordering(Ordering::Amd).unwrap(), perm);
        assert_eq!(a.ordering(Ordering::Amd).unwrap(), perm, "deterministic");
    }
    let shuffled = grid_laplacian(9)
        .permute_symmetric(&rng.permutation(81))
        .unwrap();
    assert!(is_permutation(
        &shuffled.ordering(Ordering::Amd).unwrap(),
        81
    ));
    assert!(CsrMatrix::<f64>::default()
        .ordering(Ordering::Amd)
        .unwrap()
        .is_empty());
}

#[test]
fn amd_orders_trees_and_arrowheads_without_fill() {
    // Random trees: each node joins an earlier one, then the numbering is shuffled
    let mut rng = Lcg(8);
    for n in [10, 100, 1000] {
        let parent_edges: Vec<(usize, usize)> =
            (1..n).map(|i| (i, rng.next() as usize % i)).collect();
        let tree = graph_matrix(n, &parent_edges)
            .permute_symmetric(&rng.permutation(n))
            .unwrap();
        assert_eq!(fill_in(&tree, Ordering::Amd), 0, "tree of {n}");
    }
    // Arrowheads: a hub joined to every other node. Small, the hub waits until only one leaf is
    //   left (then either may go last); large, it's a dense row, set aside to come last
    for (n, hub) in [(12, 5), (500, 123)] {
        let edges: Vec<(usize, usize)> = (0..n).filter(|&i| i != hub).map(|i| (hub, i)).collect();
        let arrow = graph_matrix(n, &edges);
        let perm = arrow.ordering(Ordering::Amd).unwrap();
        assert!(perm[n - 2..].contains(&hub), "n = {n}");
        if n > 160 {
            assert_eq!(perm.last(), Some(&hub), "n = {n}");
        }
        assert_eq!(fill_in(&arrow, Ordering::Amd), 0);
        assert!(fill_in(&arrow, Ordering::Natural) > 0);
    }
}

#[test]
fn amd_sets_dense_rows_aside() {
    // A 30 x 30 grid (900 nodes; dense means more than 300 neighbours) with two extra nodes
    //   joined to every grid node
    let k = 30;
    let mut edges = vec![];
    for r in 0..k {
        for c in 0..k {
            let i = r * k + c;
            if c + 1 < k {
                edges.push((i, i + 1));
            }
            if r + 1 < k {
                edges.push((i, i + k));
            }
            edges.extend([(i, 900), (i, 901)]);
        }
    }
    let a = graph_matrix(902, &edges);
    let perm = a.ordering(Ordering::Amd).unwrap();
    assert!(is_permutation(&perm, 902));
    assert_eq!(&perm[900..], [900, 901]);
}

#[test]
fn amd_fills_grids_less_than_rcm() {
    let mut rng = Lcg(13);
    for k in [10, 20, 40] {
        let grid = grid_laplacian(k)
            .permute_symmetric(&rng.permutation(k * k))
            .unwrap();
        let (amd, rcm) = (fill_in(&grid, Ordering::Amd), fill_in(&grid, Ordering::Rcm));
        assert!(amd < rcm, "k = {k}: AMD fill {amd}, RCM fill {rcm}");
    }
}

#[test]
fn amd_factors_solve_correctly() {
    let mut rng = Lcg(17);
    let a = grid_laplacian(15)
        .permute_symmetric(&rng.permutation(225))
        .unwrap();
    let factor = a
        .analyze_cholesky(Ordering::Amd)
        .unwrap()
        .factor(&a)
        .unwrap();
    let x: Vec<f64> = (0..225).map(|i| (i as f64 * 0.3).sin()).collect();
    assert_close(&factor.solve(&a.mul_vec(&x).unwrap()).unwrap(), &x, 1e-12);
}
