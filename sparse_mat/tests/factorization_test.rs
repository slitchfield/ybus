// Tests for Phase 3's factorization machinery: fill-reducing orderings, symbolic analysis and
//   the ordered factorizations.

mod common;

use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix, LuOptions, Ordering, RefactorOutcome};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::Complex;

use common::{reference_cholesky, reference_lu, reference_solve_spd};

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

// LU

/// A pseudo-random square matrix with about `fill_percent` of its off-diagonal positions
/// stored, and a diagonal that is `diagonal_scale` times its row's largest entry, so a small
/// scale forces pivoting off the diagonal.
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

/// The dense product `L U` of an LU factorization, against the dense `P A Q`.
fn assert_lu_reconstructs(a: &CscMatrix<f64>, lu: &sparse_mat::compressed_matrix::Lu<f64>) {
    let l = lu.l();
    let u = lu.u();
    assert!(l.iter().all(|(row, col, _)| row >= col));
    assert!(u.iter().all(|(row, col, _)| row <= col));
    assert!(l.diagonal().values().iter().all(|&d| d == 1.0));
    let paq = a
        .permute(lu.permutation(), lu.column_permutation())
        .unwrap();
    let product = (l * u).to_dense();
    for (got, want) in product.iter().zip(paq.to_dense()) {
        assert_close(got, &want, 1e-10);
    }
}

#[test]
fn lu_reconstructs_the_permuted_matrix() {
    for (n, fill, scale, seed) in [
        (1, 0, 1.0, 1),
        (12, 30, 0.1, 2),
        (40, 10, 2.0, 3),
        (60, 6, 0.01, 4),
    ] {
        let a = random_square(n, fill, scale, seed).to_csc();
        for ordering in [Ordering::Natural, Ordering::Rcm, Ordering::Amd] {
            for threshold in [0.001, 0.1, 1.0] {
                let lu = a
                    .analyze_lu(ordering)
                    .unwrap()
                    .factor(
                        &a,
                        LuOptions {
                            pivot_threshold: threshold,
                        },
                    )
                    .unwrap();
                assert_eq!(lu.column_permutation(), a.ordering(ordering).unwrap());
                assert!(is_permutation(lu.permutation(), n));
                assert_lu_reconstructs(&a, &lu);
            }
        }
    }
}

#[test]
fn lu_solves_like_the_reference_lu() {
    for (n, fill, scale, seed) in [(30, 15, 0.05, 5), (80, 4, 1.5, 6), (120, 3, 0.02, 7)] {
        let csr = random_square(n, fill, scale, seed);
        let a = csr.to_csc();
        let x: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).cos()).collect();
        let b = a.mul_vec(&x).unwrap();
        let reference = reference_lu(&csr).unwrap();
        let permuted: Vec<f64> = reference.perm.iter().map(|&row| b[row]).collect();
        let y = reference.l.solve_lower_triangular(&permuted).unwrap();
        let expected = reference.u.solve_upper_triangular(&y).unwrap();
        let solved = a.solve(&b).unwrap();
        assert_close(&solved, &expected, 1e-9);
        assert_close(&solved, &x, 1e-9);
        // Under true partial pivoting with no ordering, both choose the same pivots
        let partial = a
            .analyze_lu(Ordering::Natural)
            .unwrap()
            .factor(
                &a,
                LuOptions {
                    pivot_threshold: 1.0,
                },
            )
            .unwrap();
        assert_eq!(partial.permutation(), reference.perm, "n = {n}");
    }
}

#[test]
fn lu_threshold_decides_between_diagonal_and_largest() {
    // Column 0's candidates are 0.5 on the diagonal and 1.0 below it
    let a = CscMatrix::from_dense(&[[0.5, 1.0], [1.0, 1.0]]).unwrap();
    let symbolic = a.analyze_lu(Ordering::Natural).unwrap();
    let with = |pivot_threshold| symbolic.factor(&a, LuOptions { pivot_threshold }).unwrap();
    assert_eq!(
        with(0.1).permutation(),
        [0, 1],
        "0.5 >= 0.1 x 1.0 keeps the diagonal"
    );
    assert_eq!(with(0.5).permutation(), [0, 1]);
    assert_eq!(
        with(0.6).permutation(),
        [1, 0],
        "0.5 < 0.6 x 1.0 takes the largest"
    );
    assert_eq!(with(1.0).permutation(), [1, 0]);
    assert_eq!(LuOptions::default().pivot_threshold, 0.1);
}

#[test]
#[should_panic(expected = "the pivot threshold must be in (0, 1]")]
fn lu_rejects_a_threshold_outside_the_unit_interval() {
    let a = CscMatrix::<f64>::identity(2);
    let _ = a.analyze_lu(Ordering::Natural).unwrap().factor(
        &a,
        LuOptions {
            pivot_threshold: 0.0,
        },
    );
}

#[test]
fn lu_factors_complex_matrices() {
    // A Ybus-like pattern: a shuffled grid's Laplacian with complex admittances, plus a
    //   diagonal too small to pivot on in some rows
    let mut rng = Lcg(31);
    let k = 7;
    let grid = grid_laplacian(k)
        .permute_symmetric(&rng.permutation(k * k))
        .unwrap();
    let a: CscMatrix<Complex<f64>> = grid
        .to_csc()
        .map(|v| Complex::new(v * 0.3, -v * (1.0 + (v.abs() * 7.0) % 3.0)));
    let x: Vec<Complex<f64>> = (0..k * k)
        .map(|i| Complex::new((i as f64).sin(), (i as f64 * 0.5).cos()))
        .collect();
    let b = a.mul_vec(&x).unwrap();
    let lu = a.lu().unwrap();
    let solved = lu.solve(&b).unwrap();
    for (got, want) in solved.iter().zip(&x) {
        assert!((got - want).norm() < 1e-10, "{got} vs {want}");
    }
    // An unsymmetric complex matrix that needs pivoting
    let c = CscMatrix::from_dense(&[
        [
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 1.0),
            Complex::new(2.0, 0.0),
        ],
        [
            Complex::new(0.0, 3.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, -1.0),
        ],
        [
            Complex::new(1.0, 0.0),
            Complex::new(2.0, 2.0),
            Complex::new(0.0, 0.0),
        ],
    ])
    .unwrap();
    let lu = c.lu().unwrap();
    let rhs = [
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 1.0),
        Complex::new(2.0, -1.0),
    ];
    let solution = lu.solve(&rhs).unwrap();
    let residual = c.mul_vec(&solution).unwrap();
    for (r, b) in residual.iter().zip(&rhs) {
        assert!((r - b).norm() < 1e-12);
    }
    // det = 0(0 - 2(1 - i)) - (1 + i)(0 - (1 - i)) + 2(3i(2 + 2i) - 0) = 2 - 12 + 12i
    let det = lu.determinant();
    assert!((det - Complex::new(-10.0, 12.0)).norm() < 1e-12, "{det}");
}

#[test]
fn lu_determinant_follows_both_permutations() {
    let a = CscMatrix::from_dense(&[[4.0f64, 3.0], [6.0, 3.0]]).unwrap();
    assert!((a.lu().unwrap().determinant() - (-6.0)).abs() < 1e-12);
    // A row swap is needed: det [[0, 2], [3, 1]] = -6
    let b = CscMatrix::from_dense(&[[0.0f64, 2.0], [3.0, 1.0]]).unwrap();
    assert!((b.lu().unwrap().determinant() - (-6.0)).abs() < 1e-12);
    assert_eq!(
        CscMatrix::<f64>::identity(5).lu().unwrap().determinant(),
        1.0
    );
    assert_eq!(CscMatrix::<f64>::default().lu().unwrap().determinant(), 1.0);
    // Against the reference on matrices whose orderings and pivots permute both ways
    for (n, seed) in [(9, 8), (25, 9)] {
        let csr = random_square(n, 25, 0.05, seed);
        let reference = reference_lu(&csr).unwrap();
        let sign = if reference.swaps.is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        let expected = reference
            .u
            .diagonal()
            .values()
            .iter()
            .fold(sign, |d, &u| d * u);
        let det = csr.to_csc().lu().unwrap().determinant();
        assert!(
            ((det - expected) / expected).abs() < 1e-10,
            "{det} vs {expected}"
        );
    }
}

#[test]
fn lu_keeps_its_symbolic_pattern_and_stored_zeros() {
    // A stored zero in A: its column still has a pivot, and the position is in the pattern
    let a = CscMatrix::from_triplets(2, 2, vec![(0, 0, 2.0), (1, 0, 0.0), (1, 1, 3.0)]).unwrap();
    let lu = a
        .analyze_lu(Ordering::Natural)
        .unwrap()
        .factor(&a, LuOptions::default())
        .unwrap();
    assert_eq!(
        lu.l().num_stored(),
        3,
        "the multiplier of the stored zero is kept"
    );
    assert_eq!(lu.u().num_stored(), 2);
    assert_close(&lu.solve(&[2.0, 3.0]).unwrap(), &[1.0, 1.0], 1e-12);
}

#[test]
fn lu_errors() {
    // Singular: the error names the column of A that has no pivot
    let singular = CscMatrix::from_dense(&[[1.0, 2.0], [2.0, 4.0]]).unwrap();
    let natural = singular.analyze_lu(Ordering::Natural).unwrap();
    assert_eq!(
        natural.factor(&singular, LuOptions::default()).err(),
        Some(SparseError::SingularMatrix { pivot: 1 })
    );
    let empty_column = CscMatrix::from_dense(&[[1.0, 0.0], [1.0, 0.0]]).unwrap();
    assert_eq!(
        empty_column.lu().err(),
        Some(SparseError::SingularMatrix { pivot: 1 })
    );
    let rect = CscMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        rect.lu().err(),
        Some(SparseError::NotSquare { shape: (2, 3) })
    );
    assert_eq!(
        rect.solve(&[1.0, 2.0]).err(),
        Some(SparseError::NotSquare { shape: (2, 3) })
    );
    let identity = CscMatrix::<f64>::identity(3);
    assert_eq!(
        identity.lu().unwrap().solve(&[1.0]).err(),
        Some(SparseError::ShapeMismatch {
            left: (3, 3),
            right: (1, 1)
        })
    );
    assert!(identity.solve(&[1.0]).is_err());
    // Another pattern is refused
    let symbolic = identity.analyze_lu(Ordering::Amd).unwrap();
    assert_eq!(
        symbolic.factor(&singular, LuOptions::default()).err(),
        Some(SparseError::PatternMismatch)
    );
}

// Refactoring

/// `a` with every stored value changed by a pseudo-random factor in [0.5, 1.5), keeping the
/// pattern (stored zeros stay zero).
fn new_values(a: &CscMatrix<f64>, seed: u64) -> CscMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut b = a.clone();
    for value in b.values_mut() {
        *value *= 0.5 + (rng.next() % 1000) as f64 / 1000.0;
    }
    b
}

#[test]
fn refactor_reuses_the_pivots_for_new_values() {
    for (n, fill, scale, seed) in [(20, 20, 1.5, 1), (90, 4, 3.0, 2)] {
        let a = random_square(n, fill, scale, seed).to_csc();
        let mut lu = a.lu().unwrap();
        let (l_pattern, u_pattern) = (lu.l().indices().to_vec(), lu.u().indices().to_vec());
        let permutations = (lu.permutation().to_vec(), lu.column_permutation().to_vec());
        for round in 0..3 {
            let b = new_values(&a, 10 + round);
            assert_eq!(lu.refactor(&b), Ok(RefactorOutcome::Reused), "n = {n}");
            // Same pivots and patterns, new values that factor b
            assert_eq!(lu.l().indices(), l_pattern);
            assert_eq!(lu.u().indices(), u_pattern);
            assert_eq!(lu.permutation(), permutations.0);
            assert_eq!(lu.column_permutation(), permutations.1);
            assert_lu_reconstructs(&b, &lu);
            let x: Vec<f64> = (0..n).map(|i| 1.0 + (i % 5) as f64).collect();
            let rhs = b.mul_vec(&x).unwrap();
            assert_close(&lu.solve(&rhs).unwrap(), &x, 1e-9);
            assert_close(&lu.solve(&rhs).unwrap(), &b.solve(&rhs).unwrap(), 1e-9);
            let det = lu.determinant();
            let fresh = b.lu().unwrap().determinant();
            assert!(((det - fresh) / fresh).abs() < 1e-9, "{det} vs {fresh}");
        }
    }
}

#[test]
fn refactor_repivots_when_a_pivot_becomes_too_small() {
    let a = CscMatrix::from_dense(&[[4.0, 1.0, 0.0], [1.0, 4.0, 1.0], [0.0, 1.0, 4.0]]).unwrap();
    let mut lu = a
        .analyze_lu(Ordering::Natural)
        .unwrap()
        .factor(&a, LuOptions::default())
        .unwrap();
    assert_eq!(lu.permutation(), [0, 1, 2]);
    // The first pivot shrinks below 0.1 of its column's largest entry
    let mut b = a.clone();
    let p = b.position(0, 0).unwrap().unwrap();
    b.values_mut()[p] = 0.01;
    assert_eq!(lu.refactor(&b), Ok(RefactorOutcome::Repivoted));
    assert_eq!(lu.permutation()[0], 1, "row 1 pivots first now");
    assert_lu_reconstructs(&b, &lu);
    // Zeroing (0, 0) no longer touches a pivot, since row 1 pivots column 0 now
    b.values_mut()[p] = 0.0;
    assert_eq!(lu.refactor(&b), Ok(RefactorOutcome::Reused));
    // From the original pivots, an exactly zero pivot can't be reused
    let mut fresh = a
        .analyze_lu(Ordering::Natural)
        .unwrap()
        .factor(&a, LuOptions::default())
        .unwrap();
    assert_eq!(fresh.refactor(&b), Ok(RefactorOutcome::Repivoted));
    assert_lu_reconstructs(&b, &fresh);
    assert_close(
        &lu.solve(&[1.0, 2.0, 3.0]).unwrap(),
        &b.solve(&[1.0, 2.0, 3.0]).unwrap(),
        1e-12,
    );
    // Back to the original values, the pivots chosen for b fail in column 1 (1 against
    //   0.1 x 15), so it repivots, back to the diagonal
    assert_eq!(lu.refactor(&a), Ok(RefactorOutcome::Repivoted));
    assert_eq!(lu.permutation(), [0, 1, 2]);
    assert_lu_reconstructs(&a, &lu);
}

#[test]
fn refactor_handles_complex_values() {
    let mut rng = Lcg(41);
    let grid = grid_laplacian(5)
        .permute_symmetric(&rng.permutation(25))
        .unwrap()
        .to_csc();
    let a = grid.map(|v| Complex::new(v, -2.0 * v));
    let b = grid.map(|v| Complex::new(0.7 * v, -1.1 * v));
    let mut lu = a.lu().unwrap();
    assert_eq!(lu.refactor(&b), Ok(RefactorOutcome::Reused));
    let x: Vec<Complex<f64>> = (0..25).map(|i| Complex::new(i as f64, 1.0)).collect();
    let solved = lu.solve(&b.mul_vec(&x).unwrap()).unwrap();
    assert!(solved.iter().zip(&x).all(|(s, x)| (s - x).norm() < 1e-10));
}

#[test]
fn refactor_errors_leave_the_factorization_unchanged() {
    let a = CscMatrix::from_dense(&[[2.0, 1.0], [1.0, 3.0]]).unwrap();
    let mut lu = a.lu().unwrap();
    let before = (lu.l().values().to_vec(), lu.u().values().to_vec());
    // Another pattern
    let other = CscMatrix::from_dense(&[[2.0, 0.0], [1.0, 3.0]]).unwrap();
    assert_eq!(lu.refactor(&other), Err(SparseError::PatternMismatch));
    // The same pattern, but singular: refactoring falls back to a fresh factorization, which
    //   fails
    let singular = CscMatrix::from_dense(&[[1.0, 2.0], [2.0, 4.0]]).unwrap();
    assert!(matches!(
        lu.refactor(&singular),
        Err(SparseError::SingularMatrix { .. })
    ));
    assert_eq!((lu.l().values().to_vec(), lu.u().values().to_vec()), before);
    assert_close(&lu.solve(&[3.0, 4.0]).unwrap(), &[1.0, 1.0], 1e-12);
}

// Several right-hand sides

fn right_hand_sides(n: usize, count: usize) -> Vec<Vec<f64>> {
    (0..count)
        .map(|k| {
            (0..n)
                .map(|i| ((i * 7 + k * 13) as f64 * 0.21).sin())
                .collect()
        })
        .collect()
}

#[test]
fn cholesky_solves_in_place_and_many_like_solve() {
    let mut rng = Lcg(51);
    let a = grid_laplacian(9)
        .permute_symmetric(&rng.permutation(81))
        .unwrap();
    let factor = a.cholesky().unwrap();
    let rhs = right_hand_sides(81, 5);
    let solved: Vec<Vec<f64>> = rhs.iter().map(|b| factor.solve(b).unwrap()).collect();
    assert_eq!(factor.solve_many(&rhs).unwrap(), solved, "bit for bit");
    for (b, x) in rhs.iter().zip(&solved) {
        let mut in_place = b.clone();
        factor.solve_in_place(&mut in_place).unwrap();
        assert_eq!(&in_place, x, "bit for bit");
        assert_close(&a.mul_vec(x).unwrap(), b, 1e-12);
    }
    assert!(factor.solve_many(&[]).unwrap().is_empty());

    // A wrong length is refused, and leaves an in-place argument untouched
    let mut short = vec![1.0; 80];
    assert_eq!(
        factor.solve_in_place(&mut short),
        Err(SparseError::ShapeMismatch {
            left: (81, 81),
            right: (80, 1)
        })
    );
    assert_eq!(short, vec![1.0; 80]);
    let mut mixed = rhs.clone();
    mixed.insert(2, vec![0.0; 3]);
    assert_eq!(
        factor.solve_many(&mixed).err(),
        Some(SparseError::ShapeMismatch {
            left: (81, 81),
            right: (3, 1)
        })
    );
}

#[test]
fn lu_solves_in_place_and_many_like_solve() {
    let a = random_square(50, 8, 0.05, 61).to_csc();
    let lu = a.lu().unwrap();
    let rhs = right_hand_sides(50, 4);
    let solved: Vec<Vec<f64>> = rhs.iter().map(|b| lu.solve(b).unwrap()).collect();
    assert_eq!(lu.solve_many(&rhs).unwrap(), solved, "bit for bit");
    for (b, x) in rhs.iter().zip(&solved) {
        let mut in_place = b.clone();
        lu.solve_in_place(&mut in_place).unwrap();
        assert_eq!(&in_place, x, "bit for bit");
        assert_close(&a.mul_vec(x).unwrap(), b, 1e-10);
    }

    // Complex right-hand sides
    let c = a.map(|v| Complex::new(v, 0.5 * v));
    let complex_lu = c.lu().unwrap();
    let complex_rhs: Vec<Vec<Complex<f64>>> = rhs
        .iter()
        .map(|b| b.iter().map(|&v| Complex::new(v, -v)).collect())
        .collect();
    for (b, x) in complex_rhs
        .iter()
        .zip(complex_lu.solve_many(&complex_rhs).unwrap())
    {
        let residual = c.mul_vec(&x).unwrap();
        assert!(residual.iter().zip(b).all(|(r, b)| (r - b).norm() < 1e-10));
    }

    let mut short = vec![2.0; 7];
    assert_eq!(
        lu.solve_in_place(&mut short),
        Err(SparseError::ShapeMismatch {
            left: (50, 50),
            right: (7, 1)
        })
    );
    assert_eq!(short, vec![2.0; 7]);
    assert!(lu.solve_many(&[vec![0.0; 49]]).is_err());
}

// Hermitian Cholesky

/// A Hermitian positive definite matrix on a `k x k` grid: each edge gets `-(1 + 0.5i)` one way
/// and its conjugate the other, with a real diagonal of 5, which dominates the four neighbours'
/// magnitudes (about 1.12 each). The numbering is then shuffled.
fn hermitian_grid(k: usize, seed: u64) -> CsrMatrix<Complex<f64>> {
    let edge = Complex::new(-1.0, -0.5);
    let mut triplets = vec![];
    for r in 0..k {
        for c in 0..k {
            let i = r * k + c;
            triplets.push((i, i, Complex::new(5.0, 0.0)));
            for j in [(c + 1 < k).then_some(i + 1), (r + 1 < k).then_some(i + k)]
                .into_iter()
                .flatten()
            {
                triplets.extend([(i, j, edge), (j, i, edge.conj())]);
            }
        }
    }
    let mut rng = Lcg(seed);
    CsrMatrix::from_triplets(k * k, k * k, triplets)
        .unwrap()
        .permute_symmetric(&rng.permutation(k * k))
        .unwrap()
}

#[test]
fn cholesky_factors_hermitian_matrices() {
    let a = hermitian_grid(7, 71);
    assert!(a.is_hermitian() && !a.is_symmetric());
    for ordering in [Ordering::Natural, Ordering::Rcm, Ordering::Amd] {
        let factor = a.analyze_cholesky(ordering).unwrap().factor(&a).unwrap();
        let l = factor.l();
        // A real, positive diagonal
        assert!(l
            .diagonal()
            .values()
            .iter()
            .all(|d| d.im == 0.0 && d.re > 0.0));
        // P A Pᵀ = L Lᴴ
        let reordered = a.permute_symmetric(factor.permutation()).unwrap();
        let product = l * &l.adjoint().to_csr();
        let difference = &product - &reordered;
        assert!(difference.norm_frobenius() < 1e-12, "{ordering:?}");
    }
    let x: Vec<Complex<f64>> = (0..49)
        .map(|i| Complex::new((i as f64).sin(), (i as f64 * 0.3).cos()))
        .collect();
    let b = a.mul_vec(&x).unwrap();
    for solved in [
        a.solve_spd(&b).unwrap(),
        a.cholesky()
            .unwrap()
            .solve_many(std::slice::from_ref(&b))
            .unwrap()
            .remove(0),
    ] {
        assert!(solved.iter().zip(&x).all(|(s, x)| (s - x).norm() < 1e-12));
    }
    // It agrees with LU on the same system
    let by_lu = a.to_csc().solve(&b).unwrap();
    assert!(a
        .solve_spd(&b)
        .unwrap()
        .iter()
        .zip(&by_lu)
        .all(|(s, l)| (s - l).norm() < 1e-12));
}

#[test]
fn cholesky_rejects_complex_matrices_that_are_not_hermitian_positive_definite() {
    // Complex symmetric, not Hermitian
    let symmetric = CsrMatrix::from_dense(&[
        [Complex::new(4.0, 0.0), Complex::new(1.0, 1.0)],
        [Complex::new(1.0, 1.0), Complex::new(4.0, 0.0)],
    ])
    .unwrap();
    assert_eq!(symmetric.cholesky().err(), Some(SparseError::NotSymmetric));
    // Hermitian but indefinite: eigenvalues 1 ± |2 + 2i| = 1 ± 2.83
    let indefinite = CsrMatrix::from_dense(&[
        [Complex::new(1.0, 0.0), Complex::new(2.0, 2.0)],
        [Complex::new(2.0, -2.0), Complex::new(1.0, 0.0)],
    ])
    .unwrap();
    assert!(indefinite.is_hermitian());
    assert!(matches!(
        indefinite.cholesky().err(),
        Some(SparseError::NotPositiveDefinite { .. })
    ));
}

#[test]
fn conjugate_gradient_solves_hermitian_systems() {
    let a = hermitian_grid(8, 81);
    let x: Vec<Complex<f64>> = (0..64)
        .map(|i| Complex::new(1.0 + (i % 3) as f64, (i as f64 * 0.2).sin()))
        .collect();
    let b = a.mul_vec(&x).unwrap();
    let by_cg = a.conjugate_gradient(&b, 1e-13, 500).unwrap();
    let by_cholesky = a.solve_spd(&b).unwrap();
    for ((cg, chol), x) in by_cg.iter().zip(&by_cholesky).zip(&x) {
        assert!((cg - chol).norm() < 1e-10, "{cg} vs {chol}");
        assert!((cg - x).norm() < 1e-10);
    }
    // Both orientations agree bit for bit, as for real matrices
    assert_eq!(
        a.to_csc().conjugate_gradient(&b, 1e-13, 500).unwrap(),
        by_cg
    );
    // CG assumes a Hermitian matrix without checking. On [[0, i], [i, 0]], which is symmetric
    //   but not Hermitian, the first curvature ⟨b, A b⟩ is zero, so it stops at once
    let not_hermitian = CsrMatrix::from_dense(&[
        [Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)],
        [Complex::new(0.0, 1.0), Complex::new(0.0, 0.0)],
    ])
    .unwrap();
    assert!(matches!(
        not_hermitian.conjugate_gradient(
            &[Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            1e-12,
            10
        ),
        Err(SparseError::NoConvergence { iterations: 0 })
    ));
}
