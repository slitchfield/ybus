// The front tree's shape (SymbolicLdlt::front_tree), on matrices whose fronts are worked out by
//   hand: the natural ordering without amalgamation, so the fronts are the fundamental
//   supernodes.

use ybus_sparse::compressed_matrix::{Amalgamation, CsrMatrix, FrontTree, Ordering, SymbolicLdlt};

fn tree(n: usize, entries: &[(usize, usize, f64)]) -> FrontTree {
    let a = CsrMatrix::from_triplets(n, n, entries.to_vec()).unwrap();
    let cholesky = a.analyze_cholesky(Ordering::Natural).unwrap();
    SymbolicLdlt::with_amalgamation(cholesky, Amalgamation::NONE).front_tree()
}

#[test]
fn two_independent_blocks_are_two_roots() {
    // Each 2×2 block is one front, c = 2 and r = 0 (m = 2): it factors in (1·2)/2 + 0 = 1
    //   multiply-add, and sweeps in 2·1/2 + 0 = 1
    let block = |o: usize| {
        [
            (o, o, 4.0),
            (o, o + 1, 1.0),
            (o + 1, o, 1.0),
            (o + 1, o + 1, 4.0),
        ]
    };
    let entries: Vec<_> = block(0).into_iter().chain(block(2)).collect();
    assert_eq!(
        tree(4, &entries),
        FrontTree {
            fronts: 2,
            roots: 2,
            height: 1,
            work: 2,
            critical_work: 1,
            solve_work: 2,
            solve_critical_work: 1,
        }
    );
}

#[test]
fn a_star_is_two_leaves_under_its_hub() {
    // Columns 0 and 1 each touch only the hub, 2: two leaf fronts (c = 1, r = 1, m = 2, each
    //   factoring in 1·2/2 = 1 and sweeping in 0 + 1 = 1) under the hub's (c = 1, r = 0: 0
    //   and 0), so a leaf-to-root chain holds 1 of the 2
    let entries = [
        (0, 0, 4.0),
        (1, 1, 4.0),
        (2, 2, 4.0),
        (0, 2, 1.0),
        (2, 0, 1.0),
        (1, 2, 1.0),
        (2, 1, 1.0),
    ];
    assert_eq!(
        tree(3, &entries),
        FrontTree {
            fronts: 3,
            roots: 1,
            height: 2,
            work: 2,
            critical_work: 1,
            solve_work: 2,
            solve_critical_work: 1,
        }
    );
}
