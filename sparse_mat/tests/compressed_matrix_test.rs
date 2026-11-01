// Tests for CompressedMatrix (CsrMatrix and CscMatrix): construction, validation, the fixed
//   pattern and its stored zeros, transposing and reorienting, equality and formatting.

use sparse_mat::compressed_matrix::{CompressedMatrix, CscMatrix, CsrMatrix};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::sparse_matrix::SparseMatrix;

/// [[1, 0, 2], [0, 0, 3]]
fn sample() -> SparseMatrix<f64> {
    SparseMatrix::from_triplets(2, 3, vec![(1, 2, 3.0), (0, 0, 1.0), (0, 2, 2.0)]).unwrap()
}

#[test]
fn sparse_to_csr_and_csc_store_the_same_entries() {
    let csr = sample().to_csr();
    assert_eq!(csr.shape(), (2, 3));
    assert_eq!(csr.pointers(), [0, 2, 3]);
    assert_eq!(csr.indices(), [0, 2, 2]);
    assert_eq!(csr.values(), [1.0, 2.0, 3.0]);

    let csc = sample().to_csc();
    assert_eq!(csc.pointers(), [0, 1, 1, 3]);
    assert_eq!(csc.indices(), [0, 0, 1]);
    assert_eq!(csc.values(), [1.0, 2.0, 3.0]);

    let dense = vec![vec![1.0, 0.0, 2.0], vec![0.0, 0.0, 3.0]];
    assert_eq!(csr.to_dense(), dense);
    assert_eq!(csc.to_dense(), dense);
    assert_eq!(CsrMatrix::from_sparse(&sample()), csr);
    assert_eq!(CscMatrix::from(&sample()), csc);
}

#[test]
fn raw_parts_are_validated_in_either_orientation() {
    let csr = |rows, cols, p: Vec<usize>, i: Vec<usize>, v: Vec<f64>| {
        match CsrMatrix::try_from_raw_parts(rows, cols, p, i, v).err() {
            Some(SparseError::InvalidCompressed { reason }) => reason,
            other => panic!("expected InvalidCompressed, got {other:?}"),
        }
    };
    let csc = |rows, cols, p: Vec<usize>, i: Vec<usize>, v: Vec<f64>| {
        match CscMatrix::try_from_raw_parts(rows, cols, p, i, v).err() {
            Some(SparseError::InvalidCompressed { reason }) => reason,
            other => panic!("expected InvalidCompressed, got {other:?}"),
        }
    };
    assert_eq!(
        csr(2, 2, vec![0, 1], vec![0], vec![1.0]),
        "pointers must have rows + 1 entries"
    );
    assert_eq!(
        csc(2, 3, vec![0, 1], vec![0], vec![1.0]),
        "pointers must have cols + 1 entries"
    );
    assert_eq!(
        csr(1, 2, vec![0, 1], vec![0], vec![]),
        "indices and values must have the same length"
    );
    assert_eq!(
        csc(2, 1, vec![0, 2], vec![0], vec![1.0]),
        "pointers must start at 0 and end at the entry count"
    );
    assert_eq!(
        csr(3, 2, vec![0, 2, 1, 2], vec![0, 1], vec![1.0, 2.0]),
        "pointers must be non-decreasing"
    );
    assert_eq!(
        csr(1, 2, vec![0, 1], vec![2], vec![1.0]),
        "column index out of range"
    );
    assert_eq!(
        csc(2, 1, vec![0, 1], vec![2], vec![1.0]),
        "row index out of range"
    );
    assert_eq!(
        csr(1, 3, vec![0, 2], vec![1, 1], vec![1.0, 2.0]),
        "column indices must be strictly increasing within a row"
    );
    assert_eq!(
        csc(3, 1, vec![0, 2], vec![2, 1], vec![1.0, 2.0]),
        "row indices must be strictly increasing within a column"
    );
    // rows + 1 doesn't fit in a usize, so no pointer array can have the right length
    assert_eq!(
        csr(usize::MAX, 1, vec![0], vec![], vec![]),
        "pointers must have rows + 1 entries"
    );
}

#[test]
fn raw_parts_keep_zeros_and_round_trip() {
    let matrix =
        CsrMatrix::try_from_raw_parts(2, 2, vec![0, 2, 3], vec![0, 1, 1], vec![1.0, 0.0, 4.0])
            .unwrap();
    assert_eq!(matrix.num_stored(), 3);
    assert_eq!(matrix.num_nonzero(), 2);
    let (shape, pointers, indices, values) = matrix.clone().into_raw_parts();
    let rebuilt = CsrMatrix::try_from_raw_parts(shape.0, shape.1, pointers, indices, values);
    assert_eq!(rebuilt.unwrap().values(), matrix.values());
}

#[test]
fn from_triplets_sums_in_input_order_and_keeps_every_position() {
    // 1e16 + 1 rounds back to 1e16, so the order of summation decides the result
    let in_order = CsrMatrix::from_triplets(1, 1, vec![(0, 0, 1e16), (0, 0, 1.0), (0, 0, -1e16)]);
    assert_eq!(in_order.unwrap().values(), [0.0]);
    let reordered = CsrMatrix::from_triplets(1, 1, vec![(0, 0, 1e16), (0, 0, -1e16), (0, 0, 1.0)]);
    assert_eq!(reordered.unwrap().values(), [1.0]);

    // A cancellation and a zero-valued triplet are both stored, unlike in SparseMatrix
    let triplets = vec![(0, 1, 2.0), (1, 0, 5.0), (0, 1, -2.0), (1, 1, 0.0)];
    let csr = CsrMatrix::from_triplets(2, 2, triplets.clone()).unwrap();
    assert_eq!(csr.pointers(), [0, 1, 3]);
    assert_eq!(csr.indices(), [1, 0, 1]);
    assert_eq!(csr.values(), [0.0, 5.0, 0.0]);
    assert_eq!((csr.num_stored(), csr.num_nonzero()), (3, 1));
    assert_eq!(
        SparseMatrix::from_triplets(2, 2, triplets.clone())
            .unwrap()
            .num_nonzero(),
        1
    );

    let csc = CscMatrix::from_triplets(2, 2, triplets).unwrap();
    assert_eq!(csc.pointers(), [0, 1, 3]);
    assert_eq!(csc.indices(), [1, 0, 1]);
    assert_eq!(csc.values(), [5.0, 0.0, 0.0]);
    assert_eq!(csc, csr);

    assert_eq!(
        CsrMatrix::from_triplets(2, 2, vec![(0, 0, 1.0), (2, 0, 1.0)]).unwrap_err(),
        SparseError::OutOfBounds {
            row: 2,
            col: 0,
            shape: (2, 2)
        }
    );
}

#[test]
fn values_mut_keeps_the_pattern() {
    let mut matrix = sample().to_csr();
    matrix.values_mut()[1] = 0.0;
    assert_eq!(matrix.num_stored(), 3);
    assert_eq!(matrix.num_nonzero(), 2);
    assert_eq!(matrix.indices(), [0, 2, 2]);
    // Converting back drops the stored zero, since SparseMatrix never stores one
    assert_eq!(matrix.to_sparse().num_nonzero(), 2);
}

#[test]
fn transpose_swaps_the_orientation_in_constant_time() {
    let csr = sample().to_csr();
    let transposed: CscMatrix<f64> = csr.clone().transpose();
    assert_eq!(transposed.shape(), (3, 2));
    // The same arrays, read along the other axis
    assert_eq!(transposed.pointers(), csr.pointers());
    assert_eq!(transposed.indices(), csr.indices());
    assert_eq!(
        transposed.to_dense(),
        vec![vec![1.0, 0.0], vec![0.0, 0.0], vec![2.0, 3.0]]
    );
    let back: CsrMatrix<f64> = transposed.transpose();
    assert_eq!(back.pointers(), csr.pointers());
}

#[test]
fn reorienting_round_trips_with_stored_zeros() {
    let csr = CsrMatrix::from_triplets(
        3,
        2,
        vec![(2, 1, 4.0), (0, 1, 0.0), (1, 0, 3.0), (2, 0, 1.0)],
    )
    .unwrap();
    let csc = csr.to_csc();
    assert_eq!(csc.pointers(), [0, 2, 4]);
    assert_eq!(csc.indices(), [1, 2, 0, 2]);
    assert_eq!(csc.values(), [3.0, 1.0, 0.0, 4.0]);
    let back = csc.to_csr();
    assert_eq!(back.pointers(), csr.pointers());
    assert_eq!(back.indices(), csr.indices());
    assert_eq!(back.values(), csr.values());
    // Reorienting to the same orientation copies
    assert_eq!(csr.to_csr().values(), csr.values());
}

#[test]
fn equality_is_mathematical() {
    let plain = sample().to_csr();
    let with_zero = CsrMatrix::from_triplets(
        2,
        3,
        vec![(0, 0, 1.0), (0, 1, 0.0), (0, 2, 2.0), (1, 2, 3.0)],
    )
    .unwrap();
    assert_eq!(with_zero.num_stored(), 4);
    assert_eq!(plain, with_zero);
    assert_eq!(plain, sample().to_csc());
    assert_eq!(plain, sample());
    assert_eq!(sample(), plain);

    let mut different = plain.clone();
    different.values_mut()[0] = 9.0;
    assert_ne!(plain, different);
    assert_ne!(plain, CsrMatrix::<f64>::identity(2));

    let nan = CsrMatrix::from_triplets(1, 1, vec![(0, 0, f64::NAN)]).unwrap();
    assert_ne!(nan, nan.clone());
}

#[test]
fn other_constructors() {
    let dense = CsrMatrix::from_dense(&[[1.0, 0.0], [0.0, 2.0]]).unwrap();
    assert_eq!(dense.num_stored(), 2);
    assert_eq!(dense, CsrMatrix::from_diagonal([1.0, 2.0]));
    assert_eq!(CscMatrix::<f64>::identity(3).values(), [1.0, 1.0, 1.0]);
    assert_eq!(CsrMatrix::from_diagonal([0.0, 5.0]).num_stored(), 1);
    assert_eq!(
        CsrMatrix::from_dense(&[vec![1.0], vec![1.0, 2.0]]).unwrap_err(),
        SparseError::RaggedRows {
            row: 1,
            expected: 1,
            found: 2
        }
    );

    let empty = CsrMatrix::<f64>::default();
    assert_eq!((empty.shape(), empty.pointers()), ((0, 0), &[0][..]));
    let wide = CscMatrix::<i64>::from_triplets(0, 4, Vec::new()).unwrap();
    assert_eq!(wide.pointers(), [0, 0, 0, 0, 0]);
    assert!(wide.iter().next().is_none());
}

#[test]
fn integer_elements() {
    let ints = CsrMatrix::from_triplets(2, 2, vec![(0, 1, 3i64), (0, 1, 4), (1, 0, -1)]).unwrap();
    assert_eq!(ints.values(), [7, -1]);
    assert_eq!(ints.to_csc().to_dense(), vec![vec![0, 7], vec![-1, 0]]);
}

#[test]
fn debug_and_display() {
    let csr = sample().to_csr();
    assert_eq!(
        format!("{csr:?}"),
        "CsrMatrix { shape: (2, 3), entries: [(0, 0, 1.0), (0, 2, 2.0), (1, 2, 3.0)] }"
    );
    assert_eq!(
        format!("{:?}", sample().to_csc()),
        "CscMatrix { shape: (2, 3), entries: [(0, 0, 1.0), (0, 2, 2.0), (1, 2, 3.0)] }"
    );
    assert_eq!(format!("{csr}"), format!("{}", sample()));
    assert_eq!(format!("{}", sample().to_csc()), format!("{}", sample()));
}

// Products

/// Deterministic pseudo-random sparse matrix. `values` maps the generator's draws to entries:
/// small integers make products exact, so they compare against a dense reference.
fn pseudo_random_with(
    rows: usize,
    cols: usize,
    seed: u64,
    fill_percent: u64,
    values: impl Fn(u64) -> f64,
) -> CsrMatrix<f64> {
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 33
    };
    let mut triplets = vec![];
    for row in 0..rows {
        for col in 0..cols {
            if next() % 100 < fill_percent {
                triplets.push((row, col, values(next())));
            }
        }
    }
    CsrMatrix::from_triplets(rows, cols, triplets).unwrap()
}

/// Small integer-valued entries, so products are exact.
fn pseudo_random(rows: usize, cols: usize, seed: u64, fill_percent: u64) -> CsrMatrix<f64> {
    pseudo_random_with(rows, cols, seed, fill_percent, |draw| {
        (draw % 9) as f64 - 4.0
    })
}

fn dense_matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = b.first().map_or(0, Vec::len);
    a.iter()
        .map(|row| {
            (0..n)
                .map(|j| row.iter().zip(b).map(|(x, b_row)| x * b_row[j]).sum())
                .collect()
        })
        .collect()
}

#[test]
fn mul_vec_in_both_orientations() {
    let a = CsrMatrix::from_dense(&[[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]]).unwrap();
    assert_eq!(a.mul_vec(&[1.0, 10.0, 100.0]).unwrap(), vec![201.0, 30.0]);
    assert_eq!(&a * &vec![1.0, 10.0, 100.0], vec![201.0, 30.0]);
    assert_eq!(&a * &[1.0, 1.0, 1.0][..], vec![3.0, 3.0]);
    assert_eq!(
        a.to_csc().mul_vec(&[1.0, 10.0, 100.0]).unwrap(),
        vec![201.0, 30.0]
    );
}

#[test]
fn mul_vec_matches_dense() {
    let a = pseudo_random(7, 5, 11, 40);
    let x = vec![1.0, -2.0, 3.0, 0.0, 5.0];
    let column: Vec<Vec<f64>> = x.iter().map(|&v| vec![v]).collect();
    let expected: Vec<f64> = dense_matmul(&a.to_dense(), &column)
        .into_iter()
        .map(|row| row[0])
        .collect();
    assert_eq!(a.mul_vec(&x).unwrap(), expected);
    assert_eq!(a.to_csc().mul_vec(&x).unwrap(), expected);
}

#[test]
fn mul_vec_shape_mismatch() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0], [3.0, 4.0]]).unwrap();
    let expected = Some(SparseError::ShapeMismatch {
        left: (2, 2),
        right: (3, 1),
    });
    assert_eq!(a.mul_vec(&[1.0, 2.0, 3.0]).err(), expected);
    assert_eq!(a.to_csc().mul_vec(&[1.0, 2.0, 3.0]).err(), expected);
}

#[test]
fn vec_mul_is_the_transpose_product() {
    let a = pseudo_random(6, 4, 3, 50);
    let x = vec![2.0, -1.0, 0.0, 4.0, 1.0, 3.0];
    // The transpose of a CSR matrix is a CSC matrix of the same arrays
    let via_transpose = a.clone().transpose().mul_vec(&x).unwrap();
    assert_eq!(a.vec_mul(&x).unwrap(), via_transpose);
    assert_eq!(a.to_csc().vec_mul(&x).unwrap(), via_transpose);
    assert_eq!(&x * &a, via_transpose);
    assert_eq!(&x[..] * &a, via_transpose);
}

#[test]
fn vec_mul_shape_mismatch() {
    let a = CscMatrix::from_dense(&[[1.0, 2.0], [3.0, 4.0]]).unwrap();
    assert_eq!(
        a.vec_mul(&[1.0]).err(),
        Some(SparseError::ShapeMismatch {
            left: (1, 1),
            right: (2, 2)
        })
    );
}

#[test]
fn mul_dense_matches_dense() {
    let a = pseudo_random(5, 4, 21, 45);
    let b = pseudo_random(4, 3, 22, 100).to_dense();
    let expected = dense_matmul(&a.to_dense(), &b);
    assert_eq!(a.mul_dense(&b).unwrap(), expected);
    assert_eq!(a.to_csc().mul_dense(&b).unwrap(), expected);
}

#[test]
fn mul_dense_errors() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0], [3.0, 4.0]]).unwrap();
    assert_eq!(
        a.mul_dense(&[vec![1.0, 2.0]]).err(),
        Some(SparseError::ShapeMismatch {
            left: (2, 2),
            right: (1, 2)
        })
    );
    assert_eq!(
        a.mul_dense(&[vec![1.0, 2.0], vec![3.0]]).err(),
        Some(SparseError::RaggedRows {
            row: 1,
            expected: 2,
            found: 1
        })
    );
}

#[test]
fn matmul_matches_dense_in_both_orientations() {
    for (seed, (m, k, n)) in [
        (1u64, (6, 5, 4)),
        (2, (3, 8, 3)),
        (3, (10, 10, 10)),
        (4, (1, 7, 1)),
    ] {
        let a = pseudo_random(m, k, seed, 35);
        let b = pseudo_random(k, n, seed + 100, 35);
        let expected = dense_matmul(&a.to_dense(), &b.to_dense());
        let product = a.matmul(&b).unwrap();
        assert_eq!(product.shape(), (m, n));
        assert_eq!(product.to_dense(), expected);
        assert_eq!(a.to_csc().matmul(&b.to_csc()).unwrap().to_dense(), expected);
        // The operator agrees, for every operand form
        assert_eq!(&a * &b, product);
        assert_eq!(a.clone() * b.clone(), product);
        assert_eq!(a.clone() * &b, product);
        assert_eq!(&a * b.clone(), product);
    }
}

#[test]
fn matmul_identity_and_zero() {
    let a = pseudo_random(6, 6, 9, 40);
    let identity = CsrMatrix::identity(6);
    assert_eq!(a.matmul(&identity).unwrap(), a);
    assert_eq!(identity.matmul(&a).unwrap(), a);

    let zero = CsrMatrix::<f64>::from_triplets(6, 6, Vec::new()).unwrap();
    assert_eq!(a.matmul(&zero).unwrap().num_stored(), 0);
}

#[test]
fn matmul_associative_and_transpose_rule() {
    let a = pseudo_random(4, 5, 31, 40);
    let b = pseudo_random(5, 3, 32, 40);
    let c = pseudo_random(3, 6, 33, 40);
    assert_eq!(&(&a * &b) * &c, &a * &(&b * &c));
    // (AB)ᵀ = BᵀAᵀ, where the transposes are CSC, so this exercises the CSC product
    let product_transposed: CscMatrix<f64> = (&a * &b).transpose();
    assert_eq!(
        product_transposed,
        b.clone().transpose() * a.clone().transpose()
    );
}

#[test]
fn matmul_keeps_cancellations_as_stored_zeros() {
    // Entry (0, 0) of the product is 1*1 + (-1)*1 = 0, but the patterns reach it
    let a = CsrMatrix::from_dense(&[[1.0, -1.0], [0.0, 2.0]]).unwrap();
    let b = CsrMatrix::from_dense(&[[1.0, 0.0], [1.0, 5.0]]).unwrap();
    for product in [
        a.matmul(&b).unwrap(),
        a.to_csc().matmul(&b.to_csc()).unwrap().to_csr(),
    ] {
        assert_eq!(product.to_dense(), vec![vec![0.0, -5.0], vec![2.0, 10.0]]);
        assert_eq!((product.num_stored(), product.num_nonzero()), (4, 3));
    }
}

#[test]
fn matmul_shape_mismatch() {
    let a = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        a.matmul(&a).err(),
        Some(SparseError::ShapeMismatch {
            left: (2, 3),
            right: (2, 3)
        })
    );
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn mul_operator_shape_mismatch_panics() {
    let a = CscMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    let _ = &a * &a;
}

#[test]
fn matmul_huge_sparse_shape() {
    let a = CsrMatrix::from_triplets(1 << 20, 1 << 20, vec![(3, 1 << 19, 2.0)]).unwrap();
    let b = CsrMatrix::from_triplets(1 << 20, 1 << 20, vec![(1 << 19, 7, 5.0)]).unwrap();
    let product = a.matmul(&b).unwrap();
    assert_eq!(product.iter().collect::<Vec<_>>(), vec![(3, 7, 10.0)]);
}

#[test]
fn csr_and_csc_products_agree_bit_for_bit() {
    // Values like 1/7 aren't exact in binary, so a different summation order would show
    let fraction = |draw: u64| (draw % 1000) as f64 / 7.0 - 70.0;
    let a = pseudo_random_with(9, 7, 5, 45, fraction);
    let b = pseudo_random_with(7, 6, 6, 45, fraction);
    let x: Vec<f64> = (0..7).map(|i| 1.0 / (i as f64 + 3.0)).collect();
    let y: Vec<f64> = (0..9).map(|i| 1.0 / (i as f64 + 2.0)).collect();
    let (a_csc, b_csc) = (a.to_csc(), b.to_csc());

    assert_eq!(a.mul_vec(&x).unwrap(), a_csc.mul_vec(&x).unwrap());
    assert_eq!(a.vec_mul(&y).unwrap(), a_csc.vec_mul(&y).unwrap());
    let b_dense = b.to_dense();
    assert_eq!(
        a.mul_dense(&b_dense).unwrap(),
        a_csc.mul_dense(&b_dense).unwrap()
    );

    let by_rows = a.matmul(&b).unwrap().to_csc();
    let by_columns = a_csc.matmul(&b_csc).unwrap();
    assert_eq!(by_rows.pointers(), by_columns.pointers());
    assert_eq!(by_rows.indices(), by_columns.indices());
    assert_eq!(by_rows.values(), by_columns.values());
}

// Solvers

fn assert_vec_close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert!(
            (a - e).abs() <= tolerance,
            "{a} differs from {e}: {actual:?} vs {expected:?}"
        );
    }
}

/// A strictly diagonally dominant (so non-singular) pseudo-random matrix.
fn dominant(n: usize, seed: u64) -> CsrMatrix<f64> {
    let shift = (0..n).map(|i| (i, i, 4.0 * n as f64 + 1.0));
    CsrMatrix::from_triplets(n, n, pseudo_random(n, n, seed, 30).iter().chain(shift)).unwrap()
}

/// A symmetric positive definite matrix `BᵀB + n I`.
fn spd(n: usize, seed: u64) -> CsrMatrix<f64> {
    let b = pseudo_random(n, n, seed, 25);
    let gram = &b.clone().transpose().to_csr() * &b;
    let shift = (0..n).map(|i| (i, i, n as f64));
    CsrMatrix::from_triplets(n, n, gram.iter().chain(shift)).unwrap()
}

#[test]
fn triangular_solves_in_both_orientations() {
    let lower =
        CsrMatrix::from_dense(&[[2.0, 0.0, 0.0], [1.0, 4.0, 0.0], [0.0, 3.0, 5.0]]).unwrap();
    // x = [1, 2, 3] -> b = L x
    let b = lower.mul_vec(&[1.0, 2.0, 3.0]).unwrap();
    let by_rows = lower.solve_lower_triangular(&b).unwrap();
    assert_vec_close(&by_rows, &[1.0, 2.0, 3.0], 1e-12);
    // Column-oriented forward substitution subtracts in the same order
    assert_eq!(lower.to_csc().solve_lower_triangular(&b).unwrap(), by_rows);

    let upper: CscMatrix<f64> = lower.transpose();
    let b = upper.mul_vec(&[1.0, 2.0, 3.0]).unwrap();
    assert_vec_close(
        &upper.solve_upper_triangular(&b).unwrap(),
        &[1.0, 2.0, 3.0],
        1e-12,
    );
    assert_vec_close(
        &upper.to_csr().solve_upper_triangular(&b).unwrap(),
        &[1.0, 2.0, 3.0],
        1e-12,
    );
}

#[test]
fn triangular_solves_ignore_the_other_triangle() {
    // Garbage across the diagonal must not affect a solve, in either orientation
    let a = CsrMatrix::from_dense(&[[2.0, 99.0], [1.0, 4.0]]).unwrap();
    for lower in [
        a.solve_lower_triangular(&[2.0, 9.0]).unwrap(),
        a.to_csc().solve_lower_triangular(&[2.0, 9.0]).unwrap(),
    ] {
        assert_vec_close(&lower, &[1.0, 2.0], 1e-12);
    }
    for upper in [
        a.solve_upper_triangular(&[200.0, 8.0]).unwrap(),
        a.to_csc().solve_upper_triangular(&[200.0, 8.0]).unwrap(),
    ] {
        assert_vec_close(&upper, &[1.0, 2.0], 1e-12);
    }
}

#[test]
fn triangular_solve_errors() {
    let missing_diagonal = CsrMatrix::from_dense(&[[1.0, 0.0], [1.0, 0.0]]).unwrap();
    let singular = Some(SparseError::SingularMatrix { pivot: 1 });
    assert_eq!(
        missing_diagonal.solve_lower_triangular(&[1.0, 1.0]).err(),
        singular
    );
    assert_eq!(
        missing_diagonal
            .to_csc()
            .solve_lower_triangular(&[1.0, 1.0])
            .err(),
        singular
    );
    // A stored zero on the diagonal is as singular as a missing one
    let stored_zero =
        CscMatrix::try_from_raw_parts(2, 2, vec![0, 1, 2], vec![0, 1], vec![1.0, 0.0]).unwrap();
    assert_eq!(
        stored_zero.solve_upper_triangular(&[1.0, 1.0]).err(),
        singular
    );

    let a = CsrMatrix::<f64>::identity(2);
    assert_eq!(
        a.solve_upper_triangular(&[1.0]).err(),
        Some(SparseError::ShapeMismatch {
            left: (2, 2),
            right: (1, 1)
        })
    );
    let rect = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        rect.solve_lower_triangular(&[1.0, 1.0]).err(),
        Some(SparseError::NotSquare { shape: (2, 3) })
    );
}

#[test]
fn lu_reconstructs_permuted_matrix() {
    for (n, seed) in [(1, 1), (4, 2), (9, 3), (15, 4)] {
        let a = dominant(n, seed);
        let lu = a.lu().unwrap();
        // P A = L U, with row i of P A being row perm[i] of A
        let pa: Vec<Vec<f64>> = lu
            .permutation()
            .iter()
            .map(|&r| a.to_dense()[r].clone())
            .collect();
        let product = (lu.l() * lu.u()).to_dense();
        for (got, want) in product.iter().zip(&pa) {
            assert_vec_close(got, want, 1e-9);
        }
    }
}

#[test]
fn lu_solves_systems() {
    for (n, seed) in [(3, 7), (10, 8), (25, 9)] {
        let a = dominant(n, seed);
        let x: Vec<f64> = (0..n).map(|i| i as f64 - 3.0).collect();
        let b = a.mul_vec(&x).unwrap();
        assert_vec_close(&a.solve(&b).unwrap(), &x, 1e-9);
        assert_vec_close(&a.lu().unwrap().solve(&b).unwrap(), &x, 1e-9);
    }
}

#[test]
fn lu_needs_pivoting() {
    // A zero on the diagonal forces a row swap
    let a = CsrMatrix::from_dense(&[[0.0, 2.0], [3.0, 1.0]]).unwrap();
    let lu = a.lu().unwrap();
    assert_eq!(lu.permutation(), [1, 0]);
    assert_vec_close(&lu.solve(&[4.0, 5.0]).unwrap(), &[1.0, 2.0], 1e-12);

    // A matrix whose rows are all shifted so that nothing is on the diagonal
    let n = 6;
    let shifted =
        CsrMatrix::from_triplets(n, n, (0..n).map(|i| (i, (i + 1) % n, (i + 2) as f64))).unwrap();
    let x: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
    let b = shifted.mul_vec(&x).unwrap();
    assert_vec_close(&shifted.solve(&b).unwrap(), &x, 1e-12);
}

#[test]
fn lu_handles_stored_zeros() {
    // A stored zero below the pivot needs no elimination, and the factors store no zeros
    let a = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 2.0), (1, 0, 0.0), (1, 1, 3.0)]).unwrap();
    let lu = a.lu().unwrap();
    assert_eq!(lu.l().num_stored(), 2);
    assert_eq!(lu.u().num_stored(), 2);
    assert_vec_close(&lu.solve(&[2.0, 3.0]).unwrap(), &[1.0, 1.0], 1e-12);
}

#[test]
fn lu_determinant() {
    let a = CsrMatrix::from_dense(&[[4.0f64, 3.0], [6.0, 3.0]]).unwrap();
    assert!((a.lu().unwrap().determinant() - (-6.0)).abs() < 1e-12);
    // Row swap flips the sign: det [[0,2],[3,1]] = -6
    let b = CsrMatrix::from_dense(&[[0.0f64, 2.0], [3.0, 1.0]]).unwrap();
    assert!((b.lu().unwrap().determinant() - (-6.0)).abs() < 1e-12);
    assert_eq!(
        CsrMatrix::<f64>::identity(5).lu().unwrap().determinant(),
        1.0
    );
    assert_eq!(CsrMatrix::<f64>::default().lu().unwrap().determinant(), 1.0);
}

#[test]
fn lu_singular_and_shape_errors() {
    let singular = CsrMatrix::from_dense(&[[1.0, 2.0], [2.0, 4.0]]).unwrap();
    assert_eq!(
        singular.lu().err(),
        Some(SparseError::SingularMatrix { pivot: 1 })
    );
    let empty_column = CsrMatrix::from_dense(&[[1.0, 0.0], [1.0, 0.0]]).unwrap();
    assert_eq!(
        empty_column.lu().err(),
        Some(SparseError::SingularMatrix { pivot: 1 })
    );

    let rect = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        rect.lu().err(),
        Some(SparseError::NotSquare { shape: (2, 3) })
    );

    let a = CsrMatrix::<f64>::identity(3);
    assert!(a.lu().unwrap().solve(&[1.0]).is_err());
    assert!(a.solve(&[1.0]).is_err());
}

#[test]
fn solve_spd_matches_lu() {
    let a = spd(14, 5);
    let x: Vec<f64> = (0..14).map(|i| (i as f64).sin()).collect();
    let b = a.mul_vec(&x).unwrap();
    assert_vec_close(&a.solve_spd(&b).unwrap(), &x, 1e-9);
    assert_vec_close(&a.solve(&b).unwrap(), &x, 1e-9);
}

#[test]
fn conjugate_gradient_solves_in_both_orientations() {
    for (n, seed) in [(5, 1), (20, 2), (40, 3)] {
        let a = spd(n, seed);
        let x: Vec<f64> = (0..n).map(|i| 1.0 + (i % 3) as f64).collect();
        let b = a.mul_vec(&x).unwrap();
        let solution = a.conjugate_gradient(&b, 1e-12, 10 * n).unwrap();
        assert_vec_close(&solution, &x, 1e-7);
        // CSC's products agree bit for bit, so its iterates do too
        assert_eq!(
            a.to_csc().conjugate_gradient(&b, 1e-12, 10 * n).unwrap(),
            solution
        );
    }
}

#[test]
fn conjugate_gradient_edge_cases() {
    // Identity converges in one step; a zero right-hand side needs none
    let identity = CsrMatrix::<f64>::identity(4);
    let b = [1.0, 2.0, 3.0, 4.0];
    assert_eq!(identity.conjugate_gradient(&b, 1e-12, 1).unwrap(), b);
    assert_eq!(
        identity.conjugate_gradient(&[0.0; 4], 1e-12, 0).unwrap(),
        vec![0.0; 4]
    );

    // Too few iterations
    let a = spd(30, 6);
    assert_eq!(
        a.conjugate_gradient(&[1.0; 30], 1e-14, 1).err(),
        Some(SparseError::NoConvergence { iterations: 1 })
    );
    // Negative definite is detected on the first step
    let negative = CsrMatrix::from_diagonal(vec![-1.0; 3]);
    assert_eq!(
        negative.conjugate_gradient(&[1.0; 3], 1e-12, 10).err(),
        Some(SparseError::NoConvergence { iterations: 0 })
    );
    assert!(identity.conjugate_gradient(&[1.0], 1e-12, 10).is_err());
}

// Access, selection, permutation and properties

fn sample_4x6() -> CsrMatrix<f64> {
    CsrMatrix::from_triplets(
        4,
        6,
        vec![
            (0, 0, 10.0),
            (0, 1, 20.0),
            (1, 1, 30.0),
            (2, 2, 50.0),
            (1, 3, 40.0),
            (2, 3, 60.0),
            (2, 4, 70.0),
            (3, 5, 80.0),
        ],
    )
    .unwrap()
}

/// `A + Aᵀ`, built from triplets.
fn symmetrized(a: &CsrMatrix<f64>) -> CsrMatrix<f64> {
    let (rows, _) = a.shape();
    let mirrored = a.iter().map(|(row, col, value)| (col, row, value));
    CsrMatrix::from_triplets(rows, rows, a.iter().chain(mirrored)).unwrap()
}

fn inverse_of(perm: &[usize]) -> Vec<usize> {
    let mut inverse = vec![0; perm.len()];
    for (position, &index) in perm.iter().enumerate() {
        inverse[index] = position;
    }
    inverse
}

#[test]
fn row_and_column_extraction() {
    let csr = sample_4x6();
    let row = csr.row(2).unwrap();
    assert_eq!(row.len(), 6);
    assert_eq!(row.nnz(), 3);
    assert_eq!(row.indices(), [2, 3, 4]);
    assert_eq!(row.values(), [50.0, 60.0, 70.0]);
    assert_eq!(row.get(3), Some(60.0));
    assert!(row.get(0).is_none());
    assert_eq!(row.to_dense(), vec![0.0, 0.0, 50.0, 60.0, 70.0, 0.0]);
    let empty = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new())
        .unwrap()
        .row(1)
        .unwrap();
    assert!(empty.nnz() == 0 && empty.len() == 3);

    let col = csr.to_csc().col(3).unwrap();
    assert_eq!(col.len(), 4);
    assert_eq!(col.indices(), [1, 2]);
    assert_eq!(col.values(), [40.0, 60.0]);

    // Stored zeros come along
    let with_zero = CsrMatrix::from_triplets(1, 3, vec![(0, 0, 1.0), (0, 2, 0.0)]).unwrap();
    assert_eq!(with_zero.row(0).unwrap().indices(), [0, 2]);

    assert_eq!(
        csr.row(4).err(),
        Some(SparseError::OutOfRange { index: 4, len: 4 })
    );
    assert_eq!(
        csr.to_csc().col(6).err(),
        Some(SparseError::OutOfRange { index: 6, len: 6 })
    );
}

#[test]
fn dense_row_and_column_iterators() {
    let csr = sample_4x6();
    let rows: Vec<Vec<f64>> = csr.row_iter().collect();
    assert_eq!(rows, csr.to_dense());
    let columns: Vec<Vec<f64>> = csr.to_csc().col_iter().collect();
    assert_eq!(columns.len(), 6);
    assert_eq!(columns[3], vec![0.0, 40.0, 60.0, 0.0]);
}

#[test]
fn diagonal_in_both_orientations() {
    for diag in [sample_4x6().diagonal(), sample_4x6().to_csc().diagonal()] {
        assert_eq!(diag.len(), 4);
        assert_eq!(diag.indices(), [0, 1, 2]);
        assert_eq!(diag.values(), [10.0, 30.0, 50.0]);
        assert_eq!(diag.to_dense(), vec![10.0, 30.0, 50.0, 0.0]);
    }
    assert_eq!(
        CsrMatrix::<f64>::identity(3).diagonal().to_dense(),
        vec![1.0; 3]
    );
    // Tall matrix: diagonal length is the smaller dimension
    let tall = CscMatrix::<f64>::from_triplets(5, 2, Vec::new()).unwrap();
    assert_eq!(tall.diagonal().len(), 2);
}

#[test]
fn outer_borrows_a_row_or_column() {
    let csr = sample_4x6();
    assert_eq!(
        csr.outer(2).unwrap(),
        (&[2, 3, 4][..], &[50.0, 60.0, 70.0][..])
    );
    let csc = csr.to_csc();
    assert_eq!(csc.outer(3).unwrap(), (&[1, 2][..], &[40.0, 60.0][..]));
    // The slices are the arrays' own: entry k is slot pointers()[i] + k
    for i in 0..4 {
        let (indices, values) = csr.outer(i).unwrap();
        let start = csr.pointers()[i];
        assert_eq!(indices, &csr.indices()[start..start + indices.len()]);
        assert_eq!(values.as_ptr(), csr.values()[start..].as_ptr());
    }
    // An empty row, and stored zeros
    let empty: (&[usize], &[f64]) = (&[], &[]);
    assert_eq!(
        CsrMatrix::<f64>::from_triplets(2, 3, Vec::new())
            .unwrap()
            .outer(1)
            .unwrap(),
        empty
    );
    let with_zero = CsrMatrix::from_triplets(1, 3, vec![(0, 0, 1.0), (0, 2, 0.0)]).unwrap();
    assert_eq!(with_zero.outer(0).unwrap(), (&[0, 2][..], &[1.0, 0.0][..]));

    assert_eq!(
        csr.outer(4).err(),
        Some(SparseError::OutOfRange { index: 4, len: 4 })
    );
    assert_eq!(
        csc.outer(6).err(),
        Some(SparseError::OutOfRange { index: 6, len: 6 })
    );
}

#[test]
fn outer_mut_refills_in_place() {
    let mut csr = sample_4x6();
    {
        let (indices, values) = csr.outer_mut(2).unwrap();
        assert_eq!(indices, [2, 3, 4]);
        for (value, &col) in values.iter_mut().zip(indices) {
            *value = col as f64;
        }
    }
    // Setting a value to zero keeps it stored
    csr.outer_mut(0).unwrap().1[1] = 0.0;
    assert_eq!(csr.indices(), sample_4x6().indices());
    assert_eq!(csr.row(2).unwrap().values(), [2.0, 3.0, 4.0]);
    assert_eq!(csr.num_stored(), 8);
    assert_eq!(csr.num_nonzero(), 7);

    let mut csc = sample_4x6().to_csc();
    csc.outer_mut(3).unwrap().1.copy_from_slice(&[-1.0, -2.0]);
    assert_eq!(csc.col(3).unwrap().values(), [-1.0, -2.0]);
    assert_eq!(
        csc.outer_mut(6).err(),
        Some(SparseError::OutOfRange { index: 6, len: 6 })
    );
}

#[test]
fn position_finds_each_stored_entry() {
    let csr = sample_4x6();
    let csc = csr.to_csc();
    let dense = csr.to_dense();
    for (row, dense_row) in dense.iter().enumerate() {
        for (col, &value) in dense_row.iter().enumerate() {
            let stored = value != 0.0;
            for (position, values) in [
                (csr.position(row, col).unwrap(), csr.values()),
                (csc.position(row, col).unwrap(), csc.values()),
            ] {
                assert_eq!(position.is_some(), stored, "({row}, {col})");
                if let Some(p) = position {
                    assert_eq!(values[p], value);
                }
            }
        }
    }
    // A stored zero has a position
    let with_zero = CscMatrix::from_triplets(2, 2, vec![(1, 0, 0.0)]).unwrap();
    assert_eq!(with_zero.position(1, 0), Ok(Some(0)));
    assert_eq!(with_zero.position(0, 0), Ok(None));

    for (row, col) in [(4, 0), (0, 6)] {
        let error = Some(SparseError::OutOfBounds {
            row,
            col,
            shape: (4, 6),
        });
        assert_eq!(csr.position(row, col).err(), error);
        assert_eq!(csc.position(row, col).err(), error);
    }
}

#[test]
fn diagonal_positions_are_cached_and_survive_refills() {
    let mut csr = sample_4x6();
    let expected = [Some(0), Some(2), Some(4), None];
    assert_eq!(csr.diagonal_positions(), expected);
    for (i, &position) in expected.iter().enumerate() {
        assert_eq!(csr.position(i, i), Ok(position));
    }
    // The cache stays right after values change, including to zero
    csr.values_mut()[2] = 0.0;
    assert_eq!(csr.diagonal_positions(), expected);
    assert_eq!(csr.diagonal().indices(), [0, 1, 2]);
    assert_eq!(csr.diagonal().values(), [10.0, 0.0, 50.0]);

    // A clone and a transpose keep it; reorienting finds the new slots
    let clone = csr.clone();
    assert_eq!(clone.diagonal_positions(), expected);
    let transposed = csr.clone().transpose();
    assert_eq!(transposed.diagonal_positions(), expected);
    let csc = csr.to_csc();
    let csc_expected: Vec<Option<usize>> = (0..4).map(|i| csc.position(i, i).unwrap()).collect();
    assert_eq!(csc.diagonal_positions(), csc_expected);
    assert_eq!(csc_expected, [Some(0), Some(2), Some(3), None]);
    assert_eq!(csc.to_csr().diagonal_positions(), expected);

    // Rectangular and empty matrices
    let tall = CscMatrix::from_triplets(5, 2, vec![(1, 1, 7.0), (4, 1, 1.0)]).unwrap();
    assert_eq!(tall.diagonal_positions(), [None, Some(0)]);
    assert!(CsrMatrix::<f64>::default().diagonal_positions().is_empty());
}

#[test]
fn submatrix_in_both_orientations() {
    let csr = sample_4x6();
    let expected = vec![vec![30.0, 0.0, 40.0], vec![0.0, 50.0, 60.0]];
    let block = csr.submatrix(1..3, 1..4).unwrap();
    assert_eq!(block.shape(), (2, 3));
    assert_eq!(block.to_dense(), expected);
    assert_eq!(
        csr.to_csc().submatrix(1..3, 1..4).unwrap().to_dense(),
        expected
    );
    assert_eq!(csr.submatrix(0..4, 0..6).unwrap(), csr);
    assert_eq!(csr.submatrix(2..2, 0..6).unwrap().shape(), (0, 6));
}

#[test]
fn submatrix_invalid_range() {
    let csr = sample_4x6();
    assert_eq!(
        csr.submatrix(0..5, 0..6).err(),
        Some(SparseError::InvalidRange {
            start: 0,
            end: 5,
            len: 4
        })
    );
    #[allow(clippy::reversed_empty_ranges)]
    let reversed = csr.submatrix(0..4, 4..2);
    assert_eq!(
        reversed.err(),
        Some(SparseError::InvalidRange {
            start: 4,
            end: 2,
            len: 6
        })
    );
}

#[test]
fn density_counts_the_pattern() {
    assert_eq!(CsrMatrix::<f64>::identity(4).density(), 0.25);
    assert_eq!(
        CscMatrix::from_dense(&[[1.0, 2.0], [3.0, 4.0]])
            .unwrap()
            .density(),
        1.0
    );
    let stored_zero = CsrMatrix::from_triplets(2, 2, vec![(0, 1, 0.0)]).unwrap();
    assert_eq!(stored_zero.density(), 0.25);
    // No division by zero for an empty shape
    assert_eq!(CsrMatrix::<f64>::default().density(), 0.0);
    assert_eq!(
        CsrMatrix::<f64>::from_triplets(0, 5, Vec::new())
            .unwrap()
            .density(),
        0.0
    );
}

#[test]
fn sum_and_trace() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0, 0.0], [0.0, 3.0, 4.0], [5.0, 0.0, 6.0]]).unwrap();
    for (sum, trace) in [(a.sum(), a.trace()), (a.to_csc().sum(), a.to_csc().trace())] {
        assert_eq!(sum, 21.0);
        assert_eq!(trace, Ok(10.0));
    }
    let ints = CsrMatrix::from_dense(&[[1, 2], [3, 4]]).unwrap();
    assert!(ints.sum() == 10 && ints.trace() == Ok(5));
    let rect = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(rect.trace(), Err(SparseError::NotSquare { shape: (2, 3) }));
}

#[test]
fn is_symmetric() {
    let check = |dense: [[f64; 2]; 2]| {
        let csr = CsrMatrix::from_dense(&dense).unwrap();
        assert_eq!(csr.is_symmetric(), csr.to_csc().is_symmetric());
        csr.is_symmetric()
    };
    assert!(check([[1.0, 2.0], [2.0, 3.0]]));
    assert!(!check([[1.0, 2.0], [5.0, 3.0]]));
    // One side stored, the other not
    assert!(!check([[0.0, 2.0], [0.0, 0.0]]));
    // A stored zero opposite a missing entry is still symmetric
    assert!(CsrMatrix::from_triplets(2, 2, vec![(0, 1, 0.0)])
        .unwrap()
        .is_symmetric());
    assert!(!CsrMatrix::<f64>::from_triplets(2, 3, Vec::new())
        .unwrap()
        .is_symmetric());
    assert!(CsrMatrix::<f64>::identity(5).is_symmetric());
    assert!(CsrMatrix::<f64>::default().is_symmetric());
    assert!(symmetrized(&pseudo_random(6, 6, 77, 40)).is_symmetric());
}

#[test]
fn norms_in_both_orientations() {
    let a = CsrMatrix::from_dense(&[[1.0, -2.0, 0.0], [0.0, 3.0, 4.0]]).unwrap();
    for (frobenius, one, inf) in [
        (a.norm_frobenius(), a.norm_1(), a.norm_inf()),
        (
            a.to_csc().norm_frobenius(),
            a.to_csc().norm_1(),
            a.to_csc().norm_inf(),
        ),
    ] {
        assert!((frobenius - 30.0f64.sqrt()).abs() < 1e-12);
        // Column sums: 1, 5, 4 and row sums: 3, 7
        assert_eq!((one, inf), (5.0, 7.0));
    }
    let empty = CsrMatrix::<f64>::from_triplets(3, 3, Vec::new()).unwrap();
    assert!(empty.norm_frobenius() == 0.0 && empty.norm_1() == 0.0 && empty.norm_inf() == 0.0);

    // f32 works too, and transposing swaps the 1- and infinity-norms
    let b = CsrMatrix::from_dense(&[[1.0f32, 2.0], [3.0, 4.0]]).unwrap();
    assert_eq!(b.norm_1(), b.clone().transpose().norm_inf());
    assert_eq!(b.norm_inf(), b.clone().transpose().norm_1());
}

#[test]
fn select_matches_submatrix_for_ranges() {
    let a = pseudo_random(6, 7, 71, 45);
    let rows: Vec<usize> = (1..5).collect();
    let cols: Vec<usize> = (2..7).collect();
    assert_eq!(
        a.select(&rows, &cols).unwrap(),
        a.submatrix(1..5, 2..7).unwrap()
    );
}

#[test]
fn select_reorders_skips_and_repeats() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]).unwrap();
    // Rows 2 then 0, columns 2 then 1: non-contiguous and reversed
    let expected = vec![vec![9.0, 8.0], vec![3.0, 2.0]];
    assert_eq!(a.select(&[2, 0], &[2, 1]).unwrap().to_dense(), expected);
    assert_eq!(
        a.to_csc().select(&[2, 0], &[2, 1]).unwrap().to_dense(),
        expected
    );

    let b = CsrMatrix::from_dense(&[[1.0, 0.0], [0.0, 2.0]]).unwrap();
    let repeated = vec![
        vec![0.0, 2.0, 2.0],
        vec![0.0, 2.0, 2.0],
        vec![1.0, 0.0, 0.0],
    ];
    assert_eq!(
        b.select(&[1, 1, 0], &[0, 1, 1]).unwrap().to_dense(),
        repeated
    );
    assert_eq!(
        b.to_csc()
            .select(&[1, 1, 0], &[0, 1, 1])
            .unwrap()
            .to_dense(),
        repeated
    );
}

#[test]
fn select_matches_dense_selection_and_keeps_stored_zeros() {
    let a = pseudo_random(8, 8, 72, 40);
    let dense = a.to_dense();
    // A PV/PQ-style split: an arbitrary subset in a scrambled order
    let rows = [7, 1, 4, 2];
    let cols = [3, 0, 6, 5, 1];
    let expected: Vec<Vec<f64>> = rows
        .iter()
        .map(|&r| cols.iter().map(|&c| dense[r][c]).collect())
        .collect();
    assert_eq!(a.select(&rows, &cols).unwrap().to_dense(), expected);
    assert_eq!(
        a.to_csc().select(&rows, &cols).unwrap().to_dense(),
        expected
    );

    let stored_zero = CsrMatrix::from_triplets(2, 2, vec![(1, 0, 0.0), (0, 0, 1.0)]).unwrap();
    assert_eq!(stored_zero.select(&[1], &[0]).unwrap().num_stored(), 1);
}

#[test]
fn select_empty_and_out_of_range() {
    let a = sample_4x6();
    let none = a.select(&[], &[0, 1]).unwrap();
    assert!(none.shape() == (0, 2) && none.num_stored() == 0);
    assert_eq!(
        a.select(&[0, 4], &[0]).err(),
        Some(SparseError::OutOfRange { index: 4, len: 4 })
    );
    assert_eq!(
        a.to_csc().select(&[0], &[9]).err(),
        Some(SparseError::OutOfRange { index: 9, len: 6 })
    );
}

#[test]
fn permute_matches_select_and_inverts() {
    let a = pseudo_random(5, 4, 81, 50);
    let (rows, cols) = ([3, 0, 4, 1, 2], [2, 3, 0, 1]);
    assert_eq!(
        a.permute(&rows, &cols).unwrap(),
        a.select(&rows, &cols).unwrap()
    );

    let b = pseudo_random(6, 5, 82, 45);
    let (rows, cols) = ([5, 2, 0, 4, 1, 3], [4, 0, 3, 1, 2]);
    for permuted in [
        b.permute(&rows, &cols).unwrap(),
        b.to_csc().permute(&rows, &cols).unwrap().to_csr(),
    ] {
        assert_ne!(permuted, b);
        let restored = permuted
            .permute(&inverse_of(&rows), &inverse_of(&cols))
            .unwrap();
        assert_eq!(restored, b);
    }
    // The identity permutation changes nothing
    assert_eq!(b.permute(&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 4]).unwrap(), b);
}

#[test]
fn permute_symmetric_preserves_structure() {
    let symmetric = symmetrized(&pseudo_random(7, 7, 83, 40));
    let perm = [6, 3, 0, 5, 1, 4, 2];
    let permuted = symmetric.permute_symmetric(&perm).unwrap();
    assert!(permuted.is_symmetric());
    assert_eq!(permuted.trace(), symmetric.trace());
    let mut before = symmetric.diagonal().to_dense();
    let mut after = permuted.diagonal().to_dense();
    before.sort_by(f64::total_cmp);
    after.sort_by(f64::total_cmp);
    assert_eq!(before, after);
    // P A Pᵀ written with explicit permutation matrices
    let p =
        CsrMatrix::from_triplets(7, 7, perm.iter().enumerate().map(|(i, &j)| (i, j, 1.0))).unwrap();
    assert_eq!(
        permuted,
        &(&p * &symmetric) * &p.clone().transpose().to_csr()
    );
}

#[test]
fn permute_rows_agrees_with_lu_permutation() {
    let a = dominant(6, 84);
    let lu = a.lu().unwrap();
    let columns: Vec<usize> = (0..6).collect();
    let pa = a.permute(lu.permutation(), &columns).unwrap().to_dense();
    let product = (lu.l() * lu.u()).to_dense();
    for (got, want) in product.iter().zip(&pa) {
        assert_vec_close(got, want, 1e-9);
    }
}

#[test]
fn permute_errors() {
    let a = CsrMatrix::<f64>::from_triplets(3, 2, Vec::new()).unwrap();
    let reason = |result: Result<CsrMatrix<f64>, SparseError>| match result.err() {
        Some(SparseError::InvalidPermutation { reason }) => reason,
        other => panic!("expected InvalidPermutation, got {other:?}"),
    };
    assert_eq!(
        reason(a.permute(&[0, 1], &[0, 1])),
        "length must equal the dimension being permuted"
    );
    assert_eq!(reason(a.permute(&[0, 1, 3], &[0, 1])), "index out of range");
    assert_eq!(reason(a.permute(&[0, 1, 1], &[0, 1])), "index repeated");
    assert_eq!(reason(a.permute(&[0, 1, 2], &[1, 1])), "index repeated");
    assert_eq!(
        a.permute_symmetric(&[0, 1, 2]).err(),
        Some(SparseError::NotSquare { shape: (3, 2) })
    );
    // Empty matrices have exactly one permutation
    assert_eq!(
        CsrMatrix::<f64>::default()
            .permute_symmetric(&[])
            .unwrap()
            .shape(),
        (0, 0)
    );
}

// Arithmetic, block construction and complex parts

use sparse_mat::Complex;

fn sample_2x2() -> CsrMatrix<f64> {
    CsrMatrix::from_dense(&[[1.0, 2.0], [0.0, 3.0]]).unwrap()
}

fn pseudo_random_complex(rows: usize, cols: usize, seed: u64) -> CsrMatrix<Complex<f64>> {
    let re = pseudo_random(rows, cols, seed, 40);
    let im = pseudo_random(rows, cols, seed + 1000, 40);
    CsrMatrix::from_parts(&re, &im).unwrap()
}

fn dense_kron(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut out = vec![];
    for a_row in a {
        for b_row in b {
            out.push(
                a_row
                    .iter()
                    .flat_map(|x| b_row.iter().map(move |y| x * y))
                    .collect(),
            );
        }
    }
    out
}

#[test]
fn addition_stores_the_union_and_keeps_cancellations() {
    let a = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 1.5), (1, 1, 2.0)]).unwrap();
    let b = CsrMatrix::from_triplets(2, 2, vec![(0, 0, -1.5), (0, 1, 1.0)]).unwrap();
    for sum in [&a + &b, (&a.to_csc() + &b.to_csc()).to_csr()] {
        assert_eq!(sum.to_dense(), vec![vec![0.0, 1.0], vec![0.0, 2.0]]);
        assert_eq!((sum.num_stored(), sum.num_nonzero()), (3, 2));
    }
    assert_eq!(a.try_add(&b).unwrap(), &a + &b);
    // Operands are untouched
    assert_eq!(a.values(), [1.5, 2.0]);

    let symmetric = &sample_2x2() + &sample_2x2().transpose().to_csr();
    assert_eq!(symmetric.to_dense(), vec![vec![2.0, 2.0], vec![2.0, 6.0]]);
}

#[test]
fn subtraction_and_owned_operands() {
    let expected_sum = vec![vec![2.0, 4.0], vec![0.0, 6.0]];
    assert_eq!((&sample_2x2() + &sample_2x2()).to_dense(), expected_sum);
    assert_eq!((sample_2x2() + &sample_2x2()).to_dense(), expected_sum);
    assert_eq!((&sample_2x2() + sample_2x2()).to_dense(), expected_sum);
    assert_eq!((sample_2x2() + sample_2x2()).to_dense(), expected_sum);

    for difference in [
        sample_2x2() - sample_2x2(),
        &sample_2x2() - sample_2x2(),
        sample_2x2() - &sample_2x2(),
    ] {
        assert_eq!((difference.num_stored(), difference.num_nonzero()), (3, 0));
    }
    // An entry only on the right is negated
    let right_only = CsrMatrix::from_triplets(2, 2, vec![(1, 0, 4.0)]).unwrap();
    assert_eq!(
        (&sample_2x2() - &right_only).to_dense(),
        vec![vec![1.0, 2.0], vec![-4.0, 3.0]]
    );
}

#[test]
fn add_assign_and_sub_assign() {
    let mut a = sample_2x2();
    a += &CsrMatrix::from_dense(&[[1.0, 0.0], [0.0, 0.0]]).unwrap();
    assert_eq!(a.to_dense()[0][0], 2.0);
    a += CsrMatrix::from_dense(&[[0.0, -2.0], [0.0, 0.0]]).unwrap();
    assert_eq!((a.num_stored(), a.num_nonzero()), (3, 2));
    a -= &CsrMatrix::from_dense(&[[2.0, 0.0], [0.0, 3.0]]).unwrap();
    assert_eq!(a.num_nonzero(), 0);
    a -= CsrMatrix::from_dense(&[[1.0, 0.0], [0.0, 0.0]]).unwrap();
    assert_eq!(a.to_dense()[0][0], -1.0);
}

#[test]
fn shape_mismatches() {
    let a = CsrMatrix::<f64>::from_triplets(3, 3, Vec::new()).unwrap();
    let b = CsrMatrix::<f64>::from_triplets(2, 4, Vec::new()).unwrap();
    let mismatch = Some(SparseError::ShapeMismatch {
        left: (3, 3),
        right: (2, 4),
    });
    assert_eq!(a.try_add(&b).err(), mismatch);
    assert_eq!(a.try_sub(&b).err(), mismatch);
    assert_eq!(a.hadamard(&b).err(), mismatch);
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn add_operator_shape_mismatch_panics() {
    let a = CsrMatrix::<f64>::from_triplets(3, 3, Vec::new()).unwrap();
    let b = CsrMatrix::<f64>::from_triplets(2, 2, Vec::new()).unwrap();
    let _ = &a + &b;
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn add_assign_shape_mismatch_panics() {
    let mut a = sample_2x2();
    a += &CsrMatrix::<f64>::from_triplets(1, 1, Vec::new()).unwrap();
}

#[test]
fn element_types() {
    let ints = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 3), (0, 1, 1), (0, 0, 4)]).unwrap();
    let sum = &ints + &ints;
    assert_eq!(sum.to_dense(), vec![vec![14, 2], vec![0, 0]]);
    assert!(format!("{sum}").contains("14"));

    let doubled = &CscMatrix::<f32>::identity(3) + &CscMatrix::<f32>::identity(3);
    assert_eq!(doubled.values(), [2.0f32; 3]);

    let i = Complex::new(0.0, 1.0);
    let complex =
        CsrMatrix::from_triplets(2, 2, vec![(0, 1, i), (1, 0, Complex::new(1.0, 0.0))]).unwrap();
    assert_eq!((&complex + &complex).values()[0], Complex::new(0.0, 2.0));
}

#[test]
fn negation() {
    let a = sample_2x2();
    let negated = -&a;
    assert_eq!(negated.to_dense(), vec![vec![-1.0, -2.0], vec![0.0, -3.0]]);
    assert_eq!(a, sample_2x2());
    assert_eq!(-sample_2x2(), negated);
    assert_eq!((&a + &negated).num_nonzero(), 0);
    let ints = CsrMatrix::from_dense(&[[1, 0], [0, -4]]).unwrap();
    assert_eq!((-ints).to_dense(), vec![vec![-1, 0], vec![0, 4]]);
}

#[test]
fn scalar_multiplication_keeps_the_pattern() {
    let a = sample_2x2();
    assert_eq!((&a * 2.0).to_dense(), vec![vec![2.0, 4.0], vec![0.0, 6.0]]);
    assert_eq!(
        (a.clone() * 0.5).to_dense(),
        vec![vec![0.5, 1.0], vec![0.0, 1.5]]
    );
    assert_eq!(2.0 * &a, &a * 2.0);
    assert_eq!(2.0 * a.clone(), &a * 2.0);
    assert_eq!(a.scale(3.0).values(), [3.0, 6.0, 9.0]);
    // Multiplying by zero keeps every entry, stored as zero
    let zeroed = &a * 0.0;
    assert_eq!((zeroed.num_stored(), zeroed.num_nonzero()), (3, 0));

    let mut in_place = sample_2x2();
    in_place *= 10.0;
    assert_eq!(in_place.values(), [10.0, 20.0, 30.0]);
    // Alongside the matrix product
    assert_eq!((&a * &a).to_dense(), vec![vec![1.0, 8.0], vec![0.0, 9.0]]);
}

#[test]
fn hadamard_stores_the_intersection() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0, 0.0], [0.0, 3.0, 4.0]]).unwrap();
    let b = CsrMatrix::from_dense(&[[5.0, 0.0, 6.0], [0.0, 7.0, 8.0]]).unwrap();
    let product = a.hadamard(&b).unwrap();
    assert_eq!(
        product.to_dense(),
        vec![vec![5.0, 0.0, 0.0], vec![0.0, 21.0, 32.0]]
    );
    assert_eq!(product.num_stored(), 3);
    assert_eq!(b.hadamard(&a).unwrap(), product);
    assert_eq!(a.to_csc().hadamard(&b.to_csc()).unwrap(), product);

    let disjoint = CsrMatrix::from_dense(&[[1.0, 0.0], [0.0, 0.0]])
        .unwrap()
        .hadamard(&CsrMatrix::from_dense(&[[0.0, 2.0], [0.0, 0.0]]).unwrap())
        .unwrap();
    assert_eq!(disjoint.num_stored(), 0);
}

#[test]
fn map_and_apply_keep_the_pattern() {
    let a = CsrMatrix::from_dense(&[[1.0, 0.0], [-2.0, 3.0]]).unwrap();
    let signs = a.map(|x| if x > 0.0 { 1i32 } else { 0 });
    assert_eq!(signs.to_dense(), vec![vec![1, 0], vec![0, 1]]);
    assert_eq!((signs.num_stored(), signs.num_nonzero()), (3, 2));
    // f is only applied to stored entries, so implicit zeros stay zero
    assert_eq!(
        a.map(|x| x + 1.0).to_dense(),
        vec![vec![2.0, 0.0], vec![-1.0, 4.0]]
    );

    let mut b = CscMatrix::from_dense(&[[1.0, 0.0], [-1.0, 3.0]]).unwrap();
    b.apply(|x| x * x - 1.0);
    assert_eq!(b.to_dense(), vec![vec![0.0, 0.0], vec![0.0, 8.0]]);
    assert_eq!((b.num_stored(), b.num_nonzero()), (3, 1));
}

#[test]
fn kron_in_both_orientations() {
    let a = pseudo_random(3, 2, 5, 60);
    let b = pseudo_random(2, 4, 6, 60);
    let expected = dense_kron(&a.to_dense(), &b.to_dense());
    let product = a.kron(&b).unwrap();
    assert_eq!(product.shape(), (6, 8));
    assert_eq!(product.to_dense(), expected);
    assert_eq!(product.num_stored(), a.num_stored() * b.num_stored());
    assert_eq!(a.to_csc().kron(&b.to_csc()).unwrap().to_dense(), expected);

    let block = CsrMatrix::from_dense(&[[1.0, 2.0], [3.0, 4.0]]).unwrap();
    assert_eq!(
        CsrMatrix::identity(2).kron(&block).unwrap().to_dense(),
        vec![
            vec![1.0, 2.0, 0.0, 0.0],
            vec![3.0, 4.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 2.0],
            vec![0.0, 0.0, 3.0, 4.0],
        ]
    );
}

#[test]
fn kron_mixed_product_property() {
    // (A ⊗ B)(C ⊗ D) = (AC) ⊗ (BD)
    let (a, b) = (pseudo_random(2, 3, 41, 70), pseudo_random(3, 2, 42, 70));
    let (c, d) = (pseudo_random(3, 2, 43, 70), pseudo_random(2, 3, 44, 70));
    let left = &a.kron(&b).unwrap() * &c.kron(&d).unwrap();
    let right = (&a * &c).kron(&(&b * &d)).unwrap();
    assert_eq!(left, right);
}

#[test]
fn kron_empty_operand_and_overflow() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0]]).unwrap();
    assert_eq!(a.kron(&CsrMatrix::default()).unwrap().shape(), (0, 0));
    // A CSC matrix can have huge rows over a single column, needing only two pointers
    let huge = CscMatrix::<f64>::from_triplets(1 << 33, 1, Vec::new()).unwrap();
    assert_eq!(huge.kron(&huge).err(), Some(SparseError::ShapeOverflow));
}

#[test]
fn hstack_and_vstack_in_both_orientations() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0], [3.0, 4.0]]).unwrap();
    let b = CsrMatrix::from_dense(&[[5.0], [6.0]]).unwrap();
    let c = CsrMatrix::from_dense(&[[7.0, 8.0]]).unwrap();
    let wide = vec![vec![1.0, 2.0, 5.0], vec![3.0, 4.0, 6.0]];
    let tall = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![7.0, 8.0]];
    assert_eq!(CsrMatrix::hstack([&a, &b]).unwrap().to_dense(), wide);
    assert_eq!(CsrMatrix::vstack([&a, &c]).unwrap().to_dense(), tall);
    let (a_csc, b_csc, c_csc) = (a.to_csc(), b.to_csc(), c.to_csc());
    assert_eq!(
        CscMatrix::hstack([&a_csc, &b_csc]).unwrap().to_dense(),
        wide
    );
    assert_eq!(
        CscMatrix::vstack([&a_csc, &c_csc]).unwrap().to_dense(),
        tall
    );
    // Owned blocks work too
    assert_eq!(
        CsrMatrix::hstack(vec![a.clone(), b.clone()])
            .unwrap()
            .to_dense(),
        wide
    );
}

#[test]
fn stack_reassembles_submatrices() {
    let original = pseudo_random(5, 6, 55, 45);
    for matrix in [original.clone(), original.to_csc().to_csr()] {
        let left = matrix.submatrix(0..5, 0..2).unwrap();
        let right = matrix.submatrix(0..5, 2..6).unwrap();
        assert_eq!(CsrMatrix::hstack([&left, &right]).unwrap(), original);
        let top = matrix.submatrix(0..3, 0..6).unwrap();
        let bottom = matrix.submatrix(3..5, 0..6).unwrap();
        assert_eq!(CsrMatrix::vstack([&top, &bottom]).unwrap(), original);
    }
    let csc = original.to_csc();
    let left = csc.submatrix(0..5, 0..2).unwrap();
    let right = csc.submatrix(0..5, 2..6).unwrap();
    assert_eq!(CscMatrix::hstack([&left, &right]).unwrap(), original);
}

#[test]
fn stack_errors_and_empty() {
    let a = CsrMatrix::<f64>::from_triplets(2, 2, Vec::new()).unwrap();
    let b = CsrMatrix::<f64>::from_triplets(3, 1, Vec::new()).unwrap();
    let mismatch = Some(SparseError::ShapeMismatch {
        left: (2, 2),
        right: (3, 1),
    });
    assert_eq!(CsrMatrix::hstack([&a, &b]).err(), mismatch);
    assert_eq!(CsrMatrix::vstack([&a, &b]).err(), mismatch);

    let none: [&CsrMatrix<f64>; 0] = [];
    assert_eq!(CsrMatrix::hstack(none).unwrap().shape(), (0, 0));
    assert_eq!(CsrMatrix::vstack(none).unwrap().shape(), (0, 0));
    assert_eq!(CsrMatrix::vstack([&a]).unwrap(), a);
    // Blocks with zero columns still contribute their rows
    let tall_empty = CsrMatrix::<f64>::from_triplets(2, 0, Vec::new()).unwrap();
    assert_eq!(
        CsrMatrix::vstack([&tall_empty, &tall_empty])
            .unwrap()
            .shape(),
        (4, 0)
    );
    let tall_empty = tall_empty.to_csc();
    assert_eq!(
        CscMatrix::vstack([&tall_empty, &tall_empty])
            .unwrap()
            .shape(),
        (4, 0)
    );
}

#[test]
fn scale_rows_and_cols_match_diagonal_products() {
    let a = pseudo_random(5, 4, 61, 50);
    let row_factors = [2.0, -1.0, 0.5, 3.0, -4.0];
    let col_factors = [2.0, -1.0, 0.5, 3.0];
    let by_rows = CsrMatrix::from_diagonal(row_factors).matmul(&a).unwrap();
    let by_cols = a.matmul(&CsrMatrix::from_diagonal(col_factors)).unwrap();
    assert_eq!(a.scale_rows(&row_factors).unwrap(), by_rows);
    assert_eq!(a.scale_cols(&col_factors).unwrap(), by_cols);
    assert_eq!(a.to_csc().scale_rows(&row_factors).unwrap(), by_rows);
    assert_eq!(a.to_csc().scale_cols(&col_factors).unwrap(), by_cols);
}

#[test]
fn scale_rows_and_cols_keep_the_pattern() {
    let a = CsrMatrix::from_dense(&[[1.0, 2.0], [3.0, 4.0]]).unwrap();
    let rows = a.scale_rows(&[0.0, 1.0]).unwrap();
    assert_eq!(rows.to_dense(), vec![vec![0.0, 0.0], vec![3.0, 4.0]]);
    assert_eq!(rows.num_stored(), 4);
    let cols = a.scale_cols(&[1.0, 0.0]).unwrap();
    assert_eq!(cols.to_dense(), vec![vec![1.0, 0.0], vec![3.0, 0.0]]);
    assert_eq!(a.values(), [1.0, 2.0, 3.0, 4.0]);
    let empty = CsrMatrix::<f64>::from_triplets(0, 3, Vec::new()).unwrap();
    assert_eq!(empty.scale_rows(&[]).unwrap().shape(), (0, 3));

    let rect = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        rect.scale_rows(&[1.0; 3]).err(),
        Some(SparseError::ShapeMismatch {
            left: (3, 3),
            right: (2, 3)
        })
    );
    assert_eq!(
        rect.scale_cols(&[1.0; 2]).err(),
        Some(SparseError::ShapeMismatch {
            left: (2, 3),
            right: (2, 2)
        })
    );
}

#[test]
fn scale_rows_and_cols_complex() {
    // diag(V) · A with complex voltages, as in power flow Jacobians
    let a = CsrMatrix::from_triplets(
        2,
        2,
        vec![
            (0, 0, Complex::new(1.0, 1.0)),
            (1, 1, Complex::new(2.0, 0.0)),
        ],
    )
    .unwrap();
    let v = [Complex::new(0.0, 1.0), Complex::new(1.0, -1.0)];
    assert_eq!(
        a.scale_rows(&v).unwrap().values(),
        [Complex::new(-1.0, 1.0), Complex::new(2.0, -2.0)]
    );
    assert_eq!(
        a.scale_cols(&v).unwrap(),
        a.matmul(&CsrMatrix::from_diagonal(v)).unwrap()
    );
}

#[test]
fn real_imag_and_from_parts() {
    let a = pseudo_random_complex(6, 5, 91);
    let (re, im) = (a.real(), a.imag());
    assert_eq!((re.shape(), im.shape()), ((6, 5), (6, 5)));
    assert_eq!(CsrMatrix::from_parts(&re, &im).unwrap(), a);

    let b = CsrMatrix::from_triplets(
        2,
        2,
        vec![
            (0, 0, Complex::new(1.0, 0.0)),
            (0, 1, Complex::new(0.0, 2.0)),
            (1, 0, Complex::new(3.0, -4.0)),
        ],
    )
    .unwrap();
    assert_eq!(b.real().to_dense(), vec![vec![1.0, 0.0], vec![3.0, 0.0]]);
    assert_eq!(b.imag().to_dense(), vec![vec![0.0, 2.0], vec![-4.0, 0.0]]);
    // The parts keep the pattern, so a purely imaginary entry is a stored zero in the real part
    assert_eq!((b.real().num_stored(), b.real().num_nonzero()), (3, 2));

    let wrong = CsrMatrix::<f64>::from_triplets(3, 2, Vec::new()).unwrap();
    let re = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        CsrMatrix::from_parts(&re, &wrong).err(),
        Some(SparseError::ShapeMismatch {
            left: (2, 3),
            right: (3, 2)
        })
    );
}

#[test]
fn conj_and_adjoint() {
    let a = CsrMatrix::from_triplets(
        2,
        2,
        vec![
            (0, 0, Complex::new(1.0, 2.0)),
            (0, 1, Complex::new(0.0, -1.0)),
            (1, 1, Complex::new(3.0, 0.0)),
        ],
    )
    .unwrap();
    let conjugated = a.conj();
    assert_eq!(
        conjugated.values()[..2],
        [Complex::new(1.0, -2.0), Complex::new(0.0, 1.0)]
    );
    assert_eq!(conjugated.conj(), a);

    // The adjoint of a CSR matrix is a CSC matrix, like the transpose
    let adjoint: CscMatrix<Complex<f64>> = a.adjoint();
    assert_eq!(adjoint, a.conj().transpose());
    assert_eq!(adjoint.to_dense()[1][0], Complex::new(0.0, 1.0));
    assert_eq!(adjoint.adjoint(), a);

    // (AB)ᴴ = Bᴴ Aᴴ, on a rectangular pair
    let (c, d) = (
        pseudo_random_complex(4, 3, 92),
        pseudo_random_complex(3, 5, 93),
    );
    assert_eq!((&c * &d).adjoint(), &d.adjoint() * &c.adjoint());

    // A + Aᴴ is Hermitian, with a real diagonal
    let e = pseudo_random_complex(5, 5, 94);
    let hermitian = &e + &e.adjoint().to_csr();
    assert_eq!(hermitian.adjoint(), hermitian);
    assert!(hermitian
        .diagonal()
        .values()
        .iter()
        .all(|value| value.im == 0.0));
}

#[test]
fn complex_helpers_on_real_matrices() {
    let a = pseudo_random(4, 6, 95, 50);
    assert_eq!(a.conj(), a);
    assert_eq!(a.adjoint(), a.clone().transpose());
    assert_eq!(a.real(), a);
    let imaginary = a.imag();
    assert_eq!(imaginary.shape(), (4, 6));
    assert_eq!(
        (imaginary.num_stored(), imaginary.num_nonzero()),
        (a.num_stored(), 0)
    );
}

#[test]
fn subtraction_stores_cancellations_and_leaves_operands() {
    let b = CsrMatrix::from_dense(&[[1.0, 0.0], [4.0, 3.0]]).unwrap();
    let difference = &sample_2x2() - &b;
    assert_eq!(difference.to_dense(), vec![vec![0.0, 2.0], vec![-4.0, 0.0]]);
    assert_eq!((difference.num_stored(), difference.num_nonzero()), (4, 2));
    assert_eq!(
        sample_2x2(),
        CsrMatrix::from_dense(&[[1.0, 2.0], [0.0, 3.0]]).unwrap()
    );
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn sub_operator_shape_mismatch_panics() {
    let _ = &sample_2x2() - &CsrMatrix::<f64>::from_triplets(3, 3, Vec::new()).unwrap();
}

#[test]
fn integer_scalar_multiplication() {
    let a = CsrMatrix::from_dense(&[[1i32, 2], [0, 3]]).unwrap();
    assert_eq!((&a * 4).to_dense(), vec![vec![4, 8], vec![0, 12]]);
    assert_eq!((3i32 * &a).to_dense(), vec![vec![3, 6], vec![0, 9]]);
}

// Assembly maps

/// `count` triplets at pseudo-random positions of a `rows x cols` matrix, drawn with
/// replacement so many positions repeat, with inexact values, so summation order shows in the
/// bits. `seed` picks the positions and `value_seed` the values, so two lists can share
/// positions.
fn random_triplets(
    rows: usize,
    cols: usize,
    count: usize,
    seed: u64,
    value_seed: u64,
) -> Vec<(usize, usize, f64)> {
    let lcg = |state: &mut u64| {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state >> 33
    };
    let (mut positions, mut values) = (seed, value_seed);
    (0..count)
        .map(|_| {
            let row = lcg(&mut positions) as usize % rows;
            let col = lcg(&mut positions) as usize % cols;
            (row, col, (lcg(&mut values) % 2001) as f64 / 7.3 - 137.0)
        })
        .collect()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

/// Checks that `matrix` holds exactly the arrays `from_triplets` builds from `triplets`.
fn assert_same_arrays<O: sparse_mat::compressed_matrix::Orientation>(
    matrix: &CompressedMatrix<f64, O>,
    triplets: &[(usize, usize, f64)],
) {
    let (rows, cols) = matrix.shape();
    let fresh = CompressedMatrix::<f64, O>::from_triplets(rows, cols, triplets.to_vec()).unwrap();
    assert_eq!(matrix.pointers(), fresh.pointers());
    assert_eq!(matrix.indices(), fresh.indices());
    assert_eq!(bits(matrix.values()), bits(fresh.values()));
}

fn refills_match_from_triplets<O: sparse_mat::compressed_matrix::Orientation>() {
    let (rows, cols) = (13, 9);
    let triplets = random_triplets(rows, cols, 200, 1, 2);
    let (mut matrix, map) =
        CompressedMatrix::<f64, O>::from_triplets_mapped(rows, cols, triplets.clone()).unwrap();
    assert_same_arrays(&matrix, &triplets);
    assert_eq!(map.len(), 200);
    assert_eq!(map.shape(), (rows, cols));
    assert!(matrix.num_stored() < 200, "the positions should repeat");
    for (slot, &(row, col, _)) in map.slots().zip(&triplets) {
        assert_eq!(matrix.position(row, col), Ok(Some(slot)));
    }

    // New values for the same positions, several times over
    for value_seed in 3..8 {
        let refilled = random_triplets(rows, cols, 200, 1, value_seed);
        let values: Vec<f64> = refilled.iter().map(|triplet| triplet.2).collect();
        matrix.refill(&map, &values).unwrap();
        assert_same_arrays(&matrix, &refilled);
    }
}

#[test]
fn refills_match_from_triplets_bit_for_bit() {
    refills_match_from_triplets::<sparse_mat::compressed_matrix::Csr>();
    refills_match_from_triplets::<sparse_mat::compressed_matrix::Csc>();
}

#[test]
fn refills_keep_the_pattern_and_signed_zeros() {
    let triplets = vec![(0, 0, 1.0), (1, 1, 2.0), (0, 0, 3.0), (1, 0, 5.0)];
    let (mut csr, map) = CsrMatrix::from_triplets_mapped(2, 2, triplets).unwrap();
    assert_eq!(csr.to_dense(), vec![vec![4.0, 0.0], vec![5.0, 2.0]]);
    assert_eq!(map.slots().collect::<Vec<_>>(), [0, 2, 0, 1]);

    // Values that cancel, or are zero, stay stored
    csr.refill(&map, &[1.0, 0.0, -1.0, 6.0]).unwrap();
    assert_eq!(csr.num_stored(), 3);
    assert_eq!(csr.num_nonzero(), 1);
    assert_eq!(csr.values(), [0.0, 6.0, 0.0]);

    // A position's first triplet is written, not added to zero, so -0.0 survives as it does in
    //   from_triplets (0.0 + -0.0 would be +0.0)
    csr.refill(&map, &[-0.0, -0.0, -0.0, 1.0]).unwrap();
    let fresh = CsrMatrix::from_triplets(
        2,
        2,
        vec![(0, 0, -0.0), (1, 1, -0.0), (0, 0, -0.0), (1, 0, 1.0)],
    )
    .unwrap();
    assert_eq!(bits(csr.values()), bits(fresh.values()));
    assert!(csr.values()[0].is_sign_negative());

    // A clone refills with the same map
    let mut clone = csr.clone();
    clone.refill(&map, &[1.0, 2.0, 3.0, 4.0]).unwrap();
    assert_eq!(clone.to_dense(), vec![vec![4.0, 0.0], vec![4.0, 2.0]]);
}

#[test]
fn refill_errors_leave_the_matrix_unchanged() {
    let (mut csc, map) =
        CscMatrix::from_triplets_mapped(2, 3, vec![(0, 0, 1.0), (1, 2, 2.0)]).unwrap();
    let before = csc.clone();

    let error = csc.refill(&map, &[1.0]).unwrap_err();
    assert_eq!(
        error,
        SparseError::LengthMismatch {
            expected: 2,
            found: 1
        }
    );
    assert_eq!(error.to_string(), "expected 2 elements, found 1");

    // Maps from matrices with another shape or another stored count
    let (_, other_shape) =
        CscMatrix::from_triplets_mapped(3, 2, vec![(0, 0, 1.0), (1, 1, 2.0)]).unwrap();
    let (_, other_count) =
        CscMatrix::from_triplets_mapped(2, 3, vec![(0, 0, 1.0), (1, 2, 2.0), (1, 1, 3.0)]).unwrap();
    for other in [&other_shape, &other_count] {
        let error = csc.refill(other, &[1.0, 2.0]).unwrap_err();
        assert_eq!(error, SparseError::PatternMismatch);
        assert_eq!(
            error.to_string(),
            "the assembly map was built for a different pattern"
        );
    }
    assert_eq!(bits(csc.values()), bits(before.values()));

    assert_eq!(
        CsrMatrix::from_triplets_mapped(2, 2, vec![(0, 0, 1.0), (2, 0, 1.0)]).err(),
        Some(SparseError::OutOfBounds {
            row: 2,
            col: 0,
            shape: (2, 2)
        })
    );
    let (mut empty, map) = CsrMatrix::<f64>::from_triplets_mapped(3, 3, Vec::new()).unwrap();
    assert!(map.is_empty() && empty.num_stored() == 0);
    empty.refill(&map, &[]).unwrap();
}
