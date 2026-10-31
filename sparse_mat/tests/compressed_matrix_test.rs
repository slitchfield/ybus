// Tests for CompressedMatrix (CsrMatrix and CscMatrix): construction, validation, the fixed
//   pattern and its stored zeros, transposing and reorienting, equality and formatting.

use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix};
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
