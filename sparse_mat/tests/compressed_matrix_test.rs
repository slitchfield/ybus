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
