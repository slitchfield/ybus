use sparse_mat::{sparse_error::SparseError, sparse_matrix::SparseMatrix};

#[test]
fn sparsemat_creation() {
    let _local = SparseMatrix::<f64>::new();
    let _local2 = SparseMatrix::<f64>::empty_with_shape(3, 3);
}

#[test]
fn sparsemat_creation_huge_shape() {
    // Must not reserve space proportional to rows * cols
    let mut local = SparseMatrix::empty_with_shape(1 << 32, 1 << 32);
    local.insert(5, 7, 1.0).unwrap();
    assert!(local.num_nonzero() == 1);
}

#[test]
fn sparsemat_identity_creation() {
    let local = SparseMatrix::identity(3);
    assert!(local.num_nonzero() == 3);
    assert!(local.peek_at(0, 0).unwrap() == Some(1.0));
    assert!(local.peek_at(1, 1).unwrap() == Some(1.0));
    assert!(local.peek_at(2, 2).unwrap() == Some(1.0));
}

#[test]
fn sparsemat_transpose_creation() {
    let mut local = SparseMatrix::empty_with_shape(4, 6);
    local
        .insert_triplets(vec![
            (0, 0, 10.0),
            (0, 1, 20.0),
            (1, 1, 30.0),
            (2, 2, 50.0),
            (1, 3, 40.0),
            (2, 3, 60.0),
            (2, 4, 70.0),
            (3, 5, 80.0),
        ])
        .unwrap();

    let local2 = local.create_transpose();
    assert!(local2.shape == (6, 4));
    assert!(local2.peek_at(0, 0).unwrap() == Some(10.0));
    assert!(local2.peek_at(0, 1).unwrap().is_none());
    assert!(local2.peek_at(1, 0).unwrap() == Some(20.0));
    assert!(local2.peek_at(1, 1).unwrap() == Some(30.0));
    assert!(local2.peek_at(2, 2).unwrap() == Some(50.0));
    assert!(local2.peek_at(3, 1).unwrap() == Some(40.0));
    assert!(local2.peek_at(1, 3).unwrap().is_none());
    assert!(local2.peek_at(3, 2).unwrap() == Some(60.0));
    assert!(local2.peek_at(2, 3).unwrap().is_none());
    assert!(local2.peek_at(4, 2).unwrap() == Some(70.0));
    assert!(local2.peek_at(5, 3).unwrap() == Some(80.0));
}

#[test]
fn sparsemat_rowiter() {
    let mut local = SparseMatrix::empty_with_shape(4, 6);
    local
        .insert_triplets(vec![
            (0, 0, 10.0),
            (0, 1, 20.0),
            (1, 1, 30.0),
            (2, 2, 50.0),
            (1, 3, 40.0),
            (2, 3, 60.0),
            (2, 4, 70.0),
            (3, 5, 80.0),
        ])
        .unwrap();

    local.explicitly_compress();

    for row in local.row_iter() {
        dbg!(row);
    }
}

#[test]
fn sparsemat_rowiter_uncompressed() {
    let mut local = SparseMatrix::empty_with_shape(2, 3);
    local
        .insert_triplets(vec![(0, 2, 1.0), (1, 0, 2.0)])
        .unwrap();

    let rows: Vec<Vec<f64>> = local.row_iter().collect();
    assert!(rows == vec![vec![0.0, 0.0, 1.0], vec![2.0, 0.0, 0.0]]);
}

#[test]
fn sparsemat_rowiter_stale_after_compress() {
    let mut local = SparseMatrix::empty_with_shape(2, 2);
    local.insert(0, 0, 1.0).unwrap();
    local.explicitly_compress();
    local.insert(1, 1, 4.0).unwrap();

    let rows: Vec<Vec<f64>> = local.row_iter().collect();
    assert!(rows == vec![vec![1.0, 0.0], vec![0.0, 4.0]]);
}

#[test]
fn sparsemat_insert() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    local.insert(0, 0, 1.0).unwrap();
}

#[test]
fn sparsemat_insert_oob() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    assert_eq!(
        local.insert(4, 4, 1.0),
        Err(SparseError::OutOfBounds {
            row: 4,
            col: 4,
            shape: (3, 3)
        })
    );
}

#[test]
fn sparsemat_insert_triplets() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    local
        .insert_triplets(vec![(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)])
        .unwrap();
    assert!(local.peek_at(0, 0) == Ok(Some(1.0)));
    assert!(local.peek_at(1, 1) == Ok(Some(2.0)));
    assert!(local.peek_at(2, 2) == Ok(Some(3.0)));
}

#[test]
fn sparsemat_insert_triplets_oob() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    assert!(local.num_nonzero() == 0);
    assert!(
        local.insert_triplets(vec![(0, 0, 2.0), (4, 4, 4.0)])
            == Err(SparseError::OutOfBounds {
                row: 4,
                col: 4,
                shape: (3, 3)
            })
    );
    // Check to see we didn't partially insert the first triplet
    assert!(local.num_nonzero() == 0);
}

#[test]
fn sparsemat_remove() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    local.insert(0, 0, 1.0).unwrap();

    let output = local.clear_at(0, 0);
    assert!(output == Ok(Some(1.0)));

    let output = local.clear_at(0, 0);
    assert!(output.is_ok());
    assert!(output.unwrap().is_none());

    let output = local.clear_at(1, 1);
    assert!(output.is_ok());
    assert!(output.unwrap().is_none());
}

#[test]
fn sparsemat_remove_oob() {
    let mut local = SparseMatrix::<f64>::empty_with_shape(3, 3);

    assert!(
        local.clear_at(4, 4)
            == Err(SparseError::OutOfBounds {
                row: 4,
                col: 4,
                shape: (3, 3)
            })
    );
}

#[test]
fn sparsemat_peek() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);

    local
        .insert_triplets(vec![(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)])
        .unwrap();
    assert!(local.peek_at(0, 0).unwrap() == Some(1.0));
    assert!(local.peek_at(0, 1).unwrap().is_none());
}

#[test]
fn sparsemat_peek_oob() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);

    local
        .insert_triplets(vec![(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)])
        .unwrap();
    assert!(
        local.peek_at(4, 4)
            == Err(SparseError::OutOfBounds {
                row: 4,
                col: 4,
                shape: (3, 3)
            })
    );
}

#[test]
fn sparsemat_transposeinplace() {
    let mut local = SparseMatrix::empty_with_shape(4, 4);

    local
        .insert_triplets(vec![(0, 0, 1.0), (1, 0, 2.0), (2, 2, 3.0), (2, 3, 4.0)])
        .unwrap();
    local.transpose_inplace();

    assert!(local.peek_at(0, 0).unwrap() == Some(1.0));
    assert!(local.peek_at(0, 1).unwrap() == Some(2.0));
    assert!(local.peek_at(1, 0).unwrap().is_none());
    assert!(local.peek_at(2, 2).unwrap() == Some(3.0));
    assert!(local.peek_at(3, 2).unwrap() == Some(4.0));
    assert!(local.peek_at(2, 3).unwrap().is_none());
}

#[test]
fn sparsemat_compressedrepr() {
    let mut local = SparseMatrix::empty_with_shape(4, 4);
    local
        .insert_triplets(vec![(0, 0, 5.0), (1, 1, 8.0), (3, 1, 6.0), (2, 2, 3.0)])
        .unwrap();
    assert!(local.compressed().is_none());
    local.explicitly_compress();
    assert!(local.compressed().unwrap().data == vec![5.0, 8.0, 3.0, 6.0]);
    assert!(local.compressed().unwrap().cols == vec![0, 1, 2, 1]);
    assert!(local.compressed().unwrap().rows == vec![0, 1, 2, 3, 4]);

    let mut local2 = SparseMatrix::empty_with_shape(4, 6);
    local2
        .insert_triplets(vec![
            (0, 0, 10.0),
            (0, 1, 20.0),
            (1, 1, 30.0),
            (2, 2, 50.0),
            (1, 3, 40.0),
            (2, 3, 60.0),
            (2, 4, 70.0),
            (3, 5, 80.0),
        ])
        .unwrap();
    local2.explicitly_compress();
    assert!(
        local2.compressed().unwrap().data == vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
    );
    assert!(local2.compressed().unwrap().cols == vec![0, 1, 1, 3, 2, 3, 4, 5]);
    assert!(local2.compressed().unwrap().rows == vec![0, 2, 4, 7, 8]);

    // Unwrap as OOB not a part of this test
    local2.insert(0, 2, 1.0).unwrap();
    assert!(local2.compressed().is_none());
}

#[test]
fn sparsemat_display() {
    let mut local = SparseMatrix::empty_with_shape(4, 6);
    local
        .insert_triplets(vec![
            (0, 0, 10.0),
            (0, 1, 20.0),
            (1, 1, 30.0),
            (2, 2, 50.0),
            (1, 3, 40.0),
            (2, 3, 60.0),
            (2, 4, 70.0),
            (3, 5, 80.0),
        ])
        .unwrap();
    local.explicitly_compress();
    println!("{}", local)
}

#[test]
fn sparsemat_display_uncompressed() {
    let mut local = SparseMatrix::empty_with_shape(2, 2);
    local
        .insert_triplets(vec![(0, 0, 1.0), (1, 1, 2.0)])
        .unwrap();
    let before = format!("{}", local);

    local.explicitly_compress();
    assert!(before == format!("{}", local));
    assert!(before.contains("1.00"));
}

#[test]
#[should_panic]
fn sparsemat_bad_addition() {
    let local = SparseMatrix::<f64>::empty_with_shape(3, 3);
    let local2 = SparseMatrix::empty_with_shape(2, 2);

    let _local3 = &local + &local2;
}

#[test]
fn sparsemat_good_addition() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    local
        .insert_triplets(vec![(0, 0, 10.0), (0, 1, 20.0), (1, 1, 30.0), (2, 2, 50.0)])
        .unwrap();
    let local2 = local.create_transpose();

    let local3 = &local + &local2;
    assert!(local3.peek_at(0, 0).unwrap() == Some(20.0));
    assert!(local3.peek_at(0, 1).unwrap() == Some(20.0));
    assert!(local3.peek_at(1, 0).unwrap() == Some(20.0));
    assert!(local3.peek_at(1, 1).unwrap() == Some(60.0));
    assert!(local3.peek_at(2, 2).unwrap() == Some(100.0));
}

#[test]
fn sparsemat_insert_zero_into_empty() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    // Unwrap as OOB not a part of this test
    local.insert(0, 0, 0.0).unwrap();
    assert!(local.peek_at(0, 0).unwrap().is_none());
    assert!(local.num_nonzero() == 0);
}

#[test]
fn sparsemat_insert_zero_over_existing() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    // Unwrap as OOB not a part of this test
    local.insert(0, 0, 5.0).unwrap();
    local.insert(0, 0, 0.0).unwrap();
    assert!(local.peek_at(0, 0).unwrap().is_none());
    assert!(local.num_nonzero() == 0);
}

#[test]
fn sparsemat_insert_triplets_zeros() {
    let mut local = SparseMatrix::empty_with_shape(3, 3);
    local
        .insert_triplets(vec![(1, 1, 2.0), (2, 2, 3.0)])
        .unwrap();
    local
        .insert_triplets(vec![(1, 1, 0.0), (0, 2, 0.0)])
        .unwrap();
    assert!(local.peek_at(1, 1).unwrap().is_none());
    assert!(local.peek_at(0, 2).unwrap().is_none());
    assert!(local.peek_at(2, 2).unwrap() == Some(3.0));
    assert!(local.num_nonzero() == 1);
}

#[test]
fn sparsemat_insert_zero_invalidates_compressed() {
    let mut local = SparseMatrix::empty_with_shape(2, 2);
    local
        .insert_triplets(vec![(0, 1, 4.0), (1, 1, 3.0)])
        .unwrap();
    local.explicitly_compress();
    // Unwrap as OOB not a part of this test
    local.insert(0, 1, 0.0).unwrap();
    assert!(local.compressed().is_none());

    let rows: Vec<Vec<f64>> = local.row_iter().collect();
    assert!(rows == vec![vec![0.0, 0.0], vec![0.0, 3.0]]);
}

#[test]
fn sparsemat_addition_cancellation() {
    let mut local = SparseMatrix::empty_with_shape(2, 2);
    local
        .insert_triplets(vec![(0, 0, 1.5), (1, 1, 2.0)])
        .unwrap();
    let mut local2 = SparseMatrix::empty_with_shape(2, 2);
    local2
        .insert_triplets(vec![(0, 0, -1.5), (0, 1, 1.0)])
        .unwrap();

    let local3 = &local + &local2;
    assert!(local3.peek_at(0, 0).unwrap().is_none());
    assert!(local3.peek_at(0, 1).unwrap() == Some(1.0));
    assert!(local3.peek_at(1, 1).unwrap() == Some(2.0));
    assert!(local3.num_nonzero() == 2);
}

#[test]
fn sparsemat_prune_zeroes_noop() {
    let mut local = SparseMatrix::empty_with_shape(2, 2);
    local
        .insert_triplets(vec![(0, 0, 1.0), (1, 1, 2.0)])
        .unwrap();
    local.explicitly_compress();

    local.prune_zeroes();
    assert!(local.num_nonzero() == 2);
    assert!(local.compressed().is_some());
}

#[test]
fn sparsemat_from_triplets_basic() {
    let local = SparseMatrix::from_triplets(3, 4, vec![(0, 0, 1.0), (2, 3, 7.0)]).unwrap();
    assert!(local.shape == (3, 4));
    assert!(local.peek_at(0, 0) == Ok(Some(1.0)));
    assert!(local.peek_at(2, 3) == Ok(Some(7.0)));
    assert!(local.peek_at(1, 1) == Ok(None));
    assert!(local.num_nonzero() == 2);
    assert!(local.compressed().is_none());
}

#[test]
fn sparsemat_from_triplets_sums_duplicates() {
    let local =
        SparseMatrix::from_triplets(2, 2, vec![(0, 0, 1.0), (1, 1, 4.0), (0, 0, 2.5)]).unwrap();
    assert!(local.peek_at(0, 0) == Ok(Some(3.5)));
    assert!(local.num_nonzero() == 2);
}

#[test]
fn sparsemat_from_triplets_cancellation_removed() {
    let local = SparseMatrix::from_triplets(2, 2, vec![(1, 1, 5.0), (1, 1, -5.0)]).unwrap();
    assert!(local.peek_at(1, 1) == Ok(None));
    assert!(local.num_nonzero() == 0);
}

#[test]
fn sparsemat_from_triplets_skips_zero_inputs() {
    let local = SparseMatrix::from_triplets(2, 2, vec![(0, 0, 0.0), (1, 0, 2.0)]).unwrap();
    assert!(local.peek_at(0, 0) == Ok(None));
    assert!(local.num_nonzero() == 1);
}

#[test]
fn sparsemat_from_triplets_oob() {
    let result = SparseMatrix::from_triplets(3, 3, vec![(0, 0, 1.0), (4, 4, 2.0)]);
    assert!(
        result.err()
            == Some(SparseError::OutOfBounds {
                row: 4,
                col: 4,
                shape: (3, 3)
            })
    );
}

#[test]
fn sparsemat_from_triplets_empty() {
    let local = SparseMatrix::<f64>::from_triplets(3, 5, vec![]).unwrap();
    assert!(local.shape == (3, 5));
    assert!(local.num_nonzero() == 0);
}

#[test]
fn sparsemat_from_triplets_huge_shape() {
    // Capacity must come from the triplet count, not rows * cols
    let local = SparseMatrix::from_triplets(1 << 32, 1 << 32, [(0, 0, 1.0)]).unwrap();
    assert!(local.num_nonzero() == 1);
}

#[test]
fn sparsemat_from_triplets_row_iter() {
    let local =
        SparseMatrix::from_triplets(2, 3, vec![(0, 2, 1.0), (1, 0, 2.0), (0, 2, 1.0)]).unwrap();
    let rows: Vec<Vec<f64>> = local.row_iter().collect();
    assert!(rows == vec![vec![0.0, 0.0, 2.0], vec![2.0, 0.0, 0.0]]);
}

#[test]
fn sparsemat_from_dense_basic() {
    let local = SparseMatrix::from_dense(vec![vec![1.0, 0.0, 2.0], vec![0.0, 3.0, 0.0]]).unwrap();
    assert!(local.shape == (2, 3));
    assert!(local.peek_at(0, 0) == Ok(Some(1.0)));
    assert!(local.peek_at(0, 2) == Ok(Some(2.0)));
    assert!(local.peek_at(1, 1) == Ok(Some(3.0)));
    assert!(local.num_nonzero() == 3);
    assert!(local.compressed().is_none());
}

#[test]
fn sparsemat_from_dense_skips_zeros() {
    let local = SparseMatrix::from_dense(vec![vec![0.0, -0.0], vec![0.0, 5.0]]).unwrap();
    assert!(local.peek_at(0, 0) == Ok(None));
    assert!(local.peek_at(0, 1) == Ok(None));
    assert!(local.num_nonzero() == 1);
}

#[test]
fn sparsemat_from_dense_ragged() {
    let result = SparseMatrix::from_dense(vec![vec![1.0, 2.0], vec![3.0], vec![4.0, 5.0]]);
    assert!(
        result.err()
            == Some(SparseError::RaggedRows {
                row: 1,
                expected: 2,
                found: 1
            })
    );
}

#[test]
fn sparsemat_from_dense_empty() {
    let local = SparseMatrix::from_dense(Vec::<Vec<f64>>::new()).unwrap();
    assert!(local.shape == (0, 0));
    assert!(local.num_nonzero() == 0);

    // Rows with no columns are a valid n x 0 matrix
    let local = SparseMatrix::<f64>::from_dense(vec![vec![], vec![]]).unwrap();
    assert!(local.shape == (2, 0));
}

#[test]
fn sparsemat_from_dense_round_trip() {
    let mut original = SparseMatrix::empty_with_shape(3, 4);
    original
        .insert_triplets(vec![(0, 0, 1.0), (1, 3, -2.0), (2, 1, 4.5)])
        .unwrap();

    let rebuilt = SparseMatrix::from_dense(original.row_iter()).unwrap();
    assert!(rebuilt.shape == original.shape);
    assert!(rebuilt.num_nonzero() == original.num_nonzero());
    for (r, c, v) in [(0, 0, 1.0), (1, 3, -2.0), (2, 1, 4.5)] {
        assert!(rebuilt.peek_at(r, c) == Ok(Some(v)));
    }
}

#[test]
fn sparsemat_from_dense_array_input() {
    let local = SparseMatrix::from_dense([[1.0, 0.0], [0.0, 2.0]]).unwrap();
    assert!(local.shape == (2, 2));
    assert!(local.num_nonzero() == 2);
}

#[test]
fn sparsemat_from_dense_nan_is_stored() {
    let local = SparseMatrix::from_dense(vec![vec![f64::NAN, 0.0]]).unwrap();
    assert!(local.num_nonzero() == 1);
    assert!(local.peek_at(0, 0).unwrap().unwrap().is_nan());
}

#[test]
fn sparsemat_to_dense_basic() {
    let mut local = SparseMatrix::empty_with_shape(2, 3);
    local
        .insert_triplets(vec![(0, 2, 1.0), (1, 0, 2.0)])
        .unwrap();
    assert!(local.to_dense() == vec![vec![0.0, 0.0, 1.0], vec![2.0, 0.0, 0.0]]);
}

#[test]
fn sparsemat_to_dense_matches_row_iter() {
    let local =
        SparseMatrix::from_triplets(3, 4, vec![(0, 0, 1.0), (2, 3, 4.0), (1, 1, -2.0)]).unwrap();
    let rows: Vec<Vec<f64>> = local.row_iter().collect();
    assert!(local.to_dense() == rows);
}

#[test]
fn sparsemat_to_dense_empty_shapes() {
    assert!(SparseMatrix::<f64>::new().to_dense().is_empty());
    assert!(SparseMatrix::<f64>::empty_with_shape(2, 0).to_dense() == vec![Vec::<f64>::new(); 2]);
    // No rows means the column count is lost
    assert!(SparseMatrix::<f64>::empty_with_shape(0, 3)
        .to_dense()
        .is_empty());
}

#[test]
fn sparsemat_to_dense_round_trip() {
    let dense = vec![
        vec![1.0, 0.0, 2.0],
        vec![0.0, 0.0, 0.0],
        vec![0.0, 3.5, 0.0],
    ];
    let local = SparseMatrix::from_dense(dense.clone()).unwrap();
    assert!(local.to_dense() == dense);

    let rebuilt = SparseMatrix::from_dense(local.to_dense()).unwrap();
    assert!(rebuilt.shape == local.shape);
    assert!(rebuilt.num_nonzero() == local.num_nonzero());
}

#[test]
fn sparsemat_integer_elements() {
    let mut local = SparseMatrix::<i32>::empty_with_shape(2, 2);
    local.insert(0, 0, 3).unwrap();
    local.insert(1, 1, 4).unwrap();
    local.insert(1, 1, 0).unwrap();
    assert!(local.peek_at(0, 0) == Ok(Some(3)));
    assert!(local.peek_at(1, 1) == Ok(None));

    let other = SparseMatrix::from_triplets(2, 2, vec![(0, 0, 4), (0, 1, 1)]).unwrap();
    let sum = &local + &other;
    assert!(sum.to_dense() == vec![vec![7, 1], vec![0, 0]]);
    assert!(format!("{}", sum).contains('7'));
}

#[test]
fn sparsemat_f32_elements() {
    let local = SparseMatrix::<f32>::identity(3);
    assert!(local.num_nonzero() == 3);
    let doubled = &local + &local;
    assert!(doubled.peek_at(2, 2) == Ok(Some(2.0f32)));
}

#[test]
fn sparsemat_complex_elements() {
    use num_complex::Complex;
    let i = Complex::new(0.0, 1.0);
    let mut local = SparseMatrix::empty_with_shape(2, 2);
    local.insert(0, 1, i).unwrap();
    local.insert(1, 0, Complex::new(1.0, 0.0)).unwrap();
    let sum = &local + &local;
    assert!(sum.peek_at(0, 1) == Ok(Some(Complex::new(0.0, 2.0))));
    assert!(local.create_transpose().peek_at(1, 0) == Ok(Some(i)));
}

#[test]
fn sparsemat_compressed_csc() {
    let mut local = SparseMatrix::from_triplets(
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
    .unwrap();
    assert!(local.compressed_csc().is_none());
    local.explicitly_compress_csc();

    let csc = local.compressed_csc().unwrap();
    assert!(csc.cols == vec![0, 1, 3, 4, 6, 7, 8]);
    assert!(csc.rows == vec![0, 0, 1, 2, 1, 2, 2, 3]);
    assert!(csc.data == vec![10.0, 20.0, 30.0, 50.0, 40.0, 60.0, 70.0, 80.0]);
    // CSR was not built as a side effect
    assert!(local.compressed().is_none());

    local.insert(3, 0, 1.0).unwrap();
    assert!(local.compressed_csc().is_none());
}

#[test]
fn sparsemat_csc_matches_transpose_csr() {
    let mut local = SparseMatrix::from_triplets(
        3,
        5,
        vec![
            (0, 4, 1.0),
            (2, 0, 2.0),
            (1, 1, 3.0),
            (2, 4, 4.0),
            (0, 0, 5.0),
        ],
    )
    .unwrap();
    local.explicitly_compress_csc();
    let mut transposed = local.create_transpose();
    transposed.explicitly_compress();

    let csc = local.compressed_csc().unwrap();
    let csr = transposed.compressed().unwrap();
    assert!(csc.cols == csr.rows);
    assert!(csc.rows == csr.cols);
    assert!(csc.data == csr.data);
}

#[test]
fn sparsemat_try_add_shape_mismatch() {
    let local = SparseMatrix::<f64>::empty_with_shape(3, 3);
    let local2 = SparseMatrix::<f64>::empty_with_shape(2, 4);
    assert!(
        local.try_add(&local2).err()
            == Some(SparseError::ShapeMismatch {
                left: (3, 3),
                right: (2, 4)
            })
    );
}

#[test]
fn sparsemat_try_add_matches_operator() {
    let local = SparseMatrix::from_triplets(2, 2, vec![(0, 0, 1.5), (1, 1, 2.0)]).unwrap();
    let local2 = SparseMatrix::from_triplets(2, 2, vec![(0, 0, -1.5), (0, 1, 1.0)]).unwrap();
    let sum = local.try_add(&local2).unwrap();
    assert!(sum.to_dense() == (&local + &local2).to_dense());
    assert!(sum.num_nonzero() == 2);
    // Operands are untouched
    assert!(local.peek_at(0, 0) == Ok(Some(1.5)));
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn sparsemat_add_shape_mismatch_message() {
    let local = SparseMatrix::<f64>::empty_with_shape(3, 3);
    let local2 = SparseMatrix::<f64>::empty_with_shape(2, 2);
    let _ = &local + &local2;
}

#[test]
fn sparsemat_iter_row_major_order() {
    let local = SparseMatrix::from_triplets(
        4,
        4,
        vec![
            (3, 1, 6.0),
            (0, 3, 2.0),
            (2, 2, 5.0),
            (0, 0, 1.0),
            (3, 0, 9.0),
        ],
    )
    .unwrap();
    let entries: Vec<(usize, usize, f64)> = local.iter().collect();
    assert!(
        entries
            == vec![
                (0, 0, 1.0),
                (0, 3, 2.0),
                (2, 2, 5.0),
                (3, 0, 9.0),
                (3, 1, 6.0)
            ]
    );
}

#[test]
fn sparsemat_iter_len_and_empty() {
    let local = SparseMatrix::from_triplets(3, 3, vec![(0, 0, 1.0), (2, 2, 2.0)]).unwrap();
    let mut entries = local.iter();
    assert!(entries.len() == 2);
    entries.next();
    assert!(entries.len() == 1);

    assert!(SparseMatrix::<f64>::empty_with_shape(5, 5)
        .iter()
        .next()
        .is_none());
    assert!(SparseMatrix::<f64>::new().iter().count() == 0);
}

#[test]
fn sparsemat_iter_compressed_or_not() {
    let mut local = SparseMatrix::from_triplets(3, 3, vec![(1, 2, 4.0), (0, 1, 3.0)]).unwrap();
    let before: Vec<_> = local.iter().collect();
    local.explicitly_compress();
    let after: Vec<_> = local.iter().collect();
    assert!(before == after);
    assert!(before == vec![(0, 1, 3.0), (1, 2, 4.0)]);
}

#[test]
fn sparsemat_iter_sparse_cost_on_huge_shape() {
    let mut local = SparseMatrix::<f64>::empty_with_shape(1 << 20, 1 << 20);
    local.insert(5, 7, 1.0).unwrap();
    local.insert(1 << 19, 3, 2.0).unwrap();
    let entries: Vec<_> = local.iter().collect();
    assert!(entries == vec![(5, 7, 1.0), (1 << 19, 3, 2.0)]);
}

fn sample_4x6() -> SparseMatrix {
    SparseMatrix::from_triplets(
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

#[test]
fn sparsemat_row_extraction() {
    let local = sample_4x6();
    let row = local.row(2).unwrap();
    assert!(row.len() == 6);
    assert!(row.nnz() == 3);
    assert!(row.indices() == [2, 3, 4]);
    assert!(row.values() == [50.0, 60.0, 70.0]);
    assert!(row.get(3) == Some(60.0));
    assert!(row.get(0).is_none());
    assert!(row.to_dense() == vec![0.0, 0.0, 50.0, 60.0, 70.0, 0.0]);
    // An empty row
    let empty = SparseMatrix::<f64>::empty_with_shape(2, 3).row(1).unwrap();
    assert!(empty.nnz() == 0 && empty.len() == 3);
}

#[test]
fn sparsemat_col_extraction() {
    let mut local = sample_4x6();
    let col = local.col(3).unwrap();
    assert!(col.len() == 4);
    assert!(col.indices() == [1, 2]);
    assert!(col.values() == [40.0, 60.0]);

    // Same answer from a fresh CSC cache
    local.explicitly_compress_csc();
    assert!(local.col(3).unwrap() == col);
}

#[test]
fn sparsemat_row_col_out_of_range() {
    let local = sample_4x6();
    assert!(local.row(4).err() == Some(SparseError::OutOfRange { index: 4, len: 4 }));
    assert!(local.col(6).err() == Some(SparseError::OutOfRange { index: 6, len: 6 }));
}

#[test]
fn sparsemat_col_iter() {
    let local = sample_4x6();
    let cols: Vec<Vec<f64>> = local.col_iter().collect();
    assert!(cols.len() == 6);
    assert!(cols[1] == vec![20.0, 30.0, 0.0, 0.0]);
    assert!(cols[5] == vec![0.0, 0.0, 0.0, 80.0]);
    // Columns of a matrix are the rows of its transpose
    let transposed_rows: Vec<Vec<f64>> = local.create_transpose().row_iter().collect();
    assert!(cols == transposed_rows);
}

#[test]
fn sparsemat_diagonal() {
    let local = sample_4x6();
    let diag = local.diagonal();
    assert!(diag.len() == 4);
    assert!(diag.indices() == [0, 1, 2]);
    assert!(diag.values() == [10.0, 30.0, 50.0]);
    assert!(diag.to_dense() == vec![10.0, 30.0, 50.0, 0.0]);

    let identity = SparseMatrix::<f64>::identity(3);
    assert!(identity.diagonal().to_dense() == vec![1.0, 1.0, 1.0]);
    // Tall matrix: diagonal length is the smaller dimension
    assert!(SparseMatrix::<f64>::empty_with_shape(5, 2).diagonal().len() == 2);
}

#[test]
fn sparsemat_submatrix() {
    let local = sample_4x6();
    let block = local.submatrix(1..3, 1..4).unwrap();
    assert!(block.shape == (2, 3));
    assert!(block.to_dense() == vec![vec![30.0, 0.0, 40.0], vec![0.0, 50.0, 60.0]]);

    let whole = local.submatrix(0..4, 0..6).unwrap();
    assert!(whole.to_dense() == local.to_dense());

    let empty = local.submatrix(2..2, 0..6).unwrap();
    assert!(empty.shape == (0, 6));
}

#[test]
fn sparsemat_submatrix_invalid_range() {
    let local = sample_4x6();
    assert!(
        local.submatrix(0..5, 0..6).err()
            == Some(SparseError::InvalidRange {
                start: 0,
                end: 5,
                len: 4
            })
    );
    #[allow(clippy::reversed_empty_ranges)]
    let reversed = local.submatrix(0..4, 4..2);
    assert!(
        reversed.err()
            == Some(SparseError::InvalidRange {
                start: 4,
                end: 2,
                len: 6
            })
    );
}

#[test]
fn sparsemat_index() {
    let local = sample_4x6();
    assert!(local[(0, 1)] == 20.0);
    assert!(local[(3, 5)] == 80.0);
    // Unstored entries read as zero
    assert!(local[(1, 0)] == 0.0);
    assert!(local[(3, 0)] == 0.0);
}

#[test]
fn sparsemat_index_integers() {
    let local = SparseMatrix::from_triplets(2, 2, vec![(0, 1, 7)]).unwrap();
    assert!(local[(0, 1)] == 7);
    assert!(local[(1, 1)] == 0);
}

#[test]
#[should_panic(expected = "out of bounds")]
fn sparsemat_index_oob() {
    let local = sample_4x6();
    let _ = local[(4, 0)];
}

#[test]
#[should_panic(expected = "out of bounds")]
fn sparsemat_index_oob_col() {
    let local = sample_4x6();
    let _ = local[(0, 6)];
}

#[test]
fn sparsemat_debug_format() {
    let local = SparseMatrix::from_triplets(2, 2, vec![(1, 1, 2.0), (0, 1, 1.0)]).unwrap();
    assert!(
        format!("{:?}", local)
            == "SparseMatrix { shape: (2, 2), entries: [(0, 1, 1.0), (1, 1, 2.0)] }"
    );
}

#[test]
fn sparsemat_partial_eq() {
    let local = SparseMatrix::from_triplets(2, 2, vec![(0, 0, 1.0), (1, 1, 2.0)]).unwrap();
    // Same entries inserted in a different order, with a different cache state
    let mut local2 = SparseMatrix::empty_with_shape(2, 2);
    local2.insert(1, 1, 2.0).unwrap();
    local2.insert(0, 0, 1.0).unwrap();
    local2.explicitly_compress();
    assert!(local == local2);

    local2.insert(0, 1, 5.0).unwrap();
    assert!(local != local2);

    // Same entries but different shape
    let wider = SparseMatrix::from_triplets(2, 3, vec![(0, 0, 1.0), (1, 1, 2.0)]).unwrap();
    assert!(local != wider);
}

#[test]
fn sparsemat_partial_eq_nan_not_reflexive() {
    let local = SparseMatrix::from_dense(vec![vec![f64::NAN]]).unwrap();
    assert!(local != local.clone());
}

#[test]
fn sparsemat_into_iterator() {
    let local = SparseMatrix::from_triplets(3, 3, vec![(2, 0, 3.0), (0, 1, 1.0)]).unwrap();

    let mut by_ref = vec![];
    for (row, col, val) in &local {
        by_ref.push((row, col, val));
    }
    assert!(by_ref == vec![(0, 1, 1.0), (2, 0, 3.0)]);

    let by_value: Vec<_> = local.into_iter().collect();
    assert!(by_value == by_ref);
}

#[test]
fn sparsemat_from_array() {
    let local = SparseMatrix::from([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]]);
    assert!(local.shape == (2, 3));
    assert!(local.num_nonzero() == 3);
    assert!(local[(1, 1)] == 3.0);

    // Zero-row arrays keep their column count
    let empty: SparseMatrix<f64> = SparseMatrix::from([[0.0; 4]; 0]);
    assert!(empty.shape == (0, 4));
}

#[test]
fn sparsemat_try_from_vecs() {
    let local = SparseMatrix::try_from(vec![vec![1, 0], vec![0, 2]]).unwrap();
    assert!(local.shape == (2, 2));

    let ragged = SparseMatrix::try_from(vec![vec![1, 0], vec![2]]);
    assert!(
        ragged.err()
            == Some(SparseError::RaggedRows {
                row: 1,
                expected: 2,
                found: 1
            })
    );
}

#[test]
fn sparsemat_into_dense_vecs() {
    let local = SparseMatrix::from([[1.0, 0.0], [0.0, 2.0]]);
    let by_ref: Vec<Vec<f64>> = (&local).into();
    assert!(by_ref == vec![vec![1.0, 0.0], vec![0.0, 2.0]]);
    let by_value: Vec<Vec<f64>> = local.into();
    assert!(by_value == by_ref);
}

#[test]
fn sparsemat_from_csr_basic() {
    // [[10, 20, 0], [0, 0, 0], [0, 30, 40]]
    let local = SparseMatrix::from_csr(
        3,
        3,
        vec![0, 2, 2, 4],
        vec![0, 1, 1, 2],
        vec![10.0, 20.0, 30.0, 40.0],
    )
    .unwrap();
    assert!(local.shape == (3, 3));
    assert!(local.num_nonzero() == 4);
    assert!(local[(2, 1)] == 30.0);
    assert!(local.to_dense() == vec![vec![10.0, 20.0, 0.0], vec![0.0; 3], vec![0.0, 30.0, 40.0]]);
    // The arrays were adopted as the CSR cache
    assert!(local.compressed().unwrap().rows == [0, 2, 2, 4]);
}

#[test]
fn sparsemat_from_csr_drops_zeros_and_rebuilds_cache() {
    let local =
        SparseMatrix::from_csr(2, 2, vec![0, 2, 3], vec![0, 1, 1], vec![1.0, 0.0, 2.0]).unwrap();
    assert!(local.num_nonzero() == 2);
    assert!(local.compressed().is_none());
    assert!(local.to_dense() == vec![vec![1.0, 0.0], vec![0.0, 2.0]]);
}

#[test]
fn sparsemat_from_csr_round_trip() {
    let mut original = sample_4x6();
    original.explicitly_compress();
    let view = original.compressed().unwrap();
    let rebuilt = SparseMatrix::from_csr(
        4,
        6,
        view.rows.to_vec(),
        view.cols.to_vec(),
        view.data.to_vec(),
    )
    .unwrap();
    assert!(rebuilt == original);
}

#[test]
fn sparsemat_from_csr_invalid() {
    let reason = |result: Result<SparseMatrix, SparseError>| match result.err() {
        Some(SparseError::InvalidCompressed { reason }) => reason,
        other => panic!("expected InvalidCompressed, got {other:?}"),
    };
    assert!(
        reason(SparseMatrix::from_csr(2, 2, vec![0, 1], vec![0], vec![1.0]))
            == "pointers must have rows + 1 entries"
    );
    assert!(
        reason(SparseMatrix::from_csr(1, 2, vec![0, 1], vec![0], vec![]))
            == "indices and values must have the same length"
    );
    assert!(
        reason(SparseMatrix::from_csr(1, 2, vec![0, 2], vec![0], vec![1.0]))
            == "pointers must start at 0 and end at the entry count"
    );
    assert!(
        reason(SparseMatrix::from_csr(
            3,
            2,
            vec![0, 2, 1, 2],
            vec![0, 1],
            vec![1.0, 2.0]
        )) == "pointers must be non-decreasing"
    );
    assert!(
        reason(SparseMatrix::from_csr(1, 2, vec![0, 1], vec![2], vec![1.0]))
            == "column index out of range"
    );
    assert!(
        reason(SparseMatrix::from_csr(
            1,
            3,
            vec![0, 2],
            vec![1, 1],
            vec![1.0, 2.0]
        )) == "column indices must be strictly increasing within a row"
    );
}

#[test]
fn sparsemat_from_csr_empty() {
    let local = SparseMatrix::<f64>::from_csr(0, 0, vec![0], vec![], vec![]).unwrap();
    assert!(local.shape == (0, 0));
    let wide = SparseMatrix::<f64>::from_csr(2, 5, vec![0, 0, 0], vec![], vec![]).unwrap();
    assert!(wide.shape == (2, 5) && wide.num_nonzero() == 0);
}

#[test]
fn sparsemat_from_diagonal() {
    let local = SparseMatrix::from_diagonal(vec![1.0, 0.0, 3.0]);
    assert!(local.shape == (3, 3));
    assert!(local.num_nonzero() == 2);
    assert!(local.diagonal().to_dense() == vec![1.0, 0.0, 3.0]);
    assert!(local[(0, 1)] == 0.0);

    let empty = SparseMatrix::<f64>::from_diagonal(Vec::new());
    assert!(empty.shape == (0, 0));

    let identity = SparseMatrix::from_diagonal([1.0; 4]);
    assert!(identity == SparseMatrix::identity(4));
}

fn sample_2x2() -> SparseMatrix {
    SparseMatrix::from([[1.0, 2.0], [0.0, 3.0]])
}

#[test]
fn sparsemat_subtraction() {
    let a = sample_2x2();
    let b = SparseMatrix::from([[1.0, 0.0], [4.0, 3.0]]);
    let difference = &a - &b;
    assert!(difference.to_dense() == vec![vec![0.0, 2.0], vec![-4.0, 0.0]]);
    // Cancelled entries are not stored
    assert!(difference.num_nonzero() == 2);
    // Operands untouched
    assert!(a == sample_2x2());
}

#[test]
fn sparsemat_try_sub_shape_mismatch() {
    let a = sample_2x2();
    let b = SparseMatrix::<f64>::empty_with_shape(3, 2);
    assert!(
        a.try_sub(&b).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 2),
                right: (3, 2)
            })
    );
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn sparsemat_sub_shape_mismatch_panics() {
    let _ = &sample_2x2() - &SparseMatrix::<f64>::empty_with_shape(3, 3);
}

#[test]
fn sparsemat_owned_operand_combinations() {
    let expected_sum = vec![vec![2.0, 4.0], vec![0.0, 6.0]];
    assert!((&sample_2x2() + &sample_2x2()).to_dense() == expected_sum);
    assert!((sample_2x2() + &sample_2x2()).to_dense() == expected_sum);
    assert!((&sample_2x2() + sample_2x2()).to_dense() == expected_sum);
    assert!((sample_2x2() + sample_2x2()).to_dense() == expected_sum);

    assert!((sample_2x2() - sample_2x2()).num_nonzero() == 0);
    assert!((&sample_2x2() - sample_2x2()).num_nonzero() == 0);
    assert!((sample_2x2() - &sample_2x2()).num_nonzero() == 0);
}

#[test]
fn sparsemat_owned_lhs_invalidates_cache() {
    let mut a = sample_2x2();
    a.explicitly_compress();
    let sum = a + SparseMatrix::from([[0.0, 0.0], [5.0, 0.0]]);
    assert!(sum.compressed().is_none());
    assert!(sum[(1, 0)] == 5.0);
}

#[test]
fn sparsemat_add_assign_sub_assign() {
    let mut a = sample_2x2();
    a += &SparseMatrix::from([[1.0, 0.0], [0.0, 0.0]]);
    assert!(a[(0, 0)] == 2.0);
    a += SparseMatrix::from([[0.0, -2.0], [0.0, 0.0]]);
    assert!(a[(0, 1)] == 0.0 && a.num_nonzero() == 2);
    a -= &SparseMatrix::from([[2.0, 0.0], [0.0, 3.0]]);
    assert!(a.num_nonzero() == 0);
    a -= SparseMatrix::from([[1.0, 0.0], [0.0, 0.0]]);
    assert!(a[(0, 0)] == -1.0);
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn sparsemat_add_assign_shape_mismatch_panics() {
    let mut a = sample_2x2();
    a += &SparseMatrix::<f64>::empty_with_shape(1, 1);
}

#[test]
fn sparsemat_negation() {
    let a = sample_2x2();
    let negated = -&a;
    assert!(negated.to_dense() == vec![vec![-1.0, -2.0], vec![0.0, -3.0]]);
    assert!(a == sample_2x2());
    assert!((-sample_2x2()).to_dense() == negated.to_dense());
    assert!((&a + &negated).num_nonzero() == 0);

    let ints = SparseMatrix::from([[1, 0], [0, -4]]);
    assert!((-ints).to_dense() == vec![vec![-1, 0], vec![0, 4]]);
}

#[test]
fn sparsemat_scalar_multiplication() {
    let a = sample_2x2();
    assert!((&a * 2.0).to_dense() == vec![vec![2.0, 4.0], vec![0.0, 6.0]]);
    assert!((a.clone() * 0.5).to_dense() == vec![vec![0.5, 1.0], vec![0.0, 1.5]]);
    assert!((2.0 * &a) == (&a * 2.0));
    assert!((2.0 * a.clone()) == (&a * 2.0));
    assert!(a.scale(3.0)[(1, 1)] == 9.0);

    // Multiplying by zero leaves nothing stored
    let zeroed = &a * 0.0;
    assert!(zeroed.num_nonzero() == 0 && zeroed.shape == (2, 2));

    let mut in_place = sample_2x2();
    in_place *= 10.0;
    assert!(in_place[(0, 1)] == 20.0);
    assert!(a == sample_2x2());
}

#[test]
fn sparsemat_scalar_multiplication_integers() {
    let a = SparseMatrix::from([[1i32, 2], [0, 3]]);
    assert!((&a * 4).to_dense() == vec![vec![4, 8], vec![0, 12]]);
    assert!((3i32 * &a).to_dense() == vec![vec![3, 6], vec![0, 9]]);
}

#[test]
fn sparsemat_hadamard() {
    let a = SparseMatrix::from([[1.0, 2.0, 0.0], [0.0, 3.0, 4.0]]);
    let b = SparseMatrix::from([[5.0, 0.0, 6.0], [0.0, 7.0, 8.0]]);
    let product = a.hadamard(&b).unwrap();
    assert!(product.to_dense() == vec![vec![5.0, 0.0, 0.0], vec![0.0, 21.0, 32.0]]);
    assert!(product.num_nonzero() == 3);
    // Commutative, regardless of which operand is sparser
    assert!(b.hadamard(&a).unwrap() == product);
}

#[test]
fn sparsemat_hadamard_disjoint_and_shape_mismatch() {
    let a = SparseMatrix::from([[1.0, 0.0], [0.0, 0.0]]);
    let b = SparseMatrix::from([[0.0, 2.0], [0.0, 0.0]]);
    assert!(a.hadamard(&b).unwrap().num_nonzero() == 0);

    let other = SparseMatrix::<f64>::empty_with_shape(2, 3);
    assert!(
        a.hadamard(&other).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 2),
                right: (2, 3)
            })
    );
}

#[test]
fn sparsemat_map_changes_type_and_drops_zeros() {
    let a = SparseMatrix::from([[1.0, 0.0], [-2.0, 3.0]]);
    let signs = a.map(|x| if x > 0.0 { 1i32 } else { 0 });
    assert!(signs.to_dense() == vec![vec![1, 0], vec![0, 1]]);
    assert!(signs.num_nonzero() == 2);

    // f is only applied to stored entries, so zeros stay zero
    let shifted = a.map(|x| x + 1.0);
    assert!(shifted.to_dense() == vec![vec![2.0, 0.0], vec![-1.0, 4.0]]);
}

#[test]
fn sparsemat_apply_in_place() {
    let mut a = SparseMatrix::from([[1.0, 0.0], [-1.0, 3.0]]);
    a.explicitly_compress();
    a.apply(|x| x * x - 1.0);
    // 1 -> 0 (dropped), -1 -> 0 (dropped), 3 -> 8
    assert!(a.to_dense() == vec![vec![0.0, 0.0], vec![0.0, 8.0]]);
    assert!(a.num_nonzero() == 1);
    assert!(a.compressed().is_none());
}

/// Deterministic pseudo-random sparse matrix with small integer-valued entries, so that
/// products compare exactly against the dense reference.
fn pseudo_random(rows: usize, cols: usize, seed: u64, fill_percent: u64) -> SparseMatrix {
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
                triplets.push((row, col, (next() % 9) as f64 - 4.0));
            }
        }
    }
    SparseMatrix::from_triplets(rows, cols, triplets).unwrap()
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
fn sparsemat_mul_vec() {
    let a = SparseMatrix::from([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]]);
    assert!(a.mul_vec(&[1.0, 10.0, 100.0]).unwrap() == vec![201.0, 30.0]);
    assert!((&a * &vec![1.0, 10.0, 100.0]) == vec![201.0, 30.0]);
    assert!((&a * &[1.0, 1.0, 1.0][..]) == vec![3.0, 3.0]);
}

#[test]
fn sparsemat_mul_vec_matches_dense() {
    let a = pseudo_random(7, 5, 11, 40);
    let x = vec![1.0, -2.0, 3.0, 0.0, 5.0];
    let column: Vec<Vec<f64>> = x.iter().map(|&v| vec![v]).collect();
    let expected: Vec<f64> = dense_matmul(&a.to_dense(), &column)
        .into_iter()
        .map(|row| row[0])
        .collect();
    assert!(a.mul_vec(&x).unwrap() == expected);
}

#[test]
fn sparsemat_mul_vec_shape_mismatch() {
    let a = SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    assert!(
        a.mul_vec(&[1.0, 2.0, 3.0]).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 2),
                right: (3, 1)
            })
    );
}

#[test]
fn sparsemat_vec_mul_is_transpose_product() {
    let a = pseudo_random(6, 4, 3, 50);
    let x = vec![2.0, -1.0, 0.0, 4.0, 1.0, 3.0];
    let via_transpose = a.create_transpose().mul_vec(&x).unwrap();
    assert!(a.vec_mul(&x).unwrap() == via_transpose);
    assert!((&x * &a) == via_transpose);
    assert!((&x[..] * &a) == via_transpose);

    // Same answer with a pre-built CSC cache
    let mut cached = a.clone();
    cached.explicitly_compress_csc();
    assert!(cached.vec_mul(&x).unwrap() == via_transpose);
}

#[test]
fn sparsemat_vec_mul_shape_mismatch() {
    let a = SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    assert!(
        a.vec_mul(&[1.0]).err()
            == Some(SparseError::ShapeMismatch {
                left: (1, 1),
                right: (2, 2)
            })
    );
}

#[test]
fn sparsemat_mul_dense() {
    let a = pseudo_random(5, 4, 21, 45);
    let b = pseudo_random(4, 3, 22, 100).to_dense();
    assert!(a.mul_dense(&b).unwrap() == dense_matmul(&a.to_dense(), &b));
}

#[test]
fn sparsemat_mul_dense_errors() {
    let a = SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    assert!(
        a.mul_dense(&[vec![1.0, 2.0]]).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 2),
                right: (1, 2)
            })
    );
    assert!(
        a.mul_dense(&[vec![1.0, 2.0], vec![3.0]]).err()
            == Some(SparseError::RaggedRows {
                row: 1,
                expected: 2,
                found: 1
            })
    );
}

#[test]
fn sparsemat_matmul_matches_dense() {
    for (seed, dims) in [
        (1u64, (6, 5, 4)),
        (2, (3, 8, 3)),
        (3, (10, 10, 10)),
        (4, (1, 7, 1)),
    ] {
        let (m, k, n) = dims;
        let a = pseudo_random(m, k, seed, 35);
        let b = pseudo_random(k, n, seed + 100, 35);
        let product = a.matmul(&b).unwrap();
        assert!(product.shape == (m, n));
        assert!(product.to_dense() == dense_matmul(&a.to_dense(), &b.to_dense()));
        // The operator agrees, for every operand form
        assert!(&a * &b == product);
        assert!(a.clone() * b.clone() == product);
        assert!(a.clone() * &b == product);
        assert!(&a * b.clone() == product);
    }
}

#[test]
fn sparsemat_matmul_identity_and_zero() {
    let a = pseudo_random(6, 6, 9, 40);
    let identity = SparseMatrix::identity(6);
    assert!(a.matmul(&identity).unwrap() == a);
    assert!(identity.matmul(&a).unwrap() == a);

    let zero = SparseMatrix::<f64>::empty_with_shape(6, 6);
    assert!(a.matmul(&zero).unwrap().num_nonzero() == 0);
}

#[test]
fn sparsemat_matmul_associative_and_transpose_rule() {
    let a = pseudo_random(4, 5, 31, 40);
    let b = pseudo_random(5, 3, 32, 40);
    let c = pseudo_random(3, 6, 33, 40);
    assert!((&(&a * &b) * &c) == (&a * &(&b * &c)));
    // (AB)ᵀ = BᵀAᵀ
    assert!((&a * &b).create_transpose() == &b.create_transpose() * &a.create_transpose());
}

#[test]
fn sparsemat_matmul_drops_cancellations() {
    // Row 0 of the product is 1*1 + (-1)*1 = 0
    let a = SparseMatrix::from([[1, -1], [0, 2]]);
    let b = SparseMatrix::from([[1, 0], [1, 5]]);
    let product = a.matmul(&b).unwrap();
    assert!(product.to_dense() == vec![vec![0, -5], vec![2, 10]]);
    assert!(product.num_nonzero() == 3);
}

#[test]
fn sparsemat_matmul_shape_mismatch() {
    let a = SparseMatrix::<f64>::empty_with_shape(2, 3);
    let b = SparseMatrix::<f64>::empty_with_shape(2, 3);
    assert!(
        a.matmul(&b).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 3),
                right: (2, 3)
            })
    );
}

#[test]
#[should_panic(expected = "incompatible shapes")]
fn sparsemat_mul_operator_shape_mismatch_panics() {
    let a = SparseMatrix::<f64>::empty_with_shape(2, 3);
    let _ = &a * &a;
}

#[test]
fn sparsemat_matmul_huge_sparse_shape() {
    let mut a = SparseMatrix::<f64>::empty_with_shape(1 << 20, 1 << 20);
    a.insert(3, 1 << 19, 2.0).unwrap();
    let mut b = SparseMatrix::<f64>::empty_with_shape(1 << 20, 1 << 20);
    b.insert(1 << 19, 7, 5.0).unwrap();
    let product = a.matmul(&b).unwrap();
    assert!(product.num_nonzero() == 1 && product[(3, 7)] == 10.0);
}

#[test]
fn sparsemat_scalar_mul_still_works_alongside_matmul() {
    let a = SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    assert!((&a * 2.0).to_dense() == vec![vec![2.0, 4.0], vec![6.0, 8.0]]);
    assert!((&a * &a).to_dense() == vec![vec![7.0, 10.0], vec![15.0, 22.0]]);
}

#[test]
fn sparsemat_is_square_and_density() {
    assert!(SparseMatrix::<f64>::identity(3).is_square());
    assert!(!SparseMatrix::<f64>::empty_with_shape(2, 3).is_square());
    assert!(SparseMatrix::<f64>::new().is_square());

    assert!(SparseMatrix::<f64>::identity(4).density() == 0.25);
    assert!(SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]).density() == 1.0);
    assert!(SparseMatrix::<f64>::empty_with_shape(3, 3).density() == 0.0);
    // No division by zero for an empty shape
    assert!(SparseMatrix::<f64>::new().density() == 0.0);
    assert!(SparseMatrix::<f64>::empty_with_shape(0, 5).density() == 0.0);
}

#[test]
fn sparsemat_sum_and_trace() {
    let a = SparseMatrix::from([[1.0, 2.0, 0.0], [0.0, 3.0, 4.0], [5.0, 0.0, 6.0]]);
    assert!(a.sum() == 21.0);
    assert!(a.trace() == Ok(10.0));
    assert!(SparseMatrix::<f64>::empty_with_shape(2, 2).sum() == 0.0);

    let ints = SparseMatrix::from([[1, 2], [3, 4]]);
    assert!(ints.sum() == 10 && ints.trace() == Ok(5));
}

#[test]
fn sparsemat_trace_not_square() {
    let a = SparseMatrix::<f64>::empty_with_shape(2, 3);
    assert!(a.trace() == Err(SparseError::NotSquare { shape: (2, 3) }));
}

#[test]
fn sparsemat_is_symmetric() {
    assert!(SparseMatrix::from([[1.0, 2.0], [2.0, 3.0]]).is_symmetric());
    assert!(!SparseMatrix::from([[1.0, 2.0], [5.0, 3.0]]).is_symmetric());
    // One side stored, the other not
    assert!(!SparseMatrix::from([[0.0, 2.0], [0.0, 0.0]]).is_symmetric());
    assert!(!SparseMatrix::<f64>::empty_with_shape(2, 3).is_symmetric());
    assert!(SparseMatrix::<f64>::identity(5).is_symmetric());
    assert!(SparseMatrix::<f64>::new().is_symmetric());

    // A + Aᵀ is always symmetric
    let a = pseudo_random(6, 6, 77, 40);
    assert!((&a + &a.create_transpose()).is_symmetric());
}

#[test]
fn sparsemat_norms() {
    let a = SparseMatrix::from([[1.0, -2.0, 0.0], [0.0, 3.0, 4.0]]);
    assert!((a.norm_frobenius() - 30.0f64.sqrt()).abs() < 1e-12);
    // Column sums: 1, 5, 4 and row sums: 3, 7
    assert!(a.norm_1() == 5.0);
    assert!(a.norm_inf() == 7.0);

    let empty = SparseMatrix::<f64>::empty_with_shape(3, 3);
    assert!(empty.norm_frobenius() == 0.0 && empty.norm_1() == 0.0 && empty.norm_inf() == 0.0);

    // f32 works too, and the norms of a transpose swap the 1- and infinity-norms
    let b = SparseMatrix::from([[1.0f32, 2.0], [3.0, 4.0]]);
    assert!(b.norm_1() == b.create_transpose().norm_inf());
    assert!(b.norm_inf() == b.create_transpose().norm_1());
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
fn sparsemat_kron_matches_dense() {
    let a = pseudo_random(3, 2, 5, 60);
    let b = pseudo_random(2, 4, 6, 60);
    let product = a.kron(&b).unwrap();
    assert!(product.shape == (6, 8));
    assert!(product.to_dense() == dense_kron(&a.to_dense(), &b.to_dense()));
    assert!(product.num_nonzero() == a.num_nonzero() * b.num_nonzero());
}

#[test]
fn sparsemat_kron_identity_is_block_diagonal() {
    let b = SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    let block_diag = SparseMatrix::identity(2).kron(&b).unwrap();
    assert!(
        block_diag.to_dense()
            == vec![
                vec![1.0, 2.0, 0.0, 0.0],
                vec![3.0, 4.0, 0.0, 0.0],
                vec![0.0, 0.0, 1.0, 2.0],
                vec![0.0, 0.0, 3.0, 4.0],
            ]
    );
}

#[test]
fn sparsemat_kron_mixed_product_property() {
    // (A ⊗ B)(C ⊗ D) = (AC) ⊗ (BD)
    let (a, b) = (pseudo_random(2, 3, 41, 70), pseudo_random(3, 2, 42, 70));
    let (c, d) = (pseudo_random(3, 2, 43, 70), pseudo_random(2, 3, 44, 70));
    let left = &a.kron(&b).unwrap() * &c.kron(&d).unwrap();
    let right = (&a * &c).kron(&(&b * &d)).unwrap();
    assert!(left == right);
}

#[test]
fn sparsemat_kron_empty_operand_and_overflow() {
    let a = SparseMatrix::from([[1.0, 2.0]]);
    let empty = SparseMatrix::<f64>::new();
    let product = a.kron(&empty).unwrap();
    assert!(product.shape == (0, 0));

    let huge = SparseMatrix::<f64>::empty_with_shape(usize::MAX, 1);
    assert!(huge.kron(&huge) == Err(SparseError::ShapeOverflow));
}

#[test]
fn sparsemat_hstack_vstack() {
    let a = SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    let b = SparseMatrix::from([[5.0], [6.0]]);
    let wide = SparseMatrix::hstack([&a, &b]).unwrap();
    assert!(wide.to_dense() == vec![vec![1.0, 2.0, 5.0], vec![3.0, 4.0, 6.0]]);

    let c = SparseMatrix::from([[7.0, 8.0]]);
    let tall = SparseMatrix::vstack([&a, &c]).unwrap();
    assert!(tall.to_dense() == vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![7.0, 8.0]]);

    // Owned blocks work too
    let owned = SparseMatrix::hstack(vec![a.clone(), b.clone()]).unwrap();
    assert!(owned == wide);
}

#[test]
fn sparsemat_stack_reassembles_submatrices() {
    let original = pseudo_random(5, 6, 55, 45);
    let left = original.submatrix(0..5, 0..2).unwrap();
    let right = original.submatrix(0..5, 2..6).unwrap();
    assert!(SparseMatrix::hstack([&left, &right]).unwrap() == original);

    let top = original.submatrix(0..3, 0..6).unwrap();
    let bottom = original.submatrix(3..5, 0..6).unwrap();
    assert!(SparseMatrix::vstack([&top, &bottom]).unwrap() == original);
}

#[test]
fn sparsemat_stack_errors_and_empty() {
    let a = SparseMatrix::<f64>::empty_with_shape(2, 2);
    let b = SparseMatrix::<f64>::empty_with_shape(3, 1);
    assert!(
        SparseMatrix::hstack([&a, &b]).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 2),
                right: (3, 1)
            })
    );
    assert!(
        SparseMatrix::vstack([&a, &b]).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 2),
                right: (3, 1)
            })
    );

    let none: [&SparseMatrix; 0] = [];
    assert!(SparseMatrix::hstack(none).unwrap().shape == (0, 0));
    assert!(SparseMatrix::vstack(none).unwrap().shape == (0, 0));

    // A single block is returned as is
    assert!(SparseMatrix::vstack([&a]).unwrap() == a);
    // Blocks with zero columns still contribute their rows
    let tall_empty = SparseMatrix::<f64>::empty_with_shape(2, 0);
    assert!(
        SparseMatrix::vstack([&tall_empty, &tall_empty])
            .unwrap()
            .shape
            == (4, 0)
    );
}

fn assert_vec_close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert!(actual.len() == expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert!(
            (a - e).abs() <= tolerance,
            "{a} differs from {e}: {actual:?} vs {expected:?}"
        );
    }
}

/// A strictly diagonally dominant (so non-singular) pseudo-random matrix.
fn dominant(n: usize, seed: u64) -> SparseMatrix {
    let shift = SparseMatrix::from_diagonal(vec![4.0 * n as f64 + 1.0; n]);
    &pseudo_random(n, n, seed, 30) + &shift
}

/// A symmetric positive definite matrix `BᵀB + n I`.
fn spd(n: usize, seed: u64) -> SparseMatrix {
    let b = pseudo_random(n, n, seed, 25);
    let gram = &b.create_transpose() * &b;
    &gram + &SparseMatrix::from_diagonal(vec![n as f64; n])
}

#[test]
fn sparsemat_triangular_solves() {
    let lower = SparseMatrix::from([[2.0, 0.0, 0.0], [1.0, 4.0, 0.0], [0.0, 3.0, 5.0]]);
    // x = [1, 2, 3] -> b = L x
    let b = lower.mul_vec(&[1.0, 2.0, 3.0]).unwrap();
    assert_vec_close(
        &lower.solve_lower_triangular(&b).unwrap(),
        &[1.0, 2.0, 3.0],
        1e-12,
    );

    let upper = lower.create_transpose();
    let b = upper.mul_vec(&[1.0, 2.0, 3.0]).unwrap();
    assert_vec_close(
        &upper.solve_upper_triangular(&b).unwrap(),
        &[1.0, 2.0, 3.0],
        1e-12,
    );
}

#[test]
fn sparsemat_triangular_solves_ignore_other_triangle() {
    // Garbage above the diagonal must not affect a lower solve
    let a = SparseMatrix::from([[2.0, 99.0], [1.0, 4.0]]);
    assert_vec_close(
        &a.solve_lower_triangular(&[2.0, 9.0]).unwrap(),
        &[1.0, 2.0],
        1e-12,
    );
    assert_vec_close(
        &a.solve_upper_triangular(&[200.0, 8.0]).unwrap(),
        &[1.0, 2.0],
        1e-12,
    );
}

#[test]
fn sparsemat_triangular_solve_errors() {
    let missing_diagonal = SparseMatrix::from([[1.0, 0.0], [1.0, 0.0]]);
    assert!(
        missing_diagonal.solve_lower_triangular(&[1.0, 1.0]).err()
            == Some(SparseError::SingularMatrix { pivot: 1 })
    );
    let a = SparseMatrix::<f64>::identity(2);
    assert!(
        a.solve_upper_triangular(&[1.0]).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 2),
                right: (1, 1)
            })
    );
    let rect = SparseMatrix::<f64>::empty_with_shape(2, 3);
    assert!(
        rect.solve_lower_triangular(&[1.0, 1.0]).err()
            == Some(SparseError::NotSquare { shape: (2, 3) })
    );
}

#[test]
fn sparsemat_lu_reconstructs_permuted_matrix() {
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
fn sparsemat_lu_solves_systems() {
    for (n, seed) in [(3, 7), (10, 8), (25, 9)] {
        let a = dominant(n, seed);
        let x: Vec<f64> = (0..n).map(|i| i as f64 - 3.0).collect();
        let b = a.mul_vec(&x).unwrap();
        assert_vec_close(&a.solve(&b).unwrap(), &x, 1e-9);
        assert_vec_close(&a.lu().unwrap().solve(&b).unwrap(), &x, 1e-9);
    }
}

#[test]
fn sparsemat_lu_needs_pivoting() {
    // A zero on the diagonal forces a row swap
    let a = SparseMatrix::from([[0.0, 2.0], [3.0, 1.0]]);
    let lu = a.lu().unwrap();
    assert!(lu.permutation() == [1, 0]);
    assert_vec_close(&lu.solve(&[4.0, 5.0]).unwrap(), &[1.0, 2.0], 1e-12);

    // A matrix whose rows are all shifted so that nothing is on the diagonal
    let n = 6;
    let shifted =
        SparseMatrix::from_triplets(n, n, (0..n).map(|i| (i, (i + 1) % n, (i + 2) as f64)))
            .unwrap();
    let x: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
    let b = shifted.mul_vec(&x).unwrap();
    assert_vec_close(&shifted.solve(&b).unwrap(), &x, 1e-12);
}

#[test]
fn sparsemat_lu_determinant() {
    let a = SparseMatrix::from([[4.0f64, 3.0], [6.0, 3.0]]);
    assert!((a.lu().unwrap().determinant() - (-6.0)).abs() < 1e-12);
    // Row swap flips the sign: det [[0,2],[3,1]] = -6
    let b = SparseMatrix::from([[0.0f64, 2.0], [3.0, 1.0]]);
    assert!((b.lu().unwrap().determinant() - (-6.0)).abs() < 1e-12);
    assert!(SparseMatrix::<f64>::identity(5).lu().unwrap().determinant() == 1.0);
    assert!(SparseMatrix::<f64>::new().lu().unwrap().determinant() == 1.0);
}

#[test]
fn sparsemat_lu_singular_and_shape_errors() {
    let singular = SparseMatrix::from([[1.0, 2.0], [2.0, 4.0]]);
    assert!(singular.lu().err() == Some(SparseError::SingularMatrix { pivot: 1 }));
    let empty_column = SparseMatrix::from([[1.0, 0.0], [1.0, 0.0]]);
    assert!(empty_column.lu().err() == Some(SparseError::SingularMatrix { pivot: 1 }));

    let rect = SparseMatrix::<f64>::empty_with_shape(2, 3);
    assert!(rect.lu().err() == Some(SparseError::NotSquare { shape: (2, 3) }));

    let a = SparseMatrix::<f64>::identity(3);
    assert!(a.lu().unwrap().solve(&[1.0]).is_err());
    assert!(a.solve(&[1.0]).is_err());
}

#[test]
fn sparsemat_cholesky_reconstructs() {
    for (n, seed) in [(1, 1), (5, 2), (12, 3), (20, 4)] {
        let a = spd(n, seed);
        let l = a.cholesky().unwrap();
        // Lower triangular with a positive diagonal
        assert!(l.iter().all(|(row, col, _)| col <= row));
        assert!(l.diagonal().values().iter().all(|&d| d > 0.0));
        let reconstructed = (&l * &l.create_transpose()).to_dense();
        for (got, want) in reconstructed.iter().zip(a.to_dense()) {
            assert_vec_close(got, &want, 1e-8);
        }
    }
}

#[test]
fn sparsemat_cholesky_keeps_diagonal_matrix_sparse() {
    let a = SparseMatrix::from_diagonal(vec![4.0, 9.0, 16.0]);
    let l = a.cholesky().unwrap();
    assert!(l.num_nonzero() == 3);
    assert!(l.diagonal().to_dense() == vec![2.0, 3.0, 4.0]);
}

#[test]
fn sparsemat_cholesky_errors() {
    let not_symmetric = SparseMatrix::from([[2.0, 1.0], [0.0, 2.0]]);
    assert!(not_symmetric.cholesky().err() == Some(SparseError::NotSymmetric));
    let rect = SparseMatrix::<f64>::empty_with_shape(2, 3);
    assert!(rect.cholesky().err() == Some(SparseError::NotSquare { shape: (2, 3) }));
    // Symmetric but indefinite: eigenvalues 3 and -1
    let indefinite = SparseMatrix::from([[1.0, 2.0], [2.0, 1.0]]);
    assert!(indefinite.cholesky().err() == Some(SparseError::NotPositiveDefinite { pivot: 1 }));
    // Zero matrix has a non-positive first pivot
    let zero = SparseMatrix::<f64>::empty_with_shape(2, 2);
    assert!(zero.cholesky().err() == Some(SparseError::NotPositiveDefinite { pivot: 0 }));
}

#[test]
fn sparsemat_solve_spd_matches_lu() {
    let a = spd(14, 5);
    let x: Vec<f64> = (0..14).map(|i| (i as f64).sin()).collect();
    let b = a.mul_vec(&x).unwrap();
    assert_vec_close(&a.solve_spd(&b).unwrap(), &x, 1e-9);
    assert_vec_close(&a.solve(&b).unwrap(), &x, 1e-9);
}

#[test]
fn sparsemat_conjugate_gradient_solves() {
    for (n, seed) in [(5, 1), (20, 2), (40, 3)] {
        let a = spd(n, seed);
        let x: Vec<f64> = (0..n).map(|i| 1.0 + (i % 3) as f64).collect();
        let b = a.mul_vec(&x).unwrap();
        let solution = a.conjugate_gradient(&b, 1e-12, 10 * n).unwrap();
        assert_vec_close(&solution, &x, 1e-7);
    }
}

#[test]
fn sparsemat_conjugate_gradient_edge_cases() {
    // Identity converges in one step; a zero right-hand side needs none
    let identity = SparseMatrix::<f64>::identity(4);
    let b = [1.0, 2.0, 3.0, 4.0];
    assert!(identity.conjugate_gradient(&b, 1e-12, 1).unwrap() == b);
    assert!(identity.conjugate_gradient(&[0.0; 4], 1e-12, 0).unwrap() == vec![0.0; 4]);

    // Too few iterations
    let a = spd(30, 6);
    let rhs = vec![1.0; 30];
    assert!(
        a.conjugate_gradient(&rhs, 1e-14, 1).err()
            == Some(SparseError::NoConvergence { iterations: 1 })
    );
    // Negative definite is detected on the first step
    let negative = -SparseMatrix::<f64>::identity(3);
    assert!(
        negative.conjugate_gradient(&[1.0; 3], 1e-12, 10).err()
            == Some(SparseError::NoConvergence { iterations: 0 })
    );
    assert!(identity.conjugate_gradient(&[1.0], 1e-12, 10).is_err());
}

const SAMPLE_MTX: &str = "%%MatrixMarket matrix coordinate real general
% A 3x4 example
3 4 4

1 1 1.5
% a comment between entries
2 3 -2
3 4 1e2
3 1 0.25
";

#[test]
fn sparsemat_matrix_market_read() {
    let local = SparseMatrix::<f64>::from_matrix_market(SAMPLE_MTX).unwrap();
    assert!(local.shape == (3, 4));
    assert!(local.num_nonzero() == 4);
    assert!(local[(0, 0)] == 1.5);
    assert!(local[(1, 2)] == -2.0);
    assert!(local[(2, 3)] == 100.0);
    assert!(local[(2, 0)] == 0.25);
}

#[test]
fn sparsemat_matrix_market_round_trip() {
    for seed in [1u64, 2, 3] {
        let original = pseudo_random(7, 5, seed, 40);
        let text = original.to_matrix_market();
        assert!(text.starts_with("%%MatrixMarket matrix coordinate real general\n"));
        assert!(SparseMatrix::<f64>::from_matrix_market(&text).unwrap() == original);
    }
}

#[test]
fn sparsemat_matrix_market_round_trip_exact_floats() {
    let original = SparseMatrix::from_triplets(
        2,
        2,
        vec![
            (0, 0, 0.1 + 0.2),
            (1, 1, 1.0 / 3.0),
            (0, 1, 1e-300),
            (1, 0, -6.02e23),
        ],
    )
    .unwrap();
    let parsed = SparseMatrix::<f64>::from_matrix_market(&original.to_matrix_market()).unwrap();
    assert!(parsed == original);
}

#[test]
fn sparsemat_matrix_market_write_format() {
    let local = SparseMatrix::from([[0.0, 2.5], [3.0, 0.0]]);
    assert!(
        local.to_matrix_market()
            == "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 2 2.5\n2 1 3\n"
    );
    let mut bytes = vec![];
    local.write_matrix_market(&mut bytes).unwrap();
    assert!(String::from_utf8(bytes).unwrap() == local.to_matrix_market());
}

#[test]
fn sparsemat_matrix_market_integer_and_pattern() {
    let ints = "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 7\n2 2 -3\n";
    let local = SparseMatrix::<i32>::from_matrix_market(ints).unwrap();
    assert!(local.to_dense() == vec![vec![7, 0], vec![0, -3]]);
    assert!(local.to_matrix_market().contains("2 2 -3"));

    let pattern = "%%MatrixMarket matrix coordinate pattern general\n2 3 2\n1 3\n2 1\n";
    let local = SparseMatrix::<f64>::from_matrix_market(pattern).unwrap();
    assert!(local.to_dense() == vec![vec![0.0, 0.0, 1.0], vec![1.0, 0.0, 0.0]]);
}

#[test]
fn sparsemat_matrix_market_symmetric_and_skew() {
    let symmetric = "%%MatrixMarket matrix coordinate real symmetric\n3 3 3\n1 1 4\n2 1 1\n3 2 2\n";
    let local = SparseMatrix::<f64>::from_matrix_market(symmetric).unwrap();
    assert!(local.is_symmetric());
    assert!(local[(0, 1)] == 1.0 && local[(1, 2)] == 2.0 && local.num_nonzero() == 5);

    let skew = "%%MatrixMarket matrix coordinate real skew-symmetric\n2 2 1\n2 1 5\n";
    let local = SparseMatrix::<f64>::from_matrix_market(skew).unwrap();
    assert!(local[(1, 0)] == 5.0 && local[(0, 1)] == -5.0);
}

#[test]
fn sparsemat_matrix_market_sums_duplicates_and_drops_zeros() {
    let text =
        "%%MatrixMarket matrix coordinate real general\n2 2 4\n1 1 1\n1 1 2\n2 2 5\n2 2 -5\n";
    let local = SparseMatrix::<f64>::from_matrix_market(text).unwrap();
    assert!(local[(0, 0)] == 3.0 && local.num_nonzero() == 1);
}

#[test]
fn sparsemat_matrix_market_header_is_case_insensitive_and_empty_matrices() {
    let text = "%%matrixmarket MATRIX Coordinate Real General\n0 0 0\n";
    let local = SparseMatrix::<f64>::from_matrix_market(text).unwrap();
    assert!(local.shape == (0, 0));

    let empty = SparseMatrix::<f64>::empty_with_shape(3, 2);
    let parsed = SparseMatrix::<f64>::from_matrix_market(&empty.to_matrix_market()).unwrap();
    assert!(parsed.shape == (3, 2) && parsed.num_nonzero() == 0);
}

#[test]
fn sparsemat_matrix_market_read_from_reader() {
    let bytes = SAMPLE_MTX.as_bytes();
    let local = SparseMatrix::<f64>::read_matrix_market(std::io::BufReader::new(bytes)).unwrap();
    assert!(local.num_nonzero() == 4);
}

fn parse_failure(text: &str) -> (usize, String) {
    match SparseMatrix::<f64>::from_matrix_market(text).err() {
        Some(SparseError::Parse { line, message }) => (line, message),
        other => panic!("expected a parse error, got {other:?}"),
    }
}

#[test]
fn sparsemat_matrix_market_errors() {
    assert!(parse_failure("").1 == "empty input");
    assert!(parse_failure("not a header\n1 1 0\n").0 == 1);
    assert!(
        parse_failure("%%MatrixMarket matrix array real general\n1 1\n1\n")
            .1
            .contains("unsupported format")
    );
    assert!(
        parse_failure("%%MatrixMarket matrix coordinate complex general\n1 1 0\n")
            .1
            .contains("unsupported field")
    );
    assert!(
        parse_failure("%%MatrixMarket matrix coordinate real hermitian\n1 1 0\n")
            .1
            .contains("unsupported symmetry")
    );
    assert!(
        parse_failure("%%MatrixMarket matrix coordinate real general\n")
            .1
            .contains("missing size line")
    );
    assert!(
        parse_failure("%%MatrixMarket matrix coordinate real general\n2 2\n")
            .1
            .contains("size line")
    );
    assert!(
        parse_failure("%%MatrixMarket matrix coordinate real general\nx 2 1\n")
            .1
            .contains("invalid size")
    );
}

#[test]
fn sparsemat_matrix_market_entry_errors_name_the_line() {
    let header = "%%MatrixMarket matrix coordinate real general\n2 2 1\n";
    // Row index 0 (1-based files start at 1), past the end, and unparsable
    assert!(
        parse_failure(&format!("{header}0 1 1\n"))
            == (3, "row index '0' out of range 1..=2".to_string())
    );
    assert!(parse_failure(&format!("{header}1 3 1\n"))
        .1
        .contains("column index '3'"));
    assert!(parse_failure(&format!("{header}1 1 abc\n"))
        .1
        .contains("invalid value 'abc'"));
    assert!(parse_failure(&format!("{header}1 1\n"))
        .1
        .contains("expected 3 fields"));
    // Count mismatches
    assert!(parse_failure(&format!("{header}\n"))
        .1
        .contains("expected 1 entries, found 0"));
    assert!(parse_failure(&format!("{header}1 1 1\n2 2 2\n"))
        .1
        .contains("more than 1 entries"));
}

#[test]
fn sparsemat_matrix_market_io_error_is_reported() {
    struct FailingReader;
    impl std::io::Read for FailingReader {
        fn read(&mut self, _: &mut [u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("disk on fire"))
        }
    }
    let result = SparseMatrix::<f64>::read_matrix_market(std::io::BufReader::new(FailingReader));
    assert!(matches!(result, Err(SparseError::Io { message }) if message.contains("disk on fire")));

    struct FailingWriter;
    impl std::io::Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("no space"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let result = SparseMatrix::<f64>::identity(2).write_matrix_market(FailingWriter);
    assert!(matches!(result, Err(SparseError::Io { message }) if message.contains("no space")));
}

#[test]
fn sparsemat_from_csr_rows_overflow_is_an_error() {
    // rows + 1 does not fit in a usize, so no row_ptr can have the right length
    let result = SparseMatrix::<f64>::from_csr(usize::MAX, 1, vec![0], vec![], vec![]);
    assert!(
        result.err()
            == Some(SparseError::InvalidCompressed {
                reason: "pointers must have rows + 1 entries"
            })
    );
}

#[test]
#[should_panic(expected = "compressed storage")]
fn sparsemat_csr_for_usize_max_rows_panics_clearly() {
    // CSR needs rows + 1 pointers, which overflows
    let local = SparseMatrix::<f64>::empty_with_shape(usize::MAX, 1);
    let _ = local.iter().count();
}

#[test]
#[should_panic(expected = "compressed storage")]
fn sparsemat_csc_for_usize_max_cols_panics_clearly() {
    let local = SparseMatrix::<f64>::empty_with_shape(1, usize::MAX);
    let _ = local.col(0);
}

#[test]
#[should_panic(expected = "compressed storage")]
fn sparsemat_csr_too_large_to_allocate_panics_clearly() {
    // No overflow in rows + 1, but 2^61 pointers is more bytes than an allocation can hold
    let local = SparseMatrix::<f64>::empty_with_shape(1 << 61, 1);
    let _ = local.iter().count();
}

#[test]
fn sparsemat_matrix_market_huge_entry_count_is_not_preallocated() {
    // A header can claim any entry count; it must not be trusted as an allocation size
    let text = format!(
        "%%MatrixMarket matrix coordinate real general\n2 2 {}\n1 1 1\n",
        usize::MAX
    );
    let (line, message) = parse_failure(&text);
    assert!(line == 3);
    assert!(message == format!("expected {} entries, found 1", usize::MAX));
}

#[test]
fn sparsemat_scale_rows_matches_diagonal_product() {
    let a = pseudo_random(5, 4, 61, 50);
    let factors = [2.0, -1.0, 0.5, 3.0, -4.0];
    let expected = SparseMatrix::from_diagonal(factors).matmul(&a).unwrap();
    assert!(a.scale_rows(&factors).unwrap() == expected);
}

#[test]
fn sparsemat_scale_cols_matches_diagonal_product() {
    let a = pseudo_random(5, 4, 62, 50);
    let factors = [2.0, -1.0, 0.5, 3.0];
    let expected = a.matmul(&SparseMatrix::from_diagonal(factors)).unwrap();
    assert!(a.scale_cols(&factors).unwrap() == expected);
}

#[test]
fn sparsemat_scale_rows_cols_drop_zeros_and_keep_operand() {
    let a = SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    let rows = a.scale_rows(&[0.0, 1.0]).unwrap();
    assert!(rows.to_dense() == vec![vec![0.0, 0.0], vec![3.0, 4.0]]);
    assert!(rows.num_nonzero() == 2);
    let cols = a.scale_cols(&[1.0, 0.0]).unwrap();
    assert!(cols.to_dense() == vec![vec![1.0, 0.0], vec![3.0, 0.0]]);
    assert!(cols.num_nonzero() == 2);
    // The operand is untouched, and an empty matrix scales to an empty matrix
    assert!(a == SparseMatrix::from([[1.0, 2.0], [3.0, 4.0]]));
    let empty = SparseMatrix::<f64>::empty_with_shape(0, 3);
    assert!(empty.scale_rows(&[]).unwrap().shape == (0, 3));
}

#[test]
fn sparsemat_scale_rows_cols_complex() {
    use num_complex::Complex;
    // diag(V) · A with complex voltages, as in power flow Jacobians
    let a = SparseMatrix::from([
        [Complex::new(1.0, 1.0), Complex::new(0.0, 0.0)],
        [Complex::new(0.0, 0.0), Complex::new(2.0, 0.0)],
    ]);
    let v = [Complex::new(0.0, 1.0), Complex::new(1.0, -1.0)];
    let scaled = a.scale_rows(&v).unwrap();
    assert!(scaled[(0, 0)] == Complex::new(-1.0, 1.0));
    assert!(scaled[(1, 1)] == Complex::new(2.0, -2.0));
    assert!(a.scale_cols(&v).unwrap() == a.matmul(&SparseMatrix::from_diagonal(v)).unwrap());
}

#[test]
fn sparsemat_scale_rows_cols_shape_mismatch() {
    let a = SparseMatrix::<f64>::empty_with_shape(2, 3);
    assert!(
        a.scale_rows(&[1.0; 3]).err()
            == Some(SparseError::ShapeMismatch {
                left: (3, 3),
                right: (2, 3)
            })
    );
    assert!(
        a.scale_cols(&[1.0; 2]).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 3),
                right: (2, 2)
            })
    );
}

#[test]
fn sparsemat_select_matches_submatrix_for_ranges() {
    let a = pseudo_random(6, 7, 71, 45);
    let rows: Vec<usize> = (1..5).collect();
    let cols: Vec<usize> = (2..7).collect();
    assert!(a.select(&rows, &cols).unwrap() == a.submatrix(1..5, 2..7).unwrap());
}

#[test]
fn sparsemat_select_reorders_and_skips() {
    let a = SparseMatrix::from([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
    // Rows 2 then 0, columns 2 then 1: non-contiguous and reversed
    let picked = a.select(&[2, 0], &[2, 1]).unwrap();
    assert!(picked.to_dense() == vec![vec![9.0, 8.0], vec![3.0, 2.0]]);
}

#[test]
fn sparsemat_select_repeated_indices() {
    let a = SparseMatrix::from([[1.0, 0.0], [0.0, 2.0]]);
    let repeated = a.select(&[1, 1, 0], &[0, 1, 1]).unwrap();
    assert!(
        repeated.to_dense()
            == vec![
                vec![0.0, 2.0, 2.0],
                vec![0.0, 2.0, 2.0],
                vec![1.0, 0.0, 0.0]
            ]
    );
}

#[test]
fn sparsemat_select_matches_dense_selection() {
    let a = pseudo_random(8, 8, 72, 40);
    let dense = a.to_dense();
    // A PV/PQ-style split: an arbitrary subset in a scrambled order
    let rows = [7, 1, 4, 2];
    let cols = [3, 0, 6, 5, 1];
    let expected: Vec<Vec<f64>> = rows
        .iter()
        .map(|&r| cols.iter().map(|&c| dense[r][c]).collect())
        .collect();
    assert!(a.select(&rows, &cols).unwrap().to_dense() == expected);
}

#[test]
fn sparsemat_select_empty_and_out_of_range() {
    let a = sample_4x6();
    let none = a.select(&[], &[0, 1]).unwrap();
    assert!(none.shape == (0, 2) && none.num_nonzero() == 0);
    assert!(a.select(&[0, 4], &[0]).err() == Some(SparseError::OutOfRange { index: 4, len: 4 }));
    assert!(a.select(&[0], &[9]).err() == Some(SparseError::OutOfRange { index: 9, len: 6 }));
}

fn inverse_of(perm: &[usize]) -> Vec<usize> {
    let mut inverse = vec![0; perm.len()];
    for (position, &index) in perm.iter().enumerate() {
        inverse[index] = position;
    }
    inverse
}

#[test]
fn sparsemat_permute_matches_select() {
    let a = pseudo_random(5, 4, 81, 50);
    let (rows, cols) = ([3, 0, 4, 1, 2], [2, 3, 0, 1]);
    assert!(a.permute(&rows, &cols).unwrap() == a.select(&rows, &cols).unwrap());
}

#[test]
fn sparsemat_permute_then_inverse_is_identity() {
    let a = pseudo_random(6, 5, 82, 45);
    let (rows, cols) = ([5, 2, 0, 4, 1, 3], [4, 0, 3, 1, 2]);
    let permuted = a.permute(&rows, &cols).unwrap();
    assert!(permuted != a);
    let restored = permuted
        .permute(&inverse_of(&rows), &inverse_of(&cols))
        .unwrap();
    assert!(restored == a);
    // The identity permutation changes nothing
    assert!(a.permute(&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 4]).unwrap() == a);
}

#[test]
fn sparsemat_permute_symmetric_preserves_structure() {
    let base = pseudo_random(7, 7, 83, 40);
    let symmetric = &base + &base.create_transpose();
    let perm = [6, 3, 0, 5, 1, 4, 2];
    let permuted = symmetric.permute_symmetric(&perm).unwrap();
    assert!(permuted.is_symmetric());
    assert!(permuted.trace() == symmetric.trace());
    let mut before = symmetric.diagonal().to_dense();
    let mut after = permuted.diagonal().to_dense();
    before.sort_by(f64::total_cmp);
    after.sort_by(f64::total_cmp);
    assert!(before == after);
    // P A Pᵀ written with explicit permutation matrices
    let p = SparseMatrix::from_triplets(7, 7, perm.iter().enumerate().map(|(i, &j)| (i, j, 1.0)))
        .unwrap();
    assert!(permuted == &(&p * &symmetric) * &p.create_transpose());
}

#[test]
fn sparsemat_permute_rows_agrees_with_lu_permutation() {
    let a = dominant(6, 84);
    let lu = a.lu().unwrap();
    let columns: Vec<usize> = (0..6).collect();
    let pa = a.permute(lu.permutation(), &columns).unwrap();
    assert!((&pa - &(lu.l() * lu.u())).norm_frobenius() < 1e-9);
}

#[test]
fn sparsemat_permute_errors() {
    let a = SparseMatrix::<f64>::empty_with_shape(3, 2);
    let reason = |result: Result<SparseMatrix, SparseError>| match result.err() {
        Some(SparseError::InvalidPermutation { reason }) => reason,
        other => panic!("expected InvalidPermutation, got {other:?}"),
    };
    assert!(
        reason(a.permute(&[0, 1], &[0, 1])) == "length must equal the dimension being permuted"
    );
    assert!(reason(a.permute(&[0, 1, 3], &[0, 1])) == "index out of range");
    assert!(reason(a.permute(&[0, 1, 1], &[0, 1])) == "index repeated");
    assert!(reason(a.permute(&[0, 1, 2], &[1, 1])) == "index repeated");
    assert!(
        a.permute_symmetric(&[0, 1, 2]).err() == Some(SparseError::NotSquare { shape: (3, 2) })
    );
    // Empty matrices have exactly one permutation
    let empty = SparseMatrix::<f64>::new();
    assert!(empty.permute_symmetric(&[]).unwrap().shape == (0, 0));
}

fn pseudo_random_complex(
    rows: usize,
    cols: usize,
    seed: u64,
) -> SparseMatrix<num_complex::Complex<f64>> {
    let re = pseudo_random(rows, cols, seed, 40);
    let im = pseudo_random(rows, cols, seed + 1000, 40);
    SparseMatrix::from_parts(&re, &im).unwrap()
}

#[test]
fn sparsemat_from_parts_round_trips_real_and_imag() {
    let a = pseudo_random_complex(6, 5, 91);
    let (re, im) = (a.real(), a.imag());
    assert!(re.shape == (6, 5) && im.shape == (6, 5));
    assert!(SparseMatrix::from_parts(&re, &im).unwrap() == a);
}

#[test]
fn sparsemat_real_and_imag_drop_zero_parts() {
    use num_complex::Complex;
    let a = SparseMatrix::from([
        [Complex::new(1.0, 0.0), Complex::new(0.0, 2.0)],
        [Complex::new(3.0, -4.0), Complex::new(0.0, 0.0)],
    ]);
    let re = a.real();
    let im = a.imag();
    assert!(re.to_dense() == vec![vec![1.0, 0.0], vec![3.0, 0.0]]);
    assert!(im.to_dense() == vec![vec![0.0, 2.0], vec![-4.0, 0.0]]);
    // The purely imaginary entry is absent from the real part, and vice versa
    assert!(re.num_nonzero() == 2 && im.num_nonzero() == 2);
}

#[test]
fn sparsemat_conj_and_adjoint() {
    use num_complex::Complex;
    let a = SparseMatrix::from([
        [Complex::new(1.0, 2.0), Complex::new(0.0, -1.0)],
        [Complex::new(0.0, 0.0), Complex::new(3.0, 0.0)],
    ]);
    let conjugated = a.conj();
    assert!(conjugated[(0, 0)] == Complex::new(1.0, -2.0));
    assert!(conjugated[(0, 1)] == Complex::new(0.0, 1.0));
    assert!(conjugated.conj() == a);

    let adjoint = a.adjoint();
    assert!(adjoint == a.conj().create_transpose());
    assert!(adjoint[(1, 0)] == Complex::new(0.0, 1.0));
    assert!(adjoint.adjoint() == a);
}

#[test]
fn sparsemat_adjoint_reverses_products() {
    // (AB)ᴴ = Bᴴ Aᴴ, on a rectangular pair
    let a = pseudo_random_complex(4, 3, 92);
    let b = pseudo_random_complex(3, 5, 93);
    assert!((&a * &b).adjoint() == &b.adjoint() * &a.adjoint());
}

#[test]
fn sparsemat_hermitian_matrix_is_its_own_adjoint() {
    let a = pseudo_random_complex(5, 5, 94);
    let hermitian = &a + &a.adjoint();
    assert!(hermitian.adjoint() == hermitian);
    // Its diagonal is real
    assert!(hermitian
        .diagonal()
        .values()
        .iter()
        .all(|val| val.im == 0.0));
}

#[test]
fn sparsemat_complex_helpers_on_real_matrices() {
    let a = pseudo_random(4, 6, 95, 50);
    assert!(a.conj() == a);
    assert!(a.adjoint() == a.create_transpose());
    assert!(a.real() == a);
    let imaginary = a.imag();
    assert!(imaginary.shape == (4, 6) && imaginary.num_nonzero() == 0);

    let single = SparseMatrix::from([[1.0f32, 0.0], [0.0, -2.0]]);
    assert!(single.adjoint() == single.create_transpose());
}

#[test]
fn sparsemat_from_parts_shape_mismatch() {
    let re = SparseMatrix::<f64>::empty_with_shape(2, 3);
    let im = SparseMatrix::<f64>::empty_with_shape(3, 2);
    assert!(
        SparseMatrix::from_parts(&re, &im).err()
            == Some(SparseError::ShapeMismatch {
                left: (2, 3),
                right: (3, 2)
            })
    );
}
