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
#[should_panic(expected = "incompatible shapes")]
fn sparsemat_sub_shape_mismatch_panics() {
    let _ = &sample_2x2() - &SparseMatrix::<f64>::empty_with_shape(3, 3);
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
fn sparsemat_scalar_multiplication_integers() {
    let a = SparseMatrix::from([[1i32, 2], [0, 3]]);
    assert!((&a * 4).to_dense() == vec![vec![4, 8], vec![0, 12]]);
    assert!((3i32 * &a).to_dense() == vec![vec![3, 6], vec![0, 9]]);
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
