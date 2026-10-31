#![cfg(feature = "rayon")]

use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix};
use sparse_mat::sparse_error::SparseError;

/// Deterministic pseudo-random sparse matrix with small entries.
fn pseudo_random(rows: usize, cols: usize, seed: u64, fill_percent: u64) -> CsrMatrix<f64> {
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
                triplets.push((row, col, (next() % 9) as f64 - 4.5));
            }
        }
    }
    CsrMatrix::from_triplets(rows, cols, triplets).unwrap()
}

/// The same arrays, so the comparison covers the pattern and stored zeros too.
fn assert_same_arrays(a: &CsrMatrix<f64>, b: &CsrMatrix<f64>) {
    assert_eq!(a.shape(), b.shape());
    assert_eq!(a.pointers(), b.pointers());
    assert_eq!(a.indices(), b.indices());
    assert_eq!(a.values(), b.values());
}

#[test]
fn par_mul_vec_equals_sequential_exactly() {
    for seed in 1..=5 {
        let a = pseudo_random(60, 45, seed, 20);
        let x: Vec<f64> = (0..45).map(|i| (i as f64).sin()).collect();
        assert_eq!(a.par_mul_vec(&x).unwrap(), a.mul_vec(&x).unwrap());
    }
}

#[test]
fn par_vec_mul_equals_sequential_exactly() {
    for seed in 1..=5 {
        let a = pseudo_random(45, 60, seed, 20).to_csc();
        let x: Vec<f64> = (0..45).map(|i| (i as f64).cos()).collect();
        assert_eq!(a.par_vec_mul(&x).unwrap(), a.vec_mul(&x).unwrap());
    }
}

#[test]
fn par_matmul_equals_sequential_exactly() {
    for seed in 1..=4 {
        let a = pseudo_random(40, 30, seed, 15);
        let b = pseudo_random(30, 50, seed + 50, 15);
        assert_same_arrays(&a.par_matmul(&b).unwrap(), &a.matmul(&b).unwrap());
        let by_columns = a.to_csc().par_matmul(&b.to_csc()).unwrap();
        assert_same_arrays(&by_columns.to_csr(), &a.matmul(&b).unwrap());
    }
}

#[test]
fn par_matmul_keeps_cancellations_as_stored_zeros() {
    let a = CsrMatrix::from_dense(&[[1.0, -1.0], [0.0, 2.0]]).unwrap();
    let b = CsrMatrix::from_dense(&[[1.0, 0.0], [1.0, 5.0]]).unwrap();
    let product = a.par_matmul(&b).unwrap();
    assert_same_arrays(&product, &a.matmul(&b).unwrap());
    assert_eq!((product.num_stored(), product.num_nonzero()), (4, 3));
}

#[test]
fn par_products_report_shape_mismatches() {
    let a = CsrMatrix::<f64>::from_triplets(2, 3, Vec::new()).unwrap();
    assert_eq!(
        a.par_mul_vec(&[1.0]).err(),
        Some(SparseError::ShapeMismatch {
            left: (2, 3),
            right: (1, 1)
        })
    );
    assert_eq!(
        a.to_csc().par_vec_mul(&[1.0, 2.0, 3.0]).err(),
        Some(SparseError::ShapeMismatch {
            left: (1, 3),
            right: (2, 3)
        })
    );
    assert_eq!(
        a.par_matmul(&a).err(),
        Some(SparseError::ShapeMismatch {
            left: (2, 3),
            right: (2, 3)
        })
    );
}

#[test]
fn par_products_on_empty_matrices() {
    let a = CsrMatrix::<f64>::default();
    assert_eq!(a.par_mul_vec(&[]).unwrap(), Vec::<f64>::new());
    assert_eq!(
        CscMatrix::<f64>::default().par_vec_mul(&[]).unwrap(),
        Vec::<f64>::new()
    );
    assert_eq!(a.par_matmul(&a).unwrap().shape(), (0, 0));
}
