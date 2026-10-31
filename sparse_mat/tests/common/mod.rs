//! Reference implementations for the factorization tests: the algorithms `sparse_mat` used
//! before Phase 3, rewritten against the public API, so the new ones can be compared with them.

use std::collections::BTreeMap;

use sparse_mat::compressed_matrix::CsrMatrix;
use sparse_mat::sparse_error::SparseError;

/// The unordered up-looking Cholesky factorization `sparse_mat` had before Phase 3: `A = L Lᵀ`,
/// row by row, with each row's work held in a `BTreeMap` and entries that come out exactly zero
/// dropped. Errors as the library's `cholesky` reports them, in the matrix's own numbering.
pub fn reference_cholesky(a: &CsrMatrix<f64>) -> Result<CsrMatrix<f64>, SparseError> {
    let (n, cols) = a.shape();
    if n != cols {
        return Err(SparseError::NotSquare { shape: (n, cols) });
    }
    if *a != a.clone().transpose() {
        return Err(SparseError::NotSymmetric);
    }
    let (a_ptr, a_idx, a_values) = (a.pointers(), a.indices(), a.values());
    let mut diagonal = vec![0.0; n];
    let mut columns: Vec<Vec<(usize, f64)>> = vec![vec![]; n];
    let (mut ptr, mut idx, mut data) = (vec![0], vec![], vec![]);
    for i in 0..n {
        let mut work: BTreeMap<usize, f64> = BTreeMap::new();
        let mut a_ii = 0.0;
        for p in a_ptr[i]..a_ptr[i + 1] {
            let col = a_idx[p];
            if col < i {
                work.insert(col, a_values[p]);
            } else {
                if col == i {
                    a_ii = a_values[p];
                }
                break;
            }
        }
        let mut remaining = a_ii;
        while let Some((j, w_j)) = work.pop_first() {
            let l_ij = w_j / diagonal[j];
            if l_ij == 0.0 {
                continue;
            }
            for &(m, l_mj) in &columns[j] {
                let current = work.get(&m).copied().unwrap_or(0.0);
                work.insert(m, current - l_mj * l_ij);
            }
            columns[j].push((i, l_ij));
            idx.push(j);
            data.push(l_ij);
            remaining -= l_ij * l_ij;
        }
        if remaining.is_nan() || remaining <= 0.0 {
            return Err(SparseError::NotPositiveDefinite { pivot: i });
        }
        diagonal[i] = remaining.sqrt();
        idx.push(i);
        data.push(diagonal[i]);
        ptr.push(idx.len());
    }
    CsrMatrix::try_from_raw_parts(n, n, ptr, idx, data)
}

/// The solve `solve_spd` did before Phase 3: forward substitution with the unordered factor,
/// then back substitution with its transpose, reoriented to CSR.
pub fn reference_solve_spd(a: &CsrMatrix<f64>, b: &[f64]) -> Result<Vec<f64>, SparseError> {
    let l = reference_cholesky(a)?;
    let y = l.solve_lower_triangular(b)?;
    let l_transposed: CsrMatrix<f64> = l.to_csc().transpose();
    l_transposed.solve_upper_triangular(&y)
}
