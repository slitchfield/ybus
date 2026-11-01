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

/// The factors of [`reference_lu`]: `P A = L U`, where row `i` of `P A` is row `perm[i]` of `A`,
/// `L` is unit lower triangular and `U` upper triangular, both in CSR form and storing no zeros
/// (except `L`'s unit diagonal), and `swaps` counts the row swaps.
pub struct ReferenceLu {
    pub l: CsrMatrix<f64>,
    pub u: CsrMatrix<f64>,
    pub perm: Vec<usize>,
    pub swaps: usize,
}

/// The LU factorization `sparse_mat` had before Phase 3: right-looking Gaussian elimination with
/// true partial pivoting (the first row of largest magnitude below the diagonal), on rows held
/// as `BTreeMap`s, with entries that come out exactly zero dropped. Unordered, and at least
/// quadratic in the dimension. Errors as the library's LU reported them: `NotSquare`, or
/// `SingularMatrix` naming the elimination step with no nonzero pivot.
pub fn reference_lu(a: &CsrMatrix<f64>) -> Result<ReferenceLu, SparseError> {
    let (n, cols) = a.shape();
    if n != cols {
        return Err(SparseError::NotSquare { shape: (n, cols) });
    }
    let (a_ptr, a_idx, a_values) = (a.pointers(), a.indices(), a.values());
    let mut rows: Vec<BTreeMap<usize, f64>> = (0..n)
        .map(|i| {
            (a_ptr[i]..a_ptr[i + 1])
                .map(|p| (a_idx[p], a_values[p]))
                .collect()
        })
        .collect();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut swaps = 0;
    for k in 0..n {
        let magnitude = |row: &BTreeMap<usize, f64>| row.get(&k).map_or(0.0, |v| v.abs());
        let mut best = k;
        let mut best_magnitude = magnitude(&rows[k]);
        for (i, row) in rows.iter().enumerate().skip(k + 1) {
            if magnitude(row) > best_magnitude {
                best = i;
                best_magnitude = magnitude(row);
            }
        }
        if best_magnitude == 0.0 {
            return Err(SparseError::SingularMatrix { pivot: k });
        }
        if best != k {
            rows.swap(k, best);
            perm.swap(k, best);
            swaps += 1;
        }
        let pivot = rows[k][&k];
        let pivot_row: Vec<(usize, f64)> = rows[k].range(k + 1..).map(|(&j, &v)| (j, v)).collect();
        for row in rows.iter_mut().skip(k + 1) {
            let Some(&below) = row.get(&k) else {
                continue;
            };
            let multiplier = below / pivot;
            if multiplier == 0.0 {
                row.remove(&k);
                continue;
            }
            row.insert(k, multiplier);
            for &(j, u_kj) in &pivot_row {
                let updated = row.get(&j).copied().unwrap_or(0.0) - multiplier * u_kj;
                if updated == 0.0 {
                    row.remove(&j);
                } else {
                    row.insert(j, updated);
                }
            }
        }
    }
    let (mut l_ptr, mut l_idx, mut l_data) = (vec![0], vec![], vec![]);
    let (mut u_ptr, mut u_idx, mut u_data) = (vec![0], vec![], vec![]);
    for (i, row) in rows.iter().enumerate() {
        for (&j, &value) in row.range(..i).filter(|(_, &value)| value != 0.0) {
            l_idx.push(j);
            l_data.push(value);
        }
        l_idx.push(i);
        l_data.push(1.0);
        l_ptr.push(l_idx.len());
        for (&j, &value) in row.range(i..).filter(|(_, &value)| value != 0.0) {
            u_idx.push(j);
            u_data.push(value);
        }
        u_ptr.push(u_idx.len());
    }
    Ok(ReferenceLu {
        l: CsrMatrix::try_from_raw_parts(n, n, l_ptr, l_idx, l_data)?,
        u: CsrMatrix::try_from_raw_parts(n, n, u_ptr, u_idx, u_data)?,
        perm,
        swaps,
    })
}
