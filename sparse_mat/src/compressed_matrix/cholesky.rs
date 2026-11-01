//! Ordered sparse Cholesky: symbolic analysis of a pattern, then numeric factorization.
//!
//! [`CsrMatrix::analyze_cholesky`] chooses an ordering `P` and works out, from the pattern alone,
//! the elimination tree of `P A Pᵀ` and the exact pattern of its Cholesky factor `L`. The
//! factorization then only computes values, so its cost and memory are known in advance.

use super::ordering::Adjacency;
use num_traits::Float;

use super::{CompressedMatrix, CsrMatrix, Ordering};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// The symbolic analysis of a symmetric matrix's pattern for an ordered Cholesky factorization
/// `P A Pᵀ = L Lᵀ`, from [`CsrMatrix::analyze_cholesky`].
///
/// It depends only on the pattern of `A + Aᵀ` (stored zeros included) and the ordering, so one
/// analysis serves every matrix with that pattern.
#[derive(Clone, Debug)]
pub struct SymbolicCholesky {
    /// Row `i` of `P A Pᵀ` is row `perm[i]` of `A`
    perm: Vec<usize>,
    /// The parent of each node of the elimination tree of `P A Pᵀ`, or `None` for a root
    parent: Vec<Option<usize>>,
    postorder: Vec<usize>,
    /// The pattern of `L` in CSR form: row `i`'s columns, increasing, ending with `i` itself
    l_ptr: Vec<usize>,
    l_idx: Vec<usize>,
    /// The pattern of the analyzed matrix, which a matrix being factored must have
    a_ptr: Vec<usize>,
    a_idx: Vec<usize>,
}

impl<T: Scalar> CsrMatrix<T> {
    /// Analyzes this matrix's pattern for a Cholesky factorization ordered by `ordering`:
    /// chooses the permutation `P`, builds the elimination tree of `P A Pᵀ`, and finds the
    /// exact pattern of its factor `L`.
    ///
    /// Only the pattern is read, symmetrized as `A + Aᵀ`, so the matrix needn't be symmetric
    /// yet; `SymbolicCholesky::factor` checks that. Costs O(nnz(L) log) after the ordering.
    /// Returns `SparseError::NotSquare` for a non-square matrix.
    pub fn analyze_cholesky(&self, ordering: Ordering) -> Result<SymbolicCholesky, SparseError> {
        let perm = self.ordering(ordering)?;
        let n = perm.len();
        let mut inverse = vec![0; n];
        for (new, &old) in perm.iter().enumerate() {
            inverse[old] = new;
        }
        let graph = Adjacency::of(self);
        // Row i of the strict lower triangle of P (A + Aᵀ) Pᵀ, in the new numbering
        let lower = |i: usize| {
            graph
                .neighbours(perm[i])
                .iter()
                .map(|&old| inverse[old])
                .filter(move |&j| j < i)
        };

        // The elimination tree, by Liu's algorithm: for each row i and each j < i in it, climb
        //   from j towards the root, compressing the path to i as it goes
        let mut parent = vec![None; n];
        let mut ancestor: Vec<Option<usize>> = vec![None; n];
        for i in 0..n {
            for mut j in lower(i) {
                while j != i {
                    let next = ancestor[j];
                    ancestor[j] = Some(i);
                    match next {
                        None => {
                            parent[j] = Some(i);
                            break;
                        }
                        Some(next) => j = next,
                    }
                }
            }
        }

        // Row i of L is the set of nodes on the tree paths from row i's entries up to i: its
        //   row subtree. Each row is sorted, so a factorization can take its entries in
        //   increasing column order (a valid order, since a parent always follows its children)
        let mut l_ptr = Vec::with_capacity(n + 1);
        l_ptr.push(0);
        let mut l_idx = vec![];
        let mut mark = vec![usize::MAX; n];
        for i in 0..n {
            let start = l_idx.len();
            mark[i] = i;
            for mut j in lower(i) {
                while mark[j] != i {
                    mark[j] = i;
                    l_idx.push(j);
                    j = parent[j].expect("a node below i in its row subtree has a parent");
                }
            }
            l_idx[start..].sort_unstable();
            l_idx.push(i);
            l_ptr.push(l_idx.len());
        }

        Ok(SymbolicCholesky {
            postorder: postorder(&parent),
            perm,
            parent,
            l_ptr,
            l_idx,
            a_ptr: self.ptr.clone(),
            a_idx: self.idx.clone(),
        })
    }
}

impl SymbolicCholesky {
    /// The dimension of the analyzed matrix.
    pub fn dim(&self) -> usize {
        self.perm.len()
    }

    /// The ordering `P`, in the convention of `permute_symmetric`: row and column `i` of
    /// `P A Pᵀ` are row and column `permutation()[i]` of `A`.
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// The elimination tree of `P A Pᵀ`: the parent of each node, or `None` for a root (one per
    /// connected component). A parent always has a larger index than its children.
    pub fn elimination_tree(&self) -> &[Option<usize>] {
        &self.parent
    }

    /// A postorder of the elimination tree: every node comes after all of its descendants, and
    /// each subtree's nodes are contiguous.
    pub fn postorder(&self) -> &[usize] {
        &self.postorder
    }

    /// The number of entries in each column of `L`, its diagonal included.
    pub fn column_counts(&self) -> Vec<usize> {
        let mut counts = vec![0; self.dim()];
        for &j in &self.l_idx {
            counts[j] += 1;
        }
        counts
    }

    /// The number of entries `L` will store, its diagonal included: exactly what the factor
    /// holds, since it stores its whole symbolic pattern, numerical zeros included.
    pub fn predicted_entries(&self) -> usize {
        self.l_idx.len()
    }
}

impl SymbolicCholesky {
    /// Factors `matrix`, which must have the analyzed pattern, as `P A Pᵀ = L Lᵀ`.
    ///
    /// The up-looking algorithm: row `i` of `L` solves a triangular system with the rows above
    /// it, taking its entries in increasing column order. `L` stores the whole symbolic pattern,
    /// so an entry that happens to be zero stays stored, and its memory is known in advance.
    /// With `Ordering::Natural`, every operation happens in the order the previous, unordered
    /// Cholesky used, so its values are the same bit for bit.
    ///
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, `SparseError::NotSymmetric` (judged mathematically, ignoring stored zeros),
    /// or `SparseError::NotPositiveDefinite` if a pivot is not strictly positive, naming its row
    /// of `A` (in the original numbering).
    pub fn factor<T: Scalar + Float>(
        &self,
        matrix: &CsrMatrix<T>,
    ) -> Result<Cholesky<T>, SparseError> {
        let n = self.dim();
        if matrix.shape != (n, n) || matrix.ptr != self.a_ptr || matrix.idx != self.a_idx {
            return Err(SparseError::PatternMismatch);
        }
        if *matrix != matrix.clone().transpose() {
            return Err(SparseError::NotSymmetric);
        }
        let c = matrix.permute_symmetric(&self.perm)?;

        // L's columns below the diagonal, filled as the rows are computed: column j's entries
        //   from rows already done are col_rows/col_values[col_start[j]..col_fill[j]]
        let mut col_start = vec![0; n + 1];
        for &j in &self.l_idx {
            col_start[j + 1] += 1;
        }
        for j in 0..n {
            // Less the diagonal, which isn't kept here
            col_start[j + 1] += col_start[j];
        }
        for (j, start) in col_start.iter_mut().enumerate() {
            *start -= j;
        }
        let mut col_fill = col_start.clone();
        let below = self.l_idx.len() - n;
        let mut col_rows = vec![0; below];
        let mut col_values = vec![T::zero(); below];

        let mut l_values = vec![T::zero(); self.l_idx.len()];
        let mut diagonal = vec![T::zero(); n];
        let mut work = vec![T::zero(); n];
        for i in 0..n {
            // Row i of P A Pᵀ: the part left of the diagonal goes into the work vector
            let mut remaining = T::zero();
            for p in c.ptr[i]..c.ptr[i + 1] {
                let col = c.idx[p];
                if col < i {
                    work[col] = c.data[p];
                } else {
                    if col == i {
                        remaining = c.data[p];
                    }
                    break;
                }
            }
            // Every position the updates reach is in row i's pattern, so taking the pattern's
            //   entries leaves the work vector zeroed for the next row
            let diagonal_slot = self.l_ptr[i + 1] - 1;
            let row = self.l_ptr[i]..diagonal_slot;
            for (&j, l_value) in self.l_idx[row.clone()].iter().zip(&mut l_values[row]) {
                let l_ij = work[j] / diagonal[j];
                work[j] = T::zero();
                for q in col_start[j]..col_fill[j] {
                    let m = col_rows[q];
                    work[m] = work[m] - col_values[q] * l_ij;
                }
                col_rows[col_fill[j]] = i;
                col_values[col_fill[j]] = l_ij;
                col_fill[j] += 1;
                *l_value = l_ij;
                remaining = remaining - l_ij * l_ij;
            }
            if remaining.is_nan() || remaining <= T::zero() {
                return Err(SparseError::NotPositiveDefinite {
                    pivot: self.perm[i],
                });
            }
            diagonal[i] = remaining.sqrt();
            l_values[diagonal_slot] = diagonal[i];
        }

        // Lᵀ in CSR form is L in CSC form: each column's diagonal, then the entries below it
        let mut lt_ptr = Vec::with_capacity(n + 1);
        lt_ptr.push(0);
        let mut lt_idx = Vec::with_capacity(self.l_idx.len());
        let mut lt_values = Vec::with_capacity(self.l_idx.len());
        for j in 0..n {
            lt_idx.push(j);
            lt_values.push(diagonal[j]);
            lt_idx.extend_from_slice(&col_rows[col_start[j]..col_fill[j]]);
            lt_values.extend_from_slice(&col_values[col_start[j]..col_fill[j]]);
            lt_ptr.push(lt_idx.len());
        }
        Ok(Cholesky {
            l: CompressedMatrix::from_parts_unchecked(
                (n, n),
                self.l_ptr.clone(),
                self.l_idx.clone(),
                l_values,
            ),
            l_transpose: CompressedMatrix::from_parts_unchecked((n, n), lt_ptr, lt_idx, lt_values),
            perm: self.perm.clone(),
        })
    }
}

/// An ordered Cholesky factorization `P A Pᵀ = L Lᵀ` of a symmetric positive definite matrix,
/// from [`SymbolicCholesky::factor`] or [`CsrMatrix::cholesky`]. Keep it to solve many
/// right-hand sides against the same matrix.
///
/// It holds `L` in CSR form, and also `Lᵀ` in CSR form for the back substitution, so it takes
/// twice the factor's memory.
#[derive(Clone, Debug)]
pub struct Cholesky<T: Scalar = f64> {
    l: CsrMatrix<T>,
    l_transpose: CsrMatrix<T>,
    perm: Vec<usize>,
}

impl<T: Scalar + Float> Cholesky<T> {
    /// The lower triangular factor `L` of `P A Pᵀ`, with a positive diagonal. It stores its
    /// whole symbolic pattern, so some entries may be zero.
    pub fn l(&self) -> &CsrMatrix<T> {
        &self.l
    }

    /// The ordering `P`: row and column `i` of `P A Pᵀ` are row and column `permutation()[i]`
    /// of `A`.
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// Solves `A x = b`: permutes `b`, substitutes forward with `L` and back with `Lᵀ`, and
    /// permutes the result back. Returns `SparseError::ShapeMismatch` unless `b` has one entry
    /// per row of `A`.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        if b.len() != self.perm.len() {
            return Err(SparseError::ShapeMismatch {
                left: self.l.shape,
                right: (b.len(), 1),
            });
        }
        let permuted: Vec<T> = self.perm.iter().map(|&row| b[row]).collect();
        let y = self.l.solve_lower_triangular(&permuted)?;
        let z = self.l_transpose.solve_upper_triangular(&y)?;
        let mut x = vec![T::zero(); z.len()];
        for (&row, value) in self.perm.iter().zip(z) {
            x[row] = value;
        }
        Ok(x)
    }
}

impl<T: Scalar + Float> CsrMatrix<T> {
    /// Factors a symmetric positive definite matrix as `P A Pᵀ = L Lᵀ`, ordered by reverse
    /// Cuthill–McKee: `analyze_cholesky(Ordering::Rcm)` followed by `factor`. Analyze once
    /// instead to factor several matrices with the same pattern.
    ///
    /// Returns `SparseError::NotSquare`, `SparseError::NotSymmetric`, or
    /// `SparseError::NotPositiveDefinite` naming the row of `A` whose pivot failed.
    pub fn cholesky(&self) -> Result<Cholesky<T>, SparseError> {
        self.analyze_cholesky(Ordering::Rcm)?.factor(self)
    }

    /// Solves `A x = b` for a symmetric positive definite matrix by an ordered Cholesky
    /// factorization (`cholesky`), failing with its errors, or `SparseError::ShapeMismatch` if
    /// `b` doesn't have one entry per row. Factor once with `cholesky` instead when solving
    /// several right-hand sides.
    pub fn solve_spd(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.check_system(b.len())?;
        self.cholesky()?.solve(b)
    }
}

/// A postorder of the forest `parent` describes: each node's children are visited in increasing
/// order, depth first, and the trees in increasing order of their roots.
fn postorder(parent: &[Option<usize>]) -> Vec<usize> {
    let n = parent.len();
    // Children lists, built so that each list is in increasing order
    let mut first_child = vec![None; n];
    let mut next_sibling = vec![None; n];
    for node in (0..n).rev() {
        if let Some(p) = parent[node] {
            next_sibling[node] = first_child[p];
            first_child[p] = Some(node);
        }
    }
    let mut order = Vec::with_capacity(n);
    let mut stack = vec![];
    for root in (0..n).filter(|&node| parent[node].is_none()) {
        stack.push(root);
        while let Some(&node) = stack.last() {
            match first_child[node] {
                // Descend into the next unvisited child, unlinking it from the list
                Some(child) => {
                    first_child[node] = next_sibling[child];
                    stack.push(child);
                }
                None => {
                    stack.pop();
                    order.push(node);
                }
            }
        }
    }
    order
}
