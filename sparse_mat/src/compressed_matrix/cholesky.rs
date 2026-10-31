//! Ordered sparse Cholesky: symbolic analysis of a pattern, then numeric factorization.
//!
//! [`CsrMatrix::analyze_cholesky`] chooses an ordering `P` and works out, from the pattern alone,
//! the elimination tree of `P A Pᵀ` and the exact pattern of its Cholesky factor `L`. The
//! factorization then only computes values, so its cost and memory are known in advance.

use super::ordering::Adjacency;
use super::{CsrMatrix, Ordering};
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
    /// The pattern of `L` in CSR form, without its row pointers: row by row, each row's
    /// columns increasing and ending with the row itself
    l_idx: Vec<usize>,
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
        }

        Ok(SymbolicCholesky {
            postorder: postorder(&parent),
            perm,
            parent,
            l_idx,
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
