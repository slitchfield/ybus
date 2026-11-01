//! Fill-reducing orderings: permutations that make a matrix cheaper to factor.
//!
//! An ordering depends only on the pattern of `A + Aᵀ` (which positions are stored, stored zeros
//! included, ignoring the diagonal), so it is the same for a matrix in either orientation and
//! for any values. It is returned as a permutation in the convention of
//! [`CompressedMatrix::permute_symmetric`]: entry `(i, j)` of the reordered matrix is entry
//! `(perm[i], perm[j])` of the original.

use std::collections::VecDeque;

use super::{CompressedMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// Which ordering to compute, or to factor with.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Ordering {
    /// No reordering: the identity permutation.
    Natural,
    /// Reverse Cuthill–McKee: a breadth-first numbering from a pseudo-peripheral node, reversed,
    /// which keeps entries near the diagonal (small bandwidth and profile).
    Rcm,
}

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// The permutation `ordering` chooses for this matrix, in the convention of
    /// `permute_symmetric`. It depends only on the pattern of `A + Aᵀ`, so it is the same in
    /// either orientation and whatever the values.
    ///
    /// Returns `SparseError::NotSquare` for a non-square matrix.
    pub fn ordering(&self, ordering: Ordering) -> Result<Vec<usize>, SparseError> {
        if !self.is_square() {
            return Err(SparseError::NotSquare { shape: self.shape });
        }
        Ok(match ordering {
            Ordering::Natural => (0..self.shape.0).collect(),
            Ordering::Rcm => reverse_cuthill_mckee(&Adjacency::of(self)),
        })
    }
}

/// The graph of `A + Aᵀ` without self-loops: node `i`'s neighbours are
/// `neighbours[start[i]..start[i + 1]]`, sorted and without duplicates.
pub(crate) struct Adjacency {
    pub(crate) start: Vec<usize>,
    pub(crate) neighbours: Vec<usize>,
}

impl Adjacency {
    /// The graph of a square matrix's symmetrized pattern, in O(nnz + n).
    pub(crate) fn of<T: Scalar, O: Orientation>(matrix: &CompressedMatrix<T, O>) -> Adjacency {
        let n = matrix.shape.0;
        // Each off-diagonal entry (i, j) gives the edges i–j and j–i; count, then fill
        let mut start = vec![0; n + 1];
        for (row, col, _) in matrix.iter().filter(|&(row, col, _)| row != col) {
            start[row + 1] += 1;
            start[col + 1] += 1;
        }
        for i in 0..n {
            start[i + 1] += start[i];
        }
        let mut cursor = start.clone();
        let mut neighbours = vec![0; start[n]];
        for (row, col, _) in matrix.iter().filter(|&(row, col, _)| row != col) {
            neighbours[cursor[row]] = col;
            cursor[row] += 1;
            neighbours[cursor[col]] = row;
            cursor[col] += 1;
        }
        // Sort each list and drop the duplicates a symmetric pattern produces, compacting in
        //   place
        let mut write = 0;
        let mut compact_start = vec![0; n + 1];
        for i in 0..n {
            let list = &mut neighbours[start[i]..start[i + 1]];
            list.sort_unstable();
            let mut previous = None;
            for k in start[i]..start[i + 1] {
                let node = neighbours[k];
                if previous != Some(node) {
                    neighbours[write] = node;
                    write += 1;
                    previous = Some(node);
                }
            }
            compact_start[i + 1] = write;
        }
        neighbours.truncate(write);
        Adjacency {
            start: compact_start,
            neighbours,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.start.len() - 1
    }

    pub(crate) fn neighbours(&self, node: usize) -> &[usize] {
        &self.neighbours[self.start[node]..self.start[node + 1]]
    }

    pub(crate) fn degree(&self, node: usize) -> usize {
        self.start[node + 1] - self.start[node]
    }
}

/// Reverse Cuthill–McKee: for each connected component, a breadth-first numbering from a
/// pseudo-peripheral node that visits each node's unnumbered neighbours in increasing degree
/// (ties by index, so the result is deterministic); the whole numbering is then reversed.
fn reverse_cuthill_mckee(graph: &Adjacency) -> Vec<usize> {
    let n = graph.len();
    let mut order = Vec::with_capacity(n);
    let mut numbered = vec![false; n];
    // Scratch for the level structures, reused across searches
    let mut level = vec![usize::MAX; n];
    let mut queue = VecDeque::new();
    let mut candidates = vec![];
    for seed in 0..n {
        if numbered[seed] {
            continue;
        }
        let root = pseudo_peripheral(graph, seed, &mut level, &mut queue);
        // Cuthill–McKee from the root
        let first = order.len();
        numbered[root] = true;
        order.push(root);
        let mut head = first;
        while head < order.len() {
            let node = order[head];
            head += 1;
            candidates.clear();
            candidates.extend(
                graph
                    .neighbours(node)
                    .iter()
                    .copied()
                    .filter(|&next| !numbered[next]),
            );
            candidates.sort_unstable_by_key(|&next| (graph.degree(next), next));
            for &next in &candidates {
                numbered[next] = true;
                order.push(next);
            }
        }
    }
    order.reverse();
    order
}

/// A pseudo-peripheral node of `seed`'s component, by George and Liu's method: from the
/// current node, build its level structure (the breadth-first distances), move to a node of
/// smallest degree in the last level, and repeat while the number of levels grows.
fn pseudo_peripheral(
    graph: &Adjacency,
    seed: usize,
    level: &mut [usize],
    queue: &mut VecDeque<usize>,
) -> usize {
    let mut root = seed;
    let (mut depth, mut last_level) = level_structure(graph, root, level, queue);
    loop {
        let candidate = last_level
            .iter()
            .copied()
            .min_by_key(|&node| (graph.degree(node), node))
            .expect("the last level holds at least the root");
        let (candidate_depth, candidate_last) = level_structure(graph, candidate, level, queue);
        if candidate_depth > depth {
            root = candidate;
            depth = candidate_depth;
            last_level = candidate_last;
        } else {
            return root;
        }
    }
}

/// The breadth-first level structure from `root`: its number of levels, and the nodes of the
/// last one. Leaves `level` reset to `usize::MAX` for the nodes it visited.
fn level_structure(
    graph: &Adjacency,
    root: usize,
    level: &mut [usize],
    queue: &mut VecDeque<usize>,
) -> (usize, Vec<usize>) {
    let mut visited = vec![root];
    level[root] = 0;
    queue.clear();
    queue.push_back(root);
    while let Some(node) = queue.pop_front() {
        for &next in graph.neighbours(node) {
            if level[next] == usize::MAX {
                level[next] = level[node] + 1;
                visited.push(next);
                queue.push_back(next);
            }
        }
    }
    let depth = level[*visited.last().expect("the root was visited")];
    let last_level = visited
        .iter()
        .copied()
        .filter(|&node| level[node] == depth)
        .collect();
    for &node in &visited {
        level[node] = usize::MAX;
    }
    (depth + 1, last_level)
}
