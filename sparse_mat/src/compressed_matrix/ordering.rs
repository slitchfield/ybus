//! Fill-reducing orderings: permutations that make a matrix cheaper to factor.
//!
//! An ordering depends only on the pattern of `A + Aᵀ` (which positions are stored, stored zeros
//! included, ignoring the diagonal), so it is the same for a matrix in either orientation and
//! for any values. It is returned as a permutation in the convention of
//! [`CompressedMatrix::permute_symmetric`]: entry `(i, j)` of the reordered matrix is entry
//! `(perm[i], perm[j])` of the original.

use std::collections::{HashMap, VecDeque};

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
    /// Approximate minimum degree: repeatedly eliminates a node of (approximately) smallest
    /// degree in the graph of what remains, which keeps the factor's fill-in small. Nodes with
    /// many more neighbours than the rest (dense rows) are ordered last.
    Amd,
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
            Ordering::Amd => approximate_minimum_degree(&Adjacency::of(self)),
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

/// The role of a node during minimum degree elimination.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Node {
    /// Not yet eliminated, and the principal variable of its supervariable
    Variable,
    /// Merged into another variable's supervariable: eliminated with it
    Merged,
    /// Eliminated: a clique over its variables, held as the list of them
    Element,
    /// An element whose variables all belong to a later element, so it no longer matters
    Absorbed,
    /// A dense row, set aside to be ordered last
    Dense,
}

/// Lists of the nodes in each degree, for finding a node of smallest degree in O(1).
struct DegreeLists {
    head: Vec<Option<usize>>,
    next: Vec<Option<usize>>,
    previous: Vec<Option<usize>>,
}

impl DegreeLists {
    fn new(n: usize) -> DegreeLists {
        DegreeLists {
            head: vec![None; n + 1],
            next: vec![None; n],
            previous: vec![None; n],
        }
    }

    fn insert(&mut self, node: usize, degree: usize) {
        self.next[node] = self.head[degree];
        self.previous[node] = None;
        if let Some(first) = self.head[degree] {
            self.previous[first] = Some(node);
        }
        self.head[degree] = Some(node);
    }

    fn remove(&mut self, node: usize, degree: usize) {
        match self.previous[node] {
            Some(before) => self.next[before] = self.next[node],
            None => self.head[degree] = self.next[node],
        }
        if let Some(after) = self.next[node] {
            self.previous[after] = self.previous[node];
        }
    }
}

/// Approximate minimum degree ordering, after P. R. Amestoy, T. A. Davis and I. S. Duff, "An
/// approximate minimum degree ordering algorithm" (1996), and T. A. Davis, *Direct Methods for
/// Sparse Linear Systems* (2006), ch. 7.
///
/// Elimination works on the quotient graph: an eliminated node becomes an *element*, the list
/// of the variables its elimination joined into a clique, so fill-in is never stored. Each
/// variable `i` keeps its adjacent elements `E_i` and variables `A_i`. Choosing a variable `p`
/// of least approximate degree forms the element `L_p` (its variables and those of its
/// elements), absorbs `p`'s elements, and updates every variable in `L_p`:
/// - its degree is bounded above by the smallest of: the remaining weight; its old degree plus
///   `|L_p|`; and `|A_i| + |L_p| + Σ |L_e \ L_p|` over its other elements, where each
///   `|L_e \ L_p|` comes from one pass over `L_p` (the approximate degree)
/// - an element with `L_e ⊆ L_p` is absorbed into `p` (aggressive absorption)
/// - variables with identical adjacency merge into one supervariable, found by hashing, and are
///   then eliminated together
///
/// Nodes with more than `max(16, 10 √n)` neighbours are set aside at the start and ordered
/// last, as the book recommends. Ties go to the node most recently given its degree, so the
/// result is deterministic.
fn approximate_minimum_degree(graph: &Adjacency) -> Vec<usize> {
    let n = graph.len();
    let mut order = Vec::with_capacity(n);
    if n == 0 {
        return order;
    }
    let dense_threshold = (10.0 * (n as f64).sqrt()).max(16.0) as usize;
    let mut kind: Vec<Node> = (0..n)
        .map(|i| {
            if graph.degree(i) > dense_threshold && n > 2 {
                Node::Dense
            } else {
                Node::Variable
            }
        })
        .collect();

    // Weights (supervariable sizes), adjacent variables and elements, element lists and sizes
    let mut weight = vec![1usize; n];
    let mut variables: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            graph
                .neighbours(i)
                .iter()
                .copied()
                .filter(|&j| kind[j] != Node::Dense)
                .collect()
        })
        .collect();
    let mut elements: Vec<Vec<usize>> = vec![vec![]; n];
    let mut element_list: Vec<Vec<usize>> = vec![vec![]; n];
    let mut element_size = vec![0usize; n];
    let mut members: Vec<Vec<usize>> = vec![vec![]; n];

    let mut degree: Vec<usize> = (0..n).map(|i| variables[i].len()).collect();
    let mut lists = DegreeLists::new(n);
    let mut remaining = 0;
    for i in (0..n).rev() {
        if kind[i] == Node::Variable {
            lists.insert(i, degree[i]);
            remaining += 1;
        }
    }

    // Per-pivot scratch: marks for L_p, and |L_e \ L_p| for the elements touched
    let mut in_pivot = vec![usize::MAX; n];
    let mut outside = vec![0usize; n];
    let mut outside_stamp = vec![usize::MAX; n];
    let mut external = vec![0usize; n];
    let mut min_degree = 0;
    let mut step = 0;

    while remaining > 0 {
        // A variable of least approximate degree
        while lists.head[min_degree].is_none() {
            min_degree += 1;
        }
        let p = lists.head[min_degree].expect("a degree list is non-empty");
        lists.remove(p, min_degree);
        step += 1;

        // Form L_p from p's variables and the lists of p's elements, which it absorbs
        in_pivot[p] = step;
        let mut pivot_list = vec![];
        for &j in &variables[p] {
            if kind[j] == Node::Variable && in_pivot[j] != step {
                in_pivot[j] = step;
                pivot_list.push(j);
            }
        }
        for e in std::mem::take(&mut elements[p]) {
            if kind[e] != Node::Element {
                continue;
            }
            for &j in &element_list[e] {
                if kind[j] == Node::Variable && in_pivot[j] != step {
                    in_pivot[j] = step;
                    pivot_list.push(j);
                }
            }
            kind[e] = Node::Absorbed;
            element_list[e] = vec![];
        }
        variables[p] = vec![];
        kind[p] = Node::Element;
        remaining -= weight[p];
        order.push(p);
        order.append(&mut members[p]);

        // |L_e \ L_p| for each element adjacent to L_p, by subtracting L_p's weights from |L_e|
        for &i in &pivot_list {
            lists.remove(i, degree[i]);
            for &e in &elements[i] {
                if kind[e] == Node::Element {
                    if outside_stamp[e] != step {
                        outside_stamp[e] = step;
                        outside[e] = element_size[e];
                    }
                    outside[e] -= weight[i];
                }
            }
        }

        // Update each variable of L_p: its elements gain p and lose those absorbed, its
        //   variables lose L_p, and the parts of its degree bound that don't depend on merges
        let mut by_hash: HashMap<usize, Vec<usize>> = HashMap::new();
        for &i in &pivot_list {
            let mut kept = vec![p];
            let mut from_elements = 0;
            for &e in &elements[i] {
                if kind[e] != Node::Element {
                    continue;
                }
                if outside[e] == 0 {
                    // Aggressive absorption: every variable of L_e is in L_p
                    kind[e] = Node::Absorbed;
                    element_list[e] = vec![];
                } else {
                    from_elements += outside[e];
                    kept.push(e);
                }
            }
            let neighbours: Vec<usize> = variables[i]
                .iter()
                .copied()
                .filter(|&j| kind[j] == Node::Variable && in_pivot[j] != step)
                .collect();
            let from_variables: usize = neighbours.iter().map(|&j| weight[j]).sum();
            external[i] = from_variables + from_elements;
            kept.sort_unstable();
            let mut neighbours = neighbours;
            neighbours.sort_unstable();
            let hash = kept
                .iter()
                .chain(&neighbours)
                .fold(0usize, |h, &x| h.wrapping_mul(31).wrapping_add(x));
            elements[i] = kept;
            variables[i] = neighbours;
            by_hash.entry(hash).or_default().push(i);
        }

        // Supervariables: variables of L_p with the same elements and variables are
        //   indistinguishable from here on, so each merges into the first of its kind
        for group in by_hash.values().filter(|group| group.len() > 1) {
            for (a, &i) in group.iter().enumerate() {
                if kind[i] != Node::Variable {
                    continue;
                }
                for &j in &group[a + 1..] {
                    if kind[j] == Node::Variable
                        && elements[i] == elements[j]
                        && variables[i] == variables[j]
                    {
                        weight[i] += weight[j];
                        weight[j] = 0;
                        kind[j] = Node::Merged;
                        let merged = std::mem::take(&mut members[j]);
                        members[i].push(j);
                        members[i].extend(merged);
                        elements[j] = vec![];
                        variables[j] = vec![];
                    }
                }
            }
        }

        // The new element, over the principal variables left, and their approximate degrees
        pivot_list.retain(|&i| kind[i] == Node::Variable);
        let size: usize = pivot_list.iter().map(|&i| weight[i]).sum();
        element_size[p] = size;
        for &i in &pivot_list {
            let others = size - weight[i];
            let bound = (remaining - weight[i])
                .min(degree[i] + others)
                .min(external[i] + others);
            degree[i] = bound;
            lists.insert(i, bound);
            min_degree = min_degree.min(bound);
        }
        element_list[p] = pivot_list;
    }

    // Dense rows last
    order.extend((0..n).filter(|&i| kind[i] == Node::Dense));
    debug_assert_eq!(order.len(), n);
    order
}
