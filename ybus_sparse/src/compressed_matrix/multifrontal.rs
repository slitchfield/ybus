//! Sparse symmetric (Hermitian) indefinite `L D Lᴴ` by the multifrontal method, with threshold
//! Bunch–Kaufman pivoting and delayed pivots: [`Ldlt`].

use std::ops::Range;

use num_complex::ComplexFloat;
use num_traits::NumCast;

use super::bunch_kaufman::{factor_front, DBlock};
use super::cholesky::SymbolicCholesky;
use super::factorization::sealed::Substitute;
use super::factorization::Factorization;
use super::{CompressedMatrix, CsrMatrix, Inertia, Ordering};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// Options for [`SymbolicCholesky::factor_ldlt`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LdltOptions {
    /// The threshold `u` of the pivot test, in `(0, 0.5]`: every entry of `L` is at most `1/u`
    /// in magnitude. Larger is more stable and delays more pivots; smaller keeps more of the
    /// ordering's sparsity. The default, 0.01, is MA57's.
    pub threshold: f64,
}

impl Default for LdltOptions {
    fn default() -> LdltOptions {
        LdltOptions { threshold: 0.01 }
    }
}

/// One front's share of the factorization
#[derive(Clone, Debug)]
struct Front<T> {
    /// The front's variables (rows of `P A Pᵀ`) in their final order: the `eliminated` pivots
    /// first, then the rest, which went to the parent
    vars: Vec<usize>,
    eliminated: usize,
    /// `L`'s columns for the pivots, over all of `vars`, column-major (`vars.len()` rows each),
    /// below the diagonal (a 2×2 block's lower-left entry is `L`'s zero); the diagonal and what
    /// lies above it aren't used
    l: Vec<T>,
    blocks: Vec<DBlock<T>>,
}

/// A sparse symmetric (Hermitian) indefinite factorization `P A Pᵀ = Q L D Lᴴ Qᵀ` with threshold
/// Bunch–Kaufman pivoting, from [`SymbolicCholesky::factor_ldlt`] or
/// [`CsrMatrix::ldlt`]. `P` is the analysis's fill-reducing ordering and `Q` the
/// pivoting's further permutation, `L` is unit lower triangular and `D` block diagonal with 1×1
/// and 2×2 blocks. It factors any nonsingular symmetric matrix, with no regularization, and gives
/// its inertia.
///
/// **The multifrontal method** (Duff and Reid, 1983; Liu, 1992), on the Cholesky analysis's
/// elimination tree:
///
/// - **Fronts:** the columns are grouped into fundamental supernodes, chains of columns whose
///   factor columns nest, each with one dense front: its own variables, then the rows its
///   columns reach below them.
/// - **Assembly:** a front gathers `A`'s entries in its own columns, and adds its children's
///   contribution blocks, their Schur complements over the variables they didn't eliminate
///   (extend-add).
/// - **Pivots:** the front's fully summed variables, its own and any its children delayed, are
///   eliminated as far as threshold Bunch–Kaufman pivoting allows (see [`DenseLdlt`]). A variable
///   that fails the test is delayed: it goes up with the contribution block and is fully summed
///   in the parent, whose larger front may offer a pivot. At a root every variable is fully
///   summed, and a pivot always exists unless the matrix is singular.
/// - **Cost:** without delays, the fronts' variables are exactly the Cholesky factor's pattern,
///   so the analysis predicts the work; each delay enlarges the fronts above it. The counts are
///   reported ([`delayed`](Self::delayed), [`l_entries`](Self::l_entries)).
///
/// **References:**
///
/// - I. S. Duff and J. K. Reid, "The multifrontal solution of indefinite sparse symmetric linear
///   equations," *ACM Trans. Math. Softw.* 9(3):302–325, 1983, doi:10.1145/356044.356047.
/// - J. W. H. Liu, "The multifrontal method for sparse matrix solution: theory and practice,"
///   *SIAM Rev.* 34(1):82–109, 1992, doi:10.1137/1034004.
///
/// [`DenseLdlt`]: super::experimental::DenseLdlt
#[derive(Clone, Debug)]
pub struct Ldlt<T: Scalar = f64> {
    n: usize,
    /// Row `i` of `P A Pᵀ` is row `perm[i]` of `A`
    perm: Vec<usize>,
    /// In elimination order
    fronts: Vec<Front<T>>,
    delayed: usize,
    /// The analysis's front tree, for the parallel solves
    #[cfg(feature = "rayon")]
    tree: std::sync::Arc<SolveTree>,
}

/// The fundamental supernodes: `(start, end)` column ranges in increasing order, each with the
/// rows below `end` that its first column reaches
fn supernodes(symbolic: &SymbolicCholesky) -> Vec<(usize, usize, Vec<usize>)> {
    let n = symbolic.dim();
    // L's columns: row i's columns are l_idx[l_ptr[i]..l_ptr[i + 1]], ending with i
    let mut counts = vec![0usize; n];
    for &j in &symbolic.l_idx {
        counts[j] += 1;
    }
    let mut col_ptr = vec![0; n + 1];
    for j in 0..n {
        col_ptr[j + 1] = col_ptr[j] + counts[j];
    }
    let mut fill = col_ptr.clone();
    let mut col_rows = vec![0; symbolic.l_idx.len()];
    for i in 0..n {
        for &j in &symbolic.l_idx[symbolic.l_ptr[i]..symbolic.l_ptr[i + 1]] {
            col_rows[fill[j]] = i;
            fill[j] += 1;
        }
    }
    let parent = symbolic.elimination_tree();
    let mut children = vec![0usize; n];
    for p in parent.iter().flatten() {
        children[*p] += 1;
    }
    let mut result = Vec::new();
    let mut start = 0;
    for j in 0..n {
        let continues = j + 1 < n
            && parent[j] == Some(j + 1)
            && children[j + 1] == 1
            && counts[j] == counts[j + 1] + 1;
        if !continues {
            let end = j + 1;
            // Rows increase within a column (rows are visited in order)
            let below = col_rows[col_ptr[start]..col_ptr[start + 1]]
                .iter()
                .copied()
                .filter(|&r| r >= end)
                .collect();
            result.push((start, end, below));
            start = end;
        }
    }
    result
}

impl SymbolicCholesky {
    /// Factors `matrix`, which must have the analyzed pattern and be symmetric (for a complex
    /// element type, Hermitian), as `P A Pᵀ = Q L D Lᴴ Qᵀ` with threshold Bunch–Kaufman pivoting
    /// on the multifrontal method (see [`Ldlt`]). The analysis's ordering and
    /// elimination tree are used as they are; pivoting only reorders within and between fronts.
    ///
    /// This builds a [`SymbolicLdlt`] on each call. To factor several matrices with the
    /// same pattern, build one once and use [`SymbolicLdlt::factor`], whose results are
    /// the same bit for bit.
    ///
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, `SparseError::NotSymmetric` if `A ≠ Aᴴ`, and `SparseError::SingularMatrix`
    /// naming a row of `A` where no pivot could be found, the matrix being singular to working
    /// precision. Panics if `options.threshold` is outside `(0, 0.5]`.
    pub fn factor_ldlt<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        matrix: &CsrMatrix<T>,
        options: LdltOptions,
    ) -> Result<Ldlt<T>, SparseError> {
        SymbolicLdlt::new(self.clone()).factor(matrix, options)
    }
}

/// One fundamental supernode: its columns, the rows below them that its first column reaches,
/// and the supernode its parent column belongs to
#[derive(Clone, Debug)]
struct Supernode {
    start: usize,
    end: usize,
    below: Vec<usize>,
    parent: Option<usize>,
}

/// The symbolic analysis for [`Ldlt`]: the Cholesky analysis ([`SymbolicCholesky`]:
/// the ordering, the elimination tree and L's pattern), with what the multifrontal factorization
/// needs of the pattern, worked out once rather than at every factorization:
///
/// - **the fundamental supernodes,** each front's columns and rows, and its parent;
/// - **the gather map** from `A`'s stored values to `P A Pᵀ`'s, so a factorization gathers the
///   values instead of permuting the matrix;
/// - **each stored entry's mirror,** the position of `(j, i)` for `(i, j)`, so the symmetry
///   check is one pass over the values, accepting and refusing exactly the matrices
///   [`is_hermitian`](CsrMatrix::is_hermitian) does.
///
/// Together these were 7.5–11.7% of the factorization's time when they were redone each time
/// (`docs/benchmarks/ipm-profile-2026-11-08.md`). Build it with
/// [`CsrMatrix::analyze_ldlt`] or from a Cholesky analysis with
/// [`SymbolicLdlt::new`], and factor with [`SymbolicLdlt::factor`].
#[derive(Clone, Debug)]
pub struct SymbolicLdlt {
    cholesky: SymbolicCholesky,
    nodes: Vec<Supernode>,
    /// `P A Pᵀ`'s pattern in CSR form
    c_ptr: Vec<usize>,
    c_idx: Vec<usize>,
    /// For each of `P A Pᵀ`'s stored entries, the position of its value among `A`'s
    source: Vec<usize>,
    /// For each of `A`'s stored entries `(i, j)`, the position of `(j, i)`, or `usize::MAX` if
    /// `A`'s pattern doesn't store it
    mirror: Vec<usize>,
    /// Each front's children, in increasing index, the order the factorization assembles their
    /// contribution blocks in (read by the parallel factorization)
    #[cfg_attr(not(feature = "rayon"), allow(dead_code))]
    children: Vec<Vec<usize>>,
    /// Each front's subtree's factorization work (see [`FrontTree`]): its own and its
    /// descendants'
    #[cfg_attr(not(feature = "rayon"), allow(dead_code))]
    subtree_work: Vec<u64>,
    /// Each front's subtree's first front: the subtree is the run from it to the front, as the
    /// fronts follow a postorder
    #[cfg_attr(not(feature = "rayon"), allow(dead_code))]
    first: Vec<usize>,
    /// The front tree as the parallel solves walk it, shared by every factorization
    #[cfg(feature = "rayon")]
    solve_tree: std::sync::Arc<SolveTree>,
}

/// The front tree as the parallel solves walk it (feature `rayon`): each front's own columns,
/// and each subtree's fronts, children and one sweep's work. A subtree is a run of fronts, from
/// its first to its root, and of columns, from its first front's first to its root's last.
#[cfg(feature = "rayon")]
#[derive(Debug)]
struct SolveTree {
    /// Front `k`'s own columns are `start[k]..start[k + 1]`
    start: Vec<usize>,
    first: Vec<usize>,
    children: Vec<Vec<usize>>,
    roots: Vec<usize>,
    /// Each subtree's work for one sweep, its fronts' `FrontTree` solve work
    work: Vec<u64>,
    /// Whether every subtree's fronts and columns are runs, with its children's subtrees side
    ///   by side in index order, as a postorder gives. Without the default amalgamation's
    ///   postorder they may not be, and the solves stay serial
    contiguous: bool,
}

#[cfg(feature = "rayon")]
impl SolveTree {
    fn new(nodes: &[Supernode], children: &[Vec<usize>], first: &[usize]) -> SolveTree {
        let mut start: Vec<usize> = nodes.iter().map(|node| node.start).collect();
        start.push(nodes.last().map_or(0, |node| node.end));
        let mut work: Vec<u64> = nodes.iter().map(|node| front_work(node).1).collect();
        for (k, node) in nodes.iter().enumerate() {
            if let Some(parent) = node.parent {
                work[parent] += work[k];
            }
        }
        let roots: Vec<usize> = (0..nodes.len())
            .filter(|&k| nodes[k].parent.is_none())
            .collect();
        // Side by side: each of a run's subtrees starting where the last ended, and ending
        //   where the next begins
        let side_by_side = |run: &[usize], from: usize, to: usize| {
            let mut at = from;
            for &k in run {
                if first[k] != at {
                    return false;
                }
                at = k + 1;
            }
            at == to
        };
        let contiguous = nodes.windows(2).all(|pair| pair[0].end == pair[1].start)
            && side_by_side(&roots, 0, nodes.len())
            && (0..nodes.len()).all(|k| side_by_side(&children[k], first[k], k));
        SolveTree {
            start,
            first: first.to_vec(),
            children: children.to_vec(),
            roots,
            work,
            contiguous,
        }
    }

    /// Subtree `s`'s columns
    fn columns(&self, s: usize) -> Range<usize> {
        self.start[self.first[s]]..self.start[s + 1]
    }

    /// The subtrees `kids` (side by side, in index order) as the parallel solves' tasks, in
    ///   order: each subtree whose work passes `cutoff` alone, for its own children to split
    ///   further, and runs of the others, each run's work passing `cutoff` but its last
    fn groups(&self, kids: &[usize], cutoff: u64) -> Vec<Group> {
        let mut groups = Vec::new();
        let mut run: Option<(usize, usize, u64)> = None;
        for &k in kids {
            if self.work[k] >= cutoff && !self.children[k].is_empty() {
                if let Some((from, to, _)) = run.take() {
                    groups.push(Group::Run(from..to));
                }
                groups.push(Group::Subtree(k));
                continue;
            }
            let (from, _, work) = run.take().unwrap_or((self.first[k], k, 0));
            let work = work + self.work[k];
            if work >= cutoff {
                groups.push(Group::Run(from..k + 1));
            } else {
                run = Some((from, k + 1, work));
            }
        }
        if let Some((from, to, _)) = run {
            groups.push(Group::Run(from..to));
        }
        groups
    }
}

/// One of the parallel solves' tasks: a subtree to split further, or a run of whole subtrees'
/// fronts, solved serially
#[cfg(feature = "rayon")]
enum Group {
    Subtree(usize),
    Run(Range<usize>),
}

/// The entries of `y` above a subtree, its ancestors' own columns, which a backward sweep
/// reads once they're final: the nearest ancestor's first, then up the tree
#[cfg(feature = "rayon")]
struct Above<'a, T> {
    lo: usize,
    values: &'a [T],
    up: Option<&'a Above<'a, T>>,
}

#[cfg(feature = "rayon")]
impl<T: Copy> Above<'_, T> {
    fn get(&self, row: usize) -> T {
        let mut at = self;
        loop {
            if let Some(&value) = row.checked_sub(at.lo).and_then(|i| at.values.get(i)) {
                return value;
            }
            at = at.up.expect("a row outside a subtree is an ancestor's");
        }
    }
}

/// A front's factorization work and one triangular sweep's, as [`FrontTree`] counts them
fn front_work(node: &Supernode) -> (u64, u64) {
    let c = (node.end - node.start) as u64;
    let m = c + node.below.len() as u64;
    let factor = (0..c).map(|j| (m - j - 1) * (m - j) / 2).sum();
    let solve = c * c.saturating_sub(1) / 2 + c * node.below.len() as u64;
    (factor, solve)
}

/// Each front's children, in increasing index, each subtree's factorization work, and each
///   subtree's first front: a parent's index is above its children's, as the supernodes follow
///   the postordered columns
#[allow(clippy::type_complexity)]
fn tree_of(nodes: &[Supernode]) -> (Vec<Vec<usize>>, Vec<u64>, Vec<usize>) {
    let mut children = vec![Vec::new(); nodes.len()];
    let mut work: Vec<u64> = nodes.iter().map(|node| front_work(node).0).collect();
    let mut first: Vec<usize> = (0..nodes.len()).collect();
    for (k, node) in nodes.iter().enumerate() {
        if let Some(parent) = node.parent {
            debug_assert!(parent > k, "a parent front after its children");
            children[parent].push(k);
            work[parent] += work[k];
            first[parent] = first[parent].min(first[k]);
        }
    }
    (children, work, first)
}

/// The shape of a multifrontal factorization's front tree, from its symbolic analysis
/// ([`SymbolicLdlt::front_tree`]): how much work the fronts hold, and how much of it lies on the
/// tree's critical path, the heaviest chain from a leaf to a root.
///
/// The work counts multiply-adds on the symbolic fronts. For a front with `c` pivot columns
/// and `r` rows below them, so `m = c + r` rows in all:
/// - **the factorization** is `Σ_{j<c} (m − j − 1)(m − j) / 2`, the trailing lower triangle's
///   update after each pivot
/// - **one triangular sweep** is `c (c − 1) / 2 + c r`, the front's columns of `L`
///
/// A front depends on its children, so `work / critical_work` bounds the speedup any
/// parallelism over independent subtrees can give, and `solve_work / solve_critical_work` a
/// level-scheduled solve's. Pivots delayed at factorization time enlarge parents beyond these
/// symbolic sizes, so the bounds are for the factorization without delays.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FrontTree {
    /// The fronts
    pub fronts: usize,
    /// The roots, one per connected block of the matrix's graph
    pub roots: usize,
    /// The most fronts on a path from a leaf to a root
    pub height: usize,
    /// The factorization's work over all fronts
    pub work: u64,
    /// The factorization's work on the critical path
    pub critical_work: u64,
    /// One triangular sweep's work over all fronts
    pub solve_work: u64,
    /// One triangular sweep's work on its critical path
    pub solve_critical_work: u64,
}

/// How [`SymbolicLdlt`] merges supernodes into larger fronts (relaxed amalgamation,
/// Ashcraft and Grimes, "The influence of relaxed supernode partitions on the multifrontal
/// method," *ACM Trans. Math. Softw.* 15(4):291–309, 1989, doi:10.1145/76909.76910).
///
/// A child merges into its parent when its columns come just before the parent's (as a parent's
/// last child's do in a postorder), and the merged front, with `c` columns, of which a fraction
/// `z` of its dense lower trapezoid's entries would be explicit zeros, passes the rule:
///
/// `c < columns[0]`, or `c < columns[1]` and `z < zeros[0]`, or `c < columns[2]` and
/// `z < zeros[1]`, or `z < zeros[2]`.
///
/// Fewer, larger fronts cost fewer assemblies, allocations and kernel calls, for some
/// arithmetic on the zeros. The pivots are chosen within fronts, so merging changes the pivot
/// sequence, and with it the factorization's rounding.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Amalgamation {
    /// The column thresholds
    pub columns: [usize; 3],
    /// The zero fractions allowed below each
    pub zeros: [f64; 3],
}

impl Amalgamation {
    /// No merging: the fronts are the fundamental supernodes on the analysis's own ordering.
    pub const NONE: Amalgamation = Amalgamation {
        columns: [0; 3],
        zeros: [0.0; 3],
    };

    fn merges(&self, columns: usize, zero_fraction: f64) -> bool {
        columns < self.columns[0]
            || (columns < self.columns[1] && zero_fraction < self.zeros[0])
            || (columns < self.columns[2] && zero_fraction < self.zeros[1])
            || zero_fraction < self.zeros[2]
    }
}

impl Default for Amalgamation {
    /// CHOLMOD's defaults: merge into fronts of fewer than 4 columns always, of fewer than 16
    /// with under 80% zeros, of fewer than 48 with under 10%, and any with under 5%.
    fn default() -> Amalgamation {
        Amalgamation {
            columns: [4, 16, 48],
            zeros: [0.8, 0.1, 0.05],
        }
    }
}

/// `cholesky` with its ordering postordered: the same elimination tree, renumbered so that each
/// subtree's columns are contiguous and each parent comes right after its last child
fn postordered(cholesky: SymbolicCholesky) -> SymbolicCholesky {
    let order = cholesky.postorder();
    if order.iter().enumerate().all(|(k, &node)| k == node) {
        return cholesky;
    }
    let perm: Vec<usize> = order.iter().map(|&node| cholesky.perm[node]).collect();
    let n = cholesky.dim();
    let pattern: CsrMatrix<u8> = CompressedMatrix::from_parts_unchecked(
        (n, n),
        cholesky.a_ptr.clone(),
        cholesky.a_idx.clone(),
        vec![0; cholesky.a_idx.len()],
    );
    pattern.analyze_cholesky_permuted(perm)
}

/// The supernodes merged by `amalgamation`'s rule: children are taken from the last, so a chain
/// of small supernodes merges into one, and each merge only into a parent whose columns its own
/// come just before
fn amalgamate(nodes: Vec<Supernode>, amalgamation: Amalgamation) -> Vec<Supernode> {
    if amalgamation == Amalgamation::NONE {
        return nodes;
    }
    let count = nodes.len();
    // The front each supernode has merged into (itself if none), its start, and the entries of
    //   L its columns hold, explicit zeros not counted
    let mut into: Vec<usize> = (0..count).collect();
    let mut start: Vec<usize> = nodes.iter().map(|node| node.start).collect();
    let trapezoid = |columns: usize, below: usize| columns * (columns + 1) / 2 + columns * below;
    // A fundamental supernode's L block is dense, so its entries are its trapezoid's
    let mut entries: Vec<usize> = nodes
        .iter()
        .map(|node| trapezoid(node.end - node.start, node.below.len()))
        .collect();
    fn find(into: &mut [usize], mut s: usize) -> usize {
        while into[s] != s {
            into[s] = into[into[s]];
            s = into[s];
        }
        s
    }
    for s in (0..count).rev() {
        let Some(parent) = nodes[s].parent else {
            continue;
        };
        let p = find(&mut into, parent);
        if nodes[s].end != start[p] {
            continue;
        }
        let columns = nodes[p].end - start[s];
        let dense = trapezoid(columns, nodes[p].below.len());
        let actual = entries[s] + entries[p];
        let zero_fraction = (dense - actual) as f64 / dense as f64;
        if amalgamation.merges(columns, zero_fraction) {
            into[s] = p;
            start[p] = start[s];
            // The true entries, not the dense trapezoid: counting the merged front's explicit
            //   zeros as entries would hide them from every later merge's zero fraction, and
            //   merges would cascade into dense fronts
            entries[p] = actual;
        }
    }
    // The fronts left, in order of their columns, with their parents renumbered
    let mut number = vec![usize::MAX; count];
    let mut kept = Vec::new();
    for (s, slot) in number.iter_mut().enumerate() {
        if find(&mut into, s) == s {
            *slot = kept.len();
            kept.push(s);
        }
    }
    let mut merged: Vec<Supernode> = kept
        .iter()
        .map(|&s| Supernode {
            start: start[s],
            end: nodes[s].end,
            below: nodes[s].below.clone(),
            parent: nodes[s].parent,
        })
        .collect();
    for node in &mut merged {
        node.parent = node.parent.map(|p| number[find(&mut into, p)]);
    }
    // The fronts' column ranges partition the columns, and their ends haven't moved, so they're
    //   still in order of their columns, and a parent still follows its children
    debug_assert!(merged.windows(2).all(|pair| pair[0].end == pair[1].start));
    merged
}

impl SymbolicLdlt {
    /// The analysis for the multifrontal LDLᵀ on `cholesky`'s ordering, elimination tree and
    /// pattern, with the default relaxed amalgamation ([`Amalgamation::default`]): the
    /// ordering is postordered and small supernodes merged.
    pub fn new(cholesky: SymbolicCholesky) -> SymbolicLdlt {
        SymbolicLdlt::with_amalgamation(cholesky, Amalgamation::default())
    }

    /// The analysis with the given amalgamation. With [`Amalgamation::NONE`] the fronts are the
    /// fundamental supernodes on `cholesky`'s own ordering. Otherwise the ordering is first
    /// postordered (the same elimination tree, renumbered so each subtree's columns are
    /// contiguous and a parent follows its last child), so that a child can merge into its
    /// parent, and the supernodes are merged by the rule `amalgamation` sets.
    pub fn with_amalgamation(
        cholesky: SymbolicCholesky,
        amalgamation: Amalgamation,
    ) -> SymbolicLdlt {
        let cholesky = if amalgamation == Amalgamation::NONE {
            cholesky
        } else {
            postordered(cholesky)
        };
        let n = cholesky.dim();
        let raw = supernodes(&cholesky);
        // Each column's supernode, to find a supernode's parent
        let mut node_of = vec![0; n];
        for (s, (start, end, _)) in raw.iter().enumerate() {
            node_of[*start..*end].fill(s);
        }
        let parent = cholesky.elimination_tree();
        let fundamental: Vec<Supernode> = raw
            .into_iter()
            .map(|(start, end, below)| Supernode {
                start,
                end,
                below,
                parent: parent[end - 1].map(|p| node_of[p]),
            })
            .collect();
        let nodes = amalgamate(fundamental, amalgamation);
        // The gather map: the matrix whose values are their own positions, permuted as a
        //   factorization would permute A
        let (a_ptr, a_idx) = (&cholesky.a_ptr, &cholesky.a_idx);
        let positions: CsrMatrix<usize> = CompressedMatrix::from_parts_unchecked(
            (n, n),
            a_ptr.clone(),
            a_idx.clone(),
            (0..a_idx.len()).collect(),
        );
        let (_, c_ptr, c_idx, source) = positions
            .permute_symmetric(&cholesky.perm)
            .expect("the analysis's ordering is a permutation of its order")
            .into_raw_parts();
        let mut mirror = vec![usize::MAX; a_idx.len()];
        for i in 0..n {
            for p in a_ptr[i]..a_ptr[i + 1] {
                let j = a_idx[p];
                if let Ok(k) = a_idx[a_ptr[j]..a_ptr[j + 1]].binary_search(&i) {
                    mirror[p] = a_ptr[j] + k;
                }
            }
        }
        let (children, subtree_work, first) = tree_of(&nodes);
        SymbolicLdlt {
            cholesky,
            c_ptr,
            c_idx,
            source,
            mirror,
            #[cfg(feature = "rayon")]
            solve_tree: std::sync::Arc::new(SolveTree::new(&nodes, &children, &first)),
            nodes,
            subtree_work,
            first,
            children,
        }
    }

    /// The front tree's shape: its fronts' work and its critical path (see [`FrontTree`]).
    pub fn front_tree(&self) -> FrontTree {
        let n = self.nodes.len();
        let front = front_work;
        // Each front's heaviest chain down to a leaf, and its depth: a parent's index is
        //   above its children's, as the supernodes follow the postordered columns
        let mut chain = vec![(0u64, 0u64, 0usize); n];
        let mut tree = FrontTree {
            fronts: n,
            ..FrontTree::default()
        };
        for (k, node) in self.nodes.iter().enumerate() {
            let (factor, solve) = front(node);
            let below = chain[k];
            let ours = (below.0 + factor, below.1 + solve, below.2 + 1);
            tree.work += factor;
            tree.solve_work += solve;
            match node.parent {
                Some(parent) => {
                    debug_assert!(parent > k, "a parent front after its children");
                    let theirs = &mut chain[parent];
                    *theirs = (
                        theirs.0.max(ours.0),
                        theirs.1.max(ours.1),
                        theirs.2.max(ours.2),
                    );
                }
                None => {
                    tree.roots += 1;
                    tree.critical_work = tree.critical_work.max(ours.0);
                    tree.solve_critical_work = tree.solve_critical_work.max(ours.1);
                    tree.height = tree.height.max(ours.2);
                }
            }
        }
        tree
    }

    /// The Cholesky analysis it was built on.
    pub fn cholesky(&self) -> &SymbolicCholesky {
        &self.cholesky
    }

    /// Whether the values `data`, on the analyzed pattern, make a Hermitian matrix: each stored
    /// value is its mirror's conjugate, or zero where the mirror isn't stored, as
    /// [`is_hermitian`](CsrMatrix::is_hermitian) judges
    fn is_hermitian<T: Scalar + ComplexFloat>(&self, data: &[T]) -> bool {
        data.iter().zip(&self.mirror).all(|(&value, &mirror)| {
            if mirror == usize::MAX {
                value == T::zero()
            } else {
                value == data[mirror].conj()
            }
        })
    }

    /// Factors `matrix`, which must have the analyzed pattern and be symmetric (for a complex
    /// element type, Hermitian), as `P A Pᵀ = Q L D Lᴴ Qᵀ` with threshold Bunch–Kaufman pivoting
    /// on the multifrontal method (see [`Ldlt`]), with the analysis's supernodes and
    /// maps.
    ///
    /// **With the `rayon` feature,** independent subtrees of the front tree factor in parallel,
    /// those of 2²⁰ multiply-adds or more. Each front still does exactly the serial arithmetic:
    /// it assembles its children's contribution blocks in their index order. So the results
    /// are bit for bit the serial factorization's, at any thread count. (The dense kernel can
    /// also split a front's trailing update across columns, bit for bit too, but that's off:
    /// it didn't pay on the OPF's fronts.)
    ///
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, `SparseError::NotSymmetric` if `A ≠ Aᴴ`, and `SparseError::SingularMatrix`
    /// naming a row of `A` where no pivot could be found, the matrix being singular to working
    /// precision: the first such front's, in the serial order. Panics if `options.threshold` is
    /// outside `(0, 0.5]`.
    pub fn factor<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        matrix: &CsrMatrix<T>,
        options: LdltOptions,
    ) -> Result<Ldlt<T>, SparseError> {
        self.factor_with(matrix, options, Cutoffs::DEFAULT, None)
    }

    /// [`SymbolicLdlt::factor`], reusing `previous`'s storage, an earlier factorization of a
    /// matrix with this analysis's pattern, for its fronts. The result is bit for bit
    /// `factor`'s; only the allocations change. For a sequence of factorizations, such as an
    /// interior-point method's, each refactoring into the last.
    ///
    /// Returns [`SymbolicLdlt::factor`]'s errors, and panics where it does.
    pub fn refactor<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        matrix: &CsrMatrix<T>,
        options: LdltOptions,
        previous: Ldlt<T>,
    ) -> Result<Ldlt<T>, SparseError> {
        self.factor_with(matrix, options, Cutoffs::DEFAULT, Some(previous))
    }

    /// [`SymbolicLdlt::factor`] with the parallel paths' cutoffs given: the tests force every
    ///   subtree and update onto them with [`Cutoffs::ALWAYS`], and compare with
    ///   [`Cutoffs::NEVER`], the serial factorization
    pub(crate) fn factor_with<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        matrix: &CsrMatrix<T>,
        options: LdltOptions,
        cutoffs: Cutoffs,
        previous: Option<Ldlt<T>>,
    ) -> Result<Ldlt<T>, SparseError> {
        let threshold = options.threshold;
        assert!(
            threshold > 0.0 && threshold <= 0.5,
            "the threshold must be in (0, 0.5], not {threshold}"
        );
        let _: T::Real = NumCast::from(threshold).expect("the threshold is representable");
        let symbolic = &self.cholesky;
        let n = symbolic.dim();
        if matrix.shape != (n, n) || matrix.ptr != symbolic.a_ptr || matrix.idx != symbolic.a_idx {
            return Err(SparseError::PatternMismatch);
        }
        if !self.is_hermitian(&matrix.data) {
            return Err(SparseError::NotSymmetric);
        }
        // P A Pᵀ's values, gathered into its kept pattern
        let c_data: Vec<T> = self.source.iter().map(|&q| matrix.data[q]).collect();
        let context = Context {
            c_data: &c_data,
            threshold,
            node_order: cutoffs.node_order,
        };
        // The previous factorization's fronts, their storage handed out again
        let mut seed = Workspace::new(n);
        for front in previous.into_iter().flat_map(|ldlt| ldlt.fronts) {
            seed.spare_l.push(front.l);
            seed.spare_ordered.push(front.vars);
            seed.spare_dblocks.push(front.blocks);
        }
        let (fronts, delayed) = self.factor_all(&context, cutoffs, seed)?;
        Ok(Ldlt {
            n,
            perm: symbolic.perm.clone(),
            fronts,
            delayed,
            #[cfg(feature = "rayon")]
            tree: self.solve_tree.clone(),
        })
    }

    /// Every front, serially in index order
    #[cfg(not(feature = "rayon"))]
    fn factor_all<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        context: &Context<'_, T>,
        _: Cutoffs,
        mut workspace: Workspace<T>,
    ) -> Result<(Vec<Front<T>>, usize), SparseError> {
        let (fronts, outgoing, delayed) = self
            .factor_range(0..self.nodes.len(), context, &mut workspace)
            .map_err(|(_, error)| error)?;
        debug_assert!(outgoing.is_none(), "every front's parent within the range");
        Ok((fronts, delayed))
    }

    /// Every front: the subtrees whose work passes `cutoffs.subtree_work` in parallel, the
    ///   rest serially, each front from its children's blocks in their index order
    #[cfg(feature = "rayon")]
    fn factor_all<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        context: &Context<'_, T>,
        cutoffs: Cutoffs,
        seed: Workspace<T>,
    ) -> Result<(Vec<Front<T>>, usize), SparseError> {
        use rayon::prelude::*;
        // Serially when there's nothing to gain: no subtree cutoff, a one-thread pool, or a
        //   caller that asked (`serially`, inside a batch of parallel solves). Only the schedule
        //   changes; the results are the same either way
        let serial = cutoffs.subtree_work == u64::MAX
            || (cutoffs != Cutoffs::ALWAYS
                && (rayon::current_num_threads() == 1 || SERIALLY.with(|depth| depth.get() > 0)));
        if serial {
            let mut workspace = seed;
            let (fronts, _, delayed) = self
                .factor_range(0..self.nodes.len(), context, &mut workspace)
                .map_err(|(_, error)| error)?;
            return Ok((fronts, delayed));
        }
        let roots: Vec<usize> = (0..self.nodes.len())
            .filter(|&s| self.nodes[s].parent.is_none())
            .collect();
        // Workspaces shared by the tasks, each taken and given back whole: a front clears and
        //   rewrites everything it uses, so only the allocation is shared
        let pool = std::sync::Mutex::new(vec![seed]);
        let subtrees: Vec<_> = roots
            .par_iter()
            .map(|&root| self.par_subtree(root, context, cutoffs, &pool))
            .collect();
        let mut placed: Vec<Option<Front<T>>> = (0..self.nodes.len()).map(|_| None).collect();
        let mut delayed = 0;
        let mut failure: Option<(usize, SparseError)> = None;
        for subtree in subtrees {
            match subtree {
                Ok((fronts, _, d)) => {
                    delayed += d;
                    for (s, front) in fronts {
                        placed[s] = Some(front);
                    }
                }
                Err(error) => failure = first_failure(failure, error),
            }
        }
        if let Some((_, error)) = failure {
            return Err(error);
        }
        let fronts = placed
            .into_iter()
            .map(|front| front.expect("every front factored"))
            .collect();
        Ok((fronts, delayed))
    }

    /// Front `s`'s subtree: its children's subtrees in parallel where their work passes the
    ///   cutoff, then `s` from their blocks in index order. Returns its fronts with their
    ///   indices, `s`'s contribution block, and the delayed pivots, or the lowest-index failing
    ///   front's error
    #[cfg(feature = "rayon")]
    #[allow(clippy::type_complexity)]
    fn par_subtree<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        s: usize,
        context: &Context<'_, T>,
        cutoffs: Cutoffs,
        pool: &std::sync::Mutex<Vec<Workspace<T>>>,
    ) -> Result<(Vec<(usize, Front<T>)>, Option<Block<T>>, usize), (usize, SparseError)> {
        use rayon::prelude::*;
        let n = self.cholesky.dim();
        let take = || {
            let spare = pool.lock().expect("no task panicked").pop();
            spare.unwrap_or_else(|| Workspace::new(n))
        };
        let give = |workspace| pool.lock().expect("no task panicked").push(workspace);
        if self.subtree_work[s] < cutoffs.subtree_work {
            let mut workspace = take();
            let first = self.first[s];
            let result = self.factor_range(first..s + 1, context, &mut workspace);
            give(workspace);
            let (fronts, outgoing, delayed) = result?;
            let indexed = (first..s + 1).zip(fronts).collect();
            return Ok((indexed, outgoing, delayed));
        }
        let children: Vec<_> = self.children[s]
            .par_iter()
            .map(|&child| self.par_subtree(child, context, cutoffs, pool))
            .collect();
        let mut fronts = Vec::new();
        let mut contributions = Vec::new();
        let mut delayed = 0;
        let mut failure: Option<(usize, SparseError)> = None;
        for child in children {
            match child {
                Ok((child_fronts, block, d)) => {
                    fronts.extend(child_fronts);
                    contributions.extend(block);
                    delayed += d;
                }
                Err(error) => failure = first_failure(failure, error),
            }
        }
        if let Some(failure) = failure {
            return Err(failure);
        }
        let mut workspace = take();
        let result = self.factor_node(s, contributions, context, &mut workspace);
        give(workspace);
        let (front, block, d) = result.map_err(|error| (s, error))?;
        fronts.push((s, front));
        Ok((fronts, block, delayed + d))
    }

    /// The fronts in `range`, a run of whole subtrees in index order, serially, as the serial
    ///   factorization runs them: their fronts in order, the block of the last front whose parent
    ///   lies beyond the range (a subtree's root's), and their delayed pivots
    #[allow(clippy::type_complexity)]
    fn factor_range<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        range: Range<usize>,
        context: &Context<'_, T>,
        workspace: &mut Workspace<T>,
    ) -> Result<(Vec<Front<T>>, Option<Block<T>>, usize), (usize, SparseError)> {
        let first = range.start;
        // Contribution blocks waiting for their parent within the range
        let mut pending: Vec<Vec<Block<T>>> = vec![Vec::new(); range.len()];
        let mut fronts = Vec::with_capacity(range.len());
        let mut outgoing = None;
        let mut delayed = 0;
        for s in range.clone() {
            let contributions = std::mem::take(&mut pending[s - first]);
            let (front, block, d) = self
                .factor_node(s, contributions, context, workspace)
                .map_err(|error| (s, error))?;
            delayed += d;
            if let Some(block) = block {
                let parent = self.nodes[s]
                    .parent
                    .expect("a front with rows left has a parent");
                if range.contains(&parent) {
                    pending[parent - first].push(block);
                } else {
                    outgoing = Some(block);
                }
            }
            fronts.push(front);
        }
        Ok((fronts, outgoing, delayed))
    }

    /// Front `s`, from its children's contribution blocks in their index order: its factor, its
    ///   own contribution block for its parent if it leaves rows, and its delayed pivots
    #[allow(clippy::type_complexity)]
    fn factor_node<T: Scalar + ComplexFloat + Send + Sync>(
        &self,
        s: usize,
        contributions: Vec<Block<T>>,
        context: &Context<'_, T>,
        workspace: &mut Workspace<T>,
    ) -> Result<(Front<T>, Option<Block<T>>, usize), SparseError> {
        let Workspace {
            front,
            scratch,
            vars,
            spare_blocks,
            spare_vars,
            position,
            spare_l,
            spare_ordered,
            spare_dblocks,
            spare_order,
        } = workspace;
        let c_data = context.c_data;
        let node = &self.nodes[s];
        let (start, end, below) = (&node.start, &node.end, &node.below);
        // The front's variables: its own, the children's delayed ones, then the rows below
        vars.clear();
        vars.extend(*start..*end);
        for (child_vars, _) in &contributions {
            for &v in child_vars {
                if v < *start {
                    vars.push(v);
                }
            }
        }
        let fully_summed = vars.len();
        let delayed = fully_summed - (end - start);
        vars.extend_from_slice(below);
        let m = vars.len();
        for (p, &v) in vars.iter().enumerate() {
            position[v] = p;
        }
        front.clear();
        front.resize(m * m, T::zero());
        // A's lower triangle in the own columns, mirrored: row j of C holds column j's
        //   entries, and those at or below the diagonal belong to this front
        for j in *start..*end {
            let pj = position[j];
            let row = self.c_ptr[j]..self.c_ptr[j + 1];
            for (&i, &stored) in self.c_idx[row.clone()].iter().zip(&c_data[row]) {
                if i < j {
                    continue;
                }
                let pi = position[i];
                debug_assert_ne!(pi, usize::MAX, "A's entry within the factor's pattern");
                // c is Hermitian: C(i, j) = conj(C(j, i))
                let value = stored.conj();
                front[pi + pj * m] = front[pi + pj * m] + value;
                if pi != pj {
                    front[pj + pi * m] = front[pj + pi * m] + value.conj();
                }
            }
        }
        // The children's contribution blocks
        for (child_vars, block) in &contributions {
            let k = child_vars.len();
            for (b, &vb) in child_vars.iter().enumerate() {
                let pb = position[vb];
                debug_assert_ne!(pb, usize::MAX, "a child's variable is in its parent");
                for (a, &va) in child_vars.iter().enumerate() {
                    let pa = position[va];
                    front[pa + pb * m] = front[pa + pb * m] + block[a + b * k];
                }
            }
        }
        // Spent blocks' storage, handed out again: only allocation changes, since each reused
        //   block is cleared and written in full
        for (child_vars, block) in contributions {
            spare_vars.push(child_vars);
            spare_blocks.push(block);
        }

        let root = node.parent.is_none();
        debug_assert!(!root || below.is_empty());
        let factored = factor_front(
            front,
            m,
            fully_summed,
            NumCast::from(context.threshold).expect("the threshold is representable"),
            root,
            scratch,
            context.node_order,
            (
                spare_order.pop().unwrap_or_default(),
                spare_dblocks.pop().unwrap_or_default(),
            ),
        );
        for &v in vars.iter() {
            position[v] = usize::MAX;
        }
        let factored = factored.map_err(|error| match error {
            SparseError::SingularMatrix { pivot } => SparseError::SingularMatrix {
                pivot: self.cholesky.perm[vars[pivot]],
            },
            other => other,
        })?;
        let mut ordered = spare_ordered.pop().unwrap_or_default();
        ordered.clear();
        ordered.extend(factored.order.iter().map(|&p| vars[p]));
        spare_order.push(factored.order);
        let e = factored.eliminated;
        // The contribution block, everything not eliminated, goes to the parent
        let block = (e < m).then(|| {
            let rest = m - e;
            let mut block = spare_blocks.pop().unwrap_or_default();
            block.clear();
            block.resize(rest * rest, T::zero());
            for b in 0..rest {
                for a in 0..rest {
                    block[a + b * rest] = front[(e + a) + (e + b) * m];
                }
            }
            let mut rest_vars = spare_vars.pop().unwrap_or_default();
            rest_vars.clear();
            rest_vars.extend_from_slice(&ordered[e..]);
            (rest_vars, block)
        });
        // L's columns, at their exact size; the front is reused
        let mut l = spare_l.pop().unwrap_or_default();
        l.clear();
        l.extend_from_slice(&front[..m * e]);
        // Inside a 2×2 block, L's entry is zero (the kernel leaves D's there)
        let mut at = 0;
        for block in &factored.blocks {
            if block.size() == 2 {
                l[(at + 1) + at * m] = T::zero();
            }
            at += block.size();
        }
        let front = Front {
            vars: ordered,
            eliminated: e,
            l,
            blocks: factored.blocks,
        };
        Ok((front, block, delayed))
    }
}

/// A contribution block waiting for its parent: its variables, and its dense square, both
/// triangles
type Block<T> = (Vec<usize>, Vec<T>);

/// What every front's factorization reads: `P A Pᵀ`'s values, the pivot threshold, and the
/// trailing order from which the dense kernel's update splits across threads
struct Context<'a, T> {
    c_data: &'a [T],
    /// As the options give it; each front converts it to `T::Real`, as the serial
    ///   factorization did once
    threshold: f64,
    node_order: usize,
}

/// Working space reused from front to front by one thread: the front, the kernel's scratch,
/// the front's variables, spent contribution blocks' storage, and where each variable sits in
/// the current front (`usize::MAX` outside it)
struct Workspace<T> {
    front: Vec<T>,
    scratch: Vec<T>,
    vars: Vec<usize>,
    spare_blocks: Vec<Vec<T>>,
    spare_vars: Vec<Vec<usize>>,
    position: Vec<usize>,
    /// Storage for the fronts' results, from spent fronts: their `L` columns, their variables,
    /// `D`'s blocks, and the kernel's orders
    spare_l: Vec<Vec<T>>,
    spare_ordered: Vec<Vec<usize>>,
    spare_dblocks: Vec<Vec<DBlock<T>>>,
    spare_order: Vec<Vec<usize>>,
}

impl<T> Workspace<T> {
    fn new(n: usize) -> Workspace<T> {
        Workspace {
            front: Vec::new(),
            scratch: Vec::new(),
            vars: Vec::new(),
            spare_blocks: Vec::new(),
            spare_vars: Vec::new(),
            position: vec![usize::MAX; n],
            spare_l: Vec::new(),
            spare_ordered: Vec::new(),
            spare_dblocks: Vec::new(),
            spare_order: Vec::new(),
        }
    }
}

/// The work past which the factorization goes parallel (feature `rayon`): a subtree's (in
/// [`FrontTree`]'s multiply-adds), and a front's trailing order for its dense kernel's update.
/// Fixed, not a thread count, and the results don't depend on them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Cutoffs {
    pub(crate) subtree_work: u64,
    pub(crate) node_order: usize,
}

impl Cutoffs {
    /// The factorization's own, from a calibration over six OPF workloads at 14 threads
    /// (`docs/benchmarks/o5-parallel-ldlt-2026-11-20.md`): subtrees of 2²⁰ multiply-adds or
    /// more in parallel, and no column split. The split never beat the timing noise, and from
    /// trailing order 256 it made solves slower.
    pub(crate) const DEFAULT: Cutoffs = Cutoffs {
        subtree_work: 1 << 20,
        node_order: usize::MAX,
    };
    /// Every subtree and update on the parallel paths
    #[allow(dead_code)]
    pub(crate) const ALWAYS: Cutoffs = Cutoffs {
        subtree_work: 0,
        node_order: 0,
    };
    /// None: the serial factorization
    #[allow(dead_code)]
    pub(crate) const NEVER: Cutoffs = Cutoffs {
        subtree_work: u64::MAX,
        node_order: usize::MAX,
    };
}

#[cfg(feature = "rayon")]
thread_local! {
    /// How many `serially` scopes this thread is inside
    static SERIALLY: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Runs `f` with this thread's multifrontal LDLᵀ factorizations serial, though the `rayon`
/// feature is on. The results are the same either way; it's for callers that already run many
/// factorizations in parallel, such as a batch of OPF solves, where the nested parallelism only
/// competes with theirs. Without the feature it just runs `f`.
#[doc(hidden)]
pub fn serially<R>(f: impl FnOnce() -> R) -> R {
    #[cfg(feature = "rayon")]
    {
        SERIALLY.with(|depth| depth.set(depth.get() + 1));
        // Decremented however f returns
        struct Leave;
        impl Drop for Leave {
            fn drop(&mut self) {
                SERIALLY.with(|depth| depth.set(depth.get() - 1));
            }
        }
        let _leave = Leave;
        f()
    }
    #[cfg(not(feature = "rayon"))]
    f()
}

/// Runs `a` and `b`, concurrently with `rayon::join` if `parallel` and the `rayon` feature is
/// on, and in order otherwise, also on a one-thread pool and inside a [`serially`] scope. For
/// two independent computations whose results the caller combines afterwards, which are the
/// same either way.
#[doc(hidden)]
pub fn concurrently<A, B, RA, RB>(parallel: bool, a: A, b: B) -> (RA, RB)
where
    A: FnOnce() -> RA + Send,
    B: FnOnce() -> RB + Send,
    RA: Send,
    RB: Send,
{
    #[cfg(feature = "rayon")]
    if parallel && rayon::current_num_threads() > 1 && SERIALLY.with(|depth| depth.get() == 0) {
        return rayon::join(a, b);
    }
    #[cfg(not(feature = "rayon"))]
    let _ = parallel;
    (a(), b())
}

/// Of two failures, the one at the lower front index, the serial factorization's first
#[cfg(feature = "rayon")]
fn first_failure(
    kept: Option<(usize, SparseError)>,
    new: (usize, SparseError),
) -> Option<(usize, SparseError)> {
    match kept {
        Some(kept) if kept.0 < new.0 => Some(kept),
        _ => Some(new),
    }
}

impl<T: Scalar> CsrMatrix<T> {
    /// Analyzes this matrix's pattern for the multifrontal LDLᵀ ([`Ldlt`]) under
    /// `ordering`: `analyze_cholesky(ordering)`, with the supernodes and maps of
    /// [`SymbolicLdlt`].
    ///
    /// Returns `SparseError::NotSquare` for a non-square matrix, and
    /// `SparseError::UnsupportedOrdering` for `Ordering::Colamd`.
    pub fn analyze_ldlt(&self, ordering: Ordering) -> Result<SymbolicLdlt, SparseError> {
        Ok(SymbolicLdlt::new(self.analyze_cholesky(ordering)?))
    }

    /// [`analyze_ldlt`](Self::analyze_ldlt) with a permutation of the caller's, as
    /// [`analyze_cholesky_with_permutation`](Self::analyze_cholesky_with_permutation) takes it,
    /// with the default amalgamation. The amalgamation's postorder renumbers the columns but
    /// keeps every column after its descendants in the elimination tree, so an ordering that
    /// places one position after another it's coupled to keeps that order.
    ///
    /// Returns `SparseError::NotSquare` for a non-square matrix, and
    /// `SparseError::InvalidPermutation` unless `perm` is a permutation of its order.
    pub fn analyze_ldlt_with_permutation(
        &self,
        perm: Vec<usize>,
    ) -> Result<SymbolicLdlt, SparseError> {
        Ok(SymbolicLdlt::new(
            self.analyze_cholesky_with_permutation(perm)?,
        ))
    }
}

impl<T: Scalar + ComplexFloat + Send + Sync> CsrMatrix<T> {
    /// Factors a symmetric (for a complex element type, Hermitian) matrix as
    /// `P A Pᵀ = Q L D Lᴴ Qᵀ` with threshold Bunch–Kaufman pivoting, ordered by approximate
    /// minimum degree, with the default options: `analyze_cholesky(Ordering::Amd)` followed by
    /// [`factor_ldlt`](SymbolicCholesky::factor_ldlt). Analyze once to factor
    /// several matrices with the same pattern.
    ///
    /// Returns `SparseError::NotSquare`, `SparseError::NotSymmetric`, or
    /// `SparseError::SingularMatrix`.
    pub fn ldlt(&self) -> Result<Ldlt<T>, SparseError> {
        self.analyze_ldlt(Ordering::Amd)?
            .factor(self, LdltOptions::default())
    }
}

impl<T: Scalar + ComplexFloat> Ldlt<T> {
    /// `A`'s inertia, from `D`'s blocks (Sylvester's law of inertia). `A` is nonsingular, since
    /// a missing pivot fails the factorization.
    pub fn inertia(&self) -> Inertia {
        let (positive, negative) = self
            .fronts
            .iter()
            .flat_map(|front| &front.blocks)
            .map(DBlock::signs)
            .fold((0, 0), |(p, n), (bp, bn)| (p + bp, n + bn));
        Inertia { positive, negative }
    }

    /// How many times a variable was delayed from a front to its parent, over the whole
    /// factorization (a variable delayed twice counts twice): 0 when every pivot was taken
    /// where the analysis expected it.
    pub fn delayed(&self) -> usize {
        self.delayed
    }

    /// How many of `D`'s blocks are 2×2.
    pub fn two_by_two(&self) -> usize {
        self.fronts
            .iter()
            .flat_map(|front| &front.blocks)
            .filter(|block| block.size() == 2)
            .count()
    }

    /// The entries of `L` stored below its diagonal, over all fronts: what the analysis's
    /// Cholesky pattern predicts when nothing is delayed, and more when something is.
    pub fn l_entries(&self) -> usize {
        self.fronts
            .iter()
            .map(|front| {
                let (m, e) = (front.vars.len(), front.eliminated);
                e * m - e * (e + 1) / 2
            })
            .sum()
    }

    /// The ordering `P`: row and column `i` of `P A Pᵀ` are row and column `permutation()[i]`
    /// of `A`. The pivots' order within it is the fronts' business.
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }
}

/// The subtree solve work past which the solves go parallel (feature `rayon`), in
/// [`FrontTree`]'s multiply-adds: fixed, not a thread count, and the results don't depend on it.
/// 2²⁰, from a calibration on the predictor-corrector method's DC and AC solves at 8 and 14
/// threads (`docs/benchmarks/o5-parallel-solves-2026-11-21.md`): below 2¹⁶ the tasks' overhead
/// costs up to ×2.6, and from 2²² the gain on the largest case is gone.
#[cfg(feature = "rayon")]
const SOLVE_WORK: u64 = 1 << 20;

#[cfg(all(test, feature = "rayon"))]
thread_local! {
    /// The tests' cutoff for this thread's solves, in place of `SOLVE_WORK`: 0 puts every
    ///   subtree on the parallel path, even on a one-thread pool, and `u64::MAX` none
    static SOLVE_CUTOFF: std::cell::Cell<Option<u64>> = const { std::cell::Cell::new(None) };
}

impl<T: Scalar + ComplexFloat + Send + Sync> Ldlt<T> {
    /// Overwrites `y`, in `P A Pᵀ`'s numbering, with the solution of `P A Pᵀ y' = y`: forward
    /// substitution front by front in elimination order, `D`'s blocks, then back substitution in
    /// reverse. With the `rayon` feature, independent subtrees in parallel, bit for bit
    fn substitute_permuted(&self, y: &mut [T]) {
        #[cfg(feature = "rayon")]
        if let Some(cutoff) = self.parallel_cutoff() {
            return self.par_substitute(y, cutoff);
        }
        for front in &self.fronts {
            let m = front.vars.len();
            for j in 0..front.eliminated {
                let y_j = y[front.vars[j]];
                let column = &front.l[j * m..(j + 1) * m];
                for i in j + 1..m {
                    y[front.vars[i]] = y[front.vars[i]] - column[i] * y_j;
                }
            }
        }
        for front in &self.fronts {
            solve_diagonal(front, y, 0);
        }
        for front in self.fronts.iter().rev() {
            let m = front.vars.len();
            for j in (0..front.eliminated).rev() {
                let column = &front.l[j * m..(j + 1) * m];
                let mut y_j = y[front.vars[j]];
                for i in j + 1..m {
                    y_j = y_j - column[i].conj() * y[front.vars[i]];
                }
                y[front.vars[j]] = y_j;
            }
        }
    }

    /// The subtree cutoff for a parallel solve, or `None` for the serial one: as for the
    ///   factorization, on a one-thread pool, inside a `serially` scope, or when the whole tree
    ///   is below it
    #[cfg(feature = "rayon")]
    fn parallel_cutoff(&self) -> Option<u64> {
        let forced = forced_cutoff();
        let cutoff = forced.unwrap_or(SOLVE_WORK);
        let tree = &self.tree;
        let total: u64 = tree.roots.iter().map(|&root| tree.work[root]).sum();
        let idle = forced.is_none()
            && (rayon::current_num_threads() == 1 || SERIALLY.with(|depth| depth.get() > 0));
        (tree.contiguous && cutoff != u64::MAX && total >= cutoff && !idle).then_some(cutoff)
    }

    /// The parallel solve. Every entry of `y` sees the serial solve's operations in its order:
    ///
    /// - **The forward sweep:** the roots' subtrees in parallel, each splitting into its
    ///   children's in parallel, then its own front. A task updates the entries in its own
    ///   columns in place and records each update of an ancestor's entry, `(row, term)`, in
    ///   order; the parent applies its children's records in index order before its front, as
    ///   the serial sweep would have, and passes up those above it. Each front's `D` solve
    ///   follows its forward step, which is its entries' last update.
    /// - **The backward sweep:** top down, a subtree's front first, then its children's
    ///   subtrees in parallel, reading their ancestors' entries, final by then, through
    ///   [`Above`].
    #[cfg(feature = "rayon")]
    fn par_substitute(&self, y: &mut [T], cutoff: u64) {
        let groups = self.tree.groups(&self.tree.roots, cutoff);
        let records = self.forward_groups(&groups, y, 0, cutoff);
        debug_assert!(records.iter().all(Vec::is_empty), "no row above the roots");
        let top = Above {
            lo: 0,
            values: &[],
            up: None,
        };
        self.backward_groups(&groups, y, 0, cutoff, &top);
    }

    /// `groups`, side by side over `y`, which holds columns `lo..` of the whole vector, in
    ///   parallel: each one's records of its updates past its columns, in order
    #[cfg(feature = "rayon")]
    fn forward_groups(
        &self,
        groups: &[Group],
        y: &mut [T],
        lo: usize,
        cutoff: u64,
    ) -> Vec<Vec<(usize, T)>> {
        use rayon::prelude::*;
        self.split(groups, y, lo)
            .into_par_iter()
            .map(|(group, y, lo)| match *group {
                Group::Subtree(s) => self.forward_subtree(s, y, cutoff),
                Group::Run(ref fronts) => {
                    let mut records = Vec::new();
                    self.forward_run(fronts.clone(), y, lo, &mut records);
                    records
                }
            })
            .collect()
    }

    /// Subtree `s`'s forward sweep over `y`, its columns: its children's in parallel, their
    ///   records, then its front
    #[cfg(feature = "rayon")]
    fn forward_subtree(&self, s: usize, y: &mut [T], cutoff: u64) -> Vec<(usize, T)> {
        let tree = &self.tree;
        let Range { start: lo, end: hi } = tree.columns(s);
        let groups = tree.groups(&tree.children[s], cutoff);
        let below = tree.start[s] - lo;
        let records = self.forward_groups(&groups, &mut y[..below], lo, cutoff);
        let mut outgoing = Vec::new();
        for (row, term) in records.into_iter().flatten() {
            if row < hi {
                y[row - lo] = y[row - lo] - term;
            } else {
                outgoing.push((row, term));
            }
        }
        self.forward_run(s..s + 1, y, lo, &mut outgoing);
        outgoing
    }

    /// The forward sweep and `D` solves of `fronts`, whole subtrees side by side, over `y`,
    ///   which holds their columns from `lo`: the updates past them recorded in `records`
    #[cfg(feature = "rayon")]
    fn forward_run(
        &self,
        fronts: Range<usize>,
        y: &mut [T],
        lo: usize,
        records: &mut Vec<(usize, T)>,
    ) {
        let hi = lo + y.len();
        for front in &self.fronts[fronts] {
            let m = front.vars.len();
            for j in 0..front.eliminated {
                let y_j = y[front.vars[j] - lo];
                let column = &front.l[j * m..(j + 1) * m];
                for (&v, &l) in front.vars[j + 1..].iter().zip(&column[j + 1..]) {
                    if v < hi {
                        y[v - lo] = y[v - lo] - l * y_j;
                    } else {
                        records.push((v, l * y_j));
                    }
                }
            }
            solve_diagonal(front, y, lo);
        }
    }

    /// `groups`, side by side over `y` (columns `lo..`), in parallel, their ancestors' entries
    ///   in `above`
    #[cfg(feature = "rayon")]
    fn backward_groups(
        &self,
        groups: &[Group],
        y: &mut [T],
        lo: usize,
        cutoff: u64,
        above: &Above<'_, T>,
    ) {
        use rayon::prelude::*;
        self.split(groups, y, lo)
            .into_par_iter()
            .for_each(|(group, y, lo)| match *group {
                Group::Subtree(s) => self.backward_subtree(s, y, cutoff, above),
                Group::Run(ref fronts) => self.backward_run(fronts.clone(), y, lo, above),
            });
    }

    /// Subtree `s`'s backward sweep over `y`, its columns: its front, then its children's in
    ///   parallel
    #[cfg(feature = "rayon")]
    fn backward_subtree(&self, s: usize, y: &mut [T], cutoff: u64, above: &Above<'_, T>) {
        let tree = &self.tree;
        let lo = tree.columns(s).start;
        self.backward_run(s..s + 1, y, lo, above);
        let (below, own) = y.split_at_mut(tree.start[s] - lo);
        let above = Above {
            lo: tree.start[s],
            values: own,
            up: Some(above),
        };
        let groups = tree.groups(&tree.children[s], cutoff);
        self.backward_groups(&groups, below, lo, cutoff, &above);
    }

    /// The backward sweep of `fronts`, in reverse, over `y` (their columns from `lo`)
    #[cfg(feature = "rayon")]
    fn backward_run(&self, fronts: Range<usize>, y: &mut [T], lo: usize, above: &Above<'_, T>) {
        let hi = lo + y.len();
        for front in self.fronts[fronts].iter().rev() {
            let m = front.vars.len();
            for j in (0..front.eliminated).rev() {
                let column = &front.l[j * m..(j + 1) * m];
                let mut y_j = y[front.vars[j] - lo];
                for (&v, &l) in front.vars[j + 1..].iter().zip(&column[j + 1..]) {
                    let y_v = if v < hi { y[v - lo] } else { above.get(v) };
                    y_j = y_j - l.conj() * y_v;
                }
                y[front.vars[j] - lo] = y_j;
            }
        }
    }

    /// `y`, holding columns `lo..`, split into `groups`' columns, side by side, each with its
    ///   first column
    #[cfg(feature = "rayon")]
    fn split<'g, 'y>(
        &self,
        groups: &'g [Group],
        mut y: &'y mut [T],
        mut lo: usize,
    ) -> Vec<(&'g Group, &'y mut [T], usize)> {
        let tree = &self.tree;
        let mut parts = Vec::with_capacity(groups.len());
        for group in groups {
            let hi = match *group {
                Group::Subtree(s) => tree.start[s + 1],
                Group::Run(ref fronts) => tree.start[fronts.end],
            };
            let (part, rest) = std::mem::take(&mut y).split_at_mut(hi - lo);
            parts.push((group, part, lo));
            (y, lo) = (rest, hi);
        }
        parts
    }
}

/// The tests' cutoff for this thread's solves, if they set one
#[cfg(feature = "rayon")]
fn forced_cutoff() -> Option<u64> {
    #[cfg(test)]
    return SOLVE_CUTOFF.with(std::cell::Cell::get);
    #[cfg(not(test))]
    None
}

/// Front `front`'s `D` solve, on `y` holding the vector's entries from `lo`
fn solve_diagonal<T: Scalar + ComplexFloat>(front: &Front<T>, y: &mut [T], lo: usize) {
    let mut at = 0;
    for block in &front.blocks {
        match *block {
            DBlock::One(d) => {
                let v = front.vars[at] - lo;
                y[v] = y[v] / d;
            }
            DBlock::Two { a, b, c } => {
                let (v0, v1) = (front.vars[at] - lo, front.vars[at + 1] - lo);
                let det = a * c - b * b.conj();
                let (y0, y1) = (y[v0], y[v1]);
                y[v0] = (c * y0 - b.conj() * y1) / det;
                y[v1] = (a * y1 - b * y0) / det;
            }
        }
        at += block.size();
    }
}

impl<T: Scalar + ComplexFloat + Send + Sync> Factorization<T> for Ldlt<T> {
    fn dim(&self) -> usize {
        self.n
    }
}

impl<T: Scalar + ComplexFloat + Send + Sync> Substitute<T> for Ldlt<T> {
    /// `P b` into `scratch`, the substitution there in `P A Pᵀ`'s numbering, and the result
    /// back through `P`
    fn substitute(&self, x: &mut [T], scratch: &mut Vec<T>) {
        scratch.clear();
        scratch.extend(self.perm.iter().map(|&row| x[row]));
        self.substitute_permuted(scratch);
        for (&row, &value) in self.perm.iter().zip(scratch.iter()) {
            x[row] = value;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The natural ordering without amalgamation, so the fronts are the fundamental supernodes
    fn symbolic(n: usize, entries: &[(usize, usize, f64)]) -> SymbolicLdlt {
        let a = CsrMatrix::from_triplets(n, n, entries.to_vec()).unwrap();
        let cholesky = a.analyze_cholesky(Ordering::Natural).unwrap();
        SymbolicLdlt::with_amalgamation(cholesky, Amalgamation::NONE)
    }

    #[test]
    fn a_star_has_its_leaves_as_its_hubs_children() {
        // Columns 0 and 1 each touch only the hub, 2: two leaf fronts of work 1 under the hub's
        //   of work 0 (front_tree_test's star)
        let s = symbolic(
            3,
            &[
                (0, 0, 4.0),
                (1, 1, 4.0),
                (2, 2, 4.0),
                (0, 2, 1.0),
                (2, 0, 1.0),
                (1, 2, 1.0),
                (2, 1, 1.0),
            ],
        );
        assert_eq!(s.children, vec![vec![], vec![], vec![0, 1]]);
        assert_eq!(s.subtree_work, vec![1, 1, 2]);
    }

    #[cfg(feature = "rayon")]
    mod parallel {
        use super::*;

        /// A deterministic stream in [0, 1)
        fn stream(seed: u64) -> impl FnMut() -> f64 {
            let mut state = seed;
            move || {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (state >> 11) as f64 / (1u64 << 53) as f64
            }
        }

        /// A symmetric matrix of order `n`: a band of half-width `band`, and about `extra` random
        ///   off-band pairs, with diagonals of either sign (indefinite)
        fn indefinite(n: usize, band: usize, extra: usize, seed: u64) -> CsrMatrix<f64> {
            let mut random = stream(seed);
            let mut entries = vec![];
            for i in 0..n {
                entries.push((i, i, 4.0 * (random() - 0.4)));
                for j in i + 1..(i + band + 1).min(n) {
                    let v = random() - 0.5;
                    entries.push((i, j, v));
                    entries.push((j, i, v));
                }
            }
            for _ in 0..extra {
                let (i, j) = (
                    (random() * n as f64) as usize,
                    (random() * n as f64) as usize,
                );
                if i != j {
                    let v = random() - 0.5;
                    entries.push((i, j, v));
                    entries.push((j, i, v));
                }
            }
            CsrMatrix::from_triplets(n, n, entries).unwrap()
        }

        /// A saddle-point (KKT-shaped) matrix `[[H, Jᵀ], [J, 0]]`: H a grid Laplacian on `g × g`
        ///   nodes plus a random diagonal, and J one row per `stride` nodes coupling neighbours; its
        ///   zero block forces 2×2 pivots and delays
        fn saddle(g: usize, stride: usize, seed: u64) -> CsrMatrix<f64> {
            let mut random = stream(seed);
            let nh = g * g;
            let rows: Vec<usize> = (0..nh).step_by(stride).collect();
            let n = nh + rows.len();
            let mut entries = vec![];
            for r in 0..g {
                for c in 0..g {
                    let i = r * g + c;
                    entries.push((i, i, 4.0 + random()));
                    for (dr, dc) in [(0, 1), (1, 0)] {
                        let (r2, c2) = (r + dr, c + dc);
                        if r2 < g && c2 < g {
                            let j = r2 * g + c2;
                            entries.push((i, j, -1.0));
                            entries.push((j, i, -1.0));
                        }
                    }
                }
            }
            for (k, &node) in rows.iter().enumerate() {
                let row = nh + k;
                for j in [node, (node + 1) % nh] {
                    let v = 1.0 + random();
                    entries.push((row, j, v));
                    entries.push((j, row, v));
                }
            }
            CsrMatrix::from_triplets(n, n, entries).unwrap()
        }

        /// `f` on a rayon pool of `threads` threads
        fn on<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> R {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(f)
        }

        /// The parallel factorization is the serial one bit for bit, on pools of 1, 2, 7 and 28
        ///   threads, for each ordering: with its default cutoffs, and if `always`, with every
        ///   subtree and update on its parallel paths
        fn check(label: &str, a: &CsrMatrix<f64>, always: bool) {
            for ordering in [Ordering::Natural, Ordering::Amd] {
                let symbolic = a.analyze_ldlt(ordering).unwrap();
                let options = LdltOptions::default();
                let serial = format!(
                    "{:?}",
                    symbolic.factor_with(a, options, Cutoffs::NEVER, None)
                );
                for threads in [1, 2, 7, 28] {
                    if always {
                        let parallel = on(threads, || {
                            symbolic.factor_with(a, options, Cutoffs::ALWAYS, None)
                        });
                        assert_eq!(
                            format!("{parallel:?}"),
                            serial,
                            "{label}, {ordering:?}, {threads} threads"
                        );
                    }
                    let defaults = on(threads, || symbolic.factor(a, options));
                    assert_eq!(format!("{defaults:?}"), serial, "{label}, defaults");
                }
            }
        }

        #[test]
        fn the_parallel_factorization_is_the_serial_one_bit_for_bit() {
            check("indefinite", &indefinite(300, 3, 200, 1), true);
            let kkt = saddle(30, 3, 2);
            let ldlt = kkt
                .analyze_ldlt(Ordering::Amd)
                .unwrap()
                .factor(&kkt, LdltOptions::default())
                .unwrap();
            assert!(ldlt.two_by_two() > 0, "the saddle point has 2×2 pivots");
            check("saddle point", &kkt, true);
            // Large enough for fronts past the default cutoffs, which it runs with: forcing every
            //   pivot's update onto threads would only cost the test time
            check("large saddle point", &saddle(70, 2, 3), false);
        }

        /// `f` with this thread's solves at `cutoff` (see `SOLVE_CUTOFF`)
        fn at_cutoff<R>(cutoff: u64, f: impl FnOnce() -> R) -> R {
            SOLVE_CUTOFF.with(|forced| forced.set(Some(cutoff)));
            let result = f();
            SOLVE_CUTOFF.with(|forced| forced.set(None));
            result
        }

        /// Right-hand sides for `a`, random
        fn right_hand_sides<T: Scalar>(
            n: usize,
            count: usize,
            value: impl Fn(f64, f64) -> T,
        ) -> Vec<Vec<T>> {
            let mut random = stream(7);
            (0..count)
                .map(|_| {
                    (0..n)
                        .map(|_| value(random() - 0.5, random() - 0.5))
                        .collect()
                })
                .collect()
        }

        /// The parallel solves are the serial ones bit for bit, on pools of 1, 2, 7 and 28
        ///   threads, for each ordering: with every subtree on the parallel path, and with the
        ///   default cutoff. Single and many right-hand sides, and the refined solve
        fn check_solves<T>(label: &str, a: &CsrMatrix<T>, value: impl Fn(f64, f64) -> T)
        where
            T: Scalar + ComplexFloat + Send + Sync + std::fmt::Debug,
            T::Real: std::fmt::Debug,
        {
            let b = right_hand_sides(a.shape.0, 3, value);
            for ordering in [Ordering::Natural, Ordering::Amd] {
                let ldlt = a
                    .analyze_ldlt(ordering)
                    .unwrap()
                    .factor(a, LdltOptions::default())
                    .unwrap();
                let solves = || {
                    let one = ldlt.solve(&b[0]).unwrap();
                    let many = ldlt.solve_many(&b).unwrap();
                    let refined = ldlt.solve_refined(a, &b[1], 3).unwrap();
                    format!("{one:?} {many:?} {refined:?}")
                };
                let serial = at_cutoff(u64::MAX, || {
                    assert_eq!(ldlt.parallel_cutoff(), None);
                    solves()
                });
                for threads in [1, 2, 7, 28] {
                    let parallel = on(threads, || {
                        at_cutoff(0, || {
                            assert_eq!(ldlt.parallel_cutoff(), Some(0), "{label}: parallel");
                            solves()
                        })
                    });
                    assert_eq!(parallel, serial, "{label}, {ordering:?}, {threads} threads");
                    let defaults = on(threads, solves);
                    assert_eq!(defaults, serial, "{label}, {ordering:?}, defaults");
                }
            }
        }

        #[test]
        fn the_parallel_solves_are_the_serial_ones_bit_for_bit() {
            let real = |re: f64, _: f64| re;
            check_solves("indefinite", &indefinite(300, 3, 200, 1), real);
            let kkt = saddle(30, 3, 2);
            let ldlt = kkt
                .analyze_ldlt(Ordering::Amd)
                .unwrap()
                .factor(&kkt, LdltOptions::default())
                .unwrap();
            assert!(
                ldlt.two_by_two() > 0 && ldlt.delayed() > 0,
                "2×2 pivots and delays"
            );
            check_solves("saddle point", &kkt, real);
            check_solves("large saddle point", &saddle(70, 2, 3), real);
            // Two independent blocks: two roots, side by side
            let base = indefinite(150, 2, 80, 5);
            let mut entries = vec![];
            for offset in [0, 150] {
                for i in 0..150 {
                    let (cols, values) = base.outer(i).unwrap();
                    for (&j, &v) in cols.iter().zip(values) {
                        entries.push((offset + i, offset + j, v));
                    }
                }
            }
            let blocks = CsrMatrix::from_triplets(300, 300, entries).unwrap();
            assert!(
                blocks
                    .analyze_ldlt(Ordering::Amd)
                    .unwrap()
                    .front_tree()
                    .roots
                    >= 2
            );
            check_solves("two blocks", &blocks, real);
        }

        #[test]
        fn the_parallel_solves_are_the_serial_ones_for_a_hermitian_matrix() {
            use num_complex::Complex64;
            // The saddle point with an imaginary part on each off-diagonal pair, conjugate
            //   across the diagonal, so the backward sweep's conjugates matter
            let real = saddle(30, 3, 2);
            let mut random = stream(9);
            let mut entries = vec![];
            for i in 0..real.shape.0 {
                let (cols, values) = real.outer(i).unwrap();
                for (&j, &v) in cols.iter().zip(values) {
                    if i < j {
                        let im = random() - 0.5;
                        entries.push((i, j, Complex64::new(v, im)));
                        entries.push((j, i, Complex64::new(v, -im)));
                    } else if i == j {
                        entries.push((i, j, Complex64::new(v, 0.0)));
                    }
                }
            }
            let a = CsrMatrix::from_triplets(real.shape.0, real.shape.0, entries).unwrap();
            check_solves("hermitian", &a, Complex64::new);
        }

        #[test]
        fn serially_nests_and_restores_its_count_on_a_panic() {
            let depth = || SERIALLY.with(|depth| depth.get());
            assert_eq!(depth(), 0);
            serially(|| {
                assert_eq!(depth(), 1);
                serially(|| assert_eq!(depth(), 2));
                assert_eq!(depth(), 1);
            });
            let caught = std::panic::catch_unwind(|| serially(|| panic!("inside")));
            assert!(caught.is_err());
            assert_eq!(depth(), 0, "restored after the panic");
            // And the factorization inside it is the serial one
            let a = indefinite(300, 3, 200, 1);
            let symbolic = a.analyze_ldlt(Ordering::Amd).unwrap();
            let options = LdltOptions::default();
            let serial = format!(
                "{:?}",
                symbolic.factor_with(&a, options, Cutoffs::NEVER, None)
            );
            assert_eq!(
                format!("{:?}", serially(|| symbolic.factor(&a, options))),
                serial
            );
        }

        #[test]
        fn a_singular_matrix_fails_as_the_serial_factorization_does() {
            // Two independent blocks, each with a row and column stored but zero: every update
            //   into them multiplies by one of their zeros, so the columns stay exactly zero and
            //   a front in each block fails. The blocks are separate subtrees, so both failures
            //   happen, and the error must be the first in the serial order, not whichever
            //   thread finishes last
            let base = indefinite(100, 2, 50, 4);
            let mut entries = vec![];
            for offset in [0, 100] {
                for i in 0..100 {
                    let (cols, values) = base.outer(i).unwrap();
                    for (&j, &v) in cols.iter().zip(values) {
                        let zero = i == 30 || j == 30;
                        entries.push((offset + i, offset + j, if zero { 0.0 } else { v }));
                    }
                }
            }
            let a = CsrMatrix::from_triplets(200, 200, entries).unwrap();
            let serial = a.analyze_ldlt(Ordering::Amd).unwrap().factor_with(
                &a,
                LdltOptions::default(),
                Cutoffs::NEVER,
                None,
            );
            assert!(
                matches!(serial, Err(SparseError::SingularMatrix { .. })),
                "{serial:?}"
            );
            check("singular", &a, true);
        }
    }

    #[test]
    fn independent_blocks_are_two_roots() {
        let block = |o: usize| {
            [
                (o, o, 4.0),
                (o, o + 1, 1.0),
                (o + 1, o, 1.0),
                (o + 1, o + 1, 4.0),
            ]
        };
        let entries: Vec<_> = block(0).into_iter().chain(block(2)).collect();
        let s = symbolic(4, &entries);
        assert_eq!(s.children, vec![Vec::<usize>::new(), vec![]]);
        assert_eq!(s.subtree_work, vec![1, 1]);
        assert_eq!(s.subtree_work.iter().sum::<u64>(), s.front_tree().work);
    }
}
