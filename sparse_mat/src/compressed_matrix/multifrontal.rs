//! Sparse symmetric (Hermitian) indefinite `L D Lᴴ` by the multifrontal method, with threshold
//! Bunch–Kaufman pivoting and delayed pivots: [`IndefiniteLdlt`].

use num_complex::ComplexFloat;
use num_traits::NumCast;

use super::bunch_kaufman::{factor_front, DBlock};
use super::cholesky::SymbolicCholesky;
use super::{CompressedMatrix, CsrMatrix, Inertia, Ordering};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// Options for [`SymbolicCholesky::factor_indefinite`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndefiniteOptions {
    /// The threshold `u` of the pivot test, in `(0, 0.5]`: every entry of `L` is at most `1/u`
    /// in magnitude. Larger is more stable and delays more pivots; smaller keeps more of the
    /// ordering's sparsity. The default, 0.01, is MA57's.
    pub threshold: f64,
}

impl Default for IndefiniteOptions {
    fn default() -> IndefiniteOptions {
        IndefiniteOptions { threshold: 0.01 }
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
/// Bunch–Kaufman pivoting, from [`SymbolicCholesky::factor_indefinite`] or
/// [`CsrMatrix::ldlt_indefinite`]. `P` is the analysis's fill-reducing ordering and `Q` the
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
/// [`DenseLdlt`]: super::DenseLdlt
#[derive(Clone, Debug)]
pub struct IndefiniteLdlt<T: Scalar = f64> {
    n: usize,
    /// Row `i` of `P A Pᵀ` is row `perm[i]` of `A`
    perm: Vec<usize>,
    /// In elimination order
    fronts: Vec<Front<T>>,
    delayed: usize,
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
    /// on the multifrontal method (see [`IndefiniteLdlt`]). The analysis's ordering and
    /// elimination tree are used as they are; pivoting only reorders within and between fronts.
    ///
    /// This builds a [`SymbolicIndefinite`] on each call. To factor several matrices with the
    /// same pattern, build one once and use [`SymbolicIndefinite::factor`], whose results are
    /// the same bit for bit.
    ///
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, `SparseError::NotSymmetric` if `A ≠ Aᴴ`, and `SparseError::SingularMatrix`
    /// naming a row of `A` where no pivot could be found, the matrix being singular to working
    /// precision. Panics if `options.threshold` is outside `(0, 0.5]`.
    pub fn factor_indefinite<T: Scalar + ComplexFloat>(
        &self,
        matrix: &CsrMatrix<T>,
        options: IndefiniteOptions,
    ) -> Result<IndefiniteLdlt<T>, SparseError> {
        SymbolicIndefinite::new(self.clone()).factor(matrix, options)
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

/// The symbolic analysis for [`IndefiniteLdlt`]: the Cholesky analysis ([`SymbolicCholesky`]:
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
/// [`CsrMatrix::analyze_indefinite`] or from a Cholesky analysis with
/// [`SymbolicIndefinite::new`], and factor with [`SymbolicIndefinite::factor`].
#[derive(Clone, Debug)]
pub struct SymbolicIndefinite {
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
}

/// How [`SymbolicIndefinite`] merges supernodes into larger fronts (relaxed amalgamation,
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

impl SymbolicIndefinite {
    /// The analysis for the multifrontal LDLᵀ on `cholesky`'s ordering, elimination tree and
    /// pattern, with the default relaxed amalgamation ([`Amalgamation::default`]): the
    /// ordering is postordered and small supernodes merged.
    pub fn new(cholesky: SymbolicCholesky) -> SymbolicIndefinite {
        SymbolicIndefinite::with_amalgamation(cholesky, Amalgamation::default())
    }

    /// The analysis with the given amalgamation. With [`Amalgamation::NONE`] the fronts are the
    /// fundamental supernodes on `cholesky`'s own ordering. Otherwise the ordering is first
    /// postordered (the same elimination tree, renumbered so each subtree's columns are
    /// contiguous and a parent follows its last child), so that a child can merge into its
    /// parent, and the supernodes are merged by the rule `amalgamation` sets.
    pub fn with_amalgamation(
        cholesky: SymbolicCholesky,
        amalgamation: Amalgamation,
    ) -> SymbolicIndefinite {
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
        SymbolicIndefinite {
            cholesky,
            nodes,
            c_ptr,
            c_idx,
            source,
            mirror,
        }
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
    /// on the multifrontal method (see [`IndefiniteLdlt`]), with the analysis's supernodes and
    /// maps.
    ///
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, `SparseError::NotSymmetric` if `A ≠ Aᴴ`, and `SparseError::SingularMatrix`
    /// naming a row of `A` where no pivot could be found, the matrix being singular to working
    /// precision. Panics if `options.threshold` is outside `(0, 0.5]`.
    pub fn factor<T: Scalar + ComplexFloat>(
        &self,
        matrix: &CsrMatrix<T>,
        options: IndefiniteOptions,
    ) -> Result<IndefiniteLdlt<T>, SparseError> {
        let threshold = options.threshold;
        assert!(
            threshold > 0.0 && threshold <= 0.5,
            "the threshold must be in (0, 0.5], not {threshold}"
        );
        let threshold: T::Real = NumCast::from(threshold).expect("the threshold is representable");
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
        let nodes = &self.nodes;

        // Contribution blocks waiting for their parent: (variables, dense square, both triangles)
        let mut pending: Vec<Vec<(Vec<usize>, Vec<T>)>> = vec![Vec::new(); nodes.len()];
        // Reused working space: the front, the kernel's scratch, the front's variables, and
        //   spent contribution blocks' storage, handed out again (so blocks are still consumed
        //   in the order they were made, and every sum is as before)
        let mut front: Vec<T> = Vec::new();
        let mut scratch: Vec<T> = Vec::new();
        let mut vars: Vec<usize> = Vec::new();
        let mut spare_blocks: Vec<Vec<T>> = Vec::new();
        let mut spare_vars: Vec<Vec<usize>> = Vec::new();
        // Where each variable sits in the current front, or usize::MAX
        let mut position = vec![usize::MAX; n];
        let mut fronts = Vec::with_capacity(nodes.len());
        let mut delayed = 0;
        for (s, node) in nodes.iter().enumerate() {
            let (start, end, below) = (&node.start, &node.end, &node.below);
            let contributions = std::mem::take(&mut pending[s]);
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
            delayed += fully_summed - (end - start);
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
            for (child_vars, block) in contributions {
                spare_vars.push(child_vars);
                spare_blocks.push(block);
            }

            let root = node.parent.is_none();
            debug_assert!(!root || below.is_empty());
            let factored = factor_front(&mut front, m, fully_summed, threshold, root, &mut scratch)
                .map_err(|error| match error {
                    SparseError::SingularMatrix { pivot } => SparseError::SingularMatrix {
                        pivot: symbolic.perm[vars[pivot]],
                    },
                    other => other,
                })?;
            for &v in &vars {
                position[v] = usize::MAX;
            }
            let ordered: Vec<usize> = factored.order.iter().map(|&p| vars[p]).collect();
            let e = factored.eliminated;
            // The contribution block, everything not eliminated, goes to the parent
            if e < m {
                let rest = m - e;
                let mut block = spare_blocks.pop().unwrap_or_default();
                block.clear();
                block.resize(rest * rest, T::zero());
                for b in 0..rest {
                    for a in 0..rest {
                        block[a + b * rest] = front[(e + a) + (e + b) * m];
                    }
                }
                let p = node.parent.expect("a front with rows left has a parent");
                let mut rest_vars = spare_vars.pop().unwrap_or_default();
                rest_vars.clear();
                rest_vars.extend_from_slice(&ordered[e..]);
                pending[p].push((rest_vars, block));
            }
            // L's columns, at their exact size; the front is reused
            let mut l = front[..m * e].to_vec();
            // Inside a 2×2 block, L's entry is zero (the kernel leaves D's there)
            let mut at = 0;
            for block in &factored.blocks {
                if block.size() == 2 {
                    l[(at + 1) + at * m] = T::zero();
                }
                at += block.size();
            }
            fronts.push(Front {
                vars: ordered,
                eliminated: e,
                l,
                blocks: factored.blocks,
            });
        }
        Ok(IndefiniteLdlt {
            n,
            perm: symbolic.perm.clone(),
            fronts,
            delayed,
        })
    }
}

impl<T: Scalar> CsrMatrix<T> {
    /// Analyzes this matrix's pattern for the multifrontal LDLᵀ ([`IndefiniteLdlt`]) under
    /// `ordering`: `analyze_cholesky(ordering)`, with the supernodes and maps of
    /// [`SymbolicIndefinite`].
    ///
    /// Returns `SparseError::NotSquare` for a non-square matrix, and
    /// `SparseError::UnsupportedOrdering` for `Ordering::Colamd`.
    pub fn analyze_indefinite(
        &self,
        ordering: Ordering,
    ) -> Result<SymbolicIndefinite, SparseError> {
        Ok(SymbolicIndefinite::new(self.analyze_cholesky(ordering)?))
    }
}

impl<T: Scalar + ComplexFloat> CsrMatrix<T> {
    /// Factors a symmetric (for a complex element type, Hermitian) matrix as
    /// `P A Pᵀ = Q L D Lᴴ Qᵀ` with threshold Bunch–Kaufman pivoting, ordered by approximate
    /// minimum degree, with the default options: `analyze_cholesky(Ordering::Amd)` followed by
    /// [`factor_indefinite`](SymbolicCholesky::factor_indefinite). Analyze once to factor
    /// several matrices with the same pattern.
    ///
    /// Returns `SparseError::NotSquare`, `SparseError::NotSymmetric`, or
    /// `SparseError::SingularMatrix`.
    pub fn ldlt_indefinite(&self) -> Result<IndefiniteLdlt<T>, SparseError> {
        self.analyze_indefinite(Ordering::Amd)?
            .factor(self, IndefiniteOptions::default())
    }
}

impl<T: Scalar + ComplexFloat> IndefiniteLdlt<T> {
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

    /// Solves `A x = b`. Returns `SparseError::ShapeMismatch` unless `b` has one entry per row
    /// of `A`.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        self.check_rhs(b.len())?;
        let mut y: Vec<T> = self.perm.iter().map(|&row| b[row]).collect();
        self.substitute(&mut y);
        let mut x = vec![T::zero(); self.n];
        for (&row, &value) in self.perm.iter().zip(&y) {
            x[row] = value;
        }
        Ok(x)
    }

    /// Solves `A x = b` for each right-hand side in `rhs`, each bit for bit what `solve` gives.
    /// Returns `SparseError::ShapeMismatch` naming the first right-hand side of the wrong length.
    pub fn solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        for b in rhs {
            self.check_rhs(b.len())?;
        }
        rhs.iter().map(|b| self.solve(b)).collect()
    }

    fn check_rhs(&self, len: usize) -> Result<(), SparseError> {
        if len == self.n {
            Ok(())
        } else {
            Err(SparseError::ShapeMismatch {
                left: (self.n, self.n),
                right: (len, 1),
            })
        }
    }

    /// Overwrites `y`, in `P A Pᵀ`'s numbering, with the solution of `P A Pᵀ y' = y`: forward
    /// substitution front by front in elimination order, `D`'s blocks, then back substitution in
    /// reverse
    fn substitute(&self, y: &mut [T]) {
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
            let mut at = 0;
            for block in &front.blocks {
                match *block {
                    DBlock::One(d) => {
                        let v = front.vars[at];
                        y[v] = y[v] / d;
                    }
                    DBlock::Two { a, b, c } => {
                        let (v0, v1) = (front.vars[at], front.vars[at + 1]);
                        let det = a * c - b * b.conj();
                        let (y0, y1) = (y[v0], y[v1]);
                        y[v0] = (c * y0 - b.conj() * y1) / det;
                        y[v1] = (a * y1 - b * y0) / det;
                    }
                }
                at += block.size();
            }
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
}

#[cfg(feature = "rayon")]
impl<T: Scalar + ComplexFloat + Send + Sync> IndefiniteLdlt<T> {
    /// [`IndefiniteLdlt::solve_many`] with the right-hand sides solved in parallel (feature
    /// `rayon`), each by one thread with exactly `solve`'s operations, so the results are bit for
    /// bit `solve_many`'s. Returns `SparseError::ShapeMismatch` as `solve_many` does, before
    /// solving any.
    pub fn par_solve_many(&self, rhs: &[Vec<T>]) -> Result<Vec<Vec<T>>, SparseError> {
        use rayon::prelude::*;
        for b in rhs {
            self.check_rhs(b.len())?;
        }
        rhs.par_iter().map(|b| self.solve(b)).collect()
    }
}
