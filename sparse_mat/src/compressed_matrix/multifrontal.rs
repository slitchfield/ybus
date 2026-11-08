//! Sparse symmetric (Hermitian) indefinite `L D Lᴴ` by the multifrontal method, with threshold
//! Bunch–Kaufman pivoting and delayed pivots: [`IndefiniteLdlt`].

use num_complex::ComplexFloat;
use num_traits::NumCast;

use super::bunch_kaufman::{factor_front, DBlock};
use super::cholesky::SymbolicCholesky;
use super::{CsrMatrix, Inertia, Ordering};
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
    /// Returns `SparseError::PatternMismatch` if `matrix`'s shape or pattern differs from the
    /// analyzed one, `SparseError::NotSymmetric` if `A ≠ Aᴴ`, and `SparseError::SingularMatrix`
    /// naming a row of `A` where no pivot could be found, the matrix being singular to working
    /// precision. Panics if `options.threshold` is outside `(0, 0.5]`.
    pub fn factor_indefinite<T: Scalar + ComplexFloat>(
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
        let n = self.dim();
        if matrix.shape != (n, n) || matrix.ptr != self.a_ptr || matrix.idx != self.a_idx {
            return Err(SparseError::PatternMismatch);
        }
        if !matrix.is_hermitian() {
            return Err(SparseError::NotSymmetric);
        }
        let c = matrix.permute_symmetric(&self.perm)?;
        let nodes = supernodes(self);
        // Each column's supernode, to find a supernode's parent
        let mut node_of = vec![0; n];
        for (s, (start, end, _)) in nodes.iter().enumerate() {
            node_of[*start..*end].fill(s);
        }
        let parent = self.elimination_tree();

        // Contribution blocks waiting for their parent: (variables, dense square, both triangles)
        let mut pending: Vec<Vec<(Vec<usize>, Vec<T>)>> = vec![Vec::new(); nodes.len()];
        // Where each variable sits in the current front, or usize::MAX
        let mut position = vec![usize::MAX; n];
        let mut fronts = Vec::with_capacity(nodes.len());
        let mut delayed = 0;
        for (s, (start, end, below)) in nodes.iter().enumerate() {
            let contributions = std::mem::take(&mut pending[s]);
            // The front's variables: its own, the children's delayed ones, then the rows below
            let mut vars: Vec<usize> = (*start..*end).collect();
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
            let mut front = vec![T::zero(); m * m];
            // A's lower triangle in the own columns, mirrored: row j of C holds column j's
            //   entries, and those at or below the diagonal belong to this front
            for j in *start..*end {
                let pj = position[j];
                for q in c.ptr[j]..c.ptr[j + 1] {
                    let i = c.idx[q];
                    if i < j {
                        continue;
                    }
                    let pi = position[i];
                    debug_assert_ne!(pi, usize::MAX, "A's entry within the factor's pattern");
                    // c is Hermitian: C(i, j) = conj(C(j, i))
                    let value = c.data[q].conj();
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
            drop(contributions);

            let root = parent[end - 1].is_none();
            debug_assert!(!root || below.is_empty());
            let factored =
                factor_front(&mut front, m, fully_summed, threshold, root).map_err(|error| {
                    match error {
                        SparseError::SingularMatrix { pivot } => SparseError::SingularMatrix {
                            pivot: self.perm[vars[pivot]],
                        },
                        other => other,
                    }
                })?;
            for &v in &vars {
                position[v] = usize::MAX;
            }
            let ordered: Vec<usize> = factored.order.iter().map(|&p| vars[p]).collect();
            let e = factored.eliminated;
            // The contribution block, everything not eliminated, goes to the parent
            if e < m {
                let rest = m - e;
                let mut block = vec![T::zero(); rest * rest];
                for b in 0..rest {
                    for a in 0..rest {
                        block[a + b * rest] = front[(e + a) + (e + b) * m];
                    }
                }
                let p = node_of[parent[end - 1].expect("a front with rows left has a parent")];
                pending[p].push((ordered[e..].to_vec(), block));
            }
            front.truncate(m * e);
            // Inside a 2×2 block, L's entry is zero (the kernel leaves D's there)
            let mut at = 0;
            for block in &factored.blocks {
                if block.size() == 2 {
                    front[(at + 1) + at * m] = T::zero();
                }
                at += block.size();
            }
            fronts.push(Front {
                vars: ordered,
                eliminated: e,
                l: front,
                blocks: factored.blocks,
            });
        }
        Ok(IndefiniteLdlt {
            n,
            perm: self.perm.clone(),
            fronts,
            delayed,
        })
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
        self.analyze_cholesky(Ordering::Amd)?
            .factor_indefinite(self, IndefiniteOptions::default())
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
