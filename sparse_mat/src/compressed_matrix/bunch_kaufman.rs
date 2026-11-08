//! Dense symmetric (Hermitian) indefinite `L D Lᴴ` with threshold Bunch–Kaufman pivoting: the
//! kernel a multifrontal factorization runs on each front, and [`DenseLdlt`], the same kernel
//! on a whole dense matrix. See [`DenseLdlt`] for the method and its references.

use num_complex::ComplexFloat;
use num_traits::{Float, NumCast, Zero};

use crate::compressed_matrix::Inertia;
use crate::sparse_error::SparseError;
use crate::Scalar;

/// One block of `D` in a pivoted `L D Lᴴ` factorization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DBlock<T> {
    /// A 1×1 pivot: `d`, real (with zero imaginary part for a complex element type)
    One(T),
    /// A 2×2 pivot, `[[a, conj(b)], [b, c]]`, with `a` and `c` real
    Two {
        /// The upper-left entry
        a: T,
        /// The lower-left entry
        b: T,
        /// The lower-right entry
        c: T,
    },
}

impl<T: Scalar + ComplexFloat> DBlock<T> {
    /// The block's order: 1 or 2.
    pub fn size(&self) -> usize {
        match self {
            DBlock::One(_) => 1,
            DBlock::Two { .. } => 2,
        }
    }

    /// The block's eigenvalues' signs, as `(positive, negative)`. A 2×2 block with a negative
    /// determinant has one of each; with a positive one, both take the sign of its diagonal.
    pub(crate) fn signs(&self) -> (usize, usize) {
        let zero = <T::Real as Zero>::zero();
        match *self {
            DBlock::One(d) => {
                if d.re() > zero {
                    (1, 0)
                } else {
                    (0, 1)
                }
            }
            DBlock::Two { a, b, c } => {
                let det = a.re() * c.re() - magnitude_squared(b);
                if det < zero {
                    (1, 1)
                } else if a.re() > zero {
                    (2, 0)
                } else {
                    (0, 2)
                }
            }
        }
    }
}

/// `|z|²`, `re² + im²`
fn magnitude_squared<T: ComplexFloat>(z: T) -> T::Real {
    (z * z.conj()).re()
}

/// A real value as the element type
fn from_real<T: Scalar + ComplexFloat>(value: T::Real) -> T {
    <T as NumCast>::from(value).expect("a real value converts to its own element type")
}

/// What [`factor_front`] did to a front.
#[derive(Clone, Debug)]
pub(crate) struct FrontFactor<T> {
    /// The front's positions in their final order: `order[i]` is the position now at `i`.
    /// The first `eliminated` are the pivots, in pivot order; then the delayed fully summed
    /// positions; then the rest.
    pub(crate) order: Vec<usize>,
    /// How many positions were eliminated
    pub(crate) eliminated: usize,
    /// `D`'s blocks, in pivot order
    pub(crate) blocks: Vec<DBlock<T>>,
}

/// Factors as much of a dense Hermitian front as the threshold test allows.
///
/// `a` is the front, `m × m`, column-major, both triangles stored and equal (`a_ji =
/// conj(a_ij)`). Its first `k` positions are fully summed, and only they may be pivots. On
/// return, in the final order of `FrontFactor::order`:
///
/// - the strictly lower part of the first `eliminated` columns holds `L`, except that a 2×2
///   block's lower-left entry is `D`'s, not `L`'s (whose entry there is zero);
/// - the trailing `m − eliminated` square holds the Schur complement, both triangles, still
///   exactly Hermitian (only the lower triangle is computed, then mirrored).
///
/// With `complete`, every position must be fully summed (`k = m`), and a step where no
/// candidate passes means the matrix is singular (see [`DenseLdlt`]), an error.
/// Without it, the factorization stops when no candidate passes, and the rest are delayed.
///
/// `scratch` is working space, reused across calls so a pivot allocates nothing.
///
/// Returns `SparseError::SingularMatrix` naming the front position where no pivot was found
/// with `complete`, or of a pivot that is zero or NaN.
pub(crate) fn factor_front<T: Scalar + ComplexFloat>(
    a: &mut [T],
    m: usize,
    k: usize,
    threshold: T::Real,
    complete: bool,
    scratch: &mut Vec<T>,
) -> Result<FrontFactor<T>, SparseError> {
    debug_assert_eq!(a.len(), m * m);
    debug_assert!(k <= m && (!complete || k == m));
    let zero = <T::Real as Zero>::zero();
    let mut order: Vec<usize> = (0..m).collect();
    let mut blocks = Vec::new();
    let at = |i: usize, j: usize| i + j * m;
    let mut p = 0;
    while p < k {
        let mut chosen = None;
        for q in p..k {
            if let Some(pivot) = threshold_pivot(a, m, p, k, q, threshold) {
                chosen = Some(pivot);
                break;
            }
        }
        let pivot = match chosen {
            Some(pivot) => pivot,
            None if complete => return Err(SparseError::SingularMatrix { pivot: order[p] }),
            None => break,
        };
        match pivot {
            Pivot::One(q) => {
                swap(a, m, &mut order, p, q);
                let d = a[at(p, p)].re();
                if d == zero || Float::is_nan(d) {
                    return Err(SparseError::SingularMatrix { pivot: order[p] });
                }
                eliminate_one(a, m, p, d, scratch);
                blocks.push(DBlock::One(from_real(d)));
                p += 1;
            }
            Pivot::Two(q, r) => {
                swap(a, m, &mut order, p, q);
                // r may have been at p, and is now where q was
                let r = if r == p { q } else { r };
                swap(a, m, &mut order, p + 1, r);
                let (d11, d21, d22) = (a[at(p, p)].re(), a[at(p + 1, p)], a[at(p + 1, p + 1)].re());
                let det = d11 * d22 - magnitude_squared(d21);
                if det == zero || Float::is_nan(det) {
                    return Err(SparseError::SingularMatrix { pivot: order[p] });
                }
                eliminate_two(a, m, p, (d11, d21, d22, det), scratch);
                blocks.push(DBlock::Two {
                    a: from_real(d11),
                    b: d21,
                    c: from_real(d22),
                });
                p += 2;
            }
        }
    }
    Ok(FrontFactor {
        order,
        eliminated: p,
        blocks,
    })
}

/// A chosen pivot, by current positions: 1×1 at `q`, or 2×2 on `q` and `r`
#[derive(Clone, Copy, Debug)]
enum Pivot {
    One(usize),
    Two(usize, usize),
}

/// The largest magnitude in column `j` over rows `from..to`, skipping the rows in `skip`, with
/// its row (the first, on ties), or `None` if there are no such rows
fn column_max<T: ComplexFloat>(
    a: &[T],
    m: usize,
    j: usize,
    from: usize,
    to: usize,
    skip: &[usize],
) -> (T::Real, Option<usize>) {
    let mut best = <T::Real as Zero>::zero();
    let mut row = None;
    for i in from..to {
        if skip.contains(&i) {
            continue;
        }
        let value = a[i + j * m].abs();
        if row.is_none() || value > best {
            best = value;
            row = Some(i);
        }
    }
    (best, row)
}

/// The threshold test on candidate `q`, with positions `p..` still to be factored and `p..k`
/// fully summed: a 1×1 pivot, a 2×2 one with the fully summed candidate whose entry in `q`'s
/// column is largest, or neither
fn threshold_pivot<T: Scalar + ComplexFloat>(
    a: &[T],
    m: usize,
    p: usize,
    k: usize,
    q: usize,
    u: T::Real,
) -> Option<Pivot> {
    let zero = <T::Real as Zero>::zero();
    let at = |i: usize, j: usize| i + j * m;
    let (gamma_q, _) = column_max(a, m, q, p, m, &[q]);
    let a_qq = a[at(q, q)].re();
    if a_qq != zero && Float::abs(a_qq) >= u * gamma_q {
        return Some(Pivot::One(q));
    }
    let (_, r) = column_max(a, m, q, p, k, &[q]);
    let r = r?;
    let (gamma_q, _) = column_max(a, m, q, p, m, &[q, r]);
    let (gamma_r, _) = column_max(a, m, r, p, m, &[q, r]);
    let a_rr = a[at(r, r)].re();
    let b = a[at(r, q)].abs();
    let det = Float::abs(a_qq * a_rr - magnitude_squared(a[at(r, q)]));
    if det == zero || Float::is_nan(det) {
        return None;
    }
    // |P⁻¹| [γ'_q; γ'_r] ≤ [1/u; 1/u], with |P⁻¹| = [|a_rr| b; b |a_qq|] / |det|
    let bound = det / u;
    if Float::abs(a_rr) * gamma_q + b * gamma_r <= bound
        && b * gamma_q + Float::abs(a_qq) * gamma_r <= bound
    {
        Some(Pivot::Two(q, r))
    } else {
        None
    }
}

/// Swaps positions `i` and `j` symmetrically: rows (all columns, so `L`'s stored columns move
/// with them) and columns, and the order
fn swap<T: Copy>(a: &mut [T], m: usize, order: &mut [usize], i: usize, j: usize) {
    if i == j {
        return;
    }
    for col in 0..m {
        a.swap(i + col * m, j + col * m);
    }
    for row in 0..m {
        a.swap(row + i * m, row + j * m);
    }
    order.swap(i, j);
}

/// Eliminates the 1×1 pivot `d` at position `p`: `L`'s column below it is `a_ip / d`, and the
/// trailing square loses `a_ip conj(l_j)`, computed on the lower triangle and mirrored
fn eliminate_one<T: Scalar + ComplexFloat>(
    a: &mut [T],
    m: usize,
    p: usize,
    d: T::Real,
    scratch: &mut Vec<T>,
) {
    let at = |i: usize, j: usize| i + j * m;
    let d_t: T = from_real(d);
    // w, the pivot's column below it, then l = w / d, in the caller's scratch
    let k = m - p - 1;
    scratch.clear();
    scratch.extend((p + 1..m).map(|i| a[at(i, p)]));
    for ii in 0..k {
        let value = scratch[ii] / d_t;
        scratch.push(value);
    }
    let (w, l) = scratch.split_at(k);
    for (jj, j) in (p + 1..m).enumerate() {
        let conj_l = l[jj].conj();
        for (ii, i) in (j..m).enumerate().map(|(n, i)| (n + jj, i)) {
            a[at(i, j)] = a[at(i, j)] - w[ii] * conj_l;
        }
        a[at(j, j)] = from_real(a[at(j, j)].re());
        for i in j + 1..m {
            a[at(j, i)] = a[at(i, j)].conj();
        }
    }
    for (ii, i) in (p + 1..m).enumerate() {
        a[at(i, p)] = l[ii];
        a[at(p, i)] = T::zero();
    }
    a[at(p, p)] = d_t;
}

/// Eliminates the 2×2 pivot `[[d11, conj(d21)], [d21, d22]]` at positions `p` and `p + 1`:
/// `L`'s two columns below it are `[a_ip, a_i,p+1] D⁻¹`, and the trailing square loses
/// `a_ip conj(l_j,p) + a_i,p+1 conj(l_j,p+1)`, computed on the lower triangle and mirrored
fn eliminate_two<T: Scalar + ComplexFloat>(
    a: &mut [T],
    m: usize,
    p: usize,
    (d11, d21, d22, det): (T::Real, T, T::Real, T::Real),
    scratch: &mut Vec<T>,
) {
    let at = |i: usize, j: usize| i + j * m;
    let (d11_t, d22_t, det_t): (T, T, T) = (from_real(d11), from_real(d22), from_real(det));
    // x and y, the pivots' columns below them, then [l0, l1] = [x, y] D⁻¹ with
    //   D⁻¹ = [[d22, −conj(d21)], [−d21, d11]] / det, in the caller's scratch
    let k = m - p - 2;
    scratch.clear();
    scratch.extend((p + 2..m).map(|i| a[at(i, p)]));
    scratch.extend((p + 2..m).map(|i| a[at(i, p + 1)]));
    for ii in 0..k {
        let (x, y) = (scratch[ii], scratch[k + ii]);
        scratch.push((x * d22_t - y * d21) / det_t);
    }
    for ii in 0..k {
        let (x, y) = (scratch[ii], scratch[k + ii]);
        scratch.push((y * d11_t - x * d21.conj()) / det_t);
    }
    let (xy, l) = scratch.split_at(2 * k);
    let (x, y) = xy.split_at(k);
    let (l0, l1) = l.split_at(k);
    for (jj, j) in (p + 2..m).enumerate() {
        let (c0, c1) = (l0[jj].conj(), l1[jj].conj());
        for i in j..m {
            let ii = i - p - 2;
            a[at(i, j)] = a[at(i, j)] - x[ii] * c0 - y[ii] * c1;
        }
        a[at(j, j)] = from_real(a[at(j, j)].re());
        for i in j + 1..m {
            a[at(j, i)] = a[at(i, j)].conj();
        }
    }
    for (ii, i) in (p + 2..m).enumerate() {
        a[at(i, p)] = l0[ii];
        a[at(i, p + 1)] = l1[ii];
        a[at(p, i)] = T::zero();
        a[at(p + 1, i)] = T::zero();
    }
    a[at(p, p)] = d11_t;
    a[at(p + 1, p)] = d21;
    a[at(p, p + 1)] = d21.conj();
    a[at(p + 1, p + 1)] = d22_t;
}

/// A dense symmetric (Hermitian) indefinite factorization `P A Pᵀ = L D Lᴴ` with threshold
/// Bunch–Kaufman pivoting, from [`DenseLdlt::factor`]: `L` unit lower triangular, `D` block
/// diagonal with 1×1 and 2×2 blocks.
///
/// `D` is block diagonal, with 1×1 and 2×2 blocks ([`DBlock`]), and a symmetric permutation `P`
/// chooses them as the factorization goes, so `P A Pᵀ = L D Lᴴ` with `L` unit lower
/// triangular. A 2×2 pivot handles what no 1×1 pivot can, such as a zero diagonal next to a
/// large off-diagonal entry (`[[0, 1], [1, 0]]`), without any regularization.
///
/// **The threshold test** (Duff and Reid, 1983; MA57, Duff 2004), with a threshold `u` in
/// `(0, 0.5]`. A candidate `q` is judged against `γ_q`, its column's largest off-diagonal
/// magnitude in the part still to be factored:
///
/// - **a 1×1 pivot** if `|a_qq| ≥ u γ_q` (and `a_qq ≠ 0`);
/// - **otherwise a 2×2 pivot** on `q` and `r`, `r` the candidate whose entry `|a_rq|` is
///   largest, if `|P⁻¹| [γ'_q; γ'_r] ≤ [1/u; 1/u]` entrywise. Here `P` is the 2×2 block, and
///   `γ'` are the two columns' largest magnitudes outside it;
/// - **otherwise** the next candidate is tried.
///
/// Either test bounds every entry of `L` by `1/u`. The candidates are tried in order, starting
/// afresh after each pivot, so a matrix whose diagonal passes keeps the order it came in, and
/// with it a fill-reducing ordering's sparsity.
///
/// **A pivot always exists when every position is fully summed** and the matrix is nonsingular,
/// because `u ≤ 1/2` (Duff and Reid's argument). Let `μ = |a_rq|` be the largest off-diagonal
/// magnitude left. If neither `q` nor `r` passes as a 1×1 pivot, both diagonals are below `u μ`,
/// so `|det P| ≥ μ²(1 − u²)` and the 2×2 test's left side is at most `(1 + u)/(1 − u²) =
/// 1/(1 − u) ≤ 1/u`. So a whole matrix factors without the unbounded `L` that Bunch and Kaufman's
/// own rule (1977) allows (Ashcraft, Grimes and Lewis, 1998), and failing to find a pivot means a
/// singular matrix. In a front of a sparse factorization, a candidate that fails is instead
/// delayed to the parent front, whose larger matrix may offer a better pivot (Duff and Reid).
///
/// **References:**
///
/// - J. R. Bunch and L. Kaufman, "Some stable methods for calculating inertia and solving
///   symmetric linear systems," *Math. Comp.* 31(137):163–179, 1977,
///   doi:10.1090/S0025-5718-1977-0428694-0: the 1×1 and 2×2 pivots, and the inertia from
///   them.
/// - I. S. Duff and J. K. Reid, "The multifrontal solution of indefinite sparse symmetric linear
///   equations," *ACM Trans. Math. Softw.* 9(3):302–325, 1983, doi:10.1145/356044.356047.
/// - C. Ashcraft, R. G. Grimes and J. G. Lewis, "Accurate symmetric indefinite linear equation
///   solvers," *SIAM J. Matrix Anal. Appl.* 20(2):513–561, 1998, doi:10.1137/S0895479896296921,
///   on why `L`'s entries need bounding, which Bunch and Kaufman's rule alone doesn't do.
/// - I. S. Duff, "MA57—a code for the solution of sparse symmetric definite and indefinite
///   systems," *ACM Trans. Math. Softw.* 30(2):118–144, 2004, doi:10.1145/992200.992202.
///
/// It's the kernel a sparse multifrontal factorization runs on each dense front, applied here to
/// a whole matrix, so every position is fully summed and none is delayed. Its cost is `O(n³)`.
#[derive(Clone, Debug)]
pub struct DenseLdlt<T: Scalar = f64> {
    n: usize,
    /// `L`, column-major, unit diagonal and zeros above it stored
    l: Vec<T>,
    blocks: Vec<DBlock<T>>,
    perm: Vec<usize>,
}

impl<T: Scalar + ComplexFloat> DenseLdlt<T> {
    /// The threshold MA57 uses by default, 0.01: every entry of `L` from a threshold pivot is
    /// at most 100 in magnitude.
    pub const DEFAULT_THRESHOLD: f64 = 0.01;

    /// Factors the dense Hermitian matrix `rows` (row-major, each row as long as there are rows)
    /// with threshold `threshold`, which must be in `(0, 0.5]`. Every entry of `L` is at most
    /// `1/threshold` in magnitude: a larger threshold keeps them smaller and takes more 2×2
    /// pivots, a smaller one keeps more of the diagonal order.
    ///
    /// Returns `SparseError::RaggedRows` or `SparseError::NotSquare` for rows of the wrong
    /// lengths, `SparseError::NotSymmetric` unless `A = Aᴴ` exactly (a complex diagonal entry
    /// must be real), and `SparseError::SingularMatrix` naming a row of `A` where no pivot could
    /// be found, the matrix being singular to working precision.
    /// Panics if `threshold` is outside `(0, 0.5]`.
    pub fn factor(rows: &[Vec<T>], threshold: f64) -> Result<DenseLdlt<T>, SparseError> {
        assert!(
            threshold > 0.0 && threshold <= 0.5,
            "the threshold must be in (0, 0.5], not {threshold}"
        );
        let n = rows.len();
        for (i, row) in rows.iter().enumerate() {
            if row.len() != n {
                return Err(if i == 0 {
                    SparseError::NotSquare {
                        shape: (n, row.len()),
                    }
                } else if row.len() != rows[0].len() {
                    SparseError::RaggedRows {
                        row: i,
                        expected: rows[0].len(),
                        found: row.len(),
                    }
                } else {
                    SparseError::NotSquare {
                        shape: (n, row.len()),
                    }
                });
            }
        }
        for (i, row) in rows.iter().enumerate() {
            for (j, &value) in row.iter().enumerate().take(i + 1) {
                if value != rows[j][i].conj() {
                    return Err(SparseError::NotSymmetric);
                }
            }
        }
        let mut a = vec![T::zero(); n * n];
        for (i, row) in rows.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                a[i + j * n] = value;
            }
        }
        let threshold: T::Real = NumCast::from(threshold).expect("the threshold is representable");
        let front = factor_front(&mut a, n, n, threshold, true, &mut Vec::new())?;
        debug_assert_eq!(front.eliminated, n);
        // L: the strictly lower part, with each 2×2 block's lower-left entry (D's) cleared, and
        //   the unit diagonal
        let mut start = 0;
        for block in &front.blocks {
            if block.size() == 2 {
                a[start + 1 + start * n] = T::zero();
            }
            start += block.size();
        }
        for j in 0..n {
            for i in 0..j {
                a[i + j * n] = T::zero();
            }
            a[j + j * n] = T::one();
        }
        Ok(DenseLdlt {
            n,
            l: a,
            blocks: front.blocks,
            perm: front.order,
        })
    }

    /// `L`, unit lower triangular, as dense rows of `P A Pᵀ`.
    pub fn l(&self) -> Vec<Vec<T>> {
        (0..self.n)
            .map(|i| (0..self.n).map(|j| self.l[i + j * self.n]).collect())
            .collect()
    }

    /// `D`'s blocks, in pivot order: block `b` covers the rows after those of the blocks before
    /// it.
    pub fn d_blocks(&self) -> &[DBlock<T>] {
        &self.blocks
    }

    /// The ordering `P`: row and column `i` of `P A Pᵀ` are row and column `permutation()[i]`
    /// of `A`.
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// `A`'s inertia, from `D`'s blocks (Sylvester's law of inertia). `A` is nonsingular, since
    /// a zero pivot fails the factorization.
    pub fn inertia(&self) -> Inertia {
        let (positive, negative) = self
            .blocks
            .iter()
            .map(DBlock::signs)
            .fold((0, 0), |(p, n), (bp, bn)| (p + bp, n + bn));
        Inertia { positive, negative }
    }

    /// Solves `A x = b`: permutes `b`, substitutes forward with `L`, solves with each block of
    /// `D`, substitutes back with `Lᴴ`, and permutes back. Returns `SparseError::ShapeMismatch`
    /// unless `b` has one entry per row of `A`.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, SparseError> {
        let n = self.n;
        if b.len() != n {
            return Err(SparseError::ShapeMismatch {
                left: (n, n),
                right: (b.len(), 1),
            });
        }
        let mut y: Vec<T> = self.perm.iter().map(|&row| b[row]).collect();
        for j in 0..n {
            let y_j = y[j];
            for (i, y_i) in y.iter_mut().enumerate().skip(j + 1) {
                *y_i = *y_i - self.l[i + j * n] * y_j;
            }
        }
        let mut start = 0;
        for block in &self.blocks {
            match *block {
                DBlock::One(d) => y[start] = y[start] / d,
                DBlock::Two { a, b, c } => {
                    let det: T = from_real(a.re() * c.re() - magnitude_squared(b));
                    let (y0, y1) = (y[start], y[start + 1]);
                    y[start] = (c * y0 - b.conj() * y1) / det;
                    y[start + 1] = (a * y1 - b * y0) / det;
                }
            }
            start += block.size();
        }
        for j in (0..n).rev() {
            let mut y_j = y[j];
            for (i, &y_i) in y.iter().enumerate().skip(j + 1) {
                y_j = y_j - self.l[i + j * n].conj() * y_i;
            }
            y[j] = y_j;
        }
        let mut x = vec![T::zero(); n];
        for (&row, &value) in self.perm.iter().zip(&y) {
            x[row] = value;
        }
        Ok(x)
    }
}
