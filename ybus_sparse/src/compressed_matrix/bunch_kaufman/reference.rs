//! TEMPORARY, test-only: the kernel as it was before the lower-triangle-only and blocked
//! versions (O4 target 4), kept to check them against bit for bit. Removed with the record;
//! the original is in history at `67d936b`'s `bunch_kaufman.rs`.

use num_complex::ComplexFloat;
use num_traits::{Float, Zero};

use super::{from_real, magnitude_squared, DBlock, FrontFactor, Pivot};
use crate::sparse_error::SparseError;
use crate::Scalar;

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
