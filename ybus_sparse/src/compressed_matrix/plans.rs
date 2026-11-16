//! Kept plans for products and sums: work out a result's pattern once, then compute its values
//! from new operand values with no hashing, sorting or allocation, bit for bit as `matmul` and
//! `try_add` compute them.
//!
//! Patterns follow from the operation, never from the values (see `ops.rs`), so a product or
//! sum of operands whose patterns stay fixed always has the same pattern. Code that forms the
//! same products and sums over and over, such as an optimal power flow's Hessian at every
//! iteration, can plan each one once. A plan records the floating-point operations the
//! operation performs, in the order it performs them, so replaying it gives exactly the same
//! values:
//!
//! - **A product** (`ProductPlan`) records each term `a[p] · b[q]` with the result slot it adds
//!   into, in `matmul`'s traversal order: each row of the left factor, its entries in order,
//!   and the matching row of the right factor in order. Each result entry starts at zero and
//!   adds its terms in that order, as `matmul` does.
//! - **A sum** (`SumPlan`) records each result slot's sources: `a + b` where both operands store
//!   the position, or a copy of the one that does, as `try_add` does.

use std::marker::PhantomData;

use super::{CompressedMatrix, Csr, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// The pattern of a plan's operand or result, which a plan checks its inputs against.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Operand {
    shape: (usize, usize),
    num_stored: usize,
}

impl Operand {
    fn of<T: Scalar, O: Orientation>(matrix: &CompressedMatrix<T, O>) -> Self {
        Operand {
            shape: matrix.shape,
            num_stored: matrix.num_stored(),
        }
    }

    fn check<T: Scalar, O: Orientation>(
        &self,
        matrix: &CompressedMatrix<T, O>,
    ) -> Result<(), SparseError> {
        if Operand::of(matrix) == *self {
            Ok(())
        } else {
            Err(SparseError::PatternMismatch)
        }
    }
}

/// Checks that `values` has `expected` entries.
fn check_len<T>(values: &[T], expected: usize) -> Result<(), SparseError> {
    if values.len() == expected {
        Ok(())
    } else {
        Err(SparseError::LengthMismatch {
            expected,
            found: values.len(),
        })
    }
}

/// The product `A B` of two matrices with fixed patterns, planned by
/// [`CompressedMatrix::plan_product`] so that it can be recomputed from new values of `A` and
/// `B` in O(terms), with no hashing, sorting or allocation.
///
/// The values are bit for bit those of [`CompressedMatrix::matmul`] on the same operands, and
/// the pattern is the same, every reachable position included, even where the terms cancel.
///
/// A plan is tied to its operands' shapes and stored counts, which it checks, and to the
/// orientation `O`.
///
/// ```
/// use ybus_sparse::compressed_matrix::CsrMatrix;
///
/// let a = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 1.0), (0, 1, 2.0), (1, 1, 3.0)])?;
/// let mut b = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 4.0), (1, 0, 5.0), (1, 1, 6.0)])?;
/// let plan = a.plan_product(&b)?;
/// assert_eq!(plan.multiply(&a, &b)?, a.matmul(&b)?);
///
/// // New values for B, with the same pattern
/// b.values_mut().copy_from_slice(&[-1.0, 0.5, 2.0]);
/// assert_eq!(plan.multiply(&a, &b)?, a.matmul(&b)?);
/// # Ok::<(), ybus_sparse::sparse_error::SparseError>(())
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProductPlan<O: Orientation = Csr> {
    a: Operand,
    b: Operand,
    shape: (usize, usize),
    ptr: Vec<usize>,
    idx: Vec<usize>,
    /// Each term as `(result slot, left slot, right slot)`, in `matmul`'s order, where the left
    /// and right factors are the row-compressed arrays `matmul` multiplies: `A` and `B` for CSR,
    /// and `B`'s and `A`'s for CSC (see `swapped`).
    terms: Vec<(usize, usize, usize)>,
    /// Whether the left factor is `B` (CSC), so that each term is `b[l] · a[r]`
    swapped: bool,
    orientation: PhantomData<O>,
}

impl<O: Orientation> ProductPlan<O> {
    /// The number of entries the product stores.
    pub fn num_stored(&self) -> usize {
        self.idx.len()
    }

    /// The number of terms `a[p] · b[q]` a product adds up, the cost of `multiply_values`.
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// The product's shape.
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    /// Writes the product's values into `out`, from the values of `A` and `B`, bit for bit as
    /// `matmul` computes them: each entry starts at zero and adds its terms in `matmul`'s order.
    ///
    /// Returns `SparseError::LengthMismatch` unless `a_values`, `b_values` and `out` have one
    /// value per stored entry of `A`, `B` and the product. On an error, `out` is unchanged.
    pub fn multiply_values<T: Scalar>(
        &self,
        a_values: &[T],
        b_values: &[T],
        out: &mut [T],
    ) -> Result<(), SparseError> {
        check_len(a_values, self.a.num_stored)?;
        check_len(b_values, self.b.num_stored)?;
        check_len(out, self.num_stored())?;
        out.fill(T::zero());
        if self.swapped {
            for &(slot, left, right) in &self.terms {
                out[slot] = out[slot] + b_values[left] * a_values[right];
            }
        } else {
            for &(slot, left, right) in &self.terms {
                out[slot] = out[slot] + a_values[left] * b_values[right];
            }
        }
        Ok(())
    }

    /// The product `A B` as a new matrix, bit for bit `a.matmul(b)`.
    ///
    /// Returns `SparseError::PatternMismatch` unless `a` and `b` have the shapes and stored
    /// counts the plan was made for.
    pub fn multiply<T: Scalar>(
        &self,
        a: &CompressedMatrix<T, O>,
        b: &CompressedMatrix<T, O>,
    ) -> Result<CompressedMatrix<T, O>, SparseError> {
        self.a.check(a)?;
        self.b.check(b)?;
        let mut out = vec![T::zero(); self.num_stored()];
        self.multiply_values(&a.data, &b.data, &mut out)?;
        Ok(self.matrix(out))
    }

    /// The product's pattern holding `values`, one per stored entry.
    fn matrix<T: Scalar>(&self, values: Vec<T>) -> CompressedMatrix<T, O> {
        CompressedMatrix::from_parts_unchecked(
            self.shape,
            self.ptr.clone(),
            self.idx.clone(),
            values,
        )
    }
}

/// Where a sum's result slot gets its value: `a + b`, or a copy of one operand's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Source {
    Both(usize, usize),
    Left(usize),
    Right(usize),
}

/// The sum `A + B` of two matrices with fixed patterns, planned by
/// [`CompressedMatrix::plan_sum`] so that it can be recomputed from new values in O(entries),
/// with no merging or allocation.
///
/// The result stores the union of the two patterns, and its values are bit for bit those of
/// [`CompressedMatrix::try_add`]: `a + b` where both store the position, and a copy of the one
/// that does elsewhere.
///
/// ```
/// use ybus_sparse::compressed_matrix::CsrMatrix;
///
/// let a = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 1.0), (1, 1, 3.0)])?;
/// let b = CsrMatrix::from_triplets(2, 2, vec![(0, 0, 4.0), (0, 1, 5.0)])?;
/// let plan = a.plan_sum(&b)?;
/// assert_eq!(plan.add(&a, &b)?, a.try_add(&b)?);
/// # Ok::<(), ybus_sparse::sparse_error::SparseError>(())
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SumPlan<O: Orientation = Csr> {
    a: Operand,
    b: Operand,
    shape: (usize, usize),
    ptr: Vec<usize>,
    idx: Vec<usize>,
    sources: Vec<Source>,
    orientation: PhantomData<O>,
}

impl<O: Orientation> SumPlan<O> {
    /// The number of entries the sum stores.
    pub fn num_stored(&self) -> usize {
        self.idx.len()
    }

    /// The sum's shape.
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    /// Writes the sum's values into `out`, from the values of `A` and `B`, bit for bit as
    /// `try_add` computes them.
    ///
    /// Returns `SparseError::LengthMismatch` unless `a_values`, `b_values` and `out` have one
    /// value per stored entry of `A`, `B` and the sum. On an error, `out` is unchanged.
    pub fn add_values<T: Scalar>(
        &self,
        a_values: &[T],
        b_values: &[T],
        out: &mut [T],
    ) -> Result<(), SparseError> {
        check_len(a_values, self.a.num_stored)?;
        check_len(b_values, self.b.num_stored)?;
        check_len(out, self.num_stored())?;
        for (out, source) in out.iter_mut().zip(&self.sources) {
            *out = match *source {
                Source::Both(p, q) => a_values[p] + b_values[q],
                Source::Left(p) => a_values[p],
                Source::Right(q) => b_values[q],
            };
        }
        Ok(())
    }

    /// The sum `A + B` as a new matrix, bit for bit `a.try_add(b)`.
    ///
    /// Returns `SparseError::PatternMismatch` unless `a` and `b` have the shapes and stored
    /// counts the plan was made for.
    pub fn add<T: Scalar>(
        &self,
        a: &CompressedMatrix<T, O>,
        b: &CompressedMatrix<T, O>,
    ) -> Result<CompressedMatrix<T, O>, SparseError> {
        self.a.check(a)?;
        self.b.check(b)?;
        let mut out = vec![T::zero(); self.num_stored()];
        self.add_values(&a.data, &b.data, &mut out)?;
        Ok(CompressedMatrix::from_parts_unchecked(
            self.shape,
            self.ptr.clone(),
            self.idx.clone(),
            out,
        ))
    }
}

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// Plans the product `self · other` for operands with these patterns: see
    /// [`ProductPlan`]. Only the patterns are used, so `other` may hold values of another type,
    /// and the plan multiplies any values with the same patterns.
    ///
    /// Costs about as much as one `matmul`, and the plan holds three indices per term.
    /// Returns `SparseError::ShapeMismatch` unless `self` has as many columns as `other` has
    /// rows.
    pub fn plan_product<U: Scalar>(
        &self,
        other: &CompressedMatrix<U, O>,
    ) -> Result<ProductPlan<O>, SparseError> {
        if self.shape.1 != other.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: other.shape,
            });
        }
        // matmul multiplies the row-compressed arrays of the left and right factors: A and B
        //   for CSR, and Bᵀ and Aᵀ (whose CSR arrays are B's and A's CSC arrays) for CSC
        let swapped = !O::ROW_MAJOR;
        let (left_ptr, left_idx, right_ptr, right_idx, right_minor) = if swapped {
            (
                &other.ptr,
                &other.idx,
                &self.ptr,
                &self.idx,
                self.minor_dim(),
            )
        } else {
            (
                &self.ptr,
                &self.idx,
                &other.ptr,
                &other.idx,
                other.minor_dim(),
            )
        };
        let mut ptr = Vec::with_capacity(left_ptr.len());
        ptr.push(0);
        let mut idx = Vec::new();
        let mut terms = Vec::new();
        // For each minor index reached in the current row, its slot in the result
        let mut slot_of = vec![usize::MAX; right_minor];
        let mut row_terms: Vec<(usize, usize, usize)> = Vec::new();
        for major in 0..left_ptr.len() - 1 {
            // The terms in matmul's order: the left row's entries, then the right rows'
            row_terms.clear();
            let start = idx.len();
            let left_start = left_ptr[major];
            for (offset, &k) in left_idx[left_start..left_ptr[major + 1]].iter().enumerate() {
                let right_start = right_ptr[k];
                let right_row = &right_idx[right_start..right_ptr[k + 1]];
                for (right_offset, &minor) in right_row.iter().enumerate() {
                    if slot_of[minor] == usize::MAX {
                        slot_of[minor] = 0;
                        idx.push(minor);
                    }
                    row_terms.push((minor, left_start + offset, right_start + right_offset));
                }
            }
            // matmul sorts each row by column; the slots follow that order
            idx[start..].sort_unstable();
            for (offset, &minor) in idx[start..].iter().enumerate() {
                slot_of[minor] = start + offset;
            }
            for &(minor, p, q) in &row_terms {
                terms.push((slot_of[minor], p, q));
            }
            for &minor in &idx[start..] {
                slot_of[minor] = usize::MAX;
            }
            ptr.push(idx.len());
        }
        Ok(ProductPlan {
            a: Operand::of(self),
            b: Operand::of(other),
            shape: (self.shape.0, other.shape.1),
            ptr,
            idx,
            terms,
            swapped,
            orientation: PhantomData,
        })
    }

    /// Plans the sum `self + other` for operands with these patterns: see [`SumPlan`]. Only the
    /// patterns are used, so `other` may hold values of another type.
    ///
    /// Returns `SparseError::ShapeMismatch` if the shapes differ.
    pub fn plan_sum<U: Scalar>(
        &self,
        other: &CompressedMatrix<U, O>,
    ) -> Result<SumPlan<O>, SparseError> {
        if self.shape != other.shape {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: other.shape,
            });
        }
        let mut ptr = Vec::with_capacity(self.ptr.len());
        ptr.push(0);
        let mut idx = Vec::with_capacity(self.idx.len().max(other.idx.len()));
        let mut sources = Vec::with_capacity(idx.capacity());
        // The same walk as ops.rs's merge, recording sources instead of values
        for major in 0..self.major_dim() {
            let (mut p, p_end) = (self.ptr[major], self.ptr[major + 1]);
            let (mut q, q_end) = (other.ptr[major], other.ptr[major + 1]);
            while p < p_end || q < q_end {
                let take_left = q == q_end || (p < p_end && self.idx[p] <= other.idx[q]);
                let take_right = p == p_end || (q < q_end && other.idx[q] <= self.idx[p]);
                match (take_left, take_right) {
                    (true, true) => {
                        idx.push(self.idx[p]);
                        sources.push(Source::Both(p, q));
                        p += 1;
                        q += 1;
                    }
                    (true, false) => {
                        idx.push(self.idx[p]);
                        sources.push(Source::Left(p));
                        p += 1;
                    }
                    _ => {
                        idx.push(other.idx[q]);
                        sources.push(Source::Right(q));
                        q += 1;
                    }
                }
            }
            ptr.push(idx.len());
        }
        Ok(SumPlan {
            a: Operand::of(self),
            b: Operand::of(other),
            shape: self.shape,
            ptr,
            idx,
            sources,
            orientation: PhantomData,
        })
    }
}
