//! Arithmetic on `CompressedMatrix`: checked methods returning `Result`, the operator impls
//! that panic on incompatible shapes, and complex parts.
//!
//! Patterns follow from the operation, never from the values: a sum or difference stores the
//! union of its operands' patterns (cancellations included), an element-wise product their
//! intersection, and scaling, negation and `map` keep the pattern, even where a result is zero.
//! Binary operations take operands of the same orientation.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_complex::{Complex, ComplexFloat};

use super::{CompressedMatrix, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// Merges two matrices of the same shape and orientation over the union of their patterns:
/// `both(a, b)` where both store an entry, `only_left(a)` or `only_right(b)` where one does.
fn merge<A: Scalar, B: Scalar, C: Scalar, O: Orientation>(
    left: &CompressedMatrix<A, O>,
    right: &CompressedMatrix<B, O>,
    both: impl Fn(A, B) -> C,
    only_left: impl Fn(A) -> C,
    only_right: impl Fn(B) -> C,
) -> Result<CompressedMatrix<C, O>, SparseError> {
    if left.shape != right.shape {
        return Err(SparseError::ShapeMismatch {
            left: left.shape,
            right: right.shape,
        });
    }
    let mut ptr = Vec::with_capacity(left.ptr.len());
    ptr.push(0);
    let mut idx = Vec::with_capacity(left.idx.len().max(right.idx.len()));
    let mut data = Vec::with_capacity(idx.capacity());
    for major in 0..left.major_dim() {
        let (mut p, p_end) = (left.ptr[major], left.ptr[major + 1]);
        let (mut q, q_end) = (right.ptr[major], right.ptr[major + 1]);
        while p < p_end || q < q_end {
            let take_left = q == q_end || (p < p_end && left.idx[p] <= right.idx[q]);
            let take_right = p == p_end || (q < q_end && right.idx[q] <= left.idx[p]);
            match (take_left, take_right) {
                (true, true) => {
                    idx.push(left.idx[p]);
                    data.push(both(left.data[p], right.data[q]));
                    p += 1;
                    q += 1;
                }
                (true, false) => {
                    idx.push(left.idx[p]);
                    data.push(only_left(left.data[p]));
                    p += 1;
                }
                _ => {
                    idx.push(right.idx[q]);
                    data.push(only_right(right.data[q]));
                    q += 1;
                }
            }
        }
        ptr.push(idx.len());
    }
    Ok(CompressedMatrix::from_parts_unchecked(
        left.shape, ptr, idx, data,
    ))
}

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// Element-wise sum of two matrices, returned as a new matrix storing the union of their
    /// patterns; entries that cancel to exactly zero stay stored.
    ///
    /// Returns `SparseError::ShapeMismatch` if the shapes differ; `+` panics in that case
    /// instead.
    pub fn try_add(&self, other: &CompressedMatrix<T, O>) -> Result<Self, SparseError> {
        merge(self, other, |a, b| a + b, |a| a, |b| b)
    }

    /// Element-wise difference `self - other`, returned as a new matrix storing the union of
    /// their patterns; entries that cancel to exactly zero stay stored.
    ///
    /// Returns `SparseError::ShapeMismatch` if the shapes differ; `-` panics in that case
    /// instead. Arithmetic follows `T`'s own rules, so for unsigned integers subtracting a
    /// larger value from a smaller one overflows.
    pub fn try_sub(&self, other: &CompressedMatrix<T, O>) -> Result<Self, SparseError> {
        merge(self, other, |a, b| a - b, |a| a, |b| T::zero() - b)
    }

    /// Multiplies every entry by `factor`, returned as a new matrix with the same pattern (a
    /// zero factor leaves every entry stored as zero).
    pub fn scale(&self, factor: T) -> Self {
        self.map(|value| value * factor)
    }

    /// Element-wise (Hadamard) product of two matrices, returned as a new matrix storing the
    /// intersection of their patterns, the only positions that can be non-zero.
    ///
    /// Returns `SparseError::ShapeMismatch` if the shapes differ.
    pub fn hadamard(&self, other: &CompressedMatrix<T, O>) -> Result<Self, SparseError> {
        if self.shape != other.shape {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: other.shape,
            });
        }
        let mut ptr = Vec::with_capacity(self.ptr.len());
        ptr.push(0);
        let (mut idx, mut data) = (vec![], vec![]);
        for major in 0..self.major_dim() {
            let (mut p, p_end) = (self.ptr[major], self.ptr[major + 1]);
            let (mut q, q_end) = (other.ptr[major], other.ptr[major + 1]);
            while p < p_end && q < q_end {
                match self.idx[p].cmp(&other.idx[q]) {
                    std::cmp::Ordering::Less => p += 1,
                    std::cmp::Ordering::Greater => q += 1,
                    std::cmp::Ordering::Equal => {
                        idx.push(self.idx[p]);
                        data.push(self.data[p] * other.data[q]);
                        p += 1;
                        q += 1;
                    }
                }
            }
            ptr.push(idx.len());
        }
        Ok(Self::from_parts_unchecked(self.shape, ptr, idx, data))
    }

    /// Applies `f` to every stored entry, returned as a new matrix with the same pattern. The
    /// result may have a different element type.
    ///
    /// `f` sees stored entries only, never the implicit zeros, so `map(|x| x + 1)` does not
    /// turn the zeros into ones. Results that are zero stay stored.
    pub fn map<U: Scalar>(&self, f: impl Fn(T) -> U) -> CompressedMatrix<U, O> {
        CompressedMatrix::from_parts_unchecked(
            self.shape,
            self.ptr.clone(),
            self.idx.clone(),
            self.data.iter().map(|&value| f(value)).collect(),
        )
    }

    /// Applies `f` to every stored entry in place, keeping the pattern; see `map`.
    pub fn apply(&mut self, f: impl Fn(T) -> T) {
        for value in &mut self.data {
            *value = f(*value);
        }
    }

    /// Scales each row `i` by `factors[i]`, giving the product `diag(factors) · A` as a new
    /// matrix with the same pattern, without forming the diagonal matrix.
    ///
    /// Returns `SparseError::ShapeMismatch` unless there is one factor per row, with `factors`
    /// reported as the square diagonal matrix it stands for.
    pub fn scale_rows(&self, factors: &[T]) -> Result<Self, SparseError> {
        if factors.len() != self.shape.0 {
            return Err(SparseError::ShapeMismatch {
                left: (factors.len(), factors.len()),
                right: self.shape,
            });
        }
        Ok(self.scale_entries(|row, _, value| factors[row] * value))
    }

    /// Scales each column `j` by `factors[j]`, giving the product `A · diag(factors)` as a new
    /// matrix with the same pattern, without forming the diagonal matrix.
    ///
    /// Returns `SparseError::ShapeMismatch` unless there is one factor per column, with
    /// `factors` reported as the square diagonal matrix it stands for.
    pub fn scale_cols(&self, factors: &[T]) -> Result<Self, SparseError> {
        if factors.len() != self.shape.1 {
            return Err(SparseError::ShapeMismatch {
                left: self.shape,
                right: (factors.len(), factors.len()),
            });
        }
        Ok(self.scale_entries(|_, col, value| value * factors[col]))
    }

    /// A new matrix with each entry replaced by `scaled(row, col, value)`.
    fn scale_entries(&self, scaled: impl Fn(usize, usize, T) -> T) -> Self {
        let mut data = Vec::with_capacity(self.data.len());
        for major in 0..self.major_dim() {
            for p in self.ptr[major]..self.ptr[major + 1] {
                let (row, col) = Self::row_col(major, self.idx[p]);
                data.push(scaled(row, col, self.data[p]));
            }
        }
        Self::from_parts_unchecked(self.shape, self.ptr.clone(), self.idx.clone(), data)
    }
}

/// Implements a binary operator for every combination of owned and borrowed operands.
macro_rules! impl_matrix_binary_op {
    ($trait:ident, $method:ident, $try_method:ident, $verb:literal) => {
        #[doc = concat!("Element-wise ", $verb, " of two matrices, returned as a new matrix")]
        /// storing the union of their patterns.
        ///
        /// # Panics
        ///
        /// Panics if the two matrices have different shapes. Use the matching `try_` method to
        /// get an error instead.
        impl<T: Scalar, O: Orientation> $trait<&CompressedMatrix<T, O>>
            for &CompressedMatrix<T, O>
        {
            type Output = CompressedMatrix<T, O>;

            fn $method(self, other: &CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
                self.$try_method(other)
                    .unwrap_or_else(|err| panic!("{err}"))
            }
        }

        impl<T: Scalar, O: Orientation> $trait<CompressedMatrix<T, O>> for &CompressedMatrix<T, O> {
            type Output = CompressedMatrix<T, O>;

            fn $method(self, other: CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
                self.$method(&other)
            }
        }

        impl<T: Scalar, O: Orientation> $trait<&CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
            type Output = CompressedMatrix<T, O>;

            fn $method(self, other: &CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
                (&self).$method(other)
            }
        }

        impl<T: Scalar, O: Orientation> $trait<CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
            type Output = CompressedMatrix<T, O>;

            fn $method(self, other: CompressedMatrix<T, O>) -> CompressedMatrix<T, O> {
                (&self).$method(&other)
            }
        }
    };
}

impl_matrix_binary_op!(Add, add, try_add, "sum");
impl_matrix_binary_op!(Sub, sub, try_sub, "difference");

/// Adds `other` into `self` element-wise; the pattern becomes the union of both.
///
/// # Panics
///
/// Panics if the two matrices have different shapes.
impl<T: Scalar, O: Orientation> AddAssign<&CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
    fn add_assign(&mut self, other: &CompressedMatrix<T, O>) {
        *self = &*self + other;
    }
}

/// Adds `other` into `self` element-wise; the pattern becomes the union of both.
///
/// # Panics
///
/// Panics if the two matrices have different shapes.
impl<T: Scalar, O: Orientation> AddAssign<CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
    fn add_assign(&mut self, other: CompressedMatrix<T, O>) {
        *self += &other;
    }
}

/// Subtracts `other` from `self` element-wise; the pattern becomes the union of both.
///
/// # Panics
///
/// Panics if the two matrices have different shapes.
impl<T: Scalar, O: Orientation> SubAssign<&CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
    fn sub_assign(&mut self, other: &CompressedMatrix<T, O>) {
        *self = &*self - other;
    }
}

/// Subtracts `other` from `self` element-wise; the pattern becomes the union of both.
///
/// # Panics
///
/// Panics if the two matrices have different shapes.
impl<T: Scalar, O: Orientation> SubAssign<CompressedMatrix<T, O>> for CompressedMatrix<T, O> {
    fn sub_assign(&mut self, other: CompressedMatrix<T, O>) {
        *self -= &other;
    }
}

/// Negates every entry, returned as a new matrix with the same pattern.
impl<T: Scalar + Neg<Output = T>, O: Orientation> Neg for &CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn neg(self) -> CompressedMatrix<T, O> {
        self.map(|value| -value)
    }
}

/// Negates every entry in place.
impl<T: Scalar + Neg<Output = T>, O: Orientation> Neg for CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn neg(mut self) -> CompressedMatrix<T, O> {
        self.apply(|value| -value);
        self
    }
}

/// Multiplies every entry by a scalar; see `CompressedMatrix::scale`.
impl<T: Scalar, O: Orientation> Mul<T> for &CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn mul(self, factor: T) -> CompressedMatrix<T, O> {
        self.scale(factor)
    }
}

/// Multiplies every entry by a scalar in place.
impl<T: Scalar, O: Orientation> Mul<T> for CompressedMatrix<T, O> {
    type Output = CompressedMatrix<T, O>;

    fn mul(mut self, factor: T) -> CompressedMatrix<T, O> {
        self.apply(|value| value * factor);
        self
    }
}

/// Multiplies every entry by a scalar in place.
impl<T: Scalar, O: Orientation> MulAssign<T> for CompressedMatrix<T, O> {
    fn mul_assign(&mut self, factor: T) {
        self.apply(|value| value * factor);
    }
}

/// Lets a scalar stand on the left of `*`, as in `2.0 * &matrix`. A blanket impl over `T` is
/// not allowed for a foreign type, so each primitive scalar gets its own.
macro_rules! impl_scalar_on_left {
    ($($scalar:ty),*) => {$(
        impl<O: Orientation> Mul<&CompressedMatrix<$scalar, O>> for $scalar {
            type Output = CompressedMatrix<$scalar, O>;

            fn mul(self, matrix: &CompressedMatrix<$scalar, O>) -> CompressedMatrix<$scalar, O> {
                matrix.scale(self)
            }
        }

        impl<O: Orientation> Mul<CompressedMatrix<$scalar, O>> for $scalar {
            type Output = CompressedMatrix<$scalar, O>;

            fn mul(self, matrix: CompressedMatrix<$scalar, O>) -> CompressedMatrix<$scalar, O> {
                matrix * self
            }
        }
    )*};
}

impl_scalar_on_left!(f32, f64, i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize);

/// Complex parts and conjugates. Generic over `num_complex::ComplexFloat`, which covers `f32`
/// and `f64` as well as their complex types, so the same code works for real matrices: `conj`
/// is then a copy, `adjoint` the transpose, and `imag` all stored zeros.
impl<T: Scalar + ComplexFloat, O: Orientation> CompressedMatrix<T, O> {
    /// The element-wise complex conjugate, with the same pattern.
    pub fn conj(&self) -> Self {
        self.map(|value| value.conj())
    }

    /// The conjugate transpose `Aᴴ`, in the other orientation like `transpose`, so it costs a
    /// conjugating copy of the values. For a real element type this is the transpose.
    pub fn adjoint(&self) -> CompressedMatrix<T, O::Transposed> {
        self.conj().transpose()
    }

    /// The real part of each entry, as a matrix of the underlying real type (`f64` for
    /// `Complex<f64>`), with the same pattern.
    pub fn real(&self) -> CompressedMatrix<T::Real, O> {
        self.map(|value| value.re())
    }

    /// The imaginary part of each entry, as a matrix of the underlying real type, with the
    /// same pattern (so a real matrix's imaginary part is all stored zeros).
    pub fn imag(&self) -> CompressedMatrix<T::Real, O> {
        self.map(|value| value.im())
    }
}

impl<R: Scalar, O: Orientation> CompressedMatrix<Complex<R>, O> {
    /// Builds a complex matrix from its real and imaginary parts, over the union of their
    /// patterns: entry `(i, j)` is `re[(i, j)] + i·im[(i, j)]`. The inverse of `real` and
    /// `imag`.
    ///
    /// Returns `SparseError::ShapeMismatch` unless the two parts have the same shape.
    pub fn from_parts(
        re: &CompressedMatrix<R, O>,
        im: &CompressedMatrix<R, O>,
    ) -> Result<Self, SparseError> {
        merge(
            re,
            im,
            Complex::new,
            |re| Complex::new(re, R::zero()),
            |im| Complex::new(R::zero(), im),
        )
    }
}
