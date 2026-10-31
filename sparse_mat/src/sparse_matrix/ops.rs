//! Arithmetic on `SparseMatrix`, which forwards to `CompressedMatrix` until Phase 2 removes it
//! from `SparseMatrix`. Results drop entries that are zero, as `SparseMatrix` never stores one.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar> SparseMatrix<T> {
    /// Element-wise sum; forwards to `CompressedMatrix::try_add`.
    pub fn try_add(&self, other: &SparseMatrix<T>) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .try_add(&other.to_csr())
            .map(|sum| sum.to_sparse())
    }

    /// Element-wise difference; forwards to `CompressedMatrix::try_sub`.
    pub fn try_sub(&self, other: &SparseMatrix<T>) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .try_sub(&other.to_csr())
            .map(|difference| difference.to_sparse())
    }

    /// Multiplies every entry by `factor`; forwards to `CompressedMatrix::scale`.
    pub fn scale(&self, factor: T) -> SparseMatrix<T> {
        self.to_csr().scale(factor).to_sparse()
    }

    /// Element-wise product; forwards to `CompressedMatrix::hadamard`.
    pub fn hadamard(&self, other: &SparseMatrix<T>) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .hadamard(&other.to_csr())
            .map(|product| product.to_sparse())
    }

    /// Applies `f` to every stored entry; forwards to `CompressedMatrix::map`.
    pub fn map<U: Scalar>(&self, f: impl Fn(T) -> U) -> SparseMatrix<U> {
        self.to_csr().map(f).to_sparse()
    }

    /// Applies `f` to every stored entry in place; forwards to `CompressedMatrix::apply`.
    pub fn apply(&mut self, f: impl Fn(T) -> T) {
        let mut compressed = self.to_csr();
        compressed.apply(f);
        *self = compressed.to_sparse();
    }

    /// `diag(factors) · A`; forwards to `CompressedMatrix::scale_rows`.
    pub fn scale_rows(&self, factors: &[T]) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .scale_rows(factors)
            .map(|scaled| scaled.to_sparse())
    }

    /// `A · diag(factors)`; forwards to `CompressedMatrix::scale_cols`.
    pub fn scale_cols(&self, factors: &[T]) -> Result<SparseMatrix<T>, SparseError> {
        self.to_csr()
            .scale_cols(factors)
            .map(|scaled| scaled.to_sparse())
    }
}

/// Implements a binary operator for every combination of owned and borrowed operands.
macro_rules! impl_matrix_binary_op {
    ($trait:ident, $method:ident, $try_method:ident) => {
        /// Forwards to the matching `try_` method, panicking on different shapes.
        impl<T: Scalar> $trait<&SparseMatrix<T>> for &SparseMatrix<T> {
            type Output = SparseMatrix<T>;

            fn $method(self, other: &SparseMatrix<T>) -> SparseMatrix<T> {
                self.$try_method(other)
                    .unwrap_or_else(|err| panic!("{err}"))
            }
        }

        impl<T: Scalar> $trait<SparseMatrix<T>> for &SparseMatrix<T> {
            type Output = SparseMatrix<T>;

            fn $method(self, other: SparseMatrix<T>) -> SparseMatrix<T> {
                self.$method(&other)
            }
        }

        impl<T: Scalar> $trait<&SparseMatrix<T>> for SparseMatrix<T> {
            type Output = SparseMatrix<T>;

            fn $method(self, other: &SparseMatrix<T>) -> SparseMatrix<T> {
                (&self).$method(other)
            }
        }

        impl<T: Scalar> $trait<SparseMatrix<T>> for SparseMatrix<T> {
            type Output = SparseMatrix<T>;

            fn $method(self, other: SparseMatrix<T>) -> SparseMatrix<T> {
                (&self).$method(&other)
            }
        }
    };
}

impl_matrix_binary_op!(Add, add, try_add);
impl_matrix_binary_op!(Sub, sub, try_sub);

/// Adds `other` into `self` element-wise. Panics if the shapes differ.
impl<T: Scalar> AddAssign<&SparseMatrix<T>> for SparseMatrix<T> {
    fn add_assign(&mut self, other: &SparseMatrix<T>) {
        *self = &*self + other;
    }
}

/// Adds `other` into `self` element-wise. Panics if the shapes differ.
impl<T: Scalar> AddAssign<SparseMatrix<T>> for SparseMatrix<T> {
    fn add_assign(&mut self, other: SparseMatrix<T>) {
        *self += &other;
    }
}

/// Subtracts `other` from `self` element-wise. Panics if the shapes differ.
impl<T: Scalar> SubAssign<&SparseMatrix<T>> for SparseMatrix<T> {
    fn sub_assign(&mut self, other: &SparseMatrix<T>) {
        *self = &*self - other;
    }
}

/// Subtracts `other` from `self` element-wise. Panics if the shapes differ.
impl<T: Scalar> SubAssign<SparseMatrix<T>> for SparseMatrix<T> {
    fn sub_assign(&mut self, other: SparseMatrix<T>) {
        *self -= &other;
    }
}

/// Negates every entry.
impl<T: Scalar + Neg<Output = T>> Neg for &SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn neg(self) -> SparseMatrix<T> {
        (-self.to_csr()).to_sparse()
    }
}

/// Negates every entry.
impl<T: Scalar + Neg<Output = T>> Neg for SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn neg(self) -> SparseMatrix<T> {
        -&self
    }
}

/// Multiplies every entry by a scalar; see `SparseMatrix::scale`.
impl<T: Scalar> Mul<T> for &SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn mul(self, factor: T) -> SparseMatrix<T> {
        self.scale(factor)
    }
}

/// Multiplies every entry by a scalar.
impl<T: Scalar> Mul<T> for SparseMatrix<T> {
    type Output = SparseMatrix<T>;

    fn mul(self, factor: T) -> SparseMatrix<T> {
        self.scale(factor)
    }
}

/// Multiplies every entry by a scalar in place.
impl<T: Scalar> MulAssign<T> for SparseMatrix<T> {
    fn mul_assign(&mut self, factor: T) {
        *self = self.scale(factor);
    }
}

/// Lets a scalar stand on the left of `*`, as in `2.0 * &matrix`.
macro_rules! impl_scalar_on_left {
    ($($scalar:ty),*) => {$(
        impl Mul<&SparseMatrix<$scalar>> for $scalar {
            type Output = SparseMatrix<$scalar>;

            fn mul(self, matrix: &SparseMatrix<$scalar>) -> SparseMatrix<$scalar> {
                matrix.scale(self)
            }
        }

        impl Mul<SparseMatrix<$scalar>> for $scalar {
            type Output = SparseMatrix<$scalar>;

            fn mul(self, matrix: SparseMatrix<$scalar>) -> SparseMatrix<$scalar> {
                matrix.scale(self)
            }
        }
    )*};
}

impl_scalar_on_left!(f32, f64, i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize);
