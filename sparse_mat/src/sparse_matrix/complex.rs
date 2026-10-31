//! Complex parts on `SparseMatrix`, which forward to `CompressedMatrix` until Phase 2 removes
//! them from `SparseMatrix`.

use num_complex::{Complex, ComplexFloat};

use super::SparseMatrix;
use crate::sparse_error::SparseError;
use crate::Scalar;

impl<T: Scalar + ComplexFloat> SparseMatrix<T> {
    /// The complex conjugate; forwards to `CompressedMatrix::conj`.
    pub fn conj(&self) -> SparseMatrix<T> {
        self.to_csr().conj().to_sparse()
    }

    /// The conjugate transpose; forwards to `CompressedMatrix::adjoint`.
    pub fn adjoint(&self) -> SparseMatrix<T> {
        self.to_csr().adjoint().to_sparse()
    }

    /// The real parts; forwards to `CompressedMatrix::real`.
    pub fn real(&self) -> SparseMatrix<T::Real> {
        self.to_csr().real().to_sparse()
    }

    /// The imaginary parts; forwards to `CompressedMatrix::imag`.
    pub fn imag(&self) -> SparseMatrix<T::Real> {
        self.to_csr().imag().to_sparse()
    }
}

impl<R: Scalar> SparseMatrix<Complex<R>> {
    /// A complex matrix from its parts; forwards to `CompressedMatrix::from_parts`.
    pub fn from_parts(
        re: &SparseMatrix<R>,
        im: &SparseMatrix<R>,
    ) -> Result<SparseMatrix<Complex<R>>, SparseError> {
        crate::compressed_matrix::CsrMatrix::from_parts(&re.to_csr(), &im.to_csr())
            .map(|matrix| matrix.to_sparse())
    }
}
