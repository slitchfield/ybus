//! Shape predicates on `SparseMatrix`. Reductions and norms are on `CompressedMatrix`.

use super::SparseMatrix;
use crate::Scalar;

impl<T: Scalar> SparseMatrix<T> {
    /// Whether the matrix has as many rows as columns.
    pub fn is_square(&self) -> bool {
        self.shape.0 == self.shape.1
    }
}
