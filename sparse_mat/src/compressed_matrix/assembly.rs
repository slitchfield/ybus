//! Assembly maps: build a matrix from triplets once, then refill its values from new triplet
//! values with no sorting, searching or hashing.
//!
//! Power flow assembles the same matrices over and over with new values: the admittance
//! matrix as branch parameters change, and the Newton–Raphson Jacobian at every iteration.
//! The positions of the triplets (one per branch end, say) stay the same, so the pattern, and
//! the slot each triplet sums into, can be worked out once.

use std::marker::PhantomData;

use super::{CompressedMatrix, Csr, Orientation};
use crate::sparse_error::SparseError;
use crate::Scalar;

/// Where each triplet of a triplet list landed in the matrix built from it by
/// [`CompressedMatrix::from_triplets_mapped`], so that [`CompressedMatrix::refill`] can
/// assemble new values for the same positions in O(triplets).
///
/// A map is tied to the orientation `O` of its matrix, since a slot means a different entry in
/// the other orientation, and to its pattern: refilling checks the shape and stored count. It
/// works for the matrix it was built with and for any clone of it.
///
/// ```
/// use sparse_mat::compressed_matrix::CsrMatrix;
///
/// // Two stamps at (0, 0), summed into one entry
/// let positions = [(0, 0), (0, 1), (1, 1), (0, 0)];
/// let stamp = |values: [f64; 4]| positions.iter().zip(values).map(|(&(i, j), v)| (i, j, v));
/// let (mut a, map) = CsrMatrix::from_triplets_mapped(2, 2, stamp([1.0, 2.0, 3.0, 4.0]))?;
/// assert_eq!(a.to_dense(), vec![vec![5.0, 2.0], vec![0.0, 3.0]]);
///
/// // New values for the same stamps, without rebuilding
/// a.refill(&map, &[1.0, -1.0, 2.0, 1.0])?;
/// assert_eq!(a, CsrMatrix::from_triplets(2, 2, stamp([1.0, -1.0, 2.0, 1.0]))?);
/// # Ok::<(), sparse_mat::sparse_error::SparseError>(())
/// ```
///
/// A map only refills a matrix of its orientation:
///
/// ```compile_fail
/// use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix};
///
/// let (csr, map) = CsrMatrix::from_triplets_mapped(1, 1, vec![(0, 0, 1.0)]).unwrap();
/// let mut csc: CscMatrix<f64> = csr.to_csc();
/// csc.refill(&map, &[2.0]).unwrap(); // a CSR slot isn't a CSC slot
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AssemblyMap<O: Orientation = Csr> {
    shape: (usize, usize),
    num_stored: usize,
    /// For triplet `k`: its slot in `values()`, and whether it was the first triplet of its
    /// position in input order (which a refill writes rather than adds to).
    slots: Vec<(usize, bool)>,
    orientation: PhantomData<O>,
}

impl<O: Orientation> AssemblyMap<O> {
    /// The number of triplets the map was built from, which is the number of values a refill
    /// takes.
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    /// Whether the map was built from no triplets.
    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// The slot in `values()` that triplet `k` sums into, for each `k` in input order. Several
    /// triplets share a slot when they share a position.
    pub fn slots(&self) -> impl ExactSizeIterator<Item = usize> + '_ {
        self.slots.iter().map(|&(slot, _)| slot)
    }

    /// The shape of the matrix the map was built with.
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }
}

impl<T: Scalar, O: Orientation> CompressedMatrix<T, O> {
    /// Builds the same matrix as `from_triplets` (COO semantics: duplicates summed in input
    /// order, every given position stored), together with an [`AssemblyMap`] recording which
    /// slot of `values()` each triplet summed into.
    ///
    /// Pass the map to `refill` with new values for the same triplets to reassemble the matrix
    /// without rebuilding it. Returns `SparseError::OutOfBounds` if any triplet lies outside
    /// `rows x cols`.
    pub fn from_triplets_mapped(
        rows: usize,
        cols: usize,
        // COO semantics: duplicate (row, col) entries are summed
        triplets: impl IntoIterator<Item = (usize, usize, T)>,
    ) -> Result<(Self, AssemblyMap<O>), SparseError> {
        let triplets: Vec<(usize, usize, T)> = triplets.into_iter().collect();
        let mut slots = vec![(0, false); triplets.len()];
        let matrix = Self::assemble(rows, cols, &triplets, Some(&mut slots))?;
        let map = AssemblyMap {
            shape: (rows, cols),
            num_stored: matrix.num_stored(),
            slots,
            orientation: PhantomData,
        };
        Ok((matrix, map))
    }

    /// Reassembles the values from new triplet values: `values[k]` is the new value of triplet
    /// `k` of the list `map` was built from, at the same position. The result is exactly what
    /// `from_triplets` would build from the triplets with these values (duplicates summed in
    /// the same order, so bit for bit), and the pattern is unchanged: a position whose values
    /// now sum to zero stays stored.
    ///
    /// Costs O(triplets), with no sorting, searching or allocation. Returns
    /// `SparseError::LengthMismatch` if `values` doesn't have one value per triplet, and
    /// `SparseError::PatternMismatch` if `map` wasn't built for a matrix of this shape and
    /// stored count; on an error the matrix is unchanged.
    pub fn refill(&mut self, map: &AssemblyMap<O>, values: &[T]) -> Result<(), SparseError> {
        if map.shape != self.shape || map.num_stored != self.num_stored() {
            return Err(SparseError::PatternMismatch);
        }
        if values.len() != map.slots.len() {
            return Err(SparseError::LengthMismatch {
                expected: map.slots.len(),
                found: values.len(),
            });
        }
        // Every slot has a first triplet, which comes before the slot's other triplets, so each
        //   slot is written once and then added to in input order, as from_triplets sums
        for (&(slot, first), &value) in map.slots.iter().zip(values) {
            if first {
                self.data[slot] = value;
            } else {
                self.data[slot] = self.data[slot] + value;
            }
        }
        Ok(())
    }
}
