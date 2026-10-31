//! `serde` support (feature `serde`).
//!
//! A matrix serializes as `{"shape": [rows, cols], "entries": [[row, col, value], ...]}`
//! with entries in row-major order, so the output does not depend on insertion order.
//! Deserializing validates the entries like
//! `SparseMatrix::from_triplets`: out-of-bounds entries are an error, duplicates are summed
//! and zeros are dropped.

use serde::{Deserialize, Deserializer, Serialize, Serializer};

use super::SparseMatrix;
use crate::Scalar;

#[derive(Serialize)]
#[serde(rename = "SparseMatrix")]
struct SerializeRepr<T> {
    shape: (usize, usize),
    entries: Vec<(usize, usize, T)>,
}

#[derive(Deserialize)]
#[serde(rename = "SparseMatrix")]
struct DeserializeRepr<T> {
    shape: (usize, usize),
    entries: Vec<(usize, usize, T)>,
}

impl<T: Scalar + Serialize> Serialize for SparseMatrix<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        SerializeRepr {
            shape: self.shape,
            entries: self.iter().collect(),
        }
        .serialize(serializer)
    }
}

impl<'de, T: Scalar + Deserialize<'de>> Deserialize<'de> for SparseMatrix<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let repr = DeserializeRepr::<T>::deserialize(deserializer)?;
        SparseMatrix::from_triplets(repr.shape.0, repr.shape.1, repr.entries)
            .map_err(serde::de::Error::custom)
    }
}
