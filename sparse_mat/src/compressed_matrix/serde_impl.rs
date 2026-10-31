//! `serde` support (feature `serde`).
//!
//! A matrix serializes as its arrays: `{"shape": [rows, cols], "orientation": "CSR",
//! "pointers": [...], "indices": [...], "values": [...]}`, stored zeros included.
//! Deserializing validates the arrays like `CompressedMatrix::try_from_raw_parts`, and rejects
//! a matrix of the other orientation, so a CSR matrix can't be read as a `CscMatrix`.

use serde::{Deserialize, Deserializer, Serialize, Serializer};

use super::{CompressedMatrix, Orientation};
use crate::Scalar;

#[derive(Serialize)]
#[serde(rename = "CompressedMatrix")]
struct SerializeRepr<'a, T> {
    shape: (usize, usize),
    orientation: &'static str,
    pointers: &'a [usize],
    indices: &'a [usize],
    values: &'a [T],
}

#[derive(Deserialize)]
#[serde(rename = "CompressedMatrix")]
struct DeserializeRepr<T> {
    shape: (usize, usize),
    orientation: String,
    pointers: Vec<usize>,
    indices: Vec<usize>,
    values: Vec<T>,
}

impl<T: Scalar + Serialize, O: Orientation> Serialize for CompressedMatrix<T, O> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        SerializeRepr {
            shape: self.shape,
            orientation: O::NAME,
            pointers: &self.ptr,
            indices: &self.idx,
            values: &self.data,
        }
        .serialize(serializer)
    }
}

impl<'de, T: Scalar + Deserialize<'de>, O: Orientation> Deserialize<'de>
    for CompressedMatrix<T, O>
{
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let repr = DeserializeRepr::<T>::deserialize(deserializer)?;
        if repr.orientation != O::NAME {
            return Err(serde::de::Error::custom(format!(
                "expected a {} matrix, found {}",
                O::NAME,
                repr.orientation
            )));
        }
        CompressedMatrix::try_from_raw_parts(
            repr.shape.0,
            repr.shape.1,
            repr.pointers,
            repr.indices,
            repr.values,
        )
        .map_err(serde::de::Error::custom)
    }
}
