//! Criterion benchmarks of reassembling a matrix with new values: through an assembly map
//! (`CompressedMatrix::refill`), against building it afresh with `from_triplets`. Power flow
//! does this whenever values change but the topology doesn't: the admittance matrix for each
//! new set of branch parameters, and the Newton–Raphson Jacobian at every iteration.
//!
//! Two patterns from MATPOWER's case2383wp:
//! - **B stamps** (Ybus-shaped): the four stamps of each in-service branch that
//!   `DcModel::new` sums into Bbus, 2,383 × 2,383 from 11,584 triplets.
//! - **2×2 blocks** (Jacobian-shaped): each B stamp at `(i, j)` expanded into a 2×2 block at
//!   `(i, j)`, `(i, j + n)`, `(i + n, j)` and `(i + n, j + n)`, like the polar Jacobian's four
//!   blocks over the same pattern (before the reference and PV rows and columns are removed):
//!   4,766 × 4,766 from 46,336 triplets.
//!
//! Run with `cargo bench -p powerflow --bench assembly`. Like `dc_power_flow`, it's skipped
//! without the fetched case unless `YBUS_REQUIRE_CASES=1`, and `cargo test --all-targets`
//! runs each benchmark once, as a test.

#[path = "../tests/common/mod.rs"]
mod common;

use std::collections::HashMap;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use powerflow::case::Case;
use powerflow::dc::DcModel;
use sparse_mat::compressed_matrix::CsrMatrix;

const CASE: &str = "case2383wp";

/// A triplet list's positions, its matrix's shape, and two sets of values for them.
struct Pattern {
    dim: usize,
    positions: Vec<(usize, usize)>,
    values: Vec<f64>,
    new_values: Vec<f64>,
}

impl Pattern {
    fn triplets<'a>(&'a self, values: &'a [f64]) -> impl Iterator<Item = (usize, usize, f64)> + 'a {
        self.positions
            .iter()
            .zip(values)
            .map(|(&(row, col), &value)| (row, col, value))
    }
}

/// The B stamps of the case's in-service branches, built as `DcModel::new` builds them.
fn b_stamps(case: &Case, model: &DcModel) -> Pattern {
    let internal: HashMap<usize, usize> = model
        .bus_rows()
        .iter()
        .enumerate()
        .map(|(index, &row)| (case.bus[row].bus_label, index))
        .collect();
    let (mut positions, mut values) = (vec![], vec![]);
    for &row in model.branch_rows() {
        let branch = &case.branch[row];
        let (from, to) = (internal[&branch.from_label], internal[&branch.to_label]);
        let b = 1.0 / branch.x / branch.effective_tap();
        positions.extend([(from, from), (from, to), (to, from), (to, to)]);
        values.extend([b, -b, -b, b]);
    }
    // New values as if every reactance changed by a different small amount
    let new_values = values
        .iter()
        .enumerate()
        .map(|(k, value)| value * (1.0 + (k % 7) as f64 * 1e-3))
        .collect();
    Pattern {
        dim: model.bus_rows().len(),
        positions,
        values,
        new_values,
    }
}

/// Each stamp of `stamps` expanded into a 2×2 block over the four quadrants.
fn two_by_two_blocks(stamps: &Pattern) -> Pattern {
    let n = stamps.dim;
    let expand = |values: &[f64]| -> Vec<f64> {
        values
            .iter()
            .flat_map(|&value| [value, 0.5 * value, -0.5 * value, value])
            .collect()
    };
    Pattern {
        dim: 2 * n,
        positions: stamps
            .positions
            .iter()
            .flat_map(|&(i, j)| [(i, j), (i, j + n), (i + n, j), (i + n, j + n)])
            .collect(),
        values: expand(&stamps.values),
        new_values: expand(&stamps.new_values),
    }
}

fn assembly(c: &mut Criterion) {
    let Some(path) = common::matpower_case(CASE) else {
        return; // not fetched; the skip was reported on stderr
    };
    let case = Case::from_path(&path).unwrap();
    let model = DcModel::new(&case).unwrap();
    let stamps = b_stamps(&case, &model);
    // The stamps are DcModel's: built afresh, they give its Bbus exactly
    let bbus = CsrMatrix::from_triplets(stamps.dim, stamps.dim, stamps.triplets(&stamps.values));
    assert_eq!(bbus.unwrap().values(), model.bbus().values());
    let blocks = two_by_two_blocks(&stamps);

    let mut group = c.benchmark_group(format!("assembly/{CASE}"));
    for (name, pattern) in [("b_stamps", &stamps), ("2x2_blocks", &blocks)] {
        let (mut matrix, map) = CsrMatrix::from_triplets_mapped(
            pattern.dim,
            pattern.dim,
            pattern.triplets(&pattern.values),
        )
        .unwrap();
        // A refill gives exactly what from_triplets builds from the same values
        matrix.refill(&map, &pattern.new_values).unwrap();
        let fresh = CsrMatrix::from_triplets(
            pattern.dim,
            pattern.dim,
            pattern.triplets(&pattern.new_values),
        )
        .unwrap();
        assert_eq!(matrix.values(), fresh.values());

        group.bench_function(format!("{name}/from_triplets"), |b| {
            b.iter(|| {
                CsrMatrix::from_triplets(
                    pattern.dim,
                    pattern.dim,
                    pattern.triplets(black_box(&pattern.new_values)),
                )
                .unwrap()
            })
        });
        group.bench_function(format!("{name}/refill"), |b| {
            b.iter(|| matrix.refill(&map, black_box(&pattern.new_values)).unwrap())
        });
    }
    group.finish();
}

criterion_group!(benches, assembly);
criterion_main!(benches);
