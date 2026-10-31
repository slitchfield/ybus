//! Criterion benchmarks of the DC power flow on one medium-sized MATPOWER case: the repeatable
//! baseline for Phase 2's fixed-pattern matrices, which target building the model, and Phase
//! 3's LU, which targets the solve (`docs/PowerFlowChecklist.md`). `scripts/benchmark-dc.sh`
//! surveys every case instead, once each; its recorded results are in `docs/benchmarks/`.
//!
//! Run with `cargo bench -p powerflow --bench dc_power_flow`. The case is fetched rather than
//! committed (see `tests/common`), so without it the benchmarks are skipped, unless
//! `YBUS_REQUIRE_CASES=1`. `cargo test --all-targets` runs each benchmark once, as a test.

#[path = "../tests/common/mod.rs"]
mod common;

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion, SamplingMode};
use powerflow::case::Case;
use powerflow::dc::DcModel;

/// MATPOWER's case2383wp, the Polish grid at its 1999–2000 winter peak: 2,383 buses, 2,896
/// branches (6 of them phase shifters) and 327 generators. It's a common benchmark, it parses
/// straight from MATPOWER's files, and its numbering suits an unordered factorization badly:
/// the Cholesky factor of B′ has 40 times the entries of B′'s lower triangle, against 13 times
/// for the similarly sized case_ACTIVSg2000, so Phase 3's orderings will show.
const CASE: &str = "case2383wp";

fn dc_power_flow(c: &mut Criterion) {
    let Some(path) = common::matpower_case(CASE) else {
        return; // not fetched; the skip was reported on stderr
    };
    let case = Case::from_path(&path).unwrap();
    let model = DcModel::new(&case).unwrap();

    let mut group = c.benchmark_group(CASE);
    group.bench_function("parse", |b| {
        b.iter(|| Case::from_path(black_box(&path)).unwrap())
    });
    group.bench_function("model", |b| {
        b.iter(|| DcModel::new(black_box(&case)).unwrap())
    });
    // A solve takes about a third of a second without an ordering, so it's sampled flat (the
    //   same number of iterations in every sample), as criterion recommends for long
    //   benchmarks, with fewer samples and the time 20 of them need
    group
        .sampling_mode(SamplingMode::Flat)
        .sample_size(20)
        .measurement_time(Duration::from_secs(10));
    group.bench_function("solve", |b| {
        b.iter(|| model.solve(black_box(&case)).unwrap())
    });
    group.finish();
}

criterion_group!(benches, dc_power_flow);
criterion_main!(benches);
