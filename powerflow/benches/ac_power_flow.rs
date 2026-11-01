//! Criterion benchmarks of the AC power flow on MATPOWER's case2383wp, the case the DC benchmark
//! uses (`dc_power_flow.rs`): building the model, solving by Newton's method and by
//! fast-decoupled XB and BX, and building one Newton Jacobian, to measure the Jacobian build's
//! share of a Newton solve.
//!
//! Run with `cargo bench -p powerflow --bench ac_power_flow`. Like `dc_power_flow`, it's skipped
//! without the fetched case unless `YBUS_REQUIRE_CASES=1`, and `cargo test --all-targets` runs
//! each benchmark once, as a test.

#[path = "../tests/common/mod.rs"]
mod common;

use std::f64::consts::PI;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use powerflow::ac::{AcMethod, AcModel, AcOptions};
use powerflow::case::Case;
use sparse_mat::Complex;

const CASE: &str = "case2383wp";

fn ac_power_flow(c: &mut Criterion) {
    let Some(path) = common::matpower_case(CASE) else {
        return; // not fetched; the skip was reported on stderr
    };
    let case = Case::from_path(&path).unwrap();
    let model = AcModel::new(&case).unwrap();
    let solution = model.solve(&case, AcOptions::default()).unwrap();
    // The converged voltages, in the model's internal order, for the Jacobian benchmark
    let v: Vec<Complex<f64>> = model
        .bus_rows()
        .iter()
        .map(|&row| Complex::from_polar(solution.vm_pu[row], solution.va_deg[row] * PI / 180.0))
        .collect();

    let mut group = c.benchmark_group(format!("ac/{CASE}"));
    group.bench_function("model", |b| {
        b.iter(|| AcModel::new(black_box(&case)).unwrap())
    });
    for (name, method) in [
        ("solve_nr", AcMethod::Newton),
        ("solve_fdxb", AcMethod::FastDecoupledXb),
        ("solve_fdbx", AcMethod::FastDecoupledBx),
    ] {
        let options = AcOptions::for_method(method);
        group.bench_function(name, |b| {
            b.iter(|| model.solve(black_box(&case), options).unwrap())
        });
    }
    group.bench_function("jacobian", |b| b.iter(|| model.jacobian(black_box(&v))));
    group.finish();
}

criterion_group!(benches, ac_power_flow);
criterion_main!(benches);
