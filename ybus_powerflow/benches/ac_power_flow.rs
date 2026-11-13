//! Criterion benchmarks of the AC power flow on MATPOWER's case2383wp, the case the DC benchmark
//! uses (`dc_power_flow.rs`): building the model, solving by Newton's method and by
//! fast-decoupled XB and BX, and `AcModel::jacobian`. That lays the Jacobian out and fills it
//! once; inside a Newton solve the layout is made once per solve and only the refill repeats,
//! so the solve's time is the better measure of the Jacobian's cost.
//!
//! The `jacobian_*` benchmarks time the other parts of a round of a reactive-limit loop, on the
//! Jacobian at the solution: the AMD analysis of its pattern, a full factorization under it, a
//! refactor reusing its pivots, and a solve. `superset_ordering` is the ordering a loop can
//! compute once instead: AMD on Ybus, each bus expanded to its angle and then its magnitude.
//!
//! Run with `cargo bench -p ybus_powerflow --bench ac_power_flow`. Like `dc_power_flow`, it's
//! skipped without the fetched case unless `YBUS_REQUIRE_CASES=1`, and `cargo test --all-targets`
//! runs each benchmark once, as a test.

#[path = "../tests/common/mod.rs"]
mod common;

use std::f64::consts::PI;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use ybus_powerflow::ac::{AcMethod, AcModel, AcOptions};
use ybus_powerflow::case::Case;
use ybus_sparse::compressed_matrix::{LuOptions, Ordering};
use ybus_sparse::Complex;
use ybus_sparse::Factorization;

const CASE: &str = "case2383wp";

fn ac_power_flow(c: &mut Criterion) {
    let Some(path) = common::matpower_case(CASE) else {
        return; // not fetched; the skip was reported on stderr
    };
    let case = Case::from_path(&path).unwrap();
    let model = AcModel::new(&case).unwrap();
    let solution = model.solve(&case, &AcOptions::default()).unwrap();
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
            b.iter(|| model.solve(black_box(&case), &options).unwrap())
        });
    }
    group.bench_function("jacobian", |b| b.iter(|| model.jacobian(black_box(&v))));

    let jacobian = model.jacobian(&v);
    group.bench_function("jacobian_amd", |b| {
        b.iter(|| black_box(&jacobian).analyze_lu(Ordering::Amd).unwrap())
    });
    let symbolic = jacobian.analyze_lu(Ordering::Amd).unwrap();
    group.bench_function("jacobian_factor", |b| {
        b.iter(|| {
            symbolic
                .factor(black_box(&jacobian), LuOptions::default())
                .unwrap()
        })
    });
    let mut lu = symbolic.factor(&jacobian, LuOptions::default()).unwrap();
    group.bench_function("jacobian_refactor", |b| {
        b.iter(|| lu.refactor(black_box(&jacobian)).unwrap())
    });
    let rhs: Vec<f64> = (0..jacobian.shape().0).map(|i| (i as f64).sin()).collect();
    group.bench_function("jacobian_solve", |b| {
        b.iter(|| lu.solve(black_box(&rhs)).unwrap())
    });
    group.bench_function("superset_ordering", |b| {
        b.iter(|| {
            let buses = black_box(model.ybus()).ordering(Ordering::Amd).unwrap();
            let expanded: Vec<usize> = buses
                .iter()
                .flat_map(|&bus| [2 * bus, 2 * bus + 1])
                .collect();
            expanded
        })
    });
    group.finish();
}

criterion_group!(benches, ac_power_flow);
criterion_main!(benches);
