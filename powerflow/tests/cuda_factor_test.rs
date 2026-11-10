//! The GPU's solves against the CPU's on real factorizations (feature `cuda`): each case's B′
//! (Cholesky, as the DC power flow factors it) and its Newton Jacobian at the case voltages
//! (LU), solved for many right-hand sides by `ybus_cuda`. The tests skip, saying so, on a machine
//! without CUDA unless `YBUS_REQUIRE_CUDA=1`, and without a case unless `YBUS_REQUIRE_CASES=1`.

#![cfg(feature = "cuda")]

mod common;

use powerflow::ac::AcModel;
use powerflow::case::Case;
use powerflow::dc::DcModel;
use sparse_mat::compressed_matrix::{CsrMatrix, LuOptions, Ordering};
use sparse_mat::Complex;
use sparse_mat::Factorization;
use ybus_cuda::{Gpu, GpuFactor};

const CASES: &[&str] = &["case14", "case118", "case2383wp"];

/// The GPU's solutions agree with the CPU's to within this, relative to each solution's largest
/// entry
const BAR: f64 = 1e-10;

fn gpu(test: &str) -> Option<Gpu> {
    let gpu = Gpu::new().expect("CUDA is present but failed to start");
    if gpu.is_none() {
        if std::env::var("YBUS_REQUIRE_CUDA").is_ok_and(|v| v == "1") {
            panic!("{test}: no CUDA device, and YBUS_REQUIRE_CUDA=1");
        }
        eprintln!("{test}: skipped, no CUDA device");
    }
    gpu
}

/// B′: Bbus over the PV then PQ buses, as the DC power flow solves it
fn b_prime(model: &DcModel) -> CsrMatrix<f64> {
    let unknown: Vec<usize> = model
        .pv_buses()
        .iter()
        .chain(model.pq_buses())
        .copied()
        .collect();
    let mut position = vec![None; model.bus_rows().len()];
    for (k, &bus) in unknown.iter().enumerate() {
        position[bus] = Some(k);
    }
    let triplets = model
        .bbus()
        .iter()
        .filter_map(|(i, j, v)| Some((position[i]?, position[j]?, v)));
    CsrMatrix::from_triplets(unknown.len(), unknown.len(), triplets).unwrap()
}

/// The case's bus voltages, in the model's internal order
fn case_voltages(model: &AcModel, case: &Case) -> Vec<Complex<f64>> {
    model
        .bus_rows()
        .iter()
        .map(|&row| {
            let bus = &case.bus[row];
            Complex::from_polar(bus.vm_pu, bus.va_deg.to_radians())
        })
        .collect()
}

/// `count` deterministic right-hand sides of length `n`, with entries in [-1, 1]
fn right_hand_sides(n: usize, count: usize) -> Vec<Vec<f64>> {
    (0..count)
        .map(|j| {
            (0..n)
                .map(|i| (((i * 7919 + j * 104729) % 2001) as f64) / 1000.0 - 1.0)
                .collect()
        })
        .collect()
}

fn worst_relative_difference(gpu: &[Vec<f64>], cpu: &[Vec<f64>]) -> f64 {
    gpu.iter()
        .zip(cpu)
        .map(|(g, c)| {
            let scale = c
                .iter()
                .fold(0.0f64, |m, v| m.max(v.abs()))
                .max(f64::MIN_POSITIVE);
            g.iter()
                .zip(c)
                .fold(0.0f64, |m, (g, c)| m.max((g - c).abs()))
                / scale
        })
        .fold(0.0, f64::max)
}

#[test]
fn b_prime_and_jacobian_solves_match_the_cpu() {
    let Some(gpu) = gpu("b_prime_and_jacobian_solves_match_the_cpu") else {
        return;
    };
    for &name in CASES {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();

        let b = b_prime(&DcModel::new(&case).unwrap());
        let cholesky = b
            .analyze_cholesky(Ordering::Amd)
            .unwrap()
            .factor(&b)
            .unwrap();
        let rhs = right_hand_sides(b.shape().0, 100);
        let mut factor = GpuFactor::from_cholesky(&gpu, &cholesky).unwrap();
        factor.set_block_width(32);
        let difference = worst_relative_difference(
            &factor.solve_many(&rhs).unwrap(),
            &cholesky.solve_many(&rhs).unwrap(),
        );
        assert!(difference < BAR, "{name} B′: {difference:e}");

        let model = AcModel::new(&case).unwrap();
        let jacobian = model.jacobian(&case_voltages(&model, &case));
        let lu = jacobian
            .analyze_lu(Ordering::Amd)
            .unwrap()
            .factor(&jacobian, LuOptions::default())
            .unwrap();
        let rhs = right_hand_sides(jacobian.shape().0, 100);
        let mut factor = GpuFactor::from_lu(&gpu, &lu).unwrap();
        factor.set_block_width(32);
        let difference = worst_relative_difference(
            &factor.solve_many(&rhs).unwrap(),
            &lu.solve_many(&rhs).unwrap(),
        );
        assert!(difference < BAR, "{name} Jacobian: {difference:e}");
        eprintln!("{name}: B′ and Jacobian solves agree");
    }
}
