//! Helpers shared by ybus_cuda's tests: the GPU (or a skip), deterministic test matrices and
//! right-hand sides, and the comparison with the CPU.
#![allow(dead_code)]

use sparse_mat::compressed_matrix::{CscMatrix, CsrMatrix};
use ybus_cuda::Gpu;

/// The GPU, or `None` (after printing why) where the test should skip
pub fn gpu(test: &str) -> Option<Gpu> {
    let gpu = Gpu::new().expect("CUDA is present but failed to start");
    if gpu.is_none() {
        if std::env::var("YBUS_REQUIRE_CUDA").is_ok_and(|v| v == "1") {
            panic!("{test}: no CUDA device, and YBUS_REQUIRE_CUDA=1");
        }
        eprintln!("{test}: skipped, no CUDA device");
    }
    gpu
}

pub struct Lcg(pub u64);
impl Lcg {
    pub fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    pub fn value(&mut self) -> f64 {
        (self.next() % 2001) as f64 / 1000.0 - 1.0
    }
}

/// A pseudo-random symmetric, diagonally dominant matrix with a positive diagonal
pub fn random_spd(n: usize, fill_percent: u64, seed: u64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    let mut row_sums = vec![0.0; n];
    for i in 0..n {
        for j in 0..i {
            if rng.next() % 100 < fill_percent {
                let value = -((rng.next() % 1000) as f64 + 1.0) / 997.0;
                triplets.extend([(i, j, value), (j, i, value)]);
                row_sums[i] -= value;
                row_sums[j] -= value;
            }
        }
    }
    for (i, sum) in row_sums.iter().enumerate() {
        triplets.push((i, i, sum + 1.0 + i as f64 / n as f64));
    }
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// A pseudo-random unsymmetric matrix with a small diagonal, so the LU pivots off it
pub fn random_square(n: usize, fill_percent: u64, seed: u64) -> CscMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for i in 0..n {
        for j in 0..n {
            if i != j && rng.next() % 100 < fill_percent {
                triplets.push((i, j, 3.0 * rng.value()));
            }
        }
        triplets.push((i, i, 0.1 + rng.value().abs()));
    }
    CscMatrix::from_triplets(n, n, triplets).unwrap()
}

/// A matrix of 2×2 blocks along a chain, coupled to their neighbours, as a Jacobian's angle and
/// magnitude columns are
pub fn two_by_two_blocks(blocks: usize) -> CscMatrix<f64> {
    let mut triplets = vec![];
    for b in 0..blocks {
        let (p, q) = (2 * b, 2 * b + 1);
        triplets.extend([(p, p, 4.0), (p, q, 1.0), (q, p, -0.5), (q, q, 3.0)]);
        if b + 1 < blocks {
            let (r, s) = (p + 2, q + 2);
            triplets.extend([(p, r, -1.0), (r, p, -1.0), (q, s, -0.7), (s, q, -0.8)]);
        }
    }
    CscMatrix::from_triplets(2 * blocks, 2 * blocks, triplets).unwrap()
}

pub fn right_hand_sides(n: usize, count: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = Lcg(seed);
    (0..count)
        .map(|_| (0..n).map(|_| rng.value()).collect())
        .collect()
}

/// The largest difference between the GPU's and the CPU's solutions, relative to the CPU
/// solution's largest entry
pub fn worst_relative_difference(gpu: &[Vec<f64>], cpu: &[Vec<f64>]) -> f64 {
    assert_eq!(gpu.len(), cpu.len());
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
