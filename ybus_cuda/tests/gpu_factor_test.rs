//! GPU solves against the CPU's, on synthetic factorizations: SPD (Cholesky) and unsymmetric
//! (LU, with pivoting), 2×2 blocks, a single row, blocks narrower and wider than the number of
//! right-hand sides, and repeated runs. Each test skips, saying so, on a machine without CUDA,
//! unless `YBUS_REQUIRE_CUDA=1` makes that a failure.

mod common;

use common::{
    gpu, random_spd, random_square, right_hand_sides, two_by_two_blocks, worst_relative_difference,
};
use sparse_mat::compressed_matrix::CscMatrix;
use sparse_mat::Factorization;
use ybus_cuda::{CudaError, GpuFactor};

/// The GPU's solutions agree with the CPU's to within rounding, for each block width
const BAR: f64 = 1e-12;

#[test]
fn cholesky_solves_match_the_cpu_at_every_block_width() {
    let Some(gpu) = gpu("cholesky_solves_match_the_cpu_at_every_block_width") else {
        return;
    };
    let a = random_spd(300, 3, 7);
    let cholesky = a.cholesky().unwrap();
    let rhs = right_hand_sides(300, 37, 11);
    let cpu = cholesky.solve_many(&rhs).unwrap();
    let mut factor = GpuFactor::from_cholesky(&gpu, &cholesky).unwrap();
    // Narrower than the count (several blocks, the last padded), exactly one block, and wider
    for width in [1, 5, 37, 64, 1024, 4096] {
        factor.set_block_width(width);
        let solved = factor.solve_many(&rhs).unwrap();
        let difference = worst_relative_difference(&solved, &cpu);
        assert!(difference < BAR, "width {width}: {difference:e}");
    }
}

#[test]
fn lu_solves_match_the_cpu_with_pivoting() {
    let Some(gpu) = gpu("lu_solves_match_the_cpu_with_pivoting") else {
        return;
    };
    let a = random_square(250, 4, 3);
    let lu = a.lu().unwrap();
    assert!(
        lu.permutation().iter().enumerate().any(|(i, &p)| i != p),
        "the fixture should pivot off the diagonal"
    );
    let rhs = right_hand_sides(250, 20, 5);
    let cpu = lu.solve_many(&rhs).unwrap();
    let mut factor = GpuFactor::from_lu(&gpu, &lu).unwrap();
    factor.set_block_width(8);
    let solved = factor.solve_many(&rhs).unwrap();
    let difference = worst_relative_difference(&solved, &cpu);
    assert!(difference < BAR, "{difference:e}");
}

#[test]
fn two_by_two_blocks_and_a_single_row() {
    let Some(gpu) = gpu("two_by_two_blocks_and_a_single_row") else {
        return;
    };
    let a = two_by_two_blocks(40);
    let lu = a.lu().unwrap();
    let rhs = right_hand_sides(80, 3, 1);
    let mut factor = GpuFactor::from_lu(&gpu, &lu).unwrap();
    let difference = worst_relative_difference(
        &factor.solve_many(&rhs).unwrap(),
        &lu.solve_many(&rhs).unwrap(),
    );
    assert!(difference < BAR, "{difference:e}");

    // A single row, at block widths on both sides of 64 (cuSPARSE 12.8 gets more than 64
    //   wrong, so such a factor is solved on the host), against the CPU bit for bit
    let one = CscMatrix::from_triplets(1, 1, vec![(0, 0, 3.0)])
        .unwrap()
        .lu()
        .unwrap();
    let rhs: Vec<Vec<f64>> = (0..200).map(|j| vec![j as f64 - 70.5]).collect();
    let mut factor = GpuFactor::from_lu(&gpu, &one).unwrap();
    for width in [1, 64, 65, 1024] {
        factor.set_block_width(width);
        assert_eq!(
            factor.solve_many(&rhs).unwrap(),
            one.solve_many(&rhs).unwrap()
        );
    }
}

#[test]
fn repeated_solves_give_the_same_bits() {
    let Some(gpu) = gpu("repeated_solves_give_the_same_bits") else {
        return;
    };
    let cholesky = random_spd(500, 2, 13).cholesky().unwrap();
    let rhs = right_hand_sides(500, 50, 17);
    let mut factor = GpuFactor::from_cholesky(&gpu, &cholesky).unwrap();
    factor.set_block_width(16);
    let first = factor.solve_many(&rhs).unwrap();
    // And again after a change of block width, which redoes the analysis
    factor.set_block_width(7);
    factor.set_block_width(16);
    for _ in 0..3 {
        let again = factor.solve_many(&rhs).unwrap();
        let bits =
            |s: &[Vec<f64>]| -> Vec<u64> { s.iter().flatten().map(|v| v.to_bits()).collect() };
        assert_eq!(bits(&again), bits(&first));
    }
}

#[test]
fn empty_and_mismatched_right_hand_sides() {
    let Some(gpu) = gpu("empty_and_mismatched_right_hand_sides") else {
        return;
    };
    let cholesky = random_spd(10, 20, 1).cholesky().unwrap();
    let mut factor = GpuFactor::from_cholesky(&gpu, &cholesky).unwrap();
    assert_eq!(factor.solve_many(&[]).unwrap(), Vec::<Vec<f64>>::new());
    assert!(matches!(
        factor.solve_many(&[vec![0.0; 10], vec![0.0; 9]]),
        Err(CudaError::ShapeMismatch {
            expected: 10,
            found: 9
        })
    ));
    let timings = factor.timings();
    assert_eq!(timings.blocks, 0);
    factor.solve_many(&[vec![1.0; 10]]).unwrap();
    assert_eq!(factor.timings().blocks, 1);
    factor.reset_timings();
    assert_eq!(factor.timings().blocks, 0);
}
