//! Blocks built, solved and updated on the device, against the same formulas computed on the
//! host from the CPU's solves: pair solves, rank-one updates (including a singular one),
//! operator products with entries zeroed, and per-column maxima. The LU fixtures pivot, so the
//! right-hand side's and the solution's orders differ, as they must be told apart.

mod common;

use common::{gpu, random_spd, random_square, Lcg};
use sparse_mat::compressed_matrix::{Cholesky, CsrMatrix, Lu};
use sparse_mat::Factorization;
use ybus_cuda::{ColumnMax, GpuFactor, Pair, PairUpdate};

const BAR: f64 = 1e-12;

/// Either factorization, solved on the CPU
enum Cpu {
    Cholesky(Cholesky<f64>),
    Lu(Lu<f64>),
}

impl Cpu {
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        match self {
            Cpu::Cholesky(f) => f.solve(b),
            Cpu::Lu(f) => f.solve(b),
        }
        .unwrap()
    }
}

fn pair_vector(n: usize, pair: Pair) -> Vec<f64> {
    let mut c = vec![0.0; n];
    if let Some(p) = pair.plus {
        c[p] = 1.0;
    }
    if let Some(m) = pair.minus {
        c[m] = -1.0;
    }
    c
}

/// cᵀv for c = e_plus − e_minus
fn dot(pair: Pair, v: &[f64]) -> f64 {
    pair.plus.map_or(0.0, |p| v[p]) - pair.minus.map_or(0.0, |m| v[m])
}

/// The update `PairUpdate` describes, on the host: the result and whether it was singular
fn host_update(cpu: &Cpu, base: &[f64], update: &PairUpdate) -> (Vec<f64>, bool) {
    let n = base.len();
    let y = cpu.solve(&pair_vector(n, update.pair));
    let x: Vec<f64> = base
        .iter()
        .zip(&y)
        .map(|(b, y)| b + update.alpha * y)
        .collect();
    let w = -dot(update.pair, &y) + update.inverse_weight;
    if w == 0.0 {
        return (x, true);
    }
    let z = dot(update.pair, &x) / w;
    (x.iter().zip(&y).map(|(x, y)| x + z * y).collect(), false)
}

fn max_difference(a: &[f64], b: &[f64]) -> f64 {
    let scale = b.iter().fold(1.0f64, |m, v| m.max(v.abs()));
    a.iter()
        .zip(b)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
        / scale
}

fn pairs(n: usize, count: usize, seed: u64) -> Vec<Pair> {
    let mut rng = Lcg(seed);
    (0..count)
        .map(|k| {
            let plus = (k % 5 != 3).then(|| rng.next() as usize % n);
            let mut minus = (k % 7 != 2).then(|| rng.next() as usize % n);
            if minus.is_some() && minus == plus {
                minus = None;
            }
            Pair { plus, minus }
        })
        .collect()
}

/// A random rows × n operator
fn operator_matrix(rows: usize, n: usize, seed: u64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for i in 0..rows {
        for _ in 0..3 {
            triplets.push((i, rng.next() as usize % n, rng.value()));
        }
    }
    CsrMatrix::from_triplets(rows, n, triplets).unwrap()
}

fn factors(n: usize) -> Vec<(&'static str, Cpu)> {
    vec![
        (
            "cholesky",
            Cpu::Cholesky(random_spd(n, 3, 21).cholesky().unwrap()),
        ),
        ("lu", Cpu::Lu(random_square(n, 4, 22).lu().unwrap())),
    ]
}

fn upload<'g>(gpu: &'g ybus_cuda::Gpu, cpu: &Cpu) -> GpuFactor<'g> {
    match cpu {
        Cpu::Cholesky(f) => GpuFactor::from_cholesky(gpu, f).unwrap(),
        Cpu::Lu(f) => GpuFactor::from_lu(gpu, f).unwrap(),
    }
}

#[test]
fn pair_solves_and_their_maxima_match_the_host() {
    let Some(gpu) = gpu("pair_solves_and_their_maxima_match_the_host") else {
        return;
    };
    let n = 200;
    for (name, cpu) in factors(n) {
        let mut factor = upload(&gpu, &cpu);
        factor.set_block_width(32);
        let pairs = pairs(n, 29, 5);
        let mut block = factor.pair_solves(&pairs).unwrap();
        assert_eq!(block.columns(), 29);
        assert!(block.singular().iter().all(|&s| !s));
        let solutions = block.solutions().unwrap();
        let maxima = block.maxima().unwrap();
        for ((pair, solution), max) in pairs.iter().zip(&solutions).zip(&maxima) {
            let expected = cpu.solve(&pair_vector(n, *pair));
            let difference = max_difference(solution, &expected);
            assert!(difference < BAR, "{name}: {difference:e}");
            let largest = expected.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            assert!((max.value - largest).abs() <= BAR * largest, "{name}");
            match max.row {
                Some(row) => assert!(
                    (solution[row].abs() - max.value).abs() <= BAR * largest,
                    "{name}"
                ),
                // An empty pair's solution is zero
                None => assert!(pair.plus.is_none() && pair.minus.is_none(), "{name}"),
            }
            assert!(!max.not_finite);
        }
    }
}

#[test]
fn rank_one_updates_match_the_host_including_a_singular_one() {
    let Some(gpu) = gpu("rank_one_updates_match_the_host_including_a_singular_one") else {
        return;
    };
    let n = 150;
    let mut rng = Lcg(9);
    let base: Vec<f64> = (0..n).map(|_| rng.value()).collect();
    for (name, cpu) in factors(n) {
        let mut factor = upload(&gpu, &cpu);
        factor.set_block_width(16);
        let mut updates: Vec<PairUpdate> = pairs(n, 15, 6)
            .into_iter()
            .map(|pair| PairUpdate {
                pair,
                alpha: rng.value(),
                inverse_weight: 2.0 + rng.value(),
            })
            .collect();
        // An empty pair with a zero weight: w = 0 exactly, so singular
        updates.push(PairUpdate {
            pair: Pair {
                plus: None,
                minus: None,
            },
            alpha: 0.5,
            inverse_weight: 0.0,
        });
        let mut block = factor.pair_updates(&base, &updates).unwrap();
        let solutions = block.solutions().unwrap();
        for (k, (update, solution)) in updates.iter().zip(&solutions).enumerate() {
            let (expected, singular) = host_update(&cpu, &base, update);
            assert_eq!(block.singular()[k], singular, "{name} column {k}");
            let difference = max_difference(solution, &expected);
            assert!(difference < BAR, "{name} column {k}: {difference:e}");
        }
        assert!(block.singular()[15]);
    }
}

#[test]
fn operator_products_and_maxima_match_the_host() {
    let Some(gpu) = gpu("operator_products_and_maxima_match_the_host") else {
        return;
    };
    let n = 120;
    let rows = 90;
    let mut rng = Lcg(31);
    let base: Vec<f64> = (0..n).map(|_| rng.value()).collect();
    let offset: Vec<f64> = (0..rows).map(|_| rng.value()).collect();
    let weights: Vec<f64> = (0..rows)
        .map(|i| {
            if i % 4 == 0 {
                0.0
            } else {
                1.0 + rng.value().abs()
            }
        })
        .collect();
    let m = operator_matrix(rows, n, 32);
    for (name, cpu) in factors(n) {
        let mut factor = upload(&gpu, &cpu);
        factor.set_block_width(8);
        let mut operator = factor.operator(&m).unwrap();
        operator.set_offset(&offset).unwrap();
        operator.set_weights(&weights).unwrap();
        let updates: Vec<PairUpdate> = pairs(n, 7, 33)
            .into_iter()
            .map(|pair| PairUpdate {
                pair,
                alpha: 0.25,
                inverse_weight: 3.0,
            })
            .collect();
        let zero: Vec<Option<usize>> = (0..7).map(|k| (k % 2 == 0).then_some(k * 11)).collect();
        let mut block = factor.pair_updates(&base, &updates).unwrap();
        let products = block.product(&mut operator, &zero).unwrap();
        let maxima = block.product_maxima(&mut operator, &zero).unwrap();
        for k in 0..7 {
            let (theta, _) = host_update(&cpu, &base, &updates[k]);
            let mut expected = m.mul_vec(&theta).unwrap();
            for (e, o) in expected.iter_mut().zip(&offset) {
                *e += o;
            }
            if let Some(row) = zero[k] {
                expected[row] = 0.0;
                assert_eq!(products[k][row], 0.0, "{name}: the zeroed entry");
            }
            let difference = max_difference(&products[k], &expected);
            assert!(difference < BAR, "{name} column {k}: {difference:e}");
            // The host's maximum of |value| · weight, the smallest row on a tie
            let (mut row, mut value) = (None, 0.0);
            for (i, (e, w)) in products[k].iter().zip(&weights).enumerate() {
                if e.abs() * w > value {
                    (row, value) = (Some(i), e.abs() * w);
                }
            }
            assert_eq!(
                maxima[k],
                ColumnMax {
                    row,
                    value,
                    not_finite: false
                },
                "{name} column {k}"
            );
        }
    }
}

#[test]
fn a_single_row_factor_keeps_its_blocks_narrow() {
    let Some(gpu) = gpu("a_single_row_factor_keeps_its_blocks_narrow") else {
        return;
    };
    // cuSPARSE 12.8's SpSM is wrong for a 1×1 matrix with more than 64 columns
    let lu = sparse_mat::compressed_matrix::CscMatrix::from_triplets(1, 1, vec![(0, 0, 4.0)])
        .unwrap()
        .lu()
        .unwrap();
    let mut factor = GpuFactor::from_lu(&gpu, &lu).unwrap();
    factor.set_block_width(1024);
    assert_eq!(factor.block_width(), 64);
    let pairs = vec![
        Pair {
            plus: Some(0),
            minus: None
        };
        64
    ];
    let solutions = factor.pair_solves(&pairs).unwrap().solutions().unwrap();
    assert!(solutions.iter().all(|s| s == &vec![0.25]));
}
