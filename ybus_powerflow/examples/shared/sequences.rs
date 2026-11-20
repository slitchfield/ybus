//! The sequences and samples the OPF survey examples share (`opf_sequence`, `opf_batch`): the
//! load walk and sweep, and the seeded outage samples, the same on every machine.

// Each example uses part of it
#![allow(dead_code)]

use ybus_powerflow::case::Case;
use ybus_powerflow::islands::find_bridges;

/// The walk's seed and step
pub const WALK_SEED: u64 = 0x0005_EED0_FDA7;
pub const WALK_STEP: f64 = 0.02;
pub const WALK_STEPS: usize = 25;

/// The outage sequences' seed, and their most outages
pub const OUTAGE_SEED: u64 = 0x0000_0D7A_6E5E;
pub const OUTAGES: usize = 50;

/// SplitMix64 (Steele, Lea and Flood, 2014), for a walk that is the same on every machine
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [−1, 1)
    pub fn symmetric(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }
}

/// The sequence's steps: each bus's load factor, and the mean factor for the CSV
pub fn steps(sequence: &str, buses: usize) -> Vec<(Vec<f64>, f64)> {
    match sequence {
        "sweep" => (0..=20)
            .map(|k| {
                let factor = 0.90 + 0.01 * k as f64;
                (vec![factor; buses], factor)
            })
            .collect(),
        _ => {
            let mut random = SplitMix64(WALK_SEED);
            let mut factors = vec![1.0; buses];
            let mut steps = vec![(factors.clone(), 1.0)];
            for _ in 1..WALK_STEPS {
                for factor in &mut factors {
                    *factor *= 1.0 + WALK_STEP * random.symmetric();
                }
                let mean = factors.iter().sum::<f64>() / buses as f64;
                steps.push((factors.clone(), mean));
            }
            steps
        }
    }
}

/// At most `n` of `rows`, a seeded sample (a partial Fisher–Yates shuffle by SplitMix64), in
///   case order
pub fn sample(mut rows: Vec<usize>, n: usize) -> Vec<usize> {
    let mut random = SplitMix64(OUTAGE_SEED);
    let n = n.min(rows.len());
    for i in 0..n {
        let j = i + (random.next() % (rows.len() - i) as u64) as usize;
        rows.swap(i, j);
    }
    rows.truncate(n);
    rows.sort_unstable();
    rows
}

pub fn loaded(case: &Case, factors: &[f64]) -> Case {
    let mut case = case.clone();
    for (bus, &factor) in case.bus.iter_mut().zip(factors) {
        bus.pd_mw *= factor;
        bus.qd_mvar *= factor;
    }
    case
}

/// The outages a survey samples from: the in-service branches that aren't bridges, whose
/// outage would island part of the network (`branches`), or the in-service generators
pub fn outage_rows(case: &Case, branches: bool) -> Result<Vec<usize>, String> {
    if branches {
        let bridges = find_bridges(case).map_err(|err| err.to_string())?;
        Ok((0..case.branch.len())
            .filter(|&row| case.branch[row].status && bridges.binary_search(&row).is_err())
            .collect())
    } else {
        Ok((0..case.generator.len())
            .filter(|&row| case.generator[row].status)
            .collect())
    }
}
