//! The Newton system's matrix `W = H + Jhᵀ D Jh`, kept planned across a solve's iterations.
//!
//! Both solvers form, at every iteration, the Hessian `H` plus the inequalities' Gram matrix
//! `Jhᵀ D Jh` for a diagonal scaling `D` (MIPS's `Z⁻¹ diag(μ)`, the predictor-corrector method's
//! `μ/z`). The patterns of `H` and `Jh` rarely change during a solve, so [`WeightedGram`] plans
//! the transpose, the product and the sum once ([`ProductPlan`], [`SumPlan`]), and recomputes
//! `W` from new values with the same floating-point operations in the same order: bit for bit
//! what `jh.transpose()`, `matmul` and `try_add` give. A pattern that does change is planned
//! afresh.

use ybus_sparse::compressed_matrix::{CsrMatrix, ProductPlan, SumPlan};

/// A matrix's pattern, to tell whether a plan made for it still applies
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Pattern {
    shape: (usize, usize),
    pointers: Vec<usize>,
    indices: Vec<usize>,
}

impl Pattern {
    pub(crate) fn of(m: &CsrMatrix<f64>) -> Self {
        Pattern {
            shape: m.shape(),
            pointers: m.pointers().to_vec(),
            indices: m.indices().to_vec(),
        }
    }

    /// Whether `m` has this pattern, compared without allocating
    pub(crate) fn matches(&self, m: &CsrMatrix<f64>) -> bool {
        self.shape == m.shape() && self.pointers == m.pointers() && self.indices == m.indices()
    }
}

/// The row of each stored entry of `m`
pub(crate) fn entry_rows(m: &CsrMatrix<f64>) -> Vec<usize> {
    m.pointers()
        .windows(2)
        .enumerate()
        .flat_map(|(row, bounds)| std::iter::repeat_n(row, bounds[1] - bounds[0]))
        .collect()
}

/// `Jh`'s transpose's entries from `Jh`'s: entry `k` of `jh.transpose().to_csr()` is
/// `jh.values()[sources[k]]`, found by transposing `Jh`'s slot numbers
fn transpose_sources(jh: &CsrMatrix<f64>) -> (CsrMatrix<f64>, Vec<usize>) {
    let (rows, cols) = jh.shape();
    let slots = CsrMatrix::<usize>::try_from_raw_parts(
        rows,
        cols,
        jh.pointers().to_vec(),
        jh.indices().to_vec(),
        (1..=jh.num_stored()).collect(),
    )
    .expect("Jh's own valid pattern");
    let traced = slots.transpose().to_csr();
    let sources = traced.values().iter().map(|&slot| slot - 1).collect();
    (traced.map(|_| 0.0), sources)
}

/// The plans for one pair of patterns of `H` and `Jh`
struct Plans {
    hessian: Pattern,
    jh: Pattern,
    jh_rows: Vec<usize>,
    jh_t_sources: Vec<usize>,
    product: ProductPlan,
    sum: SumPlan,
    /// `W`'s pattern, holding zeros
    w: CsrMatrix<f64>,
}

impl Plans {
    fn new(hessian: &CsrMatrix<f64>, jh: &CsrMatrix<f64>) -> Self {
        let (jh_t, jh_t_sources) = transpose_sources(jh);
        let product = jh_t.plan_product(jh).expect("nx × niq times niq × nx");
        let gram = product.multiply(&jh_t, &jh.map(|_| 0.0)).expect("planned");
        let sum = hessian.plan_sum(&gram).expect("both nx × nx");
        let w = sum.add(&hessian.map(|_| 0.0), &gram).expect("planned");
        Plans {
            hessian: Pattern::of(hessian),
            jh: Pattern::of(jh),
            jh_rows: entry_rows(jh),
            jh_t_sources,
            product,
            sum,
            w,
        }
    }
}

/// `W = H + Jhᵀ (D Jh)`, planned for the patterns it last saw (see the module documentation).
#[derive(Default)]
pub(crate) struct WeightedGram {
    plans: Option<Plans>,
}

impl WeightedGram {
    /// The row of each of `jh`'s stored entries, for scaling its rows: the scaled values the
    /// caller passes to `compute` are its own expression, entry by entry, so that they're bit
    /// for bit what its `scale_rows` gave.
    pub(crate) fn jh_rows(&mut self, hessian: &CsrMatrix<f64>, jh: &CsrMatrix<f64>) -> &[usize] {
        &self.plans_for(hessian, jh).jh_rows
    }

    /// `hessian.try_add(&jh.transpose().to_csr().matmul(&scaled_jh))`, bit for bit, where
    /// `scaled` holds `scaled_jh`'s values: `Jh`'s pattern with each row scaled.
    pub(crate) fn compute(
        &mut self,
        hessian: &CsrMatrix<f64>,
        jh: &CsrMatrix<f64>,
        scaled: &[f64],
    ) -> CsrMatrix<f64> {
        let plans = self.plans_for(hessian, jh);
        let jh_t: Vec<f64> = plans.jh_t_sources.iter().map(|&s| jh.values()[s]).collect();
        let mut gram = vec![0.0; plans.product.num_stored()];
        plans
            .product
            .multiply_values(&jh_t, scaled, &mut gram)
            .expect("the planned lengths");
        let mut w = plans.w.clone();
        plans
            .sum
            .add_values(hessian.values(), &gram, w.values_mut())
            .expect("the planned lengths");
        w
    }

    /// The plans for these patterns, made afresh if they differ from the last ones
    fn plans_for(&mut self, hessian: &CsrMatrix<f64>, jh: &CsrMatrix<f64>) -> &Plans {
        let stale = self
            .plans
            .as_ref()
            .is_none_or(|p| !p.hessian.matches(hessian) || !p.jh.matches(jh));
        if stale {
            self.plans = Some(Plans::new(hessian, jh));
        }
        self.plans.as_ref().expect("just made")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }

        /// Values that make the order of additions matter
        fn value(&mut self) -> f64 {
            match self.next() % 6 {
                0 => 0.0,
                1 => -0.0,
                2 => 1e16,
                3 => -1e16,
                _ => (self.next() % 2001) as f64 / 1000.0 - 1.0,
            }
        }
    }

    fn random(rows: usize, cols: usize, percent: u64, rng: &mut Lcg) -> CsrMatrix<f64> {
        let mut triplets = vec![];
        for i in 0..rows {
            for j in 0..cols {
                if rng.next() % 100 < percent {
                    triplets.push((i, j, rng.value()));
                }
            }
        }
        CsrMatrix::from_triplets(rows, cols, triplets).unwrap()
    }

    fn bits(m: &CsrMatrix<f64>) -> (Vec<usize>, Vec<usize>, Vec<u64>) {
        let values = m.values().iter().map(|v| v.to_bits()).collect();
        (m.pointers().to_vec(), m.indices().to_vec(), values)
    }

    #[test]
    fn the_planned_w_matches_try_add_and_matmul_bit_for_bit() {
        let mut rng = Lcg(3);
        let mut gram = WeightedGram::default();
        let (nx, niq) = (8, 12);
        let (h1, h2) = (random(nx, nx, 30, &mut rng), random(nx, nx, 30, &mut rng));
        let (j1, j2) = (random(niq, nx, 30, &mut rng), random(niq, nx, 30, &mut rng));
        // Phases of patterns, each with new values three times: plans made, kept, and made
        //   again when only Jh's pattern changes, only H's, or both
        for (phase, (h, jh)) in [(&h1, &j1), (&h1, &j2), (&h2, &j2), (&h1, &j1)]
            .into_iter()
            .enumerate()
        {
            for round in 0..3 {
                let (mut hessian, mut jh) = (h.clone(), jh.clone());
                hessian
                    .values_mut()
                    .iter_mut()
                    .for_each(|v| *v = rng.value());
                jh.values_mut().iter_mut().for_each(|v| *v = rng.value());
                let d: Vec<f64> = (0..niq).map(|_| rng.value()).collect();
                let direct = hessian
                    .try_add(
                        &jh.clone()
                            .transpose()
                            .to_csr()
                            .matmul(&jh.scale_rows(&d).unwrap())
                            .unwrap(),
                    )
                    .unwrap();
                let rows = gram.jh_rows(&hessian, &jh).to_vec();
                let scaled: Vec<f64> = jh
                    .values()
                    .iter()
                    .zip(&rows)
                    .map(|(&v, &r)| d[r] * v)
                    .collect();
                let planned = gram.compute(&hessian, &jh, &scaled);
                assert_eq!(
                    bits(&planned),
                    bits(&direct),
                    "phase {phase}, round {round}"
                );
            }
        }
    }
}
