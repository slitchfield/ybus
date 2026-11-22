//! Feasibility restoration's problem (`LineSearch::MonotoneFilter`, `docs/LineSearchChecklist.md`,
//! PR 2): Wächter and Biegler's restoration phase (2006, §3.3) in the form this method solves,
//! so that the method itself solves it.
//!
//! At the reference point `x_R`, where the line search found no acceptable step, it minimizes the
//! constraint violation, kept near `x_R`:
//!
//! - **The variables:** `x`, then `p_g, n_g ≥ 0` per nonlinear equality, `p_h ≥ 0` per nonlinear
//!   inequality, and `p_a, n_a ≥ 0` per linear row, in that order.
//! - **The objective:** `ρ (Σ p + Σ n) + (ζ/2) ‖D_R (x − x_R)‖²`, with `ρ = 1000`, `ζ = √μ̄` and
//!   `D_R = diag(min(1, 1/|x_R,i|))`.
//! - **The rows:** `g(x) − p_g + n_g = 0`, `h(x) − p_h ≤ 0`, and `l ≤ A x − p_a + n_a ≤ u`. An
//!   inequality needs only its excess `p`, since the method's own slack takes the rest.
//! - **The bounds:** `x`'s own, and `p, n ≥ 0`.

use ybus_sparse::compressed_matrix::CsrMatrix;

use crate::problem::{Constraints, NonlinearConstraints, NonlinearMultipliers, Problem};

/// The violation's weight in the restoration problem's objective, `ρ`
pub(super) const RHO: f64 = 1000.0;

/// The restoration problem for `problem`, at the reference point `x_r`
pub(super) struct Restoration<'a, P: Problem + ?Sized> {
    problem: &'a P,
    nx: usize,
    /// The nonlinear equalities and inequalities, and the linear rows
    neq: usize,
    niq: usize,
    nlin: usize,
    x_r: Vec<f64>,
    /// `ζ D_R²`, the proximity term's weights
    weights: Vec<f64>,
}

impl<'a, P: Problem + ?Sized> Restoration<'a, P> {
    /// At `x_r` and barrier parameter `mu`, for a problem with `neq` nonlinear equalities, `niq`
    ///   nonlinear inequalities and `nlin` linear rows
    pub(super) fn new(
        problem: &'a P,
        x_r: &[f64],
        mu: f64,
        (neq, niq, nlin): (usize, usize, usize),
    ) -> Self {
        let zeta = mu.sqrt();
        let weights = x_r
            .iter()
            .map(|x| {
                let d = 1.0_f64.min(1.0 / x.abs());
                zeta * d * d
            })
            .collect();
        Restoration {
            problem,
            nx: x_r.len(),
            neq,
            niq,
            nlin,
            x_r: x_r.to_vec(),
            weights,
        }
    }

    /// The restoration problem's variables before the linear rows' `p_a` and `n_a`
    fn linear_start(&self) -> usize {
        self.nx + 2 * self.neq + self.niq
    }

    /// The linear rows, bounds and starting point of the restoration problem, from the
    ///   original's `constraints` and the residuals at `x_R`: each equality's and linear row's
    ///   `p` and `n` by the paper's closed form at barrier parameter `mu` (eq. 30), so that its
    ///   row holds exactly, and each inequality's `p` as its excess
    pub(super) fn start(
        &self,
        constraints: &Constraints,
        (g, h, linear): (&[f64], &[f64], &[f64]),
        mu: f64,
    ) -> (Constraints, Vec<f64>) {
        let (nx, nlin) = (self.nx, self.nlin);
        let n = self.linear_start() + 2 * nlin;
        // The closed form: n = (μ − ρc)/(2ρ) + √(((μ − ρc)/(2ρ))² + μc/(2ρ)), p = c + n
        let pair = |c: f64| -> (f64, f64) {
            let b = (mu - RHO * c) / (2.0 * RHO);
            let neg = b + (b * b + mu * c / (2.0 * RHO)).sqrt();
            (c + neg, neg)
        };
        let mut y = self.x_r.clone();
        let (pg, ng): (Vec<f64>, Vec<f64>) = g.iter().map(|&c| pair(c)).unzip();
        y.extend(pg);
        y.extend(ng);
        y.extend(h.iter().map(|&h| h.max(0.0)));
        // A linear row's residual is how far `A x` lies outside `[l, u]`
        let (l, u) = (constraints.l.as_deref(), constraints.u.as_deref());
        let residual = |i: usize| -> f64 {
            let v = linear[i];
            let lo = l.map_or(f64::NEG_INFINITY, |l| l[i]);
            let hi = u.map_or(f64::INFINITY, |u| u[i]);
            if v < lo {
                v - lo
            } else if v > hi {
                v - hi
            } else {
                0.0
            }
        };
        let (pa, na): (Vec<f64>, Vec<f64>) = (0..nlin).map(|i| pair(residual(i))).unzip();
        y.extend(pa);
        y.extend(na);
        // A x − p_a + n_a, beside the zeros of p_g, n_g and p_h
        let a = constraints.a.as_ref().map(|a| {
            let mut triplets = vec![];
            for i in 0..nlin {
                for p in a.pointers()[i]..a.pointers()[i + 1] {
                    triplets.push((i, a.indices()[p], a.values()[p]));
                }
                triplets.push((i, self.linear_start() + i, -1.0));
                triplets.push((i, self.linear_start() + nlin + i, 1.0));
            }
            CsrMatrix::from_triplets(nlin, n, triplets).expect("in range")
        });
        let bound = |given: Option<&Vec<f64>>, absent: f64| -> Vec<f64> {
            let mut b = given.cloned().unwrap_or_else(|| vec![absent; nx]);
            b.resize(n, 0.0);
            b
        };
        let mut xmax = bound(constraints.xmax.as_ref(), f64::INFINITY);
        for v in &mut xmax[nx..] {
            *v = f64::INFINITY;
        }
        let restoration = Constraints {
            a,
            l: constraints.l.clone(),
            u: constraints.u.clone(),
            xmin: Some(bound(constraints.xmin.as_ref(), f64::NEG_INFINITY)),
            xmax: Some(xmax),
        };
        (restoration, y)
    }

    /// `x`'s part of a restoration point
    pub(super) fn x<'y>(&self, y: &'y [f64]) -> &'y [f64] {
        &y[..self.nx]
    }

    /// A matrix's entries placed in a `rows × cols` one, then the extra entries
    fn widened(
        m: &CsrMatrix<f64>,
        (rows, cols): (usize, usize),
        extra: impl Iterator<Item = (usize, usize, f64)>,
    ) -> CsrMatrix<f64> {
        let mut triplets: Vec<(usize, usize, f64)> = vec![];
        for i in 0..m.shape().0 {
            for p in m.pointers()[i]..m.pointers()[i + 1] {
                triplets.push((i, m.indices()[p], m.values()[p]));
            }
        }
        triplets.extend(extra);
        CsrMatrix::from_triplets(rows, cols, triplets).expect("in range")
    }
}

impl<P: Problem + ?Sized> Problem for Restoration<'_, P> {
    fn dimension(&self) -> usize {
        self.linear_start() + 2 * self.nlin
    }

    fn objective(&self, y: &[f64]) -> (f64, Vec<f64>) {
        let nx = self.nx;
        let mut f = 0.0;
        let mut df = vec![RHO; y.len()];
        for i in 0..nx {
            let d = y[i] - self.x_r[i];
            f += 0.5 * self.weights[i] * d * d;
            df[i] = self.weights[i] * d;
        }
        f += RHO * y[nx..].iter().sum::<f64>();
        (f, df)
    }

    fn objective_hessian(&self, y: &[f64]) -> CsrMatrix<f64> {
        let triplets: Vec<(usize, usize, f64)> = self
            .weights
            .iter()
            .enumerate()
            .map(|(i, &w)| (i, i, w))
            .collect();
        CsrMatrix::from_triplets(y.len(), y.len(), triplets).expect("in range")
    }

    fn has_nonlinear_constraints(&self) -> bool {
        self.problem.has_nonlinear_constraints()
    }

    fn constraints(&self, y: &[f64]) -> NonlinearConstraints {
        let (nx, neq, niq) = (self.nx, self.neq, self.niq);
        let n = self.dimension();
        let c = self.problem.constraints(&y[..nx]);
        let (pg, ng, ph) = (nx, nx + neq, nx + 2 * neq);
        NonlinearConstraints {
            g: (0..neq).map(|i| c.g[i] - y[pg + i] + y[ng + i]).collect(),
            h: (0..niq).map(|i| c.h[i] - y[ph + i]).collect(),
            dg: Self::widened(
                &c.dg,
                (neq, n),
                (0..neq).flat_map(|i| [(i, pg + i, -1.0), (i, ng + i, 1.0)]),
            ),
            dh: Self::widened(&c.dh, (niq, n), (0..niq).map(|i| (i, ph + i, -1.0))),
        }
    }

    fn lagrangian_hessian(
        &self,
        y: &[f64],
        multipliers: NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        // The constraints' curvature alone, from the original's Hessian without its objective,
        //   and the proximity term's, which the objective's scale multiplies
        let curvature = self
            .problem
            .lagrangian_hessian(&y[..self.nx], multipliers, 0.0);
        let n = self.dimension();
        Self::widened(
            &curvature,
            (n, n),
            self.weights
                .iter()
                .enumerate()
                .map(|(i, &w)| (i, i, cost_mult * w)),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `min x₀ + x₁` subject to `x₀² + x₁ − 1 = 0` and `x₀ − x₁ ≤ 0`
    struct Small;

    impl Problem for Small {
        fn dimension(&self) -> usize {
            2
        }
        fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
            (x[0] + x[1], vec![1.0, 1.0])
        }
        fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
            CsrMatrix::from_triplets(2, 2, vec![]).unwrap()
        }
        fn has_nonlinear_constraints(&self) -> bool {
            true
        }
        fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
            NonlinearConstraints {
                g: vec![x[0] * x[0] + x[1] - 1.0],
                h: vec![x[0] - x[1]],
                dg: CsrMatrix::from_triplets(1, 2, vec![(0, 0, 2.0 * x[0]), (0, 1, 1.0)]).unwrap(),
                dh: CsrMatrix::from_triplets(1, 2, vec![(0, 0, 1.0), (0, 1, -1.0)]).unwrap(),
            }
        }
        fn lagrangian_hessian(
            &self,
            _: &[f64],
            m: NonlinearMultipliers<'_>,
            _: f64,
        ) -> CsrMatrix<f64> {
            CsrMatrix::from_triplets(2, 2, vec![(0, 0, 2.0 * m.eq[0])]).unwrap()
        }
    }

    fn dense(m: &CsrMatrix<f64>) -> Vec<Vec<f64>> {
        let (rows, cols) = m.shape();
        let mut d = vec![vec![0.0; cols]; rows];
        for (i, row) in d.iter_mut().enumerate() {
            for p in m.pointers()[i]..m.pointers()[i + 1] {
                row[m.indices()[p]] += m.values()[p];
            }
        }
        d
    }

    #[test]
    fn its_jacobians_match_finite_differences() {
        let r = Restoration::new(&Small, &[2.0, 3.0], 0.01, (1, 1, 0));
        assert_eq!(r.dimension(), 5);
        let y = [0.7, -0.4, 0.2, 0.3, 0.5];
        let c = r.constraints(&y);
        let (dg, dh) = (dense(&c.dg), dense(&c.dh));
        let step = 1e-7;
        for j in 0..5 {
            let mut yp = y;
            yp[j] += step;
            let cp = r.constraints(&yp);
            assert!(
                ((cp.g[0] - c.g[0]) / step - dg[0][j]).abs() < 1e-5,
                "∂g/∂y{j}"
            );
            assert!(
                ((cp.h[0] - c.h[0]) / step - dh[0][j]).abs() < 1e-5,
                "∂h/∂y{j}"
            );
        }
        let (f, df) = r.objective(&y);
        for (j, &d) in df.iter().enumerate() {
            let mut yp = y;
            yp[j] += step;
            assert!(
                ((r.objective(&yp).0 - f) / step - d).abs() < 1e-4,
                "∂f/∂y{j}"
            );
        }
    }

    #[test]
    fn its_start_satisfies_each_equality_row_exactly_with_positive_p_and_n() {
        let x_r = [2.0, 3.0];
        let r = Restoration::new(&Small, &x_r, 0.01, (1, 1, 0));
        let c = Small.constraints(&x_r);
        let (restoration, y) = r.start(&Constraints::default(), (&c.g, &c.h, &[]), 0.01);
        assert_eq!(y.len(), 5);
        assert!(y[2..].iter().all(|v| *v >= 0.0), "{y:?}");
        assert!(y[2] > 0.0 && y[3] > 0.0);
        let rc = r.constraints(&y);
        assert!(
            rc.g[0].abs() < 1e-12,
            "g − p + n = 0 at the start: {}",
            rc.g[0]
        );
        assert!(rc.h[0] <= 0.0, "h − p ≤ 0 at the start: {}", rc.h[0]);
        assert_eq!(r.x(&y), &x_r);
        // p and n are bounded below by 0, and above by nothing
        let xmin = restoration.xmin.unwrap();
        assert_eq!(&xmin[2..], &[0.0, 0.0, 0.0]);
        assert!(restoration.xmax.unwrap()[2..]
            .iter()
            .all(|v| v.is_infinite()));
    }

    #[test]
    fn a_linear_rows_p_and_n_take_its_distance_outside_its_limits() {
        let x_r = [2.0, 3.0];
        let r = Restoration::new(&Small, &x_r, 0.01, (1, 1, 1));
        let c = Small.constraints(&x_r);
        // x₀ + x₁ ≤ 4, which x_R violates by 1
        let constraints = Constraints {
            a: Some(CsrMatrix::from_triplets(1, 2, vec![(0, 0, 1.0), (0, 1, 1.0)]).unwrap()),
            u: Some(vec![4.0]),
            ..Constraints::default()
        };
        let (restoration, y) = r.start(&constraints, (&c.g, &c.h, &[5.0]), 0.01);
        let a = restoration.a.unwrap();
        let ay: f64 = dense(&a)[0].iter().zip(&y).map(|(a, y)| a * y).sum();
        assert!((ay - 4.0).abs() < 1e-12, "A x − p + n on its limit: {ay}");
    }
}
