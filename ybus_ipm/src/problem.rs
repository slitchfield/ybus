//! The nonlinear programs the solvers take:
//!
//! minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`,
//!
//! where `f`, `g` and `h` come from a [`Problem`] and the linear constraints and bounds from a
//! [`Constraints`]. A [`QuadraticProblem`] supplies `f(x) = ½ xᵀ H x + cᵀ x` for quadratic and
//! linear programs.

use ybus_sparse::compressed_matrix::CsrMatrix;

/// The nonlinear constraints' values at a point and their Jacobians, one row per constraint:
/// `h` (`h(x) ≤ 0`) with `dh` (`nh × nx`), and `g` (`g(x) = 0`) with `dg` (`ng × nx`).
#[derive(Clone, Debug, PartialEq)]
pub struct NonlinearConstraints {
    /// The inequality constraints' values, `h(x)`
    pub h: Vec<f64>,
    /// The equality constraints' values, `g(x)`
    pub g: Vec<f64>,
    /// `∂h/∂x`, one row per inequality constraint
    pub dh: CsrMatrix<f64>,
    /// `∂g/∂x`, one row per equality constraint
    pub dg: CsrMatrix<f64>,
}

/// The multipliers of the nonlinear constraints, in their order: what a Lagrangian's Hessian
/// takes.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearMultipliers<'a> {
    /// One per equality constraint `g(x) = 0`
    pub eq: &'a [f64],
    /// One per inequality constraint `h(x) ≤ 0`
    pub ineq: &'a [f64],
}

/// A nonlinear program's functions: the objective, and optionally nonlinear constraints with
/// the Lagrangian's Hessian. Matrices are `nx × nx` with both triangles stored, or as noted.
pub trait Problem {
    /// The number of variables, `nx`
    fn dimension(&self) -> usize;

    /// `f(x)` and its gradient `∇f(x)`
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>);

    /// `∇²f(x)`, used when there are no nonlinear constraints; with them, the solvers take the
    /// whole Lagrangian's Hessian from [`Problem::lagrangian_hessian`] instead.
    fn objective_hessian(&self, x: &[f64]) -> CsrMatrix<f64>;

    /// [`Problem::objective_hessian`] written into `out`, which the solvers keep across
    /// iterations: on entry it holds what the previous call left (or an empty matrix), so an
    /// implementation whose pattern hasn't changed can refill its values in place instead of
    /// allocating. On return `out` must equal what `objective_hessian(x)` gives, pattern and
    /// values. The default assigns `objective_hessian(x)`.
    fn objective_hessian_into(&self, x: &[f64], out: &mut CsrMatrix<f64>) {
        *out = self.objective_hessian(x);
    }

    /// Whether the problem has nonlinear constraints, so that [`Problem::constraints`] and
    /// [`Problem::lagrangian_hessian`] are called. False by default.
    fn has_nonlinear_constraints(&self) -> bool {
        false
    }

    /// The nonlinear constraints at `x`. Called only if `has_nonlinear_constraints`.
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        let _ = x;
        unreachable!("a problem with nonlinear constraints implements constraints")
    }

    /// [`Problem::constraints`] written into `out`, kept across calls as
    /// [`Problem::objective_hessian_into`]'s is: on return it must equal what `constraints(x)`
    /// gives. The default assigns `constraints(x)`.
    fn constraints_into(&self, x: &[f64], out: &mut NonlinearConstraints) {
        *out = self.constraints(x);
    }

    /// The Lagrangian's Hessian `cost_mult · ∇²f(x) + Σ λᵢ ∇²gᵢ(x) + Σ μᵢ ∇²hᵢ(x)` at `x`,
    /// with the nonlinear constraints' multipliers. Called only if `has_nonlinear_constraints`.
    fn lagrangian_hessian(
        &self,
        x: &[f64],
        multipliers: NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        let _ = (x, multipliers, cost_mult);
        unreachable!("a problem with nonlinear constraints implements lagrangian_hessian")
    }

    /// [`Problem::lagrangian_hessian`] written into `out`, kept across calls as
    /// [`Problem::objective_hessian_into`]'s is: on return it must equal what
    /// `lagrangian_hessian(x, multipliers, cost_mult)` gives. The default assigns that.
    fn lagrangian_hessian_into(
        &self,
        x: &[f64],
        multipliers: NonlinearMultipliers<'_>,
        cost_mult: f64,
        out: &mut CsrMatrix<f64>,
    ) {
        *out = self.lagrangian_hessian(x, multipliers, cost_mult);
    }
}

impl NonlinearConstraints {
    /// No constraints over `nx` variables: an empty buffer for
    /// [`Problem::constraints_into`].
    pub fn empty(nx: usize) -> Self {
        let empty = || CsrMatrix::from_triplets(0, nx, vec![]).expect("an empty 0 × nx matrix");
        NonlinearConstraints {
            h: vec![],
            g: vec![],
            dh: empty(),
            dg: empty(),
        }
    }
}

/// Whether `a` and `b` have the same pattern (shape and stored positions), compared without
/// allocating: when an `_into` method may refill its buffer's values in place.
pub fn same_pattern(a: &CsrMatrix<f64>, b: &CsrMatrix<f64>) -> bool {
    a.shape() == b.shape() && a.pointers() == b.pointers() && a.indices() == b.indices()
}

/// The linear constraints `l ≤ A x ≤ u` and the bounds `xmin ≤ x ≤ xmax`, any of which may be
/// absent: an absent `A` means no linear constraints, an absent `l`, `u`, `xmin` or `xmax`
/// means unbounded that way, and an infinite entry the same. A row with `l = u` is an equality.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Constraints {
    /// `A`, `nA × nx`
    pub a: Option<CsrMatrix<f64>>,
    /// `l`, one per row of `A`
    pub l: Option<Vec<f64>>,
    /// `u`, one per row of `A`
    pub u: Option<Vec<f64>>,
    /// `xmin`, one per variable
    pub xmin: Option<Vec<f64>>,
    /// `xmax`, one per variable
    pub xmax: Option<Vec<f64>>,
}

/// A quadratic (or, with no `H`, linear) objective `f(x) = ½ xᵀ H x + cᵀ x`, for MIPS's
/// quadratic-program wrapper, as MATPOWER's `qps_mips` builds it.
#[derive(Clone, Debug, PartialEq)]
pub struct QuadraticProblem {
    /// `H`, symmetric with both triangles stored, or `None` for a linear program
    pub h: Option<CsrMatrix<f64>>,
    /// `c`
    pub c: Vec<f64>,
}

impl Problem for QuadraticProblem {
    fn dimension(&self) -> usize {
        self.c.len()
    }

    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        let linear: f64 = self.c.iter().zip(x).map(|(c, x)| c * x).sum();
        match &self.h {
            Some(h) => {
                let hx = h.mul_vec(x).expect("H is nx × nx");
                let quadratic: f64 = x.iter().zip(&hx).map(|(x, hx)| x * hx).sum();
                let gradient = hx.iter().zip(&self.c).map(|(hx, c)| hx + c).collect();
                (0.5 * quadratic + linear, gradient)
            }
            None => (linear, self.c.clone()),
        }
    }

    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        let n = self.c.len();
        match &self.h {
            Some(h) => h.clone(),
            None => CsrMatrix::from_triplets(n, n, vec![]).expect("an empty n × n matrix"),
        }
    }

    /// `H`'s values copied into `out` when it already has `H`'s pattern, without allocating
    fn objective_hessian_into(&self, x: &[f64], out: &mut CsrMatrix<f64>) {
        match &self.h {
            Some(h) if same_pattern(h, out) => out.values_mut().copy_from_slice(h.values()),
            None if out.shape() == (self.c.len(), self.c.len()) && out.num_stored() == 0 => {}
            _ => *out = self.objective_hessian(x),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn objective_hessian_into_equals_objective_hessian_from_any_buffer() {
        let h =
            CsrMatrix::from_triplets(2, 2, vec![(0, 0, 2.0), (0, 1, -1.0), (1, 0, -1.0)]).unwrap();
        let qp = QuadraticProblem {
            h: Some(h.clone()),
            c: vec![1.0, 0.0],
        };
        let lp = QuadraticProblem {
            h: None,
            c: vec![1.0, 0.0],
        };
        for problem in [&qp, &lp] {
            // An empty buffer, a buffer of another pattern, then the buffer the last call left
            //   with stale values
            let mut out = CsrMatrix::from_triplets(0, 0, vec![]).unwrap();
            for buffer in [
                CsrMatrix::from_triplets(0, 0, vec![]).unwrap(),
                CsrMatrix::from_triplets(2, 2, vec![(1, 1, 9.0)]).unwrap(),
            ] {
                out = buffer;
                problem.objective_hessian_into(&[0.0, 0.0], &mut out);
                assert_eq!(out, problem.objective_hessian(&[0.0, 0.0]));
                out.values_mut().iter_mut().for_each(|v| *v = 7.0);
                problem.objective_hessian_into(&[0.0, 0.0], &mut out);
                assert_eq!(out, problem.objective_hessian(&[0.0, 0.0]));
            }
            let _ = out;
        }
    }
}
