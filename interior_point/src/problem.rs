//! The nonlinear programs the solvers take:
//!
//! minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`,
//!
//! where `f`, `g` and `h` come from a [`Problem`] and the linear constraints and bounds from a
//! [`Constraints`]. A [`QuadraticProblem`] supplies `f(x) = ½ xᵀ H x + cᵀ x` for quadratic and
//! linear programs.

use sparse_mat::compressed_matrix::CsrMatrix;

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
}
