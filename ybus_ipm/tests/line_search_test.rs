// Tests the predictor-corrector method's line search (MehrotraOptions::line_search,
//   docs/LineSearchChecklist.md): its setup, and the iteration record it adds, which with
//   LineSearch::Off records only the constraint violation, leaving the method as it was.

mod common;

use common::{dense, Hs71, Nonlinear2, Nonlinear3};
use ybus_ipm::kkt::KktFactorization;
use ybus_ipm::mehrotra::{mehrotra, IterationMode, LineSearch, MehrotraIteration, MehrotraOptions};
use ybus_ipm::problem::{NonlinearConstraints, NonlinearMultipliers};
use ybus_ipm::solution::{kkt_residuals, History};
use ybus_ipm::{Constraints, InteriorPointError, Problem, Solution};
use ybus_sparse::compressed_matrix::CsrMatrix;

/// hs71's bounds, 1 ≤ x ≤ 5
fn hs71_bounds() -> Constraints {
    Constraints {
        xmin: Some(vec![1.0; 4]),
        xmax: Some(vec![5.0; 4]),
        ..Constraints::default()
    }
}

fn records(s: &Solution) -> &[MehrotraIteration] {
    match &s.history {
        History::Mehrotra(records) => records,
        History::Mips(_) => panic!("a MIPS history from the predictor-corrector method"),
    }
}

#[test]
fn a_filter_needs_ldlt() {
    let mut options = MehrotraOptions::default();
    options.line_search = LineSearch::MonotoneFilter;
    assert_eq!(options.kkt_factorization, KktFactorization::Lu);
    let lu = mehrotra(&Hs71, &hs71_bounds(), &[1.0, 5.0, 5.0, 1.0], &options);
    assert_eq!(lu.unwrap_err(), InteriorPointError::LineSearchNeedsLdlt);
    options.kkt_factorization = KktFactorization::Ldlt;
    assert!(mehrotra(&Hs71, &hs71_bounds(), &[1.0, 5.0, 5.0, 1.0], &options).is_ok());
}

#[test]
fn off_records_the_constraint_violation_and_nothing_of_a_filter() {
    let s = mehrotra(
        &Hs71,
        &hs71_bounds(),
        &[1.0, 5.0, 5.0, 1.0],
        &MehrotraOptions::default(),
    )
    .unwrap();
    assert!(s.converged());
    let history = records(&s);
    assert!(
        history[0].theta > 1.0,
        "the start violates hs71's constraints"
    );
    let last = history.last().unwrap();
    assert!(
        last.theta < 1e-5,
        "feasible at the solution: {}",
        last.theta
    );
    for record in history {
        assert_eq!(record.mode, IterationMode::Free);
        assert_eq!((record.barrier, record.backtracks), (0.0, 0));
        assert!(!record.second_order);
    }
}

/// LDLᵀ with the monotone filter
fn monotone() -> MehrotraOptions {
    let mut options = MehrotraOptions::default();
    options.kkt_factorization = KktFactorization::Ldlt;
    options.line_search = LineSearch::MonotoneFilter;
    options
}

#[test]
fn fixed_mode_lowers_its_barrier_parameter_monotonically_to_the_optimum() {
    let s = mehrotra(&Hs71, &hs71_bounds(), &[1.0, 5.0, 5.0, 1.0], &monotone()).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    assert!(
        (s.f - 17.014_017_3).abs() < 1e-6 * 17.0,
        "the optimum, not {}",
        s.f
    );
    let steps = &records(&s)[1..];
    assert!(steps.iter().all(|r| r.mode == IterationMode::Fixed));
    assert_eq!(steps[0].barrier_parameter, 0.1, "the paper's μ₀");
    assert!(
        steps
            .windows(2)
            .all(|pair| pair[1].barrier_parameter <= pair[0].barrier_parameter),
        "μ̄ never rises"
    );
    assert!(steps.last().unwrap().barrier_parameter < 1e-6);
    // φ is recorded at the step's μ̄, f − μ̄ Σ ln s, finite at an interior point
    assert!(steps.iter().all(|r| r.barrier.is_finite()));
}

/// A KKT point's residuals, its stationarity and feasibility, as `kkt_residuals` measures them
fn kkt<P: Problem>(problem: &P, constraints: &Constraints, s: &Solution) -> f64 {
    let r = kkt_residuals(problem, constraints, &s.x, &s.multipliers).unwrap();
    r.stationarity.max(r.feasibility)
}

#[test]
fn from_a_start_where_the_method_alone_stalls_the_filter_reaches_hs71s_optimum() {
    let x0 = [1.0; 4];
    let mut off = MehrotraOptions::default();
    off.kkt_factorization = KktFactorization::Ldlt;
    let alone = mehrotra(&Hs71, &hs71_bounds(), &x0, &off).unwrap();
    assert!(!alone.converged(), "without a line search it stalls from 1");
    let s = mehrotra(&Hs71, &hs71_bounds(), &x0, &monotone()).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    assert!(
        (s.f - 17.014_017_3).abs() < 1e-6 * 17.0,
        "the optimum, not {}",
        s.f
    );
    assert!(kkt(&Hs71, &hs71_bounds(), &s) < 1e-6);
    // Its steps that weren't Armijo steps grew the filter
    assert!(records(&s).iter().any(|r| r.filter_size > 0));
}

#[test]
fn from_2_the_filter_avoids_hs71s_local_maximum() {
    // Without a filter the method converges to the constrained local maximum (f = 131.26).
    //   The filter's steps reduce φ, so it reaches another KKT point instead: one of the
    //   local points, f = 27.146 at x = (1, 5, 1.449, 3.449), not the global optimum, which a
    //   line search doesn't promise
    let s = mehrotra(&Hs71, &hs71_bounds(), &[2.0; 4], &monotone()).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    assert!(s.f < 131.0, "not the local maximum: {}", s.f);
    assert!(kkt(&Hs71, &hs71_bounds(), &s) < 1e-6);
}

#[test]
fn the_nonconvex_test_problems_reach_their_optima() {
    let bounds = Constraints {
        xmin: Some(vec![0.0; 2]),
        ..Constraints::default()
    };
    let s = mehrotra(&Nonlinear2, &bounds, &[1.1, 0.0], &monotone()).unwrap();
    assert!(
        s.converged() && (s.f + 2.0).abs() < 1e-6,
        "{:?} {}",
        s.status,
        s.f
    );
    assert!(kkt(&Nonlinear2, &bounds, &s) < 1e-6);
    let none = Constraints::default();
    let s = mehrotra(&Nonlinear3, &none, &[1.0, 1.0, 0.0], &monotone()).unwrap();
    let optimum = -5.0 * 2f64.sqrt();
    assert!(
        s.converged() && (s.f - optimum).abs() < 1e-6,
        "{:?} {}",
        s.status,
        s.f
    );
    assert!(kkt(&Nonlinear3, &none, &s) < 1e-6);
}

/// The Maratos effect's example (Nocedal and Wright, *Numerical Optimization*, Example 15.4):
///   minimize `2(x₁² + x₂² − 1) − x₁` on the unit circle, whose solution is `(1, 0)`. Near it, a
///   full Newton step from a point on the circle raises both the constraint's violation and
///   the objective, so a line search without a second-order correction shortens steps that
///   would converge quadratically
struct Maratos;

impl Problem for Maratos {
    fn dimension(&self) -> usize {
        2
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        (
            2.0 * (x[0] * x[0] + x[1] * x[1] - 1.0) - x[0],
            vec![4.0 * x[0] - 1.0, 4.0 * x[1]],
        )
    }
    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        dense(&[&[4.0, 0.0], &[0.0, 4.0]])
    }
    fn has_nonlinear_constraints(&self) -> bool {
        true
    }
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        NonlinearConstraints {
            g: vec![x[0] * x[0] + x[1] * x[1] - 1.0],
            h: vec![],
            dg: dense(&[&[2.0 * x[0], 2.0 * x[1]]]),
            dh: CsrMatrix::from_triplets(0, 2, vec![]).unwrap(),
        }
    }
    fn lagrangian_hessian(
        &self,
        _: &[f64],
        m: NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        let d = 4.0 * cost_mult + 2.0 * m.eq[0];
        dense(&[&[d, 0.0], &[0.0, d]])
    }
}

#[test]
fn the_second_order_correction_keeps_full_steps_against_the_maratos_effect() {
    let none = Constraints::default();
    let angle: f64 = 0.1;
    let s = mehrotra(&Maratos, &none, &[angle.cos(), angle.sin()], &monotone()).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    assert!(
        (s.x[0] - 1.0).abs() < 1e-6 && s.x[1].abs() < 1e-6,
        "{:?}",
        s.x
    );
    let steps = &records(&s)[1..];
    assert!(
        steps.iter().any(|r| r.second_order),
        "the full step is rejected at least once, and corrected"
    );
    assert!(
        steps.iter().all(|r| r.backtracks == 0),
        "no step is shortened: {:?}",
        steps.iter().map(|r| r.backtracks).collect::<Vec<_>>()
    );
}
