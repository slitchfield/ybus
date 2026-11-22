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

/// The restoration phases of a solve's record
fn restorations(s: &Solution) -> Vec<&MehrotraIteration> {
    records(s)
        .iter()
        .filter(|r| r.mode == IterationMode::Restoration)
        .collect()
}

/// Wächter and Biegler's counterexample ("Failure of global convergence for a class of interior
///   point methods for nonlinear programming," *Math. Program.* 88(3):565–574, 2000,
///   doi:10.1007/pl00011386): `min x₀` subject to `x₀² − x₁ − 1 = 0`, here nonlinear,
///   `x₀ − x₂ = 0.5`, here linear, and `x₁, x₂ ≥ 0`. Its solution is `(1, 0, 0.5)`
struct Counterexample;

impl Problem for Counterexample {
    fn dimension(&self) -> usize {
        3
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        (x[0], vec![1.0, 0.0, 0.0])
    }
    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        CsrMatrix::from_triplets(3, 3, vec![]).unwrap()
    }
    fn has_nonlinear_constraints(&self) -> bool {
        true
    }
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        NonlinearConstraints {
            g: vec![x[0] * x[0] - x[1] - 1.0],
            h: vec![],
            dg: dense(&[&[2.0 * x[0], -1.0, 0.0]]),
            dh: CsrMatrix::from_triplets(0, 3, vec![]).unwrap(),
        }
    }
    fn lagrangian_hessian(&self, _: &[f64], m: NonlinearMultipliers<'_>, _: f64) -> CsrMatrix<f64> {
        CsrMatrix::from_triplets(3, 3, vec![(0, 0, 2.0 * m.eq[0])]).unwrap()
    }
}

fn counterexample_constraints() -> Constraints {
    Constraints {
        a: Some(dense(&[&[1.0, 0.0, -1.0]])),
        l: Some(vec![0.5]),
        u: Some(vec![0.5]),
        xmin: Some(vec![f64::NEG_INFINITY, 0.0, 0.0]),
        ..Constraints::default()
    }
}

#[test]
fn on_wachter_and_bieglers_counterexample_restoration_reaches_the_solution() {
    // From x₀ = −0.5, the method alone runs into x₁ = x₂ = 0 at a point that isn't feasible,
    //   as that paper shows such methods do. The filter's line search fails there, and
    //   restoration finds a feasible region the method then converges in
    let x0 = [-0.5, 1.0, 1.0];
    let mut off = MehrotraOptions::default();
    off.kkt_factorization = KktFactorization::Ldlt;
    let alone = mehrotra(&Counterexample, &counterexample_constraints(), &x0, &off).unwrap();
    assert!(!alone.converged(), "without restoration: {:?}", alone.x);
    let s = mehrotra(
        &Counterexample,
        &counterexample_constraints(),
        &x0,
        &monotone(),
    )
    .unwrap();
    assert!(s.converged(), "{:?} at {:?}", s.status, s.x);
    let phases = restorations(&s);
    assert!(!phases.is_empty(), "it restores");
    assert!(phases.iter().all(|r| r.restoration_iterations > 0));
    // Entering added the point's margin to the filter, so that restoration can't return to it:
    //   here one entry more than before, since it dominates none
    let history = records(&s);
    let entered = (history.iter())
        .position(|r| r.mode == IterationMode::Restoration)
        .unwrap();
    assert_eq!(
        history[entered].filter_size,
        history[entered - 1].filter_size + 1
    );
    for (x, solution) in s.x.iter().zip([1.0, 0.0, 0.5]) {
        assert!((x - solution).abs() < 1e-6, "{:?}", s.x);
    }
    assert!(kkt(&Counterexample, &counterexample_constraints(), &s) < 1e-6);
    // The solve's count includes restoration's iterations, which its own records don't show
    let own = records(&s).len() - 1 - phases.len();
    assert_eq!(
        s.iterations,
        own + phases.len()
            + phases
                .iter()
                .map(|r| r.restoration_iterations)
                .sum::<usize>()
    );
}

#[test]
fn from_further_left_restoration_stops_at_a_local_minimizer_of_the_violation() {
    // From x₀ = −2, restoration is entered more than once, and the last ends at (−1, 0, 0).
    //   With x₁, x₂ ≥ 0 the violation near there is 1.5 − x₀ − x₀² for x₀ > −1 and
    //   0.5 − x₀ below, a strict local minimizer of it, so the solve ends there as locally
    //   infeasible
    let s = mehrotra(
        &Counterexample,
        &counterexample_constraints(),
        &[-2.0, 1.0, 1.0],
        &monotone(),
    )
    .unwrap();
    assert_eq!(s.status, ybus_ipm::Status::NumericalFailure);
    let last = records(&s).last().unwrap();
    assert_eq!(last.mode, IterationMode::Restoration);
    assert!(restorations(&s).len() > 1);
    assert!((s.x[0] + 1.0).abs() < 0.1, "{:?}", s.x);
}

/// `min x₀² + x₁²` subject to `x₀² + x₁² + 1 = 0` and `x₀ ≥ 0`, which nothing satisfies
struct Infeasible;

impl Problem for Infeasible {
    fn dimension(&self) -> usize {
        2
    }
    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        (x[0] * x[0] + x[1] * x[1], vec![2.0 * x[0], 2.0 * x[1]])
    }
    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        dense(&[&[2.0, 0.0], &[0.0, 2.0]])
    }
    fn has_nonlinear_constraints(&self) -> bool {
        true
    }
    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        NonlinearConstraints {
            g: vec![x[0] * x[0] + x[1] * x[1] + 1.0],
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
        let d = 2.0 * cost_mult + 2.0 * m.eq[0];
        dense(&[&[d, 0.0], &[0.0, d]])
    }
}

#[test]
fn an_infeasible_problem_ends_in_restoration() {
    let constraints = Constraints {
        xmin: Some(vec![0.0, f64::NEG_INFINITY]),
        ..Constraints::default()
    };
    let s = mehrotra(&Infeasible, &constraints, &[1.0, 1.0], &monotone()).unwrap();
    assert_eq!(s.status, ybus_ipm::Status::NumericalFailure);
    assert_eq!(records(&s).last().unwrap().mode, IterationMode::Restoration);
    // Near the violation's minimizer, the origin, where it's 1
    assert!(s.x.iter().all(|x| x.abs() < 0.1), "{:?}", s.x);
}

#[test]
fn problems_whose_line_search_succeeds_never_restore() {
    for x0 in [[1.0, 5.0, 5.0, 1.0], [1.0; 4], [2.0; 4]] {
        let s = mehrotra(&Hs71, &hs71_bounds(), &x0, &monotone()).unwrap();
        assert!(restorations(&s).is_empty(), "from {x0:?}");
        assert_eq!(s.iterations, records(&s).len() - 1);
    }
}
