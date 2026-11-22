// Tests the predictor-corrector method's line search (MehrotraOptions::line_search,
//   docs/LineSearchChecklist.md): its setup, and the iteration record it adds, which with
//   LineSearch::Off records only the constraint violation, leaving the method as it was.

mod common;

use common::Hs71;
use ybus_ipm::kkt::KktFactorization;
use ybus_ipm::mehrotra::{mehrotra, IterationMode, LineSearch, MehrotraIteration, MehrotraOptions};
use ybus_ipm::solution::History;
use ybus_ipm::{Constraints, InteriorPointError, Solution};

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
