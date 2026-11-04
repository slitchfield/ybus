//! The DC optimal power flow on small cases whose optimum can be worked out by hand: a binding
//! flow limit, a binding angle difference limit, piecewise linear costs, a dispatchable load,
//! and two islands with a reference bus each.

mod common;

use common::dc_fixture::{fixture, B, G, L};
use powerflow::case::{Case, GenCost};
use powerflow::opf::dc::dc_opf;
use powerflow::powerflow_error::PowerflowError;

fn linear(cost_per_mwh: f64) -> GenCost {
    GenCost::Polynomial {
        startup: 0.0,
        shutdown: 0.0,
        coeffs: vec![cost_per_mwh, 0.0],
    }
}

/// Asserts that `actual` is within `tolerance` of `expected`, entry by entry.
fn near(actual: &[f64], expected: &[f64], tolerance: f64, what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: lengths");
    for (k, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!((a - e).abs() <= tolerance, "{what}[{k}]: {a} vs {e}");
    }
}

/// MIPS stops at a relative tolerance of 1e-6 (5e-6 for feasibility), so hand-derived optima
/// are met to about that, scaled by the quantities' size
const MW: f64 = 1e-3;
const PRICE: f64 = 1e-3;

/// Two buses joined by a 0.1 p.u. line, a load at bus 2, a cheap generator at bus 1 and an
/// expensive one at bus 2.
fn two_bus(load_mw: f64, cheap: f64, expensive: f64) -> Case {
    let mut case = fixture(
        &[B(1, 3, 0.0, 0.0, 0.0), B(2, 1, load_mw, 0.0, 0.0)],
        &[G(1, 0.0, true), G(2, 0.0, true)],
        &[L(1, 2, 0.1, 0.0, 0.0, true)],
    );
    case.gencost = Some(vec![linear(cheap), linear(expensive)]);
    case
}

#[test]
fn a_binding_flow_limit_splits_the_prices() {
    // 100 MW at bus 2; the line carries at most 60 MW of the 10 $/MWh power, and the 30 $/MWh
    //   generator makes up the rest. The limit is worth the price difference, 20 $/MWh
    let mut case = two_bus(100.0, 10.0, 30.0);
    case.branch[0].rate_a_mva = 60.0;
    let s = dc_opf(&case).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    near(&s.pg_mw, &[60.0, 40.0], MW, "pg");
    near(&s.pf_mw, &[60.0], MW, "pf");
    near(&s.lam_p, &[10.0, 30.0], PRICE, "lam_p");
    near(&s.mu_sf, &[20.0], PRICE, "mu_sf");
    near(&s.mu_st, &[0.0], PRICE, "mu_st");
    assert!(
        (s.objective - (600.0 + 1200.0)).abs() < 1e-2,
        "{}",
        s.objective
    );
    // θ₂ = −60 MW / (100 MVA / 0.1 p.u.) = −0.06 rad
    near(&s.va_deg, &[0.0, -0.06f64.to_degrees()], 1e-5, "va");
}

#[test]
fn a_binding_angle_limit_acts_like_a_flow_limit() {
    // θ₁ − θ₂ ≤ 0.06 rad limits the line to 60 MW, as above. Its multiplier is 20 $/MWh times
    //   1000 MW/rad (100 MVA / 0.1 p.u.), reported per degree
    let mut case = two_bus(100.0, 10.0, 30.0);
    case.branch[0].angmin_deg = -360.0;
    case.branch[0].angmax_deg = 0.06f64.to_degrees();
    let s = dc_opf(&case).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    near(&s.pg_mw, &[60.0, 40.0], MW, "pg");
    near(&s.lam_p, &[10.0, 30.0], PRICE, "lam_p");
    near(
        &s.mu_angmax,
        &[20.0 * 1000.0 * std::f64::consts::PI / 180.0],
        1e-2,
        "mu_angmax",
    );
    near(&s.mu_angmin, &[0.0], PRICE, "mu_angmin");
    near(&s.mu_sf, &[0.0], PRICE, "mu_sf");
}

#[test]
fn piecewise_linear_costs_dispatch_segment_by_segment() {
    // Generator 1 costs 10 $/MWh up to 50 MW and 20 $/MWh beyond; generator 2 costs 15 $/MWh.
    //   For 80 MW, generator 1 runs its cheap segment and generator 2 the rest, at 15 $/MWh
    let mut case = two_bus(80.0, 0.0, 15.0);
    case.gencost.as_mut().unwrap()[0] = GenCost::PiecewiseLinear {
        startup: 0.0,
        shutdown: 0.0,
        points: vec![(0.0, 0.0), (50.0, 500.0), (100.0, 1500.0)],
    };
    let s = dc_opf(&case).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    near(&s.pg_mw, &[50.0, 30.0], MW, "pg");
    near(&s.lam_p, &[15.0, 15.0], PRICE, "lam_p");
    assert!(
        (s.objective - (500.0 + 450.0)).abs() < 1e-2,
        "{}",
        s.objective
    );
    // MIPS's own variables end with the helper, at generator 1's cost
    assert!(
        (s.solver.x[4] - 500.0).abs() < 1e-2,
        "y = {}",
        s.solver.x[4]
    );
}

#[test]
fn a_dispatchable_load_is_a_negative_generator() {
    // A load at bus 2 worth 25 $/MWh for up to 50 MW (a generator from −50 to 0 MW whose cost
    //   is 25 $/MWh of output), supplied by 30 MW at 10 $/MWh or more at 30 $/MWh. It takes
    //   the 30 cheap MW, and sets the price at its own value
    let mut case = fixture(
        &[B(1, 3, 0.0, 0.0, 0.0), B(2, 1, 0.0, 0.0, 0.0)],
        &[G(1, 0.0, true), G(1, 0.0, true), G(2, 0.0, true)],
        &[L(1, 2, 0.1, 0.0, 0.0, true)],
    );
    case.generator[0].pmax_mw = 30.0;
    case.generator[2].pmin_mw = -50.0;
    case.generator[2].pmax_mw = 0.0;
    case.gencost = Some(vec![linear(10.0), linear(30.0), linear(25.0)]);
    let s = dc_opf(&case).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    near(&s.pg_mw, &[30.0, 0.0, -30.0], MW, "pg");
    near(&s.lam_p, &[25.0, 25.0], PRICE, "lam_p");
    near(&s.mu_pmax, &[15.0, 0.0, 0.0], PRICE, "mu_pmax");
    near(&s.mu_pmin, &[0.0, 5.0, 0.0], PRICE, "mu_pmin");
    assert!(
        (s.objective - (300.0 - 750.0)).abs() < 1e-2,
        "{}",
        s.objective
    );
}

#[test]
fn each_island_balances_around_its_own_reference() {
    // Two islands, each with a reference bus: 40 MW at 10 $/MWh, and 20 MW at 20 $/MWh with
    //   its reference angle at 10°. An out-of-service generator and branch report zeros, and
    //   an isolated bus keeps its angle
    let mut case = fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 40.0, 0.0, 0.0),
            B(3, 3, 0.0, 0.0, 10.0),
            B(4, 1, 20.0, 0.0, 0.0),
            B(5, 4, 0.0, 0.0, 7.0),
        ],
        &[G(1, 0.0, true), G(3, 0.0, true), G(2, 0.0, false)],
        &[
            L(1, 2, 0.1, 0.0, 0.0, true),
            L(3, 4, 0.1, 0.0, 0.0, true),
            L(2, 3, 0.1, 0.0, 0.0, false),
        ],
    );
    case.gencost = Some(vec![linear(10.0), linear(20.0), linear(1.0)]);
    let s = dc_opf(&case).unwrap();
    assert!(s.converged(), "{:?}", s.status);
    near(&s.pg_mw, &[40.0, 20.0, 0.0], MW, "pg");
    near(&s.pf_mw, &[40.0, 20.0, 0.0], MW, "pf");
    near(&s.lam_p, &[10.0, 10.0, 20.0, 20.0, 0.0], PRICE, "lam_p");
    assert_eq!(s.va_deg[0], 0.0);
    assert!((s.va_deg[2] - 10.0).abs() < 1e-12, "{}", s.va_deg[2]);
    assert_eq!(s.va_deg[4], 7.0);
    assert!((s.objective - 800.0).abs() < 1e-2, "{}", s.objective);
}

#[test]
fn a_case_without_a_ref_bus_or_costs_is_refused() {
    let mut case = two_bus(10.0, 10.0, 30.0);
    case.bus[0].bus_type = powerflow::case::BusType::PV;
    assert_eq!(dc_opf(&case).unwrap_err(), PowerflowError::NoReferenceBus);
    let mut case = two_bus(10.0, 10.0, 30.0);
    case.gencost = None;
    assert_eq!(
        dc_opf(&case).unwrap_err(),
        PowerflowError::MissingField { field: "gencost" }
    );
}
