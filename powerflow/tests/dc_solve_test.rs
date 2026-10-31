// Tests for the DC power flow solve on hand-built fixtures, covering the paths the MATPOWER
//   reference cases don't: phase shifters, shunt conductance, isolated buses, out-of-service
//   elements, several reference buses, the LU fallback and the errors. Expected values are
//   worked out by hand.

mod common;
use common::dc_fixture::{assert_physical, fixture, B, G, L};
use powerflow::case::Case;
use powerflow::dc::{dc_power_flow, DcModel, DcSolution, DcSolver};
use powerflow::powerflow_error::PowerflowError;

use std::f64::consts::PI;

/// Hand-computed values involve one or two divisions, so they agree to within a few ulps
const TOLERANCE: f64 = 1e-12;

fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} vs {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert!(
            (a - e).abs() <= TOLERANCE * e.abs().max(1.0),
            "{actual:?} vs {expected:?}"
        );
    }
}

fn degrees(radians: f64) -> f64 {
    radians * 180.0 / PI
}

/// The triangle of `dc_model_test`: susceptances 10 (1–2), 5 (2–3) and 4 (1–3); 50 MW at bus
/// 2 and 30 MW at bus 3, served by the generator at reference bus 1. Its file Pg is 80, the
/// balanced value.
fn triangle() -> (Vec<B>, Vec<G>, Vec<L>) {
    (
        vec![
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 50.0, 0.0, 0.0),
            B(3, 1, 30.0, 0.0, 0.0),
        ],
        vec![G(1, 80.0, true)],
        vec![
            L(1, 2, 0.1, 0.0, 0.0, true),
            L(2, 3, 0.2, 0.0, 0.0, true),
            L(1, 3, 0.25, 0.0, 0.0, true),
        ],
    )
}

fn case(tables: &(Vec<B>, Vec<G>, Vec<L>)) -> Case {
    fixture(&tables.0, &tables.1, &tables.2)
}

fn run(tables: &(Vec<B>, Vec<G>, Vec<L>)) -> Result<DcSolution, PowerflowError> {
    dc_power_flow(&case(tables))
}

#[test]
fn triangle_by_hand() {
    let tables = triangle();
    let solution = run(&tables).unwrap();
    assert_close(
        &solution.va_deg,
        &[0.0, degrees(-6.0 / 110.0), degrees(-7.0 / 110.0)],
    );
    assert_close(
        &solution.pf_mw,
        &[6000.0 / 110.0, 500.0 / 110.0, 2800.0 / 110.0],
    );
    assert_close(&solution.pg_mw, &[80.0]);
    assert_eq!(solution.solver, DcSolver::Cholesky);
    assert_physical(&case(&tables), &solution);
}

#[test]
fn slack_goes_to_the_first_generator_at_the_reference_bus() {
    // Both generators at bus 1 start at 10 MW; only the first absorbs the 60 MW mismatch
    let mut tables = triangle();
    tables.1 = vec![G(1, 10.0, true), G(1, 10.0, true)];
    let solution = run(&tables).unwrap();
    assert_close(&solution.pg_mw, &[70.0, 10.0]);
    // An out-of-service first generator doesn't count
    tables.1.insert(0, G(1, 500.0, false));
    let solution = run(&tables).unwrap();
    assert_close(&solution.pg_mw, &[0.0, 70.0, 10.0]);
    assert_physical(&case(&tables), &solution);
}

#[test]
fn out_of_service_elements_report_zero() {
    // A generator at bus 2 that is off, and branch 2–3 off
    let mut tables = triangle();
    tables.1.push(G(2, 40.0, false));
    tables.2[1].5 = false;
    let solution = run(&tables).unwrap();
    assert_eq!(solution.pg_mw[1], 0.0);
    assert_eq!(solution.pf_mw[1], 0.0);
    // Radial now: each load is served directly from bus 1
    assert_close(&solution.pf_mw, &[50.0, 0.0, 30.0]);
    assert_physical(&case(&tables), &solution);
}

#[test]
fn isolated_bus_keeps_its_angle_and_its_load_is_not_served() {
    // As in MATPOWER (checked with rundcpf on case14): the isolated bus's angle is unchanged,
    //   its generator and branch report 0, and generation covers only the connected load
    let mut tables = triangle();
    tables.0.push(B(99, 4, 3.0, 0.0, 7.5));
    tables.1.push(G(99, 20.0, true));
    tables.2.push(L(99, 2, 0.1, 0.0, 0.0, true));
    let solution = run(&tables).unwrap();
    let plain = run(&triangle()).unwrap();
    assert_eq!(solution.va_deg[3], 7.5);
    assert_eq!(solution.va_deg[..3], plain.va_deg[..]);
    assert_eq!(solution.pf_mw[3], 0.0);
    assert_close(&solution.pg_mw, &[80.0, 0.0]);
    assert_physical(&case(&tables), &solution);
}

#[test]
fn shunt_conductance_is_load() {
    // 20 MW of Gs at bus 3 has the same effect as 20 MW more Pd
    let mut with_gs = triangle();
    with_gs.0[2].3 = 20.0;
    let mut with_pd = triangle();
    with_pd.0[2].2 = 50.0;
    let (gs, pd) = (run(&with_gs).unwrap(), run(&with_pd).unwrap());
    assert_eq!(gs.va_deg, pd.va_deg);
    assert_eq!(gs.pf_mw, pd.pf_mw);
    assert_close(&gs.pg_mw, &[100.0]);
    assert_physical(&case(&with_gs), &gs);
}

#[test]
fn phase_shifter_on_a_radial_branch_shifts_the_angle_not_the_flow() {
    // Bus 2's 50 MW must flow through the one branch whatever the shift; the shift of φ
    //   lowers bus 2's angle by exactly φ
    let tables = (
        vec![B(1, 3, 0.0, 0.0, 0.0), B(2, 1, 50.0, 0.0, 0.0)],
        vec![G(1, 0.0, true)],
        vec![L(1, 2, 0.1, 0.0, 10.0, true)],
    );
    let solution = run(&tables).unwrap();
    assert_close(&solution.pf_mw, &[50.0]);
    assert_close(&solution.va_deg, &[0.0, degrees(-0.05) - 10.0]);
    assert_close(&solution.pg_mw, &[50.0]);
}

#[test]
fn phase_shifter_in_a_loop_redistributes_flow() {
    let mut tables = triangle();
    tables.2[0] = L(1, 2, 0.1, 0.0, -5.0, true);
    let solution = run(&tables).unwrap();
    // An advance on 1–2 pushes more power through it than without the shifter
    assert!(solution.pf_mw[0] > 6000.0 / 110.0 + 1.0);
    assert_physical(&case(&tables), &solution);
}

#[test]
fn tap_ratio_scales_the_flow_equation() {
    let mut tables = triangle();
    tables.2[2] = L(1, 3, 0.25, 0.8, 0.0, true);
    let solution = run(&tables).unwrap();
    assert_physical(&case(&tables), &solution);
}

#[test]
fn reference_angle_offsets_every_angle() {
    let mut tables = triangle();
    tables.0[0].4 = 5.0;
    let shifted = run(&tables).unwrap();
    let plain = run(&triangle()).unwrap();
    let expected: Vec<f64> = plain.va_deg.iter().map(|va| va + 5.0).collect();
    assert_close(&shifted.va_deg, &expected);
    assert_close(&shifted.pf_mw, &plain.pf_mw);
}

#[test]
fn several_reference_buses_each_take_their_own_slack() {
    // Two islands, each with its own reference bus and generator: 1–2 and 3–4
    let tables = (
        vec![
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 40.0, 0.0, 0.0),
            B(3, 3, 0.0, 0.0, 2.0),
            B(4, 1, 25.0, 0.0, 0.0),
        ],
        vec![G(1, 0.0, true), G(3, 0.0, true)],
        vec![L(1, 2, 0.1, 0.0, 0.0, true), L(3, 4, 0.2, 0.0, 0.0, true)],
    );
    let solution = run(&tables).unwrap();
    assert_close(&solution.pg_mw, &[40.0, 25.0]);
    assert_close(
        &solution.va_deg,
        &[0.0, degrees(-0.04), 2.0, 2.0 + degrees(-0.25 * 0.2)],
    );
    assert_physical(&case(&tables), &solution);
}

#[test]
fn negative_reactance_falls_back_to_lu() {
    // A series capacitor (x = −0.05, b = −20) on 2–3 makes B′ = [[−10, 20], [20, −16]],
    //   which is indefinite, so Cholesky fails and LU solves it
    let mut tables = triangle();
    tables.2[1] = L(2, 3, -0.05, 0.0, 0.0, true);
    let solution = run(&tables).unwrap();
    assert_eq!(solution.solver, DcSolver::Lu);
    // θ = B′⁻¹ P′ with P′ = (−0.5, −0.3): det B′ = 160 − 400 = −240
    let theta2 = (-16.0 * -0.5 - 20.0 * -0.3) / -240.0;
    let theta3 = (-20.0 * -0.5 + -10.0 * -0.3) / -240.0;
    assert_close(&solution.va_deg, &[0.0, degrees(theta2), degrees(theta3)]);
    assert_physical(&case(&tables), &solution);
}

#[test]
fn cancelling_reactances_are_a_singular_network() {
    // Bus 3 hangs off bus 2 by x = 0.2 and x = −0.2 in parallel: connected, but with zero
    //   net susceptance, so its angle is undetermined
    let tables = (
        vec![
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 10.0, 0.0, 0.0),
            B(3, 1, 10.0, 0.0, 0.0),
        ],
        vec![G(1, 0.0, true)],
        vec![
            L(1, 2, 0.1, 0.0, 0.0, true),
            L(2, 3, 0.2, 0.0, 0.0, true),
            L(2, 3, -0.2, 0.0, 0.0, true),
        ],
    );
    assert_eq!(
        run(&tables).unwrap_err(),
        PowerflowError::SingularNetwork { bus_label: 3 }
    );
}

#[test]
fn only_reference_buses_need_no_solve() {
    let tables = (vec![B(1, 3, 12.0, 3.0, 4.0)], vec![G(1, 0.0, true)], vec![]);
    let solution = run(&tables).unwrap();
    assert_eq!(solution.solver, DcSolver::NoSolve);
    assert_eq!(solution.va_deg, [4.0]);
    assert_close(&solution.pg_mw, &[15.0]);
}

#[test]
fn a_model_solves_new_injections() {
    let tables = triangle();
    let base = case(&tables);
    let model = DcModel::new(&base).unwrap();
    let mut changed = base.clone();
    changed.bus[1].pd_mw = 70.0;
    changed.generator[0].pg_mw = 0.0;
    assert_eq!(
        model.solve(&changed).unwrap(),
        dc_power_flow(&changed).unwrap()
    );
    assert_eq!(model.solve(&base).unwrap(), dc_power_flow(&base).unwrap());
}

#[test]
fn a_model_rejects_a_different_case() {
    let model = DcModel::new(&case(&triangle())).unwrap();
    let mut other = triangle();
    other.1.push(G(2, 1.0, true));
    assert_eq!(
        model.solve(&case(&other)).unwrap_err(),
        PowerflowError::CaseMismatch {
            expected: (3, 1, 3),
            found: (3, 2, 3)
        }
    );
}
