// Tests for the AC power flow's Newton solve and its results: convergence control, the
//   starting point, dropped elements, reactive power sharing between generators, and the
//   physics of the solution (power balance at every bus).

mod common;
use common::ac_fixture::{ac_fixture, AcBranch as L, AcBus as B, AcGen as G};
use common::reference::{ac_deviation, committed_dir, read_ac};
use powerflow::ac::{ac_power_flow, AcModel, AcOptions, AcSolution};
use powerflow::case::{BusType, Case};
use powerflow::powerflow_error::PowerflowError;
use sparse_mat::Complex;

use std::collections::HashMap;
use std::f64::consts::PI;

/// A four-bus network: a reference bus, a PV bus, two PQ buses (one with a shunt), lines with
/// resistance and charging, and a tapped transformer.
fn network() -> Case {
    ac_fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            B(2, 2, 20.0, 5.0, 0.0, 0.0, 1.0, 0.0),
            B(3, 1, 45.0, 15.0, 0.0, 10.0, 1.0, 0.0),
            B(4, 1, 40.0, 12.0, 0.0, 0.0, 1.0, 0.0),
        ],
        &[
            G(1, 0.0, 0.0, 100.0, -100.0, 1.04, true),
            G(2, 40.0, 0.0, 60.0, -40.0, 1.02, true),
        ],
        &[
            L(1, 2, 0.02, 0.06, 0.03, 0.0, 0.0, true),
            L(2, 3, 0.03, 0.09, 0.02, 0.0, 0.0, true),
            L(1, 3, 0.01, 0.05, 0.02, 0.0, 0.0, true),
            L(3, 4, 0.005, 0.08, 0.0, 0.97, 0.0, true),
        ],
    )
}

/// Checks the physics of `solution`, independently of how it was computed: each in-service
/// branch's flows are its π-model's at the solved voltages, every connected bus balances
/// generation, load, shunt and flows, and dropped elements report zero.
fn assert_physical(case: &Case, solution: &AcSolution) {
    let index: HashMap<usize, usize> = case.bus_indices();
    let voltage = |label: usize| {
        let row = index[&label];
        Complex::from_polar(solution.vm_pu[row], solution.va_deg[row] * PI / 180.0)
    };
    let connected = |label: usize| case.bus[index[&label]].bus_type != BusType::Isolated;
    let base = case.base_mva;
    let mut balance = vec![Complex::new(0.0, 0.0); case.bus.len()];
    for (row, branch) in case.branch.iter().enumerate() {
        let live = branch.status && connected(branch.from_label) && connected(branch.to_label);
        let flows = [
            solution.pf_mw[row],
            solution.qf_mvar[row],
            solution.pt_mw[row],
            solution.qt_mvar[row],
        ];
        if !live {
            assert_eq!(flows, [0.0; 4], "branch row {row}");
            continue;
        }
        let ys = Complex::new(1.0, 0.0) / Complex::new(branch.r, branch.x);
        let tap = Complex::from_polar(branch.effective_tap(), branch.shift_deg * PI / 180.0);
        let (vf, vt) = (voltage(branch.from_label), voltage(branch.to_label));
        let half_charging = Complex::new(0.0, branch.b / 2.0);
        let i_f = (ys + half_charging) * vf / (tap * tap.conj()) - ys * vt / tap.conj();
        let i_t = -ys * vf / tap + (ys + half_charging) * vt;
        let (sf, st) = (vf * i_f.conj() * base, vt * i_t.conj() * base);
        for (got, want) in flows.iter().zip([sf.re, sf.im, st.re, st.im]) {
            assert!(
                (got - want).abs() <= 1e-9 * want.abs().max(1.0),
                "row {row}: {got} vs {want}"
            );
        }
        balance[index[&branch.from_label]] -= sf;
        balance[index[&branch.to_label]] -= st;
    }
    for (row, generator) in case.generator.iter().enumerate() {
        if !(generator.status && connected(generator.bus_label)) {
            assert_eq!((solution.pg_mw[row], solution.qg_mvar[row]), (0.0, 0.0));
            continue;
        }
        balance[index[&generator.bus_label]] +=
            Complex::new(solution.pg_mw[row], solution.qg_mvar[row]);
    }
    for (row, bus) in case.bus.iter().enumerate() {
        if bus.bus_type == BusType::Isolated {
            assert_eq!(
                (solution.vm_pu[row], solution.va_deg[row]),
                (bus.vm_pu, bus.va_deg)
            );
            continue;
        }
        let vm2 = solution.vm_pu[row] * solution.vm_pu[row];
        let net = balance[row]
            - Complex::new(bus.pd_mw, bus.qd_mvar)
            - Complex::new(bus.gs_mw, -bus.bs_mvar) * vm2;
        assert!(
            net.norm() < 1e-5,
            "bus {}: imbalance {net} MVA",
            bus.bus_label
        );
    }
}

#[test]
fn solutions_balance_at_every_bus() {
    let case = network();
    let solution = ac_power_flow(&case).unwrap();
    assert_physical(&case, &solution);
    // Voltage-controlled buses hold their setpoints, and the reference angle stays put
    assert!((solution.vm_pu[0] - 1.04).abs() < 1e-12);
    assert!((solution.vm_pu[1] - 1.02).abs() < 1e-12);
    assert_eq!(solution.va_deg[0], 0.0);
    // The PV generator keeps its scheduled real output
    assert_eq!(solution.pg_mw[1], 40.0);
    assert!(solution.iterations > 0 && solution.mismatch_pu < 1e-8);
}

#[test]
fn dropped_elements_report_zero_and_isolated_buses_keep_their_voltage() {
    let mut case = network();
    case.bus.push(case.bus[3].clone());
    let extra = case.bus.len() - 1;
    case.bus[extra].bus_label = 5;
    case.bus[extra].bus_type = BusType::Isolated;
    case.bus[extra].vm_pu = 0.95;
    case.bus[extra].va_deg = -7.0;
    let mut tie = case.branch[3].clone();
    tie.to_label = 5;
    case.branch.push(tie);
    case.branch[1].status = false;
    let mut idle = case.generator[1].clone();
    idle.status = false;
    case.generator.push(idle);
    let solution = ac_power_flow(&case).unwrap();
    assert_physical(&case, &solution);
}

#[test]
fn convergence_control() {
    let case = network();
    let model = AcModel::new(&case).unwrap();
    // One iteration isn't enough from this start; the error reports the mismatch left
    // AcOptions is non-exhaustive, so a caller changes a default rather than writing a literal
    let mut one = AcOptions::default();
    one.max_iterations = 1;
    match model.solve(&case, one) {
        Err(PowerflowError::NoConvergence {
            iterations,
            mismatch_pu,
        }) => {
            assert_eq!(iterations, 1);
            assert!(mismatch_pu > 1e-8 && mismatch_pu.is_finite());
        }
        other => panic!("expected NoConvergence, got {other:?}"),
    }
    // A looser tolerance stops sooner
    let mut loose = AcOptions::default();
    loose.tolerance = 1e-3;
    let full = model.solve(&case, AcOptions::default()).unwrap();
    let rough = model.solve(&case, loose).unwrap();
    assert!(rough.iterations < full.iterations);
    assert!(rough.mismatch_pu < 1e-3);
    assert_eq!(AcOptions::default().tolerance, 1e-8);
    assert_eq!(AcOptions::default().max_iterations, 10);

    // Starting from the solution needs no iteration, and no factorization
    let mut solved = case.clone();
    for (bus, (vm, va)) in solved
        .bus
        .iter_mut()
        .zip(full.vm_pu.iter().zip(&full.va_deg))
    {
        bus.vm_pu = *vm;
        bus.va_deg = *va;
    }
    let again = model.solve(&solved, AcOptions::default()).unwrap();
    assert_eq!((again.iterations, again.full_factorizations), (0, 0));

    // Another case's tables are refused
    let mut shorter = case.clone();
    shorter.branch.pop();
    assert!(matches!(
        model.solve(&shorter, AcOptions::default()),
        Err(PowerflowError::CaseMismatch { .. })
    ));
}

#[test]
fn a_flat_start_reaches_the_same_solution() {
    let Some(path) = common::matpower_case("case118") else {
        return; // not fetched; the skip was reported on stderr
    };
    // Flat everywhere but the reference bus, whose angle (30° in case118) fixes all the others
    let mut case = Case::from_path(&path).unwrap();
    for bus in case
        .bus
        .iter_mut()
        .filter(|bus| bus.bus_type != BusType::Ref)
    {
        bus.vm_pu = 1.0;
        bus.va_deg = 0.0;
    }
    let solution = ac_power_flow(&case).unwrap();
    let reference = read_ac(&committed_dir(), "case118").unwrap();
    let deviation = ac_deviation(&case, &solution, &reference).unwrap();
    // Converged from elsewhere, so it agrees with MATPOWER only to the tolerance: mismatches
    //   under 1e-8 p.u. are under 1e-6 MW or MVAr on case118's 100 MVA base (measured: 1.3e-8)
    assert!(deviation.largest() < 1e-6, "{deviation:?}");
    assert!(solution.mismatch_pu < 1e-8);
}

#[test]
fn generators_at_one_bus_share_reactive_power_by_their_ranges() {
    // Two generators at the PV bus 2, with ranges [-10, 30] and [-30, 90] MVAr: pfsoln gives
    //   each its minimum plus the same fraction of its range, so their outputs stay in the
    //   ratio of their ranges' positions
    let mut case = network();
    case.generator[1].pg_mw = 20.0;
    case.generator[1].qmax_mvar = 30.0;
    case.generator[1].qmin_mvar = -10.0;
    let mut second = case.generator[1].clone();
    second.qmax_mvar = 90.0;
    second.qmin_mvar = -30.0;
    case.generator.push(second);
    let solution = ac_power_flow(&case).unwrap();
    assert_physical(&case, &solution);
    let (q1, q2) = (solution.qg_mvar[1], solution.qg_mvar[2]);
    let total = q1 + q2;
    let fraction = (total - (-40.0)) / (120.0 - (-40.0) + f64::EPSILON);
    assert!((q1 - (-10.0 + fraction * 40.0)).abs() < 1e-9, "{q1}");
    assert!((q2 - (-30.0 + fraction * 120.0)).abs() < 1e-9, "{q2}");

    // With an infinite limit, the bound is each generator's share plus its finite limits'
    //   magnitudes, summed over the bus
    case.generator[2].qmax_mvar = f64::INFINITY;
    let solution = ac_power_flow(&case).unwrap();
    assert_physical(&case, &solution);
    let total = solution.qg_mvar[1] + solution.qg_mvar[2];
    let bound = (total / 2.0).abs() + 30.0 + 10.0 + (total / 2.0).abs() + 30.0;
    let (min, max) = (-10.0 - 30.0, 30.0 + bound);
    let fraction = (total - min) / (max - min + f64::EPSILON);
    assert!((solution.qg_mvar[2] - (-30.0 + fraction * (bound + 30.0))).abs() < 1e-9);
}

// Fast-decoupled

use powerflow::ac::AcMethod;

/// `network()` with a phase-shifting transformer, which makes B′ unsymmetric.
fn shifted_network() -> Case {
    let mut case = network();
    case.branch[3].shift_deg = -4.0;
    case
}

#[test]
fn fast_decoupled_reaches_newtons_solution() {
    for case in [network(), shifted_network()] {
        let model = AcModel::new(&case).unwrap();
        let newton = model.solve(&case, AcOptions::default()).unwrap();
        for method in [AcMethod::FastDecoupledXb, AcMethod::FastDecoupledBx] {
            let solution = model.solve(&case, AcOptions::for_method(method)).unwrap();
            assert_physical(&case, &solution);
            assert_eq!(solution.full_factorizations, 2);
            // Each converged to within the tolerance (1e-8 p.u. on a 100 MVA base), from the
            //   same start but by a different path
            for (fd, nr) in solution.vm_pu.iter().zip(&newton.vm_pu) {
                assert!((fd - nr).abs() < 1e-6, "{method:?}: {fd} vs {nr}");
            }
            for (fd, nr) in solution.pf_mw.iter().zip(&newton.pf_mw) {
                assert!((fd - nr).abs() < 1e-4, "{method:?}: {fd} vs {nr}");
            }
        }
    }
}

#[test]
fn fast_decoupled_options_and_errors() {
    assert_eq!(AcOptions::default().method, AcMethod::Newton);
    let xb = AcOptions::for_method(AcMethod::FastDecoupledXb);
    assert_eq!((xb.tolerance, xb.max_iterations), (1e-8, 30));
    assert_eq!(
        AcOptions::for_method(AcMethod::Newton),
        AcOptions::default()
    );

    let case = network();
    let model = AcModel::new(&case).unwrap();
    let mut one = AcOptions::for_method(AcMethod::FastDecoupledBx);
    one.max_iterations = 1;
    match model.solve(&case, one) {
        Err(PowerflowError::NoConvergence {
            iterations,
            mismatch_pu,
        }) => assert!(iterations == 1 && mismatch_pu > 1e-8),
        other => panic!("expected NoConvergence, got {other:?}"),
    }

    // A branch with resistance but no reactance is fine for Newton's method, but its B′ or B″
    //   stamp would be infinite once resistance is ignored
    let mut resistive = case.clone();
    resistive.branch[0].x = 0.0;
    let model = AcModel::new(&resistive).unwrap();
    assert!(model.solve(&resistive, AcOptions::default()).is_ok());
    for method in [AcMethod::FastDecoupledXb, AcMethod::FastDecoupledBx] {
        assert_eq!(
            model.solve(&resistive, AcOptions::for_method(method)).err(),
            Some(PowerflowError::ZeroReactance { row: 1 })
        );
    }
}

// Reactive power limits

use powerflow::ac::ReactiveLimits;

fn with_limits(mode: ReactiveLimits) -> AcOptions {
    let mut options = AcOptions::default();
    options.reactive_limits = mode;
    options
}

/// `network()` with a second PV bus (bus 4, which had only load), so two generators can reach
/// their limits.
fn two_pv_network() -> Case {
    let mut case = network();
    case.bus[3].bus_type = BusType::PV;
    let mut second = case.generator[1].clone();
    second.bus_label = 4;
    second.pg_mw = 25.0;
    second.vg_pu = 1.03;
    case.generator.push(second);
    case
}

/// `network()` with the reference generator's limits wide (±300 MVAr), so only the PV
/// generator at bus 2 reaches its limits, and that generator's setpoint `vg`.
fn limits_network(vg: f64) -> Case {
    let mut case = network();
    case.generator[0].qmax_mvar = 300.0;
    case.generator[0].qmin_mvar = -300.0;
    case.generator[1].vg_pu = vg;
    case
}

#[test]
fn reactive_limits_are_ignored_by_default() {
    let mut case = limits_network(1.06);
    let unlimited = ac_power_flow(&case).unwrap();
    case.generator[1].qmax_mvar = unlimited.qg_mvar[1] - 5.0;
    let solution = ac_power_flow(&case).unwrap();
    assert_eq!(
        solution.qg_mvar[1], unlimited.qg_mvar[1],
        "not enforced by default"
    );
    assert!(solution.limited_generators.is_empty());
    assert_eq!(AcOptions::default().reactive_limits, ReactiveLimits::Ignore);
}

#[test]
fn a_generator_past_its_upper_limit_is_held_there() {
    let mut case = limits_network(1.06);
    let unlimited = ac_power_flow(&case).unwrap();
    // The fixture's PV generator must produce enough that 5 MVAr less is still above Qmin
    assert!(
        unlimited.qg_mvar[1] - 5.0 > case.generator[1].qmin_mvar,
        "{}",
        unlimited.qg_mvar[1]
    );
    case.generator[1].qmax_mvar = unlimited.qg_mvar[1] - 5.0;
    let model = AcModel::new(&case).unwrap();
    for mode in [ReactiveLimits::AllAtOnce, ReactiveLimits::OneAtATime] {
        let solution = model.solve(&case, with_limits(mode)).unwrap();
        assert_physical(&case, &solution);
        assert_eq!(solution.limited_generators, [1]);
        assert_eq!(solution.qg_mvar[1], case.generator[1].qmax_mvar);
        // The bus no longer holds its setpoint: less reactive power, lower voltage
        assert!(solution.vm_pu[1] < 1.06 - 1e-6, "{}", solution.vm_pu[1]);
        assert!(solution.iterations > unlimited.iterations, "two solves");
    }
}

#[test]
fn a_generator_past_its_lower_limit_is_held_there() {
    let mut case = limits_network(0.97);
    let unlimited = ac_power_flow(&case).unwrap();
    assert!(
        unlimited.qg_mvar[1] + 3.0 < case.generator[1].qmax_mvar,
        "{}",
        unlimited.qg_mvar[1]
    );
    case.generator[1].qmin_mvar = unlimited.qg_mvar[1] + 3.0;
    let solution = AcModel::new(&case)
        .unwrap()
        .solve(&case, with_limits(ReactiveLimits::AllAtOnce))
        .unwrap();
    assert_physical(&case, &solution);
    assert_eq!(solution.limited_generators, [1]);
    assert_eq!(solution.qg_mvar[1], case.generator[1].qmin_mvar);
    assert!(solution.vm_pu[1] > 0.97 + 1e-6, "{}", solution.vm_pu[1]);
}

#[test]
fn generators_sharing_a_bus_are_held_at_their_limits_together() {
    // Two generators at bus 2. pfsoln shares the bus's output in proportion to their ranges,
    //   so lowering one generator's limit lowers its share too; to force violations, the bus's
    //   output must exceed both upper limits together
    let mut case = limits_network(1.06);
    let mut second = case.generator[1].clone();
    second.pg_mw = 0.0;
    case.generator.push(second);
    let unlimited = ac_power_flow(&case).unwrap();
    let total = unlimited.qg_mvar[1] + unlimited.qg_mvar[2];
    for row in [1, 2] {
        case.generator[row].qmax_mvar = total / 2.0 - 3.0;
        assert!(case.generator[row].qmax_mvar > case.generator[row].qmin_mvar);
    }
    let solution = AcModel::new(&case)
        .unwrap()
        .solve(&case, with_limits(ReactiveLimits::AllAtOnce))
        .unwrap();
    assert_physical(&case, &solution);
    assert_eq!(solution.limited_generators, [1, 2]);
    for row in [1, 2] {
        assert_eq!(solution.qg_mvar[row], case.generator[row].qmax_mvar);
    }
}

#[test]
fn inverted_limits_are_refused_rather_than_looping() {
    // With Qmin above Qmax, MATPOWER's loop alternates between the two limits forever
    let mut case = limits_network(1.06);
    case.generator[1].qmin_mvar = 10.0;
    case.generator[1].qmax_mvar = 5.0;
    let model = AcModel::new(&case).unwrap();
    assert_eq!(
        model
            .solve(&case, with_limits(ReactiveLimits::AllAtOnce))
            .err(),
        Some(PowerflowError::InvertedReactiveLimits { row: 2 })
    );
    // Ignoring limits, they don't matter
    assert!(model.solve(&case, AcOptions::default()).is_ok());
}

#[test]
fn a_reference_generator_at_its_limit_hands_the_reference_to_a_pv_bus() {
    let mut case = network();
    case.bus[0].va_deg = 10.0;
    let unlimited = ac_power_flow(&case).unwrap();
    case.generator[0].qmax_mvar = unlimited.qg_mvar[0] - 5.0;
    let solution = AcModel::new(&case)
        .unwrap()
        .solve(&case, with_limits(ReactiveLimits::AllAtOnce))
        .unwrap();
    assert_physical(&case, &solution);
    assert_eq!(solution.limited_generators, [0]);
    // Bus 2 became the reference and took up the slack, while the angles were shifted back so
    //   the original reference keeps its case angle
    assert!(
        (solution.va_deg[0] - 10.0).abs() < 1e-9,
        "{}",
        solution.va_deg[0]
    );
    assert!(
        (solution.vm_pu[1] - 1.02).abs() < 1e-12,
        "bus 2 holds its voltage"
    );
    assert_ne!(
        solution.pg_mw[1], 40.0,
        "bus 2's generator is the slack now"
    );
}

#[test]
fn several_references_cannot_lose_one_to_a_limit() {
    let mut case = two_pv_network();
    case.bus[3].bus_type = BusType::Ref;
    let unlimited = ac_power_flow(&case).unwrap();
    case.generator[2].qmax_mvar = unlimited.qg_mvar[2] - 5.0;
    assert_eq!(
        AcModel::new(&case)
            .unwrap()
            .solve(&case, with_limits(ReactiveLimits::AllAtOnce))
            .err(),
        Some(PowerflowError::ReferenceAtReactiveLimit { bus_label: 4 })
    );
}

#[test]
fn every_generator_past_the_same_limit_is_infeasible() {
    let mut case = limits_network(1.06);
    let unlimited = ac_power_flow(&case).unwrap();
    assert!(unlimited.qg_mvar[1] - 5.0 > case.generator[1].qmin_mvar);
    case.generator[0].qmax_mvar = unlimited.qg_mvar[0] - 5.0;
    case.generator[1].qmax_mvar = unlimited.qg_mvar[1] - 5.0;
    assert_eq!(
        AcModel::new(&case)
            .unwrap()
            .solve(&case, with_limits(ReactiveLimits::AllAtOnce))
            .err(),
        Some(PowerflowError::ReactiveLimitsInfeasible { generators: 2 })
    );
}

#[test]
fn one_at_a_time_converts_the_largest_violation_first() {
    let mut case = two_pv_network();
    let unlimited = ac_power_flow(&case).unwrap();
    // Both PV generators exceed their limits, generator 2 (at bus 4) by more
    case.generator[1].qmax_mvar = unlimited.qg_mvar[1] - 1.0;
    case.generator[2].qmax_mvar = unlimited.qg_mvar[2] - 6.0;
    let model = AcModel::new(&case).unwrap();
    let all = model
        .solve(&case, with_limits(ReactiveLimits::AllAtOnce))
        .unwrap();
    let one = model
        .solve(&case, with_limits(ReactiveLimits::OneAtATime))
        .unwrap();
    assert_physical(&case, &all);
    assert_physical(&case, &one);
    assert_eq!(one.limited_generators[0], 2, "the largest violation first");
    assert!(one.iterations >= all.iterations);
    for solution in [&all, &one] {
        for &row in &solution.limited_generators {
            assert_eq!(solution.qg_mvar[row], case.generator[row].qmax_mvar);
        }
    }
}
