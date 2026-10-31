//! Small hand-built cases for the DC power flow tests, and a check of a solution's physics.
//!
//! Each fixture table lists only the columns DC power flow reads; the rest are filled with
//! neutral values, and the result is parsed from MATPOWER text, so fixtures go through the same
//! reader as real cases.

use std::collections::HashMap;
use std::f64::consts::PI;

use powerflow::case::{BusType, Case};
use powerflow::dc::DcSolution;

/// A bus: label, MATPOWER type code (1 PQ, 2 PV, 3 REF, 4 isolated), Pd (MW), Gs (MW at
/// 1 p.u.) and Va (degrees).
pub struct B(pub usize, pub u8, pub f64, pub f64, pub f64);

/// A generator: bus label, Pg (MW), in service.
pub struct G(pub usize, pub f64, pub bool);

/// A branch: from label, to label, x (p.u.), tap ratio (0 for a line), shift (degrees), in
/// service.
pub struct L(pub usize, pub usize, pub f64, pub f64, pub f64, pub bool);

/// A case on a 100 MVA base from the given tables. Panics if it doesn't parse, since a fixture
/// that doesn't is a bug in the test.
pub fn fixture(buses: &[B], gens: &[G], branches: &[L]) -> Case {
    let mut text = String::from("function mpc = fixture\nmpc.version = '2';\nmpc.baseMVA = 100;\n");
    text.push_str("mpc.bus = [\n");
    for B(label, kind, pd, gs, va) in buses {
        text.push_str(&format!(
            "{label} {kind} {pd:?} 0 {gs:?} 0 1 1 {va:?} 135 1 1.1 0.9;\n"
        ));
    }
    text.push_str("];\nmpc.gen = [\n");
    for G(bus, pg, status) in gens {
        let status = u8::from(*status);
        text.push_str(&format!(
            "{bus} {pg:?} 0 100 -100 1 100 {status} 1000 0 0 0 0 0 0 0 0 0 0 0 0;\n"
        ));
    }
    text.push_str("];\nmpc.branch = [\n");
    for L(from, to, x, tap, shift, status) in branches {
        let status = u8::from(*status);
        text.push_str(&format!(
            "{from} {to} 0 {x:?} 0 0 0 0 {tap:?} {shift:?} {status} -360 360;\n"
        ));
    }
    text.push_str("];\n");
    Case::parse(&text).unwrap_or_else(|err| panic!("fixture doesn't parse: {err}\n{text}"))
}

/// Checks the physics of a solution, independently of how it was computed: each in-service
/// branch's flow is b · (θf − θt − shift), and real power balances at every connected bus.
pub fn assert_physical(case: &Case, solution: &DcSolution) {
    let index: HashMap<usize, usize> = case.bus_indices();
    let va = |label: usize| solution.va_deg[index[&label]] * PI / 180.0;
    let connected = |label: usize| case.bus[index[&label]].bus_type != BusType::Isolated;
    let mut net_mw = vec![0.0; case.bus.len()];
    for (row, branch) in case.branch.iter().enumerate() {
        let live = branch.status && connected(branch.from_label) && connected(branch.to_label);
        if !live {
            assert_eq!(solution.pf_mw[row], 0.0);
            continue;
        }
        let b = 1.0 / (branch.x * branch.effective_tap());
        let expected = b
            * (va(branch.from_label) - va(branch.to_label) - branch.shift_deg * PI / 180.0)
            * case.base_mva;
        assert!(
            (solution.pf_mw[row] - expected).abs() <= 1e-9 * expected.abs().max(1.0),
            "branch row {row}: {} vs {expected}",
            solution.pf_mw[row]
        );
        net_mw[index[&branch.from_label]] -= solution.pf_mw[row];
        net_mw[index[&branch.to_label]] += solution.pf_mw[row];
    }
    for (row, generator) in case.generator.iter().enumerate() {
        if !(generator.status && connected(generator.bus_label)) {
            assert_eq!(solution.pg_mw[row], 0.0);
            continue;
        }
        net_mw[index[&generator.bus_label]] += solution.pg_mw[row];
    }
    for (row, bus) in case.bus.iter().enumerate() {
        if bus.bus_type == BusType::Isolated {
            continue;
        }
        let mismatch = net_mw[row] - bus.pd_mw - bus.gs_mw;
        assert!(
            mismatch.abs() <= 1e-9,
            "bus {}: mismatch {mismatch} MW",
            bus.bus_label
        );
    }
}
