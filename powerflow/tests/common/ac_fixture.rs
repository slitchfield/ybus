//! Small hand-built cases for the AC power flow tests. Each table lists the columns the AC
//! power flow reads; the rest are filled with neutral values, and the result is parsed from
//! MATPOWER text, so fixtures go through the same reader as real cases.

use powerflow::case::Case;

/// A bus: label, MATPOWER type code (1 PQ, 2 PV, 3 REF, 4 isolated), Pd (MW), Qd (MVAr), Gs (MW
/// at 1 p.u.), Bs (MVAr at 1 p.u.), Vm (p.u.) and Va (degrees).
pub struct AcBus(
    pub usize,
    pub u8,
    pub f64,
    pub f64,
    pub f64,
    pub f64,
    pub f64,
    pub f64,
);

/// A generator: bus label, Pg (MW), Qg (MVAr), Qmax and Qmin (MVAr), Vg (p.u.), in service.
pub struct AcGen(
    pub usize,
    pub f64,
    pub f64,
    pub f64,
    pub f64,
    pub f64,
    pub bool,
);

/// A branch: from label, to label, r, x and b (p.u.), tap ratio (0 for a line), shift (degrees),
/// in service.
pub struct AcBranch(
    pub usize,
    pub usize,
    pub f64,
    pub f64,
    pub f64,
    pub f64,
    pub f64,
    pub bool,
);

/// A case on a 100 MVA base from the given tables. Panics if it doesn't parse, since a fixture
/// that doesn't is a bug in the test.
pub fn ac_fixture(buses: &[AcBus], gens: &[AcGen], branches: &[AcBranch]) -> Case {
    let mut text = String::from("function mpc = fixture\nmpc.version = '2';\nmpc.baseMVA = 100;\n");
    text.push_str("mpc.bus = [\n");
    for AcBus(label, kind, pd, qd, gs, bs, vm, va) in buses {
        text.push_str(&format!(
            "{label} {kind} {pd:?} {qd:?} {gs:?} {bs:?} 1 {vm:?} {va:?} 135 1 1.1 0.9;\n"
        ));
    }
    text.push_str("];\nmpc.gen = [\n");
    for AcGen(bus, pg, qg, qmax, qmin, vg, status) in gens {
        let status = u8::from(*status);
        text.push_str(&format!(
            "{bus} {pg:?} {qg:?} {qmax:?} {qmin:?} {vg:?} 100 {status} 1000 0 0 0 0 0 0 0 0 0 0 0 0;\n"
        ));
    }
    text.push_str("];\nmpc.branch = [\n");
    for AcBranch(from, to, r, x, b, tap, shift, status) in branches {
        let status = u8::from(*status);
        text.push_str(&format!(
            "{from} {to} {r:?} {x:?} {b:?} 0 0 0 {tap:?} {shift:?} {status} -360 360;\n"
        ));
    }
    text.push_str("];\n");
    Case::parse(&text).unwrap_or_else(|err| panic!("fixture doesn't parse: {err}\n{text}"))
}
