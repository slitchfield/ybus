// A branch's magnetizing admittance (Branch::g_mag and b_mag): its stamps add to Yff alone,
//   undivided by the tap; in service, the AC power flow solves as if it were a shunt at the
//   from bus, by every method; and a branch outage takes it out with the branch, where a shunt
//   folded into the bus would stay behind.
mod common;

use common::ac_fixture::{ac_fixture, AcBranch as L, AcBus as B, AcGen as G};
use ybus_powerflow::ac::{ac_power_flow_with_options, AcMethod, AcModel, AcOptions};
use ybus_powerflow::case::{Case, GenCost};
use ybus_powerflow::opf::ac::{ac_opf_with_options, AcOpfOptions};
use ybus_powerflow::opf::batch::{ac_opf_outages, BatchStart, OpfOutage};
use ybus_sparse::Complex;

use std::f64::consts::PI;

/// The magnetizing admittance used throughout (p.u.): a small loss and an inductive draw
const G_MAG: f64 = 0.004;
const B_MAG: f64 = -0.03;

/// The transformer's row in `network`
const TRANSFORMER: usize = 3;

fn close(a: Complex<f64>, b: Complex<f64>) -> bool {
    (a - b).norm() <= 1e-12 * b.norm().max(1.0)
}

/// Three buses in a triangle of lines, with a tapped, phase-shifting transformer from bus 1 to
/// bus 3 in parallel with the line between them, which carries the magnetizing admittance.
/// Generators at buses 1 and 2 have quadratic costs.
fn network() -> Case {
    let mut case = ac_fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            B(2, 2, 20.0, 10.0, 0.0, 0.0, 1.0, 0.0),
            B(3, 1, 90.0, 30.0, 0.0, 5.0, 1.0, 0.0),
        ],
        &[
            G(1, 0.0, 0.0, 300.0, -300.0, 1.03, true),
            G(2, 60.0, 0.0, 300.0, -300.0, 1.01, true),
        ],
        &[
            L(1, 2, 0.01, 0.08, 0.04, 0.0, 0.0, true),
            L(2, 3, 0.02, 0.1, 0.03, 0.0, 0.0, true),
            L(1, 3, 0.015, 0.09, 0.02, 0.0, 0.0, true),
            L(1, 3, 0.005, 0.12, 0.0, 0.98, 2.0, true),
        ],
    );
    case.branch[TRANSFORMER].g_mag = G_MAG;
    case.branch[TRANSFORMER].b_mag = B_MAG;
    let quadratic = |a, b| GenCost::Polynomial {
        startup: 0.0,
        shutdown: 0.0,
        coeffs: vec![a, b, 0.0],
    };
    case.gencost = Some(vec![quadratic(0.02, 20.0), quadratic(0.03, 25.0)]);
    case
}

/// `network` with the magnetizing admittance moved into bus 1's shunt instead, as MATPOWER's
/// `psse2mpc` folds it
fn folded() -> Case {
    let mut case = network();
    case.branch[TRANSFORMER].g_mag = 0.0;
    case.branch[TRANSFORMER].b_mag = 0.0;
    case.bus[0].gs_mw += G_MAG * case.base_mva;
    case.bus[0].bs_mvar += B_MAG * case.base_mva;
    case
}

#[test]
fn magnetizing_admittance_adds_to_yff_undivided_by_the_tap() {
    let case = network();
    let mut without = network();
    without.branch[TRANSFORMER].g_mag = 0.0;
    without.branch[TRANSFORMER].b_mag = 0.0;
    let (model, plain) = (
        AcModel::new(&case).unwrap(),
        AcModel::new(&without).unwrap(),
    );
    let mag = Complex::new(G_MAG, B_MAG);

    // The transformer's own stamps by hand: the π-model, then the magnetizing admittance
    let ys = Complex::new(1.0, 0.0) / Complex::new(0.005, 0.12);
    let tau = Complex::from_polar(0.98, 2.0 * PI / 180.0);
    let yff = ys / (0.98 * 0.98) + mag;
    let (yf, yt) = (model.yf().to_dense(), model.yt().to_dense());
    assert!(
        close(yf[TRANSFORMER][0], yff),
        "{} vs {yff}",
        yf[TRANSFORMER][0]
    );
    assert!(close(yf[TRANSFORMER][2], -ys / tau.conj()));
    assert!(close(yt[TRANSFORMER][0], -ys / tau));
    assert!(close(yt[TRANSFORMER][2], ys));

    // Ybus differs from the plain model's only at (1, 1), by the magnetizing admittance
    let (ybus, plain_ybus) = (model.ybus().to_dense(), plain.ybus().to_dense());
    for i in 0..3 {
        for j in 0..3 {
            let want = plain_ybus[i][j] + if (i, j) == (0, 0) { mag } else { 0.0.into() };
            assert!(
                close(ybus[i][j], want),
                "({i}, {j}): {} vs {want}",
                ybus[i][j]
            );
        }
    }
}

#[test]
fn in_service_it_solves_as_a_from_bus_shunt_by_every_method() {
    let (case, folded) = (network(), folded());
    for method in [
        AcMethod::Newton,
        AcMethod::FastDecoupledXb,
        AcMethod::FastDecoupledBx,
    ] {
        let mut options = AcOptions::default();
        options.method = method;
        options.tolerance = 1e-12;
        let got = ac_power_flow_with_options(&case, &options).unwrap();
        let want = ac_power_flow_with_options(&folded, &options).unwrap();
        for (a, b) in got.vm_pu.iter().zip(&want.vm_pu) {
            assert!((a - b).abs() <= 1e-10, "{method:?}: Vm {a} vs {b}");
        }
        for (a, b) in got.va_deg.iter().zip(&want.va_deg) {
            assert!((a - b).abs() <= 1e-8, "{method:?}: Va {a} vs {b}");
        }

        // The transformer's from-end flow includes what the magnetizing admittance draws,
        //   S = |Vf|² · conj(Gm + jBm), which the folded case has the bus shunt draw instead
        let vf2 = got.vm_pu[0] * got.vm_pu[0];
        let base = case.base_mva;
        let (pf, qf) = (got.pf_mw[TRANSFORMER], got.qf_mvar[TRANSFORMER]);
        let (pf_folded, qf_folded) = (want.pf_mw[TRANSFORMER], want.qf_mvar[TRANSFORMER]);
        assert!(
            (pf - (pf_folded + G_MAG * vf2 * base)).abs() <= 1e-7,
            "{method:?}: Pf"
        );
        assert!(
            (qf - (qf_folded - B_MAG * vf2 * base)).abs() <= 1e-7,
            "{method:?}: Qf"
        );
    }
}

#[test]
fn an_outage_takes_it_out_with_the_branch() {
    let options = AcOpfOptions::default();
    let outages = ac_opf_outages(
        &network(),
        &[OpfOutage::branch(TRANSFORMER)],
        &options,
        BatchStart::Cold,
    )
    .unwrap();
    let outage = outages.outages[0].as_ref().unwrap();

    // The same as an OPF of the case without the transformer at all
    let mut removed = network();
    removed.branch.remove(TRANSFORMER);
    let want = ac_opf_with_options(&removed, &options).unwrap();
    let relative = (outage.objective - want.objective).abs() / want.objective.abs();
    assert!(
        relative <= 1e-9,
        "{} vs {}",
        outage.objective,
        want.objective
    );
    for (a, b) in outage.vm_pu.iter().zip(&want.vm_pu) {
        assert!((a - b).abs() <= 1e-8, "Vm {a} vs {b}");
    }

    // Folded into the bus, the admittance would outlive the outage and change the result
    let mut folded_removed = folded();
    folded_removed.branch.remove(TRANSFORMER);
    let stale = ac_opf_with_options(&folded_removed, &options).unwrap();
    let relative = (stale.objective - want.objective).abs() / want.objective.abs();
    assert!(relative > 1e-4, "{} vs {}", stale.objective, want.objective);
}
