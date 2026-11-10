//! The AC optimal power flow on small cases. An AC optimum can't be worked out by hand, so each
//! solution is checked against what any optimum must satisfy: the power balance at every bus
//! (from the reported flows), every limit, complementary slackness (a positive multiplier only
//! on a binding limit), marginal costs equal to the bus prices at generators inside their
//! limits, and the objective equal to the costs of the dispatch. Then each case checks its
//! feature: a binding apparent-power limit, binding voltage limits, reactive costs, a
//! dispatchable load at its power factor, and the refusals. Each solving case runs on both
//! solvers, MIPS (the default) and the predictor-corrector method.

mod common;

use common::ac_fixture::{ac_fixture, AcBranch, AcBus, AcGen};
use powerflow::case::{BusType, Case, GenCost};
use powerflow::opf::ac::{ac_opf, ac_opf_with_options, AcOpfOptions, AcOpfSolution};
use powerflow::opf::cost::total_cost;
use powerflow::opf::{InteriorPointError, OpfSolver};
use powerflow::powerflow_error::PowerflowError;

/// Both solvers, at MATPOWER's AC OPF settings
fn solvers() -> [(&'static str, AcOpfOptions); 2] {
    [
        ("MIPS", AcOpfOptions::default()),
        ("Mehrotra", AcOpfOptions::mehrotra()),
    ]
}

/// MIPS stops at a feasibility tolerance of 5e-6 p.u. and gradient and complementarity
/// tolerances of 1e-6, so the checks allow about that, in MW and $/MWh
const MW: f64 = 1e-3;
const PRICE: f64 = 1e-3;
/// Complementarity, a multiplier times its limit's slack ($/hr). An interior-point method ends
/// with these small rather than zero: a slack limit can keep a small multiplier (0.0012 $/MVA-hr
/// at 0.29 MVA from its limit, in the flow-limit case below), which MIPS only zeroes below its
/// `mu_threshold`.
const COMPLEMENTARITY: f64 = 1e-2;

fn quadratic(c2: f64, c1: f64) -> GenCost {
    GenCost::Polynomial {
        startup: 0.0,
        shutdown: 0.0,
        coeffs: vec![c2, c1, 0.0],
    }
}

/// The checks every optimum must pass; see the module documentation.
fn assert_optimal(case: &Case, s: &AcOpfSolution) {
    assert!(
        s.converged(),
        "{:?} after {} iterations",
        s.status,
        s.iterations
    );
    let labels: Vec<usize> = case.bus.iter().map(|bus| bus.bus_label).collect();
    let index = |label: usize| labels.iter().position(|&l| l == label).unwrap();

    // The power balance: generation less load equals the power leaving through branches and
    //   shunts
    let mut balance: Vec<(f64, f64)> = case
        .bus
        .iter()
        .zip(&s.vm_pu)
        .map(|(bus, vm)| {
            let (p, q) = (-bus.pd_mw, -bus.qd_mvar);
            (p - bus.gs_mw * vm * vm, q + bus.bs_mvar * vm * vm)
        })
        .collect();
    for (k, g) in case.generator.iter().enumerate() {
        let b = index(g.bus_label);
        balance[b].0 += s.pg_mw[k];
        balance[b].1 += s.qg_mvar[k];
    }
    for (k, branch) in case.branch.iter().enumerate() {
        let (f, t) = (index(branch.from_label), index(branch.to_label));
        balance[f].0 -= s.pf_mw[k];
        balance[f].1 -= s.qf_mvar[k];
        balance[t].0 -= s.pt_mw[k];
        balance[t].1 -= s.qt_mvar[k];
    }
    for (i, (p, q)) in balance.iter().enumerate() {
        assert!(
            p.abs() < MW && q.abs() < MW,
            "bus {i} unbalanced by {p} MW, {q} MVAr"
        );
    }

    // The limits, and complementary slackness
    for (i, bus) in case.bus.iter().enumerate() {
        let vm = s.vm_pu[i];
        assert!(
            vm >= bus.vmin_pu - 1e-6 && vm <= bus.vmax_pu + 1e-6,
            "bus {i}: Vm {vm}"
        );
        let upper = s.mu_vmax[i] * (bus.vmax_pu - vm);
        let lower = s.mu_vmin[i] * (vm - bus.vmin_pu);
        assert!(
            upper.abs() < COMPLEMENTARITY,
            "bus {i}: mu_vmax × slack {upper}"
        );
        assert!(
            lower.abs() < COMPLEMENTARITY,
            "bus {i}: mu_vmin × slack {lower}"
        );
    }
    for (k, g) in case.generator.iter().enumerate() {
        let (p, q) = (s.pg_mw[k], s.qg_mvar[k]);
        assert!(
            p >= g.pmin_mw - MW && p <= g.pmax_mw + MW,
            "gen {k}: Pg {p}"
        );
        assert!(
            q >= g.qmin_mvar - MW && q <= g.qmax_mvar + MW,
            "gen {k}: Qg {q}"
        );
        let upper = s.mu_pmax[k] * (g.pmax_mw - p);
        let lower = s.mu_pmin[k] * (p - g.pmin_mw);
        assert!(
            upper.abs() < COMPLEMENTARITY,
            "gen {k}: mu_pmax × slack {upper}"
        );
        assert!(
            lower.abs() < COMPLEMENTARITY,
            "gen {k}: mu_pmin × slack {lower}"
        );
    }
    for (k, branch) in case.branch.iter().enumerate() {
        let rate = branch.rate_a_mva;
        let from = s.pf_mw[k].hypot(s.qf_mvar[k]);
        let to = s.pt_mw[k].hypot(s.qt_mvar[k]);
        if rate > 0.0 {
            assert!(
                from <= rate + MW && to <= rate + MW,
                "branch {k}: {from}, {to} MVA"
            );
        }
        let (sf, st) = (s.mu_sf[k] * (rate - from), s.mu_st[k] * (rate - to));
        assert!(sf.abs() < COMPLEMENTARITY, "branch {k}: mu_sf × slack {sf}");
        assert!(st.abs() < COMPLEMENTARITY, "branch {k}: mu_st × slack {st}");
    }

    // Marginal costs equal the bus prices where a generator is inside its limits
    let costs = case.gencost.as_ref().unwrap();
    let ng = case.generator.len();
    let marginal = |cost: &GenCost, x: f64| match cost {
        GenCost::Polynomial { coeffs, .. } if coeffs.len() == 3 => 2.0 * coeffs[0] * x + coeffs[1],
        GenCost::Polynomial { coeffs, .. } if coeffs.len() == 2 => coeffs[0],
        other => panic!("unexpected cost {other:?}"),
    };
    for (k, g) in case.generator.iter().enumerate() {
        let b = index(g.bus_label);
        let (p, q) = (s.pg_mw[k], s.qg_mvar[k]);
        let inside_p = p > g.pmin_mw + 1e-2 && p < g.pmax_mw - 1e-2;
        if inside_p && s.mu_pmin[k] < PRICE && s.mu_pmax[k] < PRICE {
            let mc = marginal(&costs[k], p);
            assert!(
                (mc - s.lam_p[b]).abs() < PRICE,
                "gen {k}: {mc} $/MWh vs LAM_P {}",
                s.lam_p[b]
            );
        }
        let inside_q = q > g.qmin_mvar + 1e-2 && q < g.qmax_mvar - 1e-2;
        let reactive = if costs.len() == 2 * ng {
            marginal(&costs[ng + k], q)
        } else {
            0.0
        };
        if inside_q && s.mu_qmin[k] < PRICE && s.mu_qmax[k] < PRICE {
            let lam_q = s.lam_q[b];
            assert!(
                (reactive - lam_q).abs() < PRICE,
                "gen {k}: {reactive} vs LAM_Q {lam_q}"
            );
        }
    }

    // The objective is the costs of the dispatch
    let direct = total_cost(costs, &s.pg_mw, Some(&s.qg_mvar)).unwrap();
    assert!(
        (s.objective - direct).abs() < 1e-6 * direct.abs().max(1.0),
        "objective {} vs {direct}",
        s.objective
    );
}

/// Three buses in a ring with resistive lines, loads at buses 2 and 3, a cheap generator at bus
/// 1 and a dearer one at bus 3.
fn three_bus() -> Case {
    let mut case = ac_fixture(
        &[
            AcBus(1, 3, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            AcBus(2, 1, 90.0, 30.0, 0.0, 0.0, 1.0, 0.0),
            AcBus(3, 2, 60.0, 20.0, 0.0, 10.0, 1.0, 0.0),
        ],
        &[
            AcGen(1, 0.0, 0.0, 100.0, -100.0, 1.0, true),
            AcGen(3, 0.0, 0.0, 100.0, -100.0, 1.0, true),
        ],
        &[
            AcBranch(1, 2, 0.02, 0.1, 0.02, 0.0, 0.0, true),
            AcBranch(2, 3, 0.03, 0.12, 0.02, 0.0, 0.0, true),
            AcBranch(1, 3, 0.02, 0.09, 0.01, 0.0, 0.0, true),
        ],
    );
    case.gencost = Some(vec![quadratic(0.02, 10.0), quadratic(0.05, 20.0)]);
    case
}

#[test]
fn an_unconstrained_ring_is_optimal() {
    let case = three_bus();
    for (solver, options) in solvers() {
        let s = ac_opf_with_options(&case, &options).unwrap();
        assert_optimal(&case, &s);
        // Losses make the far bus dearer than the generator bus
        assert!(s.lam_p[1] > s.lam_p[0], "{solver}: {:?}", s.lam_p);
    }
}

#[test]
fn a_binding_apparent_power_limit_raises_the_price_beyond_it() {
    let mut case = three_bus();
    case.branch[0].rate_a_mva = 40.0;
    case.branch[2].rate_a_mva = 60.0;
    for (solver, options) in solvers() {
        let s = ac_opf_with_options(&case, &options).unwrap();
        assert_optimal(&case, &s);
        let binding = (0..case.branch.len()).find(|&k| s.mu_sf[k] > PRICE || s.mu_st[k] > PRICE);
        assert!(
            binding.is_some(),
            "{solver}: no flow limit binds: {:?} {:?}",
            s.mu_sf,
            s.mu_st
        );
        let unconstrained = ac_opf_with_options(&three_bus(), &options).unwrap();
        assert!(s.objective > unconstrained.objective, "{solver}");
    }
}

#[test]
fn voltage_limits_bind_with_their_multipliers() {
    // A narrow voltage band at the loaded bus, which reactive demand pushes against
    let mut case = three_bus();
    case.bus[1].qd_mvar = 80.0;
    case.bus[1].vmin_pu = 0.99;
    case.bus[1].vmax_pu = 1.01;
    for (solver, options) in solvers() {
        let s = ac_opf_with_options(&case, &options).unwrap();
        assert_optimal(&case, &s);
        let binding = (0..case.bus.len()).any(|i| s.mu_vmin[i] > PRICE || s.mu_vmax[i] > PRICE);
        assert!(binding, "{solver}: no voltage limit binds: {:?}", s.vm_pu);
    }
}

#[test]
fn reactive_power_costs_set_the_reactive_prices() {
    let mut case = three_bus();
    case.gencost
        .as_mut()
        .unwrap()
        .extend([quadratic(0.01, 1.0), quadratic(0.02, 2.0)]);
    for (solver, options) in solvers() {
        let s = ac_opf_with_options(&case, &options).unwrap();
        assert_optimal(&case, &s);
        let without = ac_opf_with_options(&three_bus(), &options).unwrap();
        assert!(s.objective > without.objective, "{solver}");
    }
}

#[test]
fn a_dispatchable_load_keeps_its_power_factor() {
    // Up to 50 MW of load at bus 2 worth 40 $/MWh, consuming 0.4 MVAr per MW (Qmin = −20 at
    //   Pmin = −50)
    let mut case = three_bus();
    case.generator.push(case.generator[0].clone());
    let load = case.generator.last_mut().unwrap();
    load.bus_label = 2;
    (load.pg_mw, load.qg_mvar) = (0.0, 0.0);
    (load.pmin_mw, load.pmax_mw) = (-50.0, 0.0);
    (load.qmin_mvar, load.qmax_mvar) = (-20.0, 0.0);
    case.gencost.as_mut().unwrap().push(quadratic(0.0, 40.0));
    for (solver, options) in solvers() {
        let s = ac_opf_with_options(&case, &options).unwrap();
        assert_optimal(&case, &s);
        let (p, q) = (s.pg_mw[2], s.qg_mvar[2]);
        assert!(p < -1.0, "{solver}: the load takes power: {p} MW");
        assert!((q - 0.4 * p).abs() < MW, "{solver}: {q} MVAr at {p} MW");
    }
}

#[test]
fn unsupported_generators_and_cases_are_refused() {
    // A PQ capability curve
    let mut case = three_bus();
    let g = &mut case.generator[1];
    (g.pc1_mw, g.pc2_mw, g.qc1max_mvar, g.qc2max_mvar) = (10.0, 80.0, 90.0, 50.0);
    assert!(matches!(
        ac_opf(&case),
        Err(PowerflowError::UnsupportedGenerator { row: 1, .. })
    ));
    // A dispatchable load with both reactive limits nonzero, and one off its power factor
    let mut case = three_bus();
    let g = &mut case.generator[1];
    (g.pmin_mw, g.pmax_mw, g.qmin_mvar, g.qmax_mvar) = (-50.0, 0.0, -20.0, 10.0);
    assert!(matches!(
        ac_opf(&case),
        Err(PowerflowError::UnsupportedGenerator { row: 1, .. })
    ));
    let mut case = three_bus();
    let g = &mut case.generator[1];
    (g.pmin_mw, g.pmax_mw, g.qmin_mvar, g.qmax_mvar) = (-50.0, 0.0, -20.0, 0.0);
    (g.pg_mw, g.qg_mvar) = (-10.0, 0.0);
    assert!(matches!(
        ac_opf(&case),
        Err(PowerflowError::UnsupportedGenerator { row: 1, .. })
    ));
    // No reference bus, and no costs
    let mut case = three_bus();
    case.bus[0].bus_type = BusType::PV;
    assert_eq!(ac_opf(&case).unwrap_err(), PowerflowError::NoReferenceBus);
    let mut case = three_bus();
    case.gencost = None;
    assert_eq!(
        ac_opf(&case).unwrap_err(),
        PowerflowError::MissingField { field: "gencost" }
    );
}

#[test]
fn a_solver_option_out_of_range_is_an_error() {
    for (solver, mut options) in solvers() {
        match &mut options.solver {
            OpfSolver::Mips(mips) => mips.xi = 0.3,
            OpfSolver::Mehrotra(mehrotra) => mehrotra.xi = 0.3,
            _ => unreachable!(),
        }
        assert_eq!(
            ac_opf_with_options(&three_bus(), &options).unwrap_err(),
            PowerflowError::Solver {
                error: InteriorPointError::InvalidOption { name: "xi" }
            },
            "{solver}"
        );
    }
}
