// The kept OPF models (AcOpfModel, DcOpfModel): a model built from one case and solving a copy
//   with new values gives, bit for bit, what the one-shot functions give for that copy, on both
//   solvers; and a case whose topology, or for AC a pattern-deciding value, differs is refused.

mod common;

use ybus_powerflow::case::{Case, GenCost};
use ybus_powerflow::opf::ac::{ac_opf_with_options, AcOpfModel, AcOpfOptions, AcOpfSolution};
use ybus_powerflow::opf::dc::{dc_opf_with_options, DcOpfModel, DcOpfOptions, DcOpfSolution};
use ybus_powerflow::powerflow_error::PowerflowError;

/// Every result bit: the solver's point and multipliers, the objective and the iterations
fn ac_bits(s: &AcOpfSolution) -> String {
    format!(
        "{:?} {:?} {:?} {:?} {:?}",
        s.solver.x,
        s.solver.multipliers,
        s.objective.to_bits(),
        s.iterations,
        s.status
    )
}

fn dc_bits(s: &DcOpfSolution) -> String {
    format!(
        "{:?} {:?} {:?} {:?} {:?}",
        s.solver.x,
        s.solver.multipliers,
        s.objective.to_bits(),
        s.iterations,
        s.status
    )
}

/// The case with every load scaled by `factor`
fn scaled(case: &Case, factor: f64) -> Case {
    let mut case = case.clone();
    for bus in &mut case.bus {
        bus.pd_mw *= factor;
        bus.qd_mvar *= factor;
    }
    case
}

fn cases() -> Vec<(&'static str, Case)> {
    ["case9", "case30", "case118", "case2383wp"]
        .into_iter()
        .filter_map(|name| {
            common::matpower_case(name).map(|path| (name, Case::from_path(&path).unwrap()))
        })
        .collect()
}

#[test]
fn an_ac_model_solves_a_changed_case_as_the_one_shot_function_does() {
    for (name, case) in cases() {
        let model = AcOpfModel::new(&case).unwrap();
        for options in [AcOpfOptions::default(), AcOpfOptions::mehrotra()] {
            for factor in [1.0, 0.95, 1.04] {
                let changed = scaled(&case, factor);
                let kept = model.solve(&changed, &options).unwrap();
                let one_shot = ac_opf_with_options(&changed, &options).unwrap();
                assert_eq!(ac_bits(&kept), ac_bits(&one_shot), "{name} ×{factor}");
            }
        }
    }
}

#[test]
fn a_dc_model_solves_a_changed_case_as_the_one_shot_function_does() {
    for (name, case) in cases() {
        let model = DcOpfModel::new(&case).unwrap();
        for options in [DcOpfOptions::default(), DcOpfOptions::mehrotra()] {
            for factor in [1.0, 0.95, 1.04] {
                let changed = scaled(&case, factor);
                let kept = model.solve(&changed, &options).unwrap();
                let one_shot = dc_opf_with_options(&changed, &options).unwrap();
                assert_eq!(dc_bits(&kept), dc_bits(&one_shot), "{name} ×{factor}");
            }
        }
    }
}

#[test]
fn a_model_refuses_a_case_it_cant_solve() {
    let Some(path) = common::matpower_case("case30") else {
        return;
    };
    let case = Case::from_path(&path).unwrap();
    let (ac, dc) = (
        AcOpfModel::new(&case).unwrap(),
        DcOpfModel::new(&case).unwrap(),
    );
    let options = (AcOpfOptions::default(), DcOpfOptions::default());
    let field = |e: PowerflowError| match e {
        PowerflowError::TopologyMismatch { table, row, field } => (table, row, field),
        other => panic!("expected a TopologyMismatch, not {other:?}"),
    };

    // A branch out of service: a topology change, refused by both
    let mut outage = case.clone();
    outage.branch[3].status = false;
    assert_eq!(
        field(ac.solve(&outage, &options.0).unwrap_err()),
        ("branch", 4, "status")
    );
    assert_eq!(
        field(dc.solve(&outage, &options.1).unwrap_err()),
        ("branch", 4, "status")
    );

    // A flow limit removed (rate_a 0): the AC model's limited set changes. The DC model keeps
    //   nothing that depends on it, and solves the case as the one-shot function does
    let row = case
        .branch
        .iter()
        .position(|b| b.rate_a_mva != 0.0)
        .expect("a limited branch");
    let mut unlimited = case.clone();
    unlimited.branch[row].rate_a_mva = 0.0;
    assert_eq!(
        field(ac.solve(&unlimited, &options.0).unwrap_err()),
        ("branch", row + 1, "rate_a_mva")
    );
    assert_eq!(
        dc_bits(&dc.solve(&unlimited, &options.1).unwrap()),
        dc_bits(&dc_opf_with_options(&unlimited, &options.1).unwrap())
    );

    // A quadratic cost coefficient set to zero changes the AC Hessian's pattern
    let mut linear = case.clone();
    let costs = linear.gencost.as_mut().expect("costs");
    let gen = costs
        .iter()
        .position(|c| matches!(c, GenCost::Polynomial { coeffs, .. } if coeffs.len() == 3 && coeffs[0] != 0.0))
        .expect("a quadratic cost");
    if let GenCost::Polynomial { coeffs, .. } = &mut costs[gen] {
        coeffs[0] = 0.0;
    }
    assert_eq!(
        field(ac.solve(&linear, &options.0).unwrap_err()),
        ("gencost", gen + 1, "quadratic")
    );
}
