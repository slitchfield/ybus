//! Generator costs: their values, as MATPOWER's `totcost` computes them, and their optimal power
//! flow form, as `opf_setup` and `makeAy` build it, on hand-computed fixtures and on case files.

mod common;

use common::dc_fixture::{fixture, B, G, L};
use powerflow::case::{Case, GenCost};
use powerflow::opf::cost::{total_cost, CostVariable, OpfCosts, PiecewiseCost};
use powerflow::powerflow_error::PowerflowError;

fn polynomial(coeffs: &[f64]) -> GenCost {
    GenCost::Polynomial {
        startup: 100.0,
        shutdown: 50.0,
        coeffs: coeffs.to_vec(),
    }
}

fn piecewise(points: &[(f64, f64)]) -> GenCost {
    GenCost::PiecewiseLinear {
        startup: 100.0,
        shutdown: 50.0,
        points: points.to_vec(),
    }
}

/// A 100 MVA case with `gens` generators, all in service at bus 1, and the given costs
fn with_costs(gens: usize, costs: Vec<GenCost>) -> Case {
    let generators: Vec<G> = (0..gens).map(|_| G(1, 0.0, true)).collect();
    let mut case = fixture(
        &[B(1, 3, 0.0, 0.0, 0.0), B(2, 1, 50.0, 0.0, 0.0)],
        &generators,
        &[L(1, 2, 0.1, 0.0, 0.0, true)],
    );
    case.gencost = Some(costs);
    case
}

#[test]
fn polynomials_evaluate_with_the_constant_last_and_no_startup_cost() {
    assert_eq!(
        polynomial(&[0.01, 40.0, 5.0]).evaluate(100.0),
        100.0 + 4000.0 + 5.0
    );
    assert_eq!(polynomial(&[2.0, 0.0, 0.0, 1.0]).evaluate(3.0), 55.0);
    assert_eq!(polynomial(&[7.0]).evaluate(3.0), 7.0);
    assert_eq!(polynomial(&[]).evaluate(3.0), 0.0);
}

#[test]
fn piecewise_linear_costs_interpolate_and_extend_their_end_segments() {
    let cost = piecewise(&[(10.0, 100.0), (20.0, 300.0), (40.0, 900.0)]);
    assert_eq!(cost.evaluate(15.0), 200.0);
    assert_eq!(cost.evaluate(20.0), 300.0);
    assert_eq!(cost.evaluate(30.0), 600.0);
    assert_eq!(cost.evaluate(40.0), 900.0);
    // Past the ends, along the first segment (slope 20) and the last (slope 30)
    assert_eq!(cost.evaluate(0.0), -100.0);
    assert_eq!(cost.evaluate(50.0), 1200.0);
    assert_eq!(piecewise(&[(10.0, 100.0)]).evaluate(50.0), 100.0);
    assert_eq!(piecewise(&[]).evaluate(50.0), 0.0);
}

#[test]
fn total_cost_adds_reactive_costs_only_when_present_and_asked_for() {
    let costs = vec![
        polynomial(&[1.0, 0.0]),
        polynomial(&[2.0, 0.0]),
        polynomial(&[10.0, 0.0]),
        piecewise(&[(0.0, 0.0), (10.0, 50.0)]),
    ];
    let (pg, qg) = ([3.0, 4.0], [1.0, 2.0]);
    assert_eq!(
        total_cost(&costs, &pg, Some(&qg)),
        Ok(3.0 + 8.0 + 10.0 + 10.0)
    );
    assert_eq!(total_cost(&costs, &pg, None), Ok(11.0));
    assert_eq!(total_cost(&costs[..2], &pg, Some(&qg)), Ok(11.0));
    assert!(matches!(
        total_cost(&costs[..3], &pg, None),
        Err(PowerflowError::CountMismatch {
            field: "gencost",
            found: 3,
            ..
        })
    ));
    assert!(matches!(
        total_cost(&costs, &pg, Some(&qg[..1])),
        Err(PowerflowError::CountMismatch {
            field: "qg_mvar",
            found: 1,
            ..
        })
    ));
}

#[test]
fn polynomial_costs_convert_to_per_unit_after_dropping_zero_leading_coefficients() {
    let case = with_costs(
        5,
        vec![
            polynomial(&[0.01, 40.0, 5.0]),
            // [0, 0, 3, 7] is linear, and [0, 5] stays linear with a zero slope
            polynomial(&[0.0, 0.0, 3.0, 7.0]),
            polynomial(&[0.0, 5.0]),
            polynomial(&[9.0]),
            // One segment becomes the line through its points: 20 $/MWh, 100 $/hr at 0 MW
            piecewise(&[(10.0, 300.0), (60.0, 1300.0)]),
        ],
    );
    let costs = OpfCosts::new(&case, &[0, 1, 2, 3, 4], true).unwrap();
    assert_eq!(
        costs.active.quadratic,
        [2.0 * 0.01 * 100.0 * 100.0, 0.0, 0.0, 0.0, 0.0]
    );
    assert_eq!(costs.active.linear, [4000.0, 300.0, 0.0, 0.0, 2000.0]);
    assert_eq!(costs.active.constant, [5.0, 7.0, 5.0, 9.0, 100.0]);
    assert!(!costs.active.is_linear());
    assert_eq!(costs.reactive, None);
    assert!(costs.piecewise.is_empty());

    // At 0.5 p.u. each: ½·200·0.25 + 4000·0.5 + 5, 150 + 7, 5, 9, 1000 + 100
    let pg = [0.5; 5];
    assert_eq!(
        costs.evaluate(&pg, &[], &[]),
        2030.0 + 157.0 + 5.0 + 9.0 + 1100.0
    );
    let mw: Vec<f64> = pg.iter().map(|p| p * 100.0).collect();
    assert_eq!(
        total_cost(case.gencost.as_ref().unwrap(), &mw, None),
        Ok(3301.0)
    );
}

#[test]
fn piecewise_linear_costs_become_segment_constraints_on_helper_variables() {
    let case = with_costs(
        3,
        vec![
            polynomial(&[0.0, 40.0, 0.0]),
            piecewise(&[(0.0, 0.0), (50.0, 500.0), (100.0, 1500.0)]),
            polynomial(&[1.0, 0.0]),
            // Reactive costs: generator 0's piecewise linear, generator 1's quadratic, and
            //   generator 2's linear
            piecewise(&[(-20.0, 40.0), (0.0, 0.0), (20.0, 60.0)]),
            polynomial(&[0.5, 0.0, 0.0]),
            polynomial(&[2.0, 1.0]),
        ],
    );
    let costs = OpfCosts::new(&case, &[0, 1, 2], true).unwrap();
    assert!(costs.active.is_linear());
    assert_eq!(costs.active.linear, [4000.0, 0.0, 100.0]);
    let reactive = costs.reactive.as_ref().unwrap();
    assert_eq!(reactive.quadratic, [0.0, 10000.0, 0.0]);
    assert_eq!(reactive.linear, [0.0, 0.0, 200.0]);
    assert_eq!(reactive.constant, [0.0, 0.0, 1.0]);
    // Real power costs first, then reactive; slopes per p.u., bounds m p − c
    assert_eq!(
        costs.piecewise,
        [
            PiecewiseCost {
                variable: CostVariable::Active(1),
                slopes: vec![1000.0, 2000.0],
                bounds: vec![0.0, 500.0],
                last_cost: 1500.0,
            },
            PiecewiseCost {
                variable: CostVariable::Reactive(0),
                slopes: vec![-200.0, 300.0],
                bounds: vec![0.0, 0.0],
                last_cost: 60.0,
            },
        ]
    );

    // A DC optimal power flow ignores reactive costs, piecewise linear ones included
    let dc = OpfCosts::new(&case, &[0, 1, 2], false).unwrap();
    assert_eq!(dc.reactive, None);
    assert_eq!(dc.piecewise.len(), 1);
}

#[test]
fn only_the_given_generators_are_costed_in_their_order() {
    let case = with_costs(
        3,
        vec![
            polynomial(&[1.0, 0.0]),
            polynomial(&[2.0, 0.0]),
            piecewise(&[(0.0, 0.0), (50.0, 500.0), (100.0, 1500.0)]),
        ],
    );
    let costs = OpfCosts::new(&case, &[2, 0], true).unwrap();
    assert_eq!(costs.active.linear, [0.0, 100.0]);
    assert_eq!(costs.piecewise[0].variable, CostVariable::Active(0));
}

#[test]
fn unsupported_and_missing_costs_are_refused() {
    let refused = |cost: GenCost| {
        let case = with_costs(2, vec![polynomial(&[1.0, 0.0]), cost]);
        OpfCosts::new(&case, &[0, 1], true).unwrap_err()
    };
    assert!(matches!(
        refused(polynomial(&[0.0, 1.0, 0.0, 0.0, 0.0])),
        PowerflowError::UnsupportedCost { row: 1, .. }
    ));
    assert!(matches!(
        refused(piecewise(&[(0.0, 0.0)])),
        PowerflowError::UnsupportedCost { row: 1, .. }
    ));
    assert!(matches!(
        refused(piecewise(&[(5.0, 0.0), (5.0, 10.0)])),
        PowerflowError::UnsupportedCost { row: 1, .. }
    ));
    assert!(matches!(
        refused(piecewise(&[(0.0, 0.0), (5.0, 10.0), (5.0, 20.0)])),
        PowerflowError::UnsupportedCost { row: 1, .. }
    ));

    // A reactive cost's row counts past the generators
    let case = with_costs(
        1,
        vec![polynomial(&[1.0, 0.0]), polynomial(&[1.0, 0.0, 0.0, 0.0])],
    );
    assert_eq!(
        OpfCosts::new(&case, &[0], true).unwrap_err(),
        PowerflowError::UnsupportedCost {
            row: 1,
            reason: "polynomial costs above quadratic aren't supported"
        }
    );
    assert!(OpfCosts::new(&case, &[0], false).is_ok());

    let case = with_costs(1, vec![polynomial(&[1.0, 0.0])]);
    assert!(matches!(
        OpfCosts::new(&case, &[1], true),
        Err(PowerflowError::RowOutOfRange {
            table: "gen",
            row: 1,
            len: 1
        })
    ));
    let mut case = case;
    case.gencost = None;
    assert_eq!(
        OpfCosts::new(&case, &[0], true),
        Err(PowerflowError::MissingField { field: "gencost" })
    );
}

/// At each case's own dispatch, the optimal power flow's form of the costs matches the costs
/// evaluated directly: the committed cases' quadratic costs, case30pwl's and case_RTS_GMLC's
/// piecewise linear costs, and case9Q's and case30Q's reactive costs.
///
/// Each helper variable is set at the largest of its segments' lines, which is the cost when the
/// cost is convex and above it where it isn't. One cost isn't: case_RTS_GMLC's `gencost` row 73
/// (0-based), whose rounded data has a middle slope of 8.10345 between two of 8.10352, comes out
/// 9.1e-5 $/hr above at its dispatch. So the polynomial costs are checked with each helper variable at its direct
/// value, and each helper variable separately.
#[test]
fn case_costs_in_opf_form_match_their_direct_values() {
    let mut nonconvex = vec![];
    for name in [
        "case14",
        "case30",
        "case118",
        "case30pwl",
        "case_RTS_GMLC",
        "case9Q",
        "case30Q",
    ] {
        let Some(path) = common::matpower_case(name) else {
            continue;
        };
        let case = Case::from_path(&path).unwrap();
        let rows: Vec<usize> = (0..case.generator.len())
            .filter(|&row| case.generator[row].status)
            .collect();
        let costs = OpfCosts::new(&case, &rows, true).unwrap();

        let pg_mw: Vec<f64> = rows.iter().map(|&row| case.generator[row].pg_mw).collect();
        let qg_mvar: Vec<f64> = rows
            .iter()
            .map(|&row| case.generator[row].qg_mvar)
            .collect();
        let base = case.base_mva;
        let pg: Vec<f64> = pg_mw.iter().map(|p| p / base).collect();
        let qg: Vec<f64> = qg_mvar.iter().map(|q| q / base).collect();

        let all = case.gencost.as_ref().unwrap();
        let ng = case.generator.len();
        let mut direct_y = vec![];
        for cost in &costs.piecewise {
            let (row, p) = match cost.variable {
                CostVariable::Active(i) => (rows[i], pg[i]),
                CostVariable::Reactive(i) => (ng + rows[i], qg[i]),
            };
            let direct = all[row].evaluate(p * base);
            let y = cost
                .slopes
                .iter()
                .zip(&cost.bounds)
                .map(|(m, b)| m * p - b)
                .fold(f64::NEG_INFINITY, f64::max);
            let convex = cost.slopes.windows(2).all(|m| m[0] <= m[1]);
            if convex {
                let relative = (y - direct).abs() / direct.abs().max(1.0);
                assert!(
                    relative < 1e-12,
                    "{name} row {row}: {y} vs {direct} ({relative:e})"
                );
            } else {
                assert!(y >= direct, "{name} row {row}: {y} below {direct}");
                nonconvex.push((name, row, y - direct));
            }
            direct_y.push(direct);
        }

        // The in-service generators' rows, reactive ones after
        let mut selected: Vec<GenCost> = rows.iter().map(|&row| all[row].clone()).collect();
        if all.len() == 2 * ng {
            selected.extend(rows.iter().map(|&row| all[ng + row].clone()));
        }
        let direct = total_cost(&selected, &pg_mw, Some(&qg_mvar)).unwrap();
        let opf = costs.evaluate(&pg, &qg, &direct_y);
        let relative = (opf - direct).abs() / direct.abs();
        assert!(relative < 1e-12, "{name}: {opf} vs {direct} ({relative:e})");
        if name == "case30pwl" || name == "case_RTS_GMLC" {
            assert!(
                !costs.piecewise.is_empty(),
                "{name} has piecewise linear costs"
            );
        }
        if name.ends_with('Q') {
            assert!(costs.reactive.is_some(), "{name} has reactive costs");
        }
    }
    for &(name, row, above) in &nonconvex {
        assert!(
            name == "case_RTS_GMLC" && row == 73 && (above - 9.1e-5).abs() < 1e-6,
            "unexpected nonconvex cost: {name} row {row}, {above:e} above"
        );
    }
}
