// Warm starts across outages (map_start, OpfStart::Mapped): a model's own solution maps onto
//   itself as WarmStart::from_solution does; a branch or generator outage loses exactly its own
//   variables and constraints, the rest carried over unchanged, both ways; a mapped warm solve of
//   an outage meets the warm-start survey's bar against the outage's cold solve; and a case with
//   other tables is refused.

mod common;

use std::collections::VecDeque;

use common::reference::relative_difference;
use ybus_powerflow::case::{BusType, Case, GenCost};
use ybus_powerflow::opf::ac::{AcOpfModel, AcOpfOptions};
use ybus_powerflow::opf::dc::{DcOpfModel, DcOpfOptions};
use ybus_powerflow::opf::{
    History, KktStatistics, Multipliers, OpfSolver, OpfStart, Solution, Status, WarmStart,
};
use ybus_powerflow::powerflow_error::PowerflowError;

fn case(name: &str) -> Option<Case> {
    common::matpower_case(name).map(|path| Case::from_path(&path).unwrap())
}

/// Whether `case`'s connected buses stay connected with branch `out` out of service as well
fn connected_without(case: &Case, out: usize) -> bool {
    let buses: Vec<usize> = (0..case.bus.len())
        .filter(|&i| case.bus[i].bus_type != BusType::Isolated)
        .collect();
    let index = |label: usize| buses.iter().position(|&i| case.bus[i].bus_label == label);
    let mut adjacent = vec![vec![]; buses.len()];
    for (row, branch) in case.branch.iter().enumerate() {
        if row == out || !branch.status {
            continue;
        }
        if let (Some(f), Some(t)) = (index(branch.from_label), index(branch.to_label)) {
            adjacent[f].push(t);
            adjacent[t].push(f);
        }
    }
    let mut seen = vec![false; buses.len()];
    let mut queue = VecDeque::from([0]);
    seen[0] = true;
    while let Some(i) = queue.pop_front() {
        for &j in &adjacent[i] {
            if !seen[j] {
                seen[j] = true;
                queue.push_back(j);
            }
        }
    }
    seen.iter().all(|&s| s)
}

/// The first in-service branch whose outage doesn't island a bus, with a flow limit if
///   `limited`, and no angle limit (`angmin_deg`, `angmax_deg` at ±360 or 0)
fn branch_outage(case: &Case, limited: bool) -> usize {
    (0..case.branch.len())
        .find(|&row| {
            let branch = &case.branch[row];
            let unlimited_angle = |a: f64| a == 0.0 || a.abs() >= 360.0;
            branch.status
                && (!limited || branch.rate_a_mva != 0.0)
                && unlimited_angle(branch.angmin_deg)
                && unlimited_angle(branch.angmax_deg)
                && connected_without(case, row)
        })
        .expect("a branch that isn't a bridge")
}

fn without_branch(case: &Case, row: usize) -> Case {
    let mut case = case.clone();
    case.branch[row].status = false;
    case
}

fn without_generator(case: &Case, row: usize) -> Case {
    let mut case = case.clone();
    case.generator[row].status = false;
    case
}

/// A solution with `warm`'s point and multipliers, as a previous solution to map from
fn as_solution(warm: &WarmStart) -> Solution {
    Solution {
        x: warm.x.clone(),
        f: 0.0,
        status: Status::Converged,
        iterations: 0,
        multipliers: warm.multipliers.clone(),
        iterate_multipliers: warm.multipliers.clone(),
        history: History::Mips(vec![]),
        analyses: 0,
        refinement_steps: 0,
        inertia_corrections: 0,
        kkt: KktStatistics::default(),
    }
}

/// Every vector of a warm start, named, for comparing
fn vectors(w: &WarmStart) -> Vec<(&'static str, &Vec<f64>)> {
    let m: &Multipliers = &w.multipliers;
    vec![
        ("x", &w.x),
        ("eq_nonlin", &m.eq_nonlin),
        ("ineq_nonlin", &m.ineq_nonlin),
        ("mu_l", &m.mu_l),
        ("mu_u", &m.mu_u),
        ("lower", &m.lower),
        ("upper", &m.upper),
    ]
}

/// The entries of `base` that the round trip through an outage changed, and how many entries
///   the outage's warm start lacks
fn changes(base: &WarmStart, outage: &WarmStart, back: &WarmStart) -> (usize, usize) {
    let mut changed = 0;
    let mut lost = 0;
    for (((name, b), (_, o)), (_, r)) in vectors(base)
        .into_iter()
        .zip(vectors(outage))
        .zip(vectors(back))
    {
        assert_eq!(b.len(), r.len(), "{name}: the round trip's length");
        changed += b
            .iter()
            .zip(r)
            .filter(|(b, r)| b.to_bits() != r.to_bits())
            .count();
        lost += b.len() - o.len();
    }
    (changed, lost)
}

#[test]
fn a_models_own_solution_maps_onto_itself() {
    for name in ["case9", "case30", "case30pwl", "case118"] {
        let Some(case) = case(name) else { continue };
        let model = DcOpfModel::new(&case).unwrap();
        let solution = model.solve(&case, &DcOpfOptions::default()).unwrap();
        let mapped = model
            .map_start(&case, &model, &case, &solution.solver)
            .unwrap();
        assert_eq!(
            mapped,
            WarmStart::from_solution(&solution.solver),
            "{name} DC"
        );
        let model = AcOpfModel::new(&case).unwrap();
        let solution = model.solve(&case, &AcOpfOptions::mehrotra()).unwrap();
        let mapped = model
            .map_start(&case, &model, &case, &solution.solver)
            .unwrap();
        assert_eq!(
            mapped,
            WarmStart::from_solution(&solution.solver),
            "{name} AC"
        );
    }
}

#[test]
fn an_outage_loses_exactly_its_own_variables_and_constraints() {
    // case30: a limited branch's outage loses its flow rows (AC: Sf and St; DC: one); a
    //   generator's loses its outputs (AC: Pg and Qg; DC: Pg) and their bounds' multipliers.
    //   case30pwl: a generator with a piecewise cost loses its helper y too, and its segments'
    //   rows. Mapped back, everything else is unchanged, bit for bit.
    for name in ["case30", "case30pwl"] {
        let Some(base) = case(name) else { continue };
        let generator = 1;
        let segments = match &base.gencost.as_ref().expect("costs")[generator] {
            GenCost::PiecewiseLinear { points, .. } => points.len() - 1,
            GenCost::Polynomial { .. } => 0,
        };
        let helpers = usize::from(segments > 0);
        let branch = branch_outage(&base, true);
        // (case, DC lost entries, AC lost entries): a column counts in x, lower and upper; a
        //   row of A in mu_l and mu_u; a nonlinear inequality (AC's Sf, St) once
        let outages = [
            (without_branch(&base, branch), 2, 2),
            (
                without_generator(&base, generator),
                3 * (1 + helpers) + 2 * segments,
                3 * (2 + helpers) + 2 * segments,
            ),
        ];
        let dc_base = DcOpfModel::new(&base).unwrap();
        let dc_solution = dc_base.solve(&base, &DcOpfOptions::default()).unwrap();
        let ac_base = AcOpfModel::new(&base).unwrap();
        let ac_solution = ac_base.solve(&base, &AcOpfOptions::mehrotra()).unwrap();
        for (k, (outage, dc_lost, ac_lost)) in outages.iter().enumerate() {
            let label = format!("{name}, outage {k}");
            let model = DcOpfModel::new(outage).unwrap();
            let warm = model
                .map_start(outage, &dc_base, &base, &dc_solution.solver)
                .unwrap();
            let back = dc_base
                .map_start(&base, &model, outage, &as_solution(&warm))
                .unwrap();
            let own = WarmStart::from_solution(&dc_solution.solver);
            let (changed, lost) = changes(&own, &warm, &back);
            assert_eq!(lost, *dc_lost, "{label} DC: lost");
            assert!(
                changed <= lost,
                "{label} DC: {changed} changed, {lost} lost"
            );
            let model = AcOpfModel::new(outage).unwrap();
            let warm = model
                .map_start(outage, &ac_base, &base, &ac_solution.solver)
                .unwrap();
            let back = ac_base
                .map_start(&base, &model, outage, &as_solution(&warm))
                .unwrap();
            let own = WarmStart::from_solution(&ac_solution.solver);
            let (changed, lost) = changes(&own, &warm, &back);
            assert_eq!(lost, *ac_lost, "{label} AC: lost");
            assert!(
                changed <= lost,
                "{label} AC: {changed} changed, {lost} lost"
            );
        }
    }
}

/// O3's bar without LAM_Q, against the cold solve, as PR C's survey reads it (the cold solve is
///   MIPS-LU's here, so it's also the residuals' reference)
const SHIFT: f64 = 1e-2;

fn shifted(solver: OpfSolver) -> OpfSolver {
    match solver {
        OpfSolver::Mips(mut options) => {
            options.warm_shift = SHIFT;
            OpfSolver::Mips(options)
        }
        OpfSolver::Mehrotra(mut options) => {
            options.warm_shift = SHIFT;
            OpfSolver::Mehrotra(options)
        }
        other => other,
    }
}

#[test]
fn a_mapped_warm_start_solves_an_outage() {
    let Some(base) = case("case118") else { return };
    let outages = [
        ("branch", without_branch(&base, branch_outage(&base, false))),
        ("generator", without_generator(&base, 1)),
    ];
    let mut dc_options = DcOpfOptions::default();
    dc_options.solver = shifted(dc_options.solver);
    let mut ac_options = AcOpfOptions::default();
    ac_options.solver = shifted(ac_options.solver);
    let dc_base = DcOpfModel::new(&base).unwrap();
    let dc_previous = dc_base.solve(&base, &dc_options).unwrap();
    let ac_base = AcOpfModel::new(&base).unwrap();
    let ac_previous = ac_base.solve(&base, &ac_options).unwrap();
    for (kind, outage) in &outages {
        let model = DcOpfModel::new(outage).unwrap();
        let mapped = model
            .map_start(outage, &dc_base, &base, &dc_previous.solver)
            .unwrap();
        let cold = model.solve(outage, &dc_options).unwrap();
        let warm = model
            .solve_from(outage, &dc_options, OpfStart::Mapped(&mapped))
            .unwrap();
        let primal = model
            .solve_from(outage, &dc_options, OpfStart::MappedPrimal(&mapped))
            .unwrap();
        eprintln!(
            "{kind} DC: cold {}, warm {}, primal {} ({:?})",
            cold.iterations, warm.iterations, primal.iterations, primal.status
        );
        assert!(cold.converged() && warm.converged(), "{kind} DC");
        let objective = (warm.objective - cold.objective).abs() / cold.objective.abs().max(1.0);
        assert!(objective <= 1e-6, "{kind} DC: objective {objective:e}");
        assert!(
            relative_difference(&warm.lam_p, &cold.lam_p) <= 1e-3,
            "{kind} DC: LAM_P"
        );
        assert!(
            warm.iterations < cold.iterations,
            "{kind} DC: no fewer iterations"
        );

        let model = AcOpfModel::new(outage).unwrap();
        let mapped = model
            .map_start(outage, &ac_base, &base, &ac_previous.solver)
            .unwrap();
        let cold = model.solve(outage, &ac_options).unwrap();
        let warm = model
            .solve_from(outage, &ac_options, OpfStart::Mapped(&mapped))
            .unwrap();
        eprintln!(
            "{kind} AC: cold {}, warm {}",
            cold.iterations, warm.iterations
        );
        assert!(cold.converged() && warm.converged(), "{kind} AC");
        let objective = (warm.objective - cold.objective).abs() / cold.objective.abs().max(1.0);
        assert!(objective <= 1e-6, "{kind} AC: objective {objective:e}");
        assert!(
            relative_difference(&warm.lam_p, &cold.lam_p) <= 1e-3,
            "{kind} AC: LAM_P"
        );
        let (a, b) = (warm.residuals, cold.residuals);
        for (what, a, b) in [
            ("stationarity", a.stationarity, b.stationarity),
            ("feasibility", a.feasibility, b.feasibility),
            ("complementarity", a.complementarity, b.complementarity),
        ] {
            assert!(
                a <= 100.0 * b || a <= 1e-6,
                "{kind} AC: {what} {a:e} (cold {b:e})"
            );
        }
        assert!(
            warm.iterations < cold.iterations,
            "{kind} AC: no fewer iterations"
        );
    }
}

#[test]
fn a_case_with_other_tables_is_refused() {
    let Some(base) = case("case9") else { return };
    let model = DcOpfModel::new(&base).unwrap();
    let solution = model.solve(&base, &DcOpfOptions::default()).unwrap();
    let mut longer = base.clone();
    longer.branch.push(longer.branch[0].clone());
    let other = DcOpfModel::new(&longer).unwrap();
    assert!(matches!(
        other.map_start(&longer, &model, &base, &solution.solver),
        Err(PowerflowError::CaseMismatch { .. })
    ));
    let mut short = WarmStart::from_solution(&solution.solver);
    short.x.pop();
    assert!(matches!(
        model.map_start(&base, &model, &base, &as_solution(&short)),
        Err(PowerflowError::Solver { .. })
    ));
}

/// The cost of a piecewise linear cost's `points` (MW, $/hr) at `p` MW
fn piecewise_cost(points: &[(f64, f64)], p: f64) -> f64 {
    let j = (1..points.len())
        .find(|&j| p <= points[j].0)
        .unwrap_or(points.len() - 1);
    let ((p0, c0), (p1, c1)) = (points[j - 1], points[j]);
    c0 + (c1 - c0) * (p - p0) / (p1 - p0)
}

/// `a` equals `b` to a relative 1e-9, as a start differs from the solution it maps only by the
///   negligible shift
fn same(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// MIPS stopped at its start (`max_it` 0), with a negligible shift: its reported results are
///   the mapped start's, by case row
fn at_start(mut solver: OpfSolver) -> OpfSolver {
    if let OpfSolver::Mips(ref mut options) = solver {
        (options.warm_shift, options.max_it) = (1e-12, 0);
    }
    solver
}

#[test]
fn the_mapped_values_land_on_their_own_rows() {
    // Reported by case row, a mapped start's angles, outputs, prices and flow multipliers are
    //   the base solution's, row by row, and for case30pwl its objective is the base's less the
    //   outaged generator's cost: so each key names its own bus, generator and branch, which a
    //   round trip can't show (a key misnamed both ways maps back to itself)
    for name in ["case30", "case30pwl"] {
        let Some(base) = case(name) else { continue };
        let (branch, generator) = (branch_outage(&base, true), 1);
        // For DC, a flow limit binding after the outaged branch, where internal indices shift:
        //   the heaviest later branch's limit set below its DC flow. AC's limits bind already,
        //   and its apparent power exceeds this.
        let mut dc_case = base.clone();
        let flows = DcOpfModel::new(&base)
            .unwrap()
            .solve(&base, &DcOpfOptions::default())
            .unwrap()
            .pf_mw;
        let heaviest = (branch + 1..base.branch.len())
            .filter(|&k| base.branch[k].status && base.branch[k].rate_a_mva != 0.0)
            .max_by(|&a, &b| flows[a].abs().total_cmp(&flows[b].abs()))
            .expect("a limited branch after the outaged one");
        dc_case.branch[heaviest].rate_a_mva = 0.95 * flows[heaviest].abs();
        let mut dc_options = DcOpfOptions::default();
        dc_options.solver = at_start(dc_options.solver);
        let mut ac_options = AcOpfOptions::default();
        ac_options.solver = at_start(ac_options.solver);
        let dc_base = DcOpfModel::new(&dc_case).unwrap();
        let dc = dc_base.solve(&dc_case, &DcOpfOptions::default()).unwrap();
        let ac_base = AcOpfModel::new(&base).unwrap();
        let ac = ac_base.solve(&base, &AcOpfOptions::default()).unwrap();
        assert!(
            dc.converged() && ac.converged(),
            "{name}: DC {:?}, AC {:?}",
            dc.status,
            ac.status
        );
        for (kind, outage, dc_outage) in [
            (
                "branch",
                without_branch(&base, branch),
                without_branch(&dc_case, branch),
            ),
            (
                "generator",
                without_generator(&base, generator),
                without_generator(&dc_case, generator),
            ),
        ] {
            let label = format!("{name}, {kind}");
            let gens = (0..base.generator.len()).filter(|&g| outage.generator[g].status);
            let branches: Vec<usize> = (0..base.branch.len())
                .filter(|&k| outage.branch[k].status)
                .collect();

            let model = DcOpfModel::new(&dc_outage).unwrap();
            let mapped = model
                .map_start(&dc_outage, &dc_base, &dc_case, &dc.solver)
                .unwrap();
            let start = model
                .solve_from(&dc_outage, &dc_options, OpfStart::Mapped(&mapped))
                .unwrap();
            for g in gens.clone() {
                assert!(same(start.pg_mw[g], dc.pg_mw[g]), "{label} DC: Pg of {g}");
            }
            for i in 0..base.bus.len() {
                assert!(same(start.va_deg[i], dc.va_deg[i]), "{label} DC: Va of {i}");
                assert!(
                    same(start.lam_p[i], dc.lam_p[i]),
                    "{label} DC: LAM_P of {i}"
                );
            }
            let mut compared = 0;
            for &k in &branches {
                if dc.mu_sf[k] > 1e-5 || dc.mu_st[k] > 1e-5 {
                    compared += 1;
                    assert!(
                        same(start.mu_sf[k], dc.mu_sf[k]),
                        "{label} DC: mu_sf of {k}"
                    );
                    assert!(
                        same(start.mu_st[k], dc.mu_st[k]),
                        "{label} DC: mu_st of {k}"
                    );
                }
            }
            assert!(compared > 0, "{label} DC: no binding flow limit to compare");

            let model = AcOpfModel::new(&outage).unwrap();
            let mapped = model
                .map_start(&outage, &ac_base, &base, &ac.solver)
                .unwrap();
            let start = model
                .solve_from(&outage, &ac_options, OpfStart::Mapped(&mapped))
                .unwrap();
            for g in gens {
                assert!(same(start.pg_mw[g], ac.pg_mw[g]), "{label} AC: Pg of {g}");
                assert!(
                    same(start.qg_mvar[g], ac.qg_mvar[g]),
                    "{label} AC: Qg of {g}"
                );
            }
            for i in 0..base.bus.len() {
                assert!(same(start.vm_pu[i], ac.vm_pu[i]), "{label} AC: Vm of {i}");
                assert!(
                    same(start.lam_p[i], ac.lam_p[i]),
                    "{label} AC: LAM_P of {i}"
                );
                assert!(
                    same(start.lam_q[i], ac.lam_q[i]),
                    "{label} AC: LAM_Q of {i}"
                );
            }
            let mut compared = 0;
            for &k in &branches {
                if ac.mu_sf[k] > 1e-5 || ac.mu_st[k] > 1e-5 {
                    compared += 1;
                    assert!(
                        same(start.mu_sf[k], ac.mu_sf[k]),
                        "{label} AC: mu_sf of {k}"
                    );
                    assert!(
                        same(start.mu_st[k], ac.mu_st[k]),
                        "{label} AC: mu_st of {k}"
                    );
                }
            }
            assert!(compared > 0, "{label} AC: no binding flow limit to compare");
            if let (Some(GenCost::PiecewiseLinear { points, .. }), "generator") =
                (base.gencost.as_ref().map(|costs| &costs[generator]), kind)
            {
                // To the solution's accuracy, at which y meets its cost: a misnamed helper would
                //   move it by a whole generator's cost
                let lost = piecewise_cost(points, ac.pg_mw[generator]);
                let expected = ac.objective - lost;
                assert!(
                    (start.objective - expected).abs() <= 1e-6 * expected.abs(),
                    "{label} AC: objective {} against {} − {lost}",
                    start.objective,
                    ac.objective
                );
            }
        }
    }
}
