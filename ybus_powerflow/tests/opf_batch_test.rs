// Batches of optimal power flows (opf::batch, solve_many): each batch bit for bit the loop it
//   replaces, cold and warm, on both models and both solvers; an outage of a branch and a
//   generator together; an outage whose case can't be solved, alone in its error; and a warm
//   batch whose base fails, solved cold. With the rayon feature, each parallel form bit for bit
//   its serial one on pools of 1, 2, 7 and 28 threads.

mod common;

use ybus_powerflow::case::Case;
use ybus_powerflow::islands::find_bridges;
use ybus_powerflow::opf::ac::{AcOpfModel, AcOpfOptions, AcOpfSolution};
use ybus_powerflow::opf::batch::{
    ac_opf_outages, dc_opf_outages, BatchStart, OpfOutage, OpfOutages,
};
use ybus_powerflow::opf::dc::{DcOpfModel, DcOpfOptions, DcOpfSolution};
use ybus_powerflow::opf::OpfStart;
use ybus_powerflow::powerflow_error::PowerflowError;

fn case(name: &str) -> Option<Case> {
    common::matpower_case(name).map(|path| Case::from_path(&path).unwrap())
}

/// Every result, bit for bit, error or solution
fn bits<S: std::fmt::Debug>(result: &Result<S, PowerflowError>) -> String {
    format!("{result:?}")
}

/// A few outages of each kind: two branches that aren't bridges, two generators, one N-2 of
///   a branch and a generator, and, if the case has one, a bridge, whose outage islands buses
fn outages(case: &Case) -> Vec<OpfOutage> {
    let bridges = find_bridges(case).unwrap();
    let branches: Vec<usize> = (0..case.branch.len())
        .filter(|&row| case.branch[row].status && bridges.binary_search(&row).is_err())
        .take(2)
        .collect();
    let mut list: Vec<OpfOutage> = branches.iter().map(|&row| OpfOutage::branch(row)).collect();
    list.push(OpfOutage::generator(0));
    list.push(OpfOutage::generator(1));
    list.push(OpfOutage::new(vec![branches[0]], vec![1]));
    if let Some(&bridge) = bridges.first() {
        list.push(OpfOutage::branch(bridge));
    }
    list
}

/// What a batch replaces: the base solved cold, and each outage on its own model, cold or
///   from the base solution mapped onto it, if it converged
fn ac_loop(
    base: &Case,
    list: &[OpfOutage],
    options: &AcOpfOptions,
    start: BatchStart,
) -> (
    AcOpfSolution,
    bool,
    Vec<Result<AcOpfSolution, PowerflowError>>,
) {
    let model = AcOpfModel::new(base).unwrap();
    let solved = model.solve(base, options).unwrap();
    let warm = start == BatchStart::Warm && solved.converged();
    let results = list
        .iter()
        .map(|outage| {
            let changed = outage.apply(base)?;
            let own = AcOpfModel::new(&changed)?;
            if warm {
                let mapped = own.map_start(&changed, &model, base, &solved.solver)?;
                own.solve_from(&changed, options, OpfStart::Mapped(&mapped))
            } else {
                own.solve(&changed, options)
            }
        })
        .collect();
    (solved, warm, results)
}

fn dc_loop(
    base: &Case,
    list: &[OpfOutage],
    options: &DcOpfOptions,
    start: BatchStart,
) -> (
    DcOpfSolution,
    bool,
    Vec<Result<DcOpfSolution, PowerflowError>>,
) {
    let model = DcOpfModel::new(base).unwrap();
    let solved = model.solve(base, options).unwrap();
    let warm = start == BatchStart::Warm && solved.converged();
    let results = list
        .iter()
        .map(|outage| {
            let changed = outage.apply(base)?;
            let own = DcOpfModel::new(&changed)?;
            if warm {
                let mapped = own.map_start(&changed, &model, base, &solved.solver)?;
                own.solve_from(&changed, options, OpfStart::Mapped(&mapped))
            } else {
                own.solve(&changed, options)
            }
        })
        .collect();
    (solved, warm, results)
}

fn same<S: std::fmt::Debug>(
    label: &str,
    batch: &OpfOutages<S>,
    (base, warm, results): &(S, bool, Vec<Result<S, PowerflowError>>),
) {
    assert_eq!(
        format!("{:?}", batch.base),
        format!("{base:?}"),
        "{label}: the base"
    );
    assert_eq!(batch.warm, *warm, "{label}: warm");
    assert_eq!(batch.outages.len(), results.len(), "{label}: the outages");
    for (k, (a, b)) in batch.outages.iter().zip(results).enumerate() {
        assert_eq!(bits(a), bits(b), "{label}: outage {k}");
    }
}

const STARTS: [BatchStart; 2] = [BatchStart::Cold, BatchStart::Warm];

#[test]
fn a_batch_of_outages_is_the_loop_it_replaces() {
    for name in ["case30", "case118"] {
        let Some(base) = case(name) else { continue };
        let list = outages(&base);
        for start in STARTS {
            for (solver, options) in [
                ("MIPS", DcOpfOptions::default()),
                ("predictor-corrector", DcOpfOptions::mehrotra()),
            ] {
                let label = format!("{name} DC {solver} {start:?}");
                let batch = dc_opf_outages(&base, &list, &options, start).unwrap();
                same(&label, &batch, &dc_loop(&base, &list, &options, start));
                assert_eq!(batch.warm, start == BatchStart::Warm, "{label}");
            }
            for (solver, options) in [
                ("MIPS", AcOpfOptions::default()),
                ("predictor-corrector", AcOpfOptions::mehrotra()),
            ] {
                let label = format!("{name} AC {solver} {start:?}");
                let batch = ac_opf_outages(&base, &list, &options, start).unwrap();
                same(&label, &batch, &ac_loop(&base, &list, &options, start));
                assert_eq!(batch.warm, start == BatchStart::Warm, "{label}");
            }
        }
    }
}

#[test]
fn an_outage_applies_its_rows_together_and_refuses_one_out_of_range() {
    let Some(base) = case("case30") else { return };
    let changed = OpfOutage::new(vec![3], vec![1]).apply(&base).unwrap();
    assert!(!changed.branch[3].status && !changed.generator[1].status);
    let unchanged = |k: usize| changed.branch[k] == base.branch[k];
    assert!((0..base.branch.len()).filter(|&k| k != 3).all(unchanged));
    assert!(matches!(
        OpfOutage::generator(base.generator.len()).apply(&base),
        Err(PowerflowError::RowOutOfRange {
            table: "generator",
            ..
        })
    ));
    // An outage out of range is that outage's error, the others solved
    let list = [OpfOutage::branch(3), OpfOutage::branch(base.branch.len())];
    let batch = dc_opf_outages(&base, &list, &DcOpfOptions::default(), BatchStart::Warm).unwrap();
    assert!(batch.outages[0].is_ok());
    assert!(matches!(
        batch.outages[1],
        Err(PowerflowError::RowOutOfRange {
            table: "branch",
            ..
        })
    ));
}

#[test]
fn a_warm_batch_whose_base_fails_solves_its_outages_cold() {
    // No generator can serve the load, so the base can't converge
    let Some(mut base) = case("case30") else {
        return;
    };
    for generator in &mut base.generator {
        generator.pmax_mw = 0.0;
    }
    let list = [OpfOutage::branch(3)];
    let options = DcOpfOptions::default();
    let batch = dc_opf_outages(&base, &list, &options, BatchStart::Warm).unwrap();
    assert!(!batch.base.converged());
    assert!(!batch.warm);
    same(
        "case30 DC, failing base",
        &batch,
        &dc_loop(&base, &list, &options, BatchStart::Cold),
    );
}

#[test]
fn solve_many_is_solve_from_for_each_case() {
    let Some(base) = case("case118") else { return };
    let cases: Vec<Case> = [0.95, 1.0, 1.03]
        .iter()
        .map(|&factor| {
            let mut case = base.clone();
            for bus in &mut case.bus {
                bus.pd_mw *= factor;
                bus.qd_mvar *= factor;
            }
            case
        })
        .collect();
    let model = AcOpfModel::new(&base).unwrap();
    let options = AcOpfOptions::mehrotra();
    let previous = model.solve(&base, &options).unwrap();
    for start in [OpfStart::Cold, OpfStart::Warm(&previous.solver)] {
        let many = model.solve_many(&cases, &options, start);
        for (case, result) in cases.iter().zip(&many) {
            assert_eq!(bits(result), bits(&model.solve_from(case, &options, start)));
        }
    }
    let model = DcOpfModel::new(&base).unwrap();
    let options = DcOpfOptions::default();
    let previous = model.solve(&base, &options).unwrap();
    for start in [OpfStart::Cold, OpfStart::Warm(&previous.solver)] {
        let many = model.solve_many(&cases, &options, start);
        for (case, result) in cases.iter().zip(&many) {
            assert_eq!(bits(result), bits(&model.solve_from(case, &options, start)));
        }
    }
}
