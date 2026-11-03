//! Times every single-branch outage of one case by rank-k updates (`DcStudy::solve_outages`),
//! and compares a sample with refactoring, for the N-1 survey of every MATPOWER case that
//! `scripts/benchmark-n1.sh` runs. The script runs this once per case, each in its own process,
//! so that it can limit time and memory and record peak memory.
//!
//! Usage: `n1_survey <case.m>`, or `n1_survey --header`.
//!
//! Prints one CSV line, with the columns in `HEADER`. The outages run in chunks, so the results
//! of a case's thousands of outages are never all held at once; the whole pass is timed by its
//! median over repeated runs, for up to `BUDGET` or `MAX_RUNS` runs. Then a sample of outages
//! (every one on cases of up to `FULL_COMPARISON_BUSES` buses, otherwise every m-th, about
//! `SAMPLE` of them) is each solved again and compared with refactoring: a new model of the
//! case with that branch out of service, solved afresh, which is also timed.

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use powerflow::case::Case;
use powerflow::dc::{dc_power_flow, DcModel, DcSolution, DcSolver};
use powerflow::powerflow_error::PowerflowError;

/// `study_ms` is factoring B′ once. `outages` counts the in-service branches, each taken out
/// alone: `solved` of them solve and `islanding` leave a bus without a reference.
/// `woodbury_ms` is the whole pass, and `woodbury_us_per_outage` its mean per outage.
/// `compared` outages were compared with refactoring: `mismatches` counts those where one
/// solves and the other doesn't, or they report different errors, and
/// `worst_relative_deviation` is the largest difference in angles, flows or generator outputs
/// relative to the largest magnitude of that quantity in the refactored result.
/// `refactor_ms_per_outage` is the median time to build and solve a new model, and `speedup`
/// its ratio to `woodbury_us_per_outage`.
const HEADER: &str = "case,buses,branches,connected_buses,in_service_branches,solver,study_ms,\
                      outages,solved,islanding,woodbury_ms,woodbury_runs,woodbury_us_per_outage,\
                      compared,mismatches,worst_relative_deviation,refactor_ms_per_outage,speedup";

/// How long to keep repeating the pass for its median time
const BUDGET: Duration = Duration::from_secs(1);
/// The most runs of the pass
const MAX_RUNS: usize = 101;
/// Outages per call to `solve_outages`
const CHUNK: usize = 256;
/// Cases up to this many buses compare every outage with refactoring
const FULL_COMPARISON_BUSES: usize = 3000;
/// About how many outages larger cases compare
const SAMPLE: usize = 300;

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.as_slice() {
        [flag] if flag == "--header" => {
            println!("{HEADER}");
            ExitCode::SUCCESS
        }
        [case] => match survey(Path::new(case)) {
            Ok(line) => {
                println!("{line}");
                ExitCode::SUCCESS
            }
            Err(err) => {
                eprintln!("{case}: {err}");
                ExitCode::FAILURE
            }
        },
        _ => {
            eprintln!("usage: n1_survey <case.m> | --header");
            ExitCode::from(2)
        }
    }
}

fn survey(path: &Path) -> Result<String, String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let case = Case::from_path(path).map_err(|err| err.to_string())?;
    let model = DcModel::new(&case).map_err(|err| err.to_string())?;
    let (study_time, _, study) = time(|| model.study());
    let study = study.map_err(|err| err.to_string())?;
    let rows = model.branch_rows().to_vec();

    // Every single outage, in chunks, counting the outcomes
    let pass = || -> Result<(usize, usize), String> {
        let (mut solved, mut islanding) = (0, 0);
        for chunk in rows.chunks(CHUNK) {
            let sets: Vec<Vec<usize>> = chunk.iter().map(|&row| vec![row]).collect();
            for result in study
                .solve_outages(&case, &sets)
                .map_err(|err| err.to_string())?
            {
                match result {
                    Ok(_) => solved += 1,
                    Err(PowerflowError::IslandWithoutReference { .. }) => islanding += 1,
                    Err(err) => return Err(err.to_string()),
                }
            }
        }
        Ok((solved, islanding))
    };
    let (woodbury_time, woodbury_runs, counts) = time(pass);
    let (solved, islanding) = counts?;

    // The sample, against refactoring
    let step = if model.bus_rows().len() <= FULL_COMPARISON_BUSES {
        1
    } else {
        rows.len().div_ceil(SAMPLE).max(1)
    };
    let (mut compared, mut mismatches, mut worst) = (0, 0, 0.0_f64);
    let mut refactor_times = vec![];
    for &row in rows.iter().step_by(step) {
        let ours = study.solve_outage(&case, &[row]);
        let mut outaged = case.clone();
        outaged.branch[row].status = false;
        let start = Instant::now();
        let expected = dc_power_flow(&outaged);
        refactor_times.push(start.elapsed());
        compared += 1;
        match (ours, expected) {
            (Ok(ours), Ok(expected)) => worst = worst.max(relative_deviation(&ours, &expected)),
            (Err(ours), Err(expected)) if ours == expected => {}
            _ => mismatches += 1,
        }
    }
    refactor_times.sort();
    let refactor = refactor_times
        .get(refactor_times.len() / 2)
        .copied()
        .unwrap_or_default();

    let per_outage_us = woodbury_time.as_secs_f64() * 1e6 / rows.len().max(1) as f64;
    let speedup = refactor.as_secs_f64() * 1e6 / per_outage_us;
    Ok(format!(
        "{name},{},{},{},{},{},{:.3},{},{solved},{islanding},{:.3},{woodbury_runs},{:.3},\
         {compared},{mismatches},{worst:e},{:.3},{:.1}",
        case.bus.len(),
        case.branch.len(),
        model.bus_rows().len(),
        rows.len(),
        solver_name(study.solver()),
        milliseconds(study_time),
        rows.len(),
        milliseconds(woodbury_time),
        per_outage_us,
        milliseconds(refactor),
        speedup,
    ))
}

/// The largest difference between two solutions, relative to the largest magnitude of each
/// quantity in `expected` (angles, flows, generator outputs), at least 1.
fn relative_deviation(actual: &DcSolution, expected: &DcSolution) -> f64 {
    let worst = |a: &[f64], e: &[f64]| {
        let scale = e.iter().fold(1.0_f64, |m, v| m.max(v.abs()));
        a.iter()
            .zip(e)
            .fold(0.0_f64, |m, (a, e)| m.max((a - e).abs()))
            / scale
    };
    worst(&actual.va_deg, &expected.va_deg)
        .max(worst(&actual.pf_mw, &expected.pf_mw))
        .max(worst(&actual.pg_mw, &expected.pg_mw))
}

fn solver_name(solver: DcSolver) -> &'static str {
    match solver {
        DcSolver::Cholesky => "Cholesky",
        DcSolver::Lu => "Lu",
        DcSolver::NoSolve => "none",
    }
}

/// Runs `step` until `BUDGET` has passed or it has run `MAX_RUNS` times, and returns its median
/// time (the upper one, for an even number of runs) and number of runs, with the last run's
/// result.
fn time<T>(mut step: impl FnMut() -> T) -> (Duration, usize, T) {
    let start = Instant::now();
    let mut times = vec![];
    loop {
        let run = Instant::now();
        let result = step();
        times.push(run.elapsed());
        if start.elapsed() >= BUDGET || times.len() >= MAX_RUNS {
            times.sort();
            return (times[times.len() / 2], times.len(), result);
        }
    }
}

fn milliseconds(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1e3
}
