//! Times the AC power flow on one case, for the survey of every MATPOWER case that
//! `scripts/benchmark-ac.sh` runs. The script runs this once per case and method, each in its own
//! process, so that it can limit time and memory and record peak memory.
//!
//! Usage: `ac_survey [--method nr|fdxb|fdbx] <case.m> [<reference-dir>]`, or
//! `ac_survey --header`. The method defaults to `nr` (Newton's method, as `ac_power_flow` uses).
//!
//! Prints one CSV line, with the columns in `HEADER`. Parsing, building the model and solving are
//! each timed by their median over repeated runs, for up to `BUDGET` or `MAX_RUNS` runs. A solve
//! that doesn't converge, or fails, is still timed, and its `outcome` says why. Given a directory
//! of MATPOWER's results from `scripts/generate-reference.sh`, it also prints MATPOWER's
//! iteration count and the largest differences from MATPOWER's results; those columns are empty
//! when the directory has no results for the case and method (MATPOWER's own run failed, or
//! wasn't generated).

#[allow(dead_code)]
#[path = "../tests/common/reference.rs"]
mod reference;

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use powerflow::ac::{AcMethod, AcModel, AcOptions};
use powerflow::case::Case;
use powerflow::powerflow_error::PowerflowError;

/// `outcome` is `converged`, `no convergence` or the error (the survey script's own `status`
/// column says whether the process finished). `iterations` and `mismatch_pu` are
/// those of the converged solution, or of the last iteration when it didn't converge.
/// `matpower_iterations` is `runpf`'s count. The deviations are the largest absolute differences
/// from MATPOWER's results: voltage magnitudes (p.u.), angles (degrees), branch flows (the
/// largest of P and Q at either end, MW or MVAr) and generator outputs (P or Q).
const HEADER: &str = "case,buses,branches,generators,connected_buses,in_service_branches,\
                      method,parse_ms,model_ms,solve_ms,solve_runs,outcome,iterations,\
                      mismatch_pu,full_factorizations,matpower_iterations,\
                      vm_deviation_pu,va_deviation_deg,flow_deviation,gen_deviation";

/// How long to keep repeating a step for its median time
const BUDGET: Duration = Duration::from_secs(1);
/// The most runs of a step
const MAX_RUNS: usize = 101;

fn main() -> ExitCode {
    let mut args: Vec<String> = env::args().skip(1).collect();
    if args.as_slice() == ["--header"] {
        println!("{HEADER}");
        return ExitCode::SUCCESS;
    }
    let mut method = AcMethod::Newton;
    if args.first().is_some_and(|flag| flag == "--method") {
        method = match args.get(1).map(String::as_str) {
            Some("nr") => AcMethod::Newton,
            Some("fdxb") => AcMethod::FastDecoupledXb,
            Some("fdbx") => AcMethod::FastDecoupledBx,
            _ => return usage(),
        };
        args.drain(..2);
    }
    let (case, references) = match args.as_slice() {
        [case] => (case, None),
        [case, references] => (case, Some(Path::new(references))),
        _ => return usage(),
    };
    match survey(Path::new(case), references, method) {
        Ok(line) => {
            println!("{line}");
            ExitCode::SUCCESS
        }
        Err(err) => {
            eprintln!("{case}: {err}");
            ExitCode::FAILURE
        }
    }
}

fn usage() -> ExitCode {
    eprintln!("usage: ac_survey [--method nr|fdxb|fdbx] <case.m> [<reference-dir>] | --header");
    ExitCode::from(2)
}

/// The reference files' name for `method`'s results, and the survey's name for it.
fn names(method: AcMethod) -> (&'static str, &'static str) {
    match method {
        AcMethod::Newton => ("acpf", "nr"),
        AcMethod::FastDecoupledXb => ("fdxb", "fdxb"),
        AcMethod::FastDecoupledBx => ("fdbx", "fdbx"),
        _ => ("", "other"),
    }
}

/// The CSV line for the case at `path`. Errors only if the case can't be read or modelled;
/// a solve that fails is reported in the line.
fn survey(path: &Path, references: Option<&Path>, method: AcMethod) -> Result<String, String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let (parse, _, case) = time(|| Case::from_path(path));
    let case = case.map_err(|err| err.to_string())?;
    let (model_time, _, model) = time(|| AcModel::new(&case));
    let model = model.map_err(|err| err.to_string())?;
    let options = AcOptions::for_method(method);
    let (solve, solve_runs, solution) = time(|| model.solve(&case, options));

    let (outcome, iterations, mismatch, factorizations) = match &solution {
        Ok(solution) => (
            "converged".to_string(),
            solution.iterations.to_string(),
            format!("{:e}", solution.mismatch_pu),
            solution.full_factorizations.to_string(),
        ),
        Err(PowerflowError::NoConvergence {
            iterations,
            mismatch_pu,
        }) => (
            "no convergence".to_string(),
            iterations.to_string(),
            format!("{mismatch_pu:e}"),
            String::new(),
        ),
        Err(err) => (
            format!("error: {}", err.to_string().replace(',', ";")),
            String::new(),
            String::new(),
            String::new(),
        ),
    };

    let (run, method_name) = names(method);
    let comparison = match references {
        Some(dir) if dir.join(format!("{name}_{run}_summary.csv")).is_file() => {
            let reference = reference::read_ac_run(dir, name, run)?;
            match &solution {
                Ok(solution) => {
                    let d = reference::ac_deviation(&case, solution, &reference)?;
                    format!(
                        "{},{:e},{:e},{:e},{:e}",
                        reference.iterations,
                        d.vm_pu,
                        d.va_deg,
                        d.pf_mw.max(d.qf_mvar).max(d.pt_mw).max(d.qt_mvar),
                        d.pg_mw.max(d.qg_mvar)
                    )
                }
                Err(_) => format!("{},,,,", reference.iterations),
            }
        }
        _ => ",,,,".to_string(),
    };
    Ok(format!(
        "{name},{},{},{},{},{},{method_name},{:.3},{:.3},{:.3},{solve_runs},{outcome},\
         {iterations},{mismatch},{factorizations},{comparison}",
        case.bus.len(),
        case.branch.len(),
        case.generator.len(),
        model.bus_rows().len(),
        model.branch_rows().len(),
        milliseconds(parse),
        milliseconds(model_time),
        milliseconds(solve),
    ))
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
