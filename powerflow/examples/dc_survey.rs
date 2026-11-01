//! Times the DC power flow on one case, for the survey of every MATPOWER case that
//! `scripts/benchmark-dc.sh` runs. The script runs this once per case, each in its own process,
//! so that it can limit time and memory and record peak memory.
//!
//! Usage: `dc_survey [--ordering natural|rcm|amd] <case.m> [<reference-dir>]`, or
//! `dc_survey --header`. The ordering is the factorization's, Cholesky or LU (default `amd`, as
//! `DcModel::solve` uses).
//!
//! Prints one CSV line, with the columns in `HEADER`. Each step (parse, model, solve) is timed
//! by its median over repeated runs, for up to `BUDGET` or `MAX_RUNS` runs, so small cases get
//! stable figures and large ones run once. After timing, it counts the factor's entries, which
//! unlike the times don't depend on the machine (from the symbolic analysis for Cholesky, by
//! factoring again for LU), and it measures the solve's backward error. Given a directory of MATPOWER's results from
//! `scripts/generate-reference.sh`, it also prints the backward error of MATPOWER's angles and
//! the largest differences from MATPOWER's results; those columns are empty when the directory
//! has no results for the case.

#[path = "../tests/common/accuracy.rs"]
mod accuracy;
#[allow(dead_code)]
#[path = "../tests/common/reference.rs"]
mod reference;

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use powerflow::case::Case;
use powerflow::dc::{DcModel, DcSolver};
use sparse_mat::compressed_matrix::{LuOptions, Ordering};

/// `ordering` is the Cholesky factorization's. `negative_branches` counts the in-service branches with a negative series susceptance,
/// which can make B′ indefinite. `bprime_dimension` and `bprime_entries` give B′'s size and
/// stored entries, and `factor_entries` counts those of the factor the solve used: L for
/// Cholesky, or L and U for LU, counting their diagonal once. `backward_error` is the
/// normwise backward error of the solve's angles (`tests/common/accuracy.rs`), and
/// `matpower_backward_error` that of MATPOWER's angles for the same system: a control, since
/// MATPOWER's solve reorders B′ to limit fill-in.
const HEADER: &str = "case,buses,branches,generators,connected_buses,in_service_branches,\
                      negative_branches,parse_ms,model_ms,solve_ms,solve_runs,solver,ordering,\
                      bprime_dimension,bprime_entries,factor_entries,\
                      backward_error,matpower_backward_error,\
                      va_deviation_deg,pf_deviation_mw,pg_deviation_mw";

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
    let mut ordering = Ordering::Amd;
    if args.first().is_some_and(|flag| flag == "--ordering") {
        ordering = match args.get(1).map(String::as_str) {
            Some("natural") => Ordering::Natural,
            Some("rcm") => Ordering::Rcm,
            Some("amd") => Ordering::Amd,
            _ => return usage(),
        };
        args.drain(..2);
    }
    let (case, references) = match args.as_slice() {
        [case] => (case, None),
        [case, references] => (case, Some(Path::new(references))),
        _ => return usage(),
    };
    match survey(Path::new(case), references, ordering) {
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
    eprintln!(
        "usage: dc_survey [--ordering natural|rcm|amd] <case.m> [<reference-dir>] | --header"
    );
    ExitCode::from(2)
}

/// The CSV line for the case at `path`.
fn survey(path: &Path, references: Option<&Path>, ordering: Ordering) -> Result<String, String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let (parse, _, case) = time(|| Case::from_path(path));
    let case = case.map_err(|err| err.to_string())?;
    let (model_time, _, model) = time(|| DcModel::new(&case));
    let model = model.map_err(|err| err.to_string())?;
    let (solve, solve_runs, solution) = time(|| model.solve_with_ordering(&case, ordering));
    let solution = solution.map_err(|err| err.to_string())?;
    let negative_branches = model
        .branch_rows()
        .iter()
        .filter(|&&row| case.branch[row].x * case.branch[row].effective_tap() < 0.0)
        .count();
    let (bprime_dimension, bprime_entries, factor_entries) =
        entries(&model, solution.solver, ordering)?;
    let backward_error = accuracy::backward_error(&case, &model, &solution.va_deg);

    let comparison = match references {
        Some(dir) if dir.join(format!("{name}_dcpf_bus.csv")).is_file() => {
            let reference = reference::read(dir, name)?;
            let deviation = reference::deviation(&case, &solution, &reference)?;
            let matpower_va_deg: Vec<f64> = reference.bus.iter().map(|bus| bus.1).collect();
            format!(
                "{:e},{:e},{:e},{:e}",
                accuracy::backward_error(&case, &model, &matpower_va_deg),
                deviation.va_deg,
                deviation.pf_mw,
                deviation.pg_mw
            )
        }
        _ => ",,,".to_string(),
    };
    Ok(format!(
        "{name},{},{},{},{},{},{negative_branches},{:.3},{:.3},{:.3},{solve_runs},{:?},{},\
         {bprime_dimension},{bprime_entries},{factor_entries},{backward_error:e},{comparison}",
        case.bus.len(),
        case.branch.len(),
        case.generator.len(),
        model.bus_rows().len(),
        model.branch_rows().len(),
        milliseconds(parse),
        milliseconds(model_time),
        milliseconds(solve),
        solution.solver,
        ordering_name(ordering),
    ))
}

/// B′'s dimension and stored entries, and the stored entries of the factor `solver` made of it.
/// B′ is built as [`DcModel::solve`] builds it: Bbus over the PV and PQ buses, in that order.
/// Both factorizations are ordered by `ordering`. A Cholesky factor stores its whole symbolic
/// pattern, which its analysis predicts exactly; an LU factor's depends on its pivots, so the
/// LU is factored again to count it.
fn entries(
    model: &DcModel,
    solver: DcSolver,
    ordering: Ordering,
) -> Result<(usize, usize, usize), String> {
    let unknown: Vec<usize> = model
        .pv_buses()
        .iter()
        .chain(model.pq_buses())
        .copied()
        .collect();
    let bprime = model
        .bbus()
        .select(&unknown, &unknown)
        .map_err(|err| err.to_string())?;
    let factor = match solver {
        DcSolver::Cholesky => bprime
            .analyze_cholesky(ordering)
            .map_err(|err| err.to_string())?
            .predicted_entries(),
        DcSolver::Lu => {
            let columns = bprime.to_csc();
            let lu = columns
                .analyze_lu(ordering)
                .and_then(|symbolic| symbolic.factor(&columns, LuOptions::default()))
                .map_err(|err| err.to_string())?;
            lu.l().num_stored() + lu.u().num_stored() - unknown.len()
        }
        DcSolver::NoSolve => 0,
    };
    Ok((unknown.len(), bprime.num_stored(), factor))
}

fn ordering_name(ordering: Ordering) -> &'static str {
    match ordering {
        Ordering::Natural => "natural",
        Ordering::Rcm => "rcm",
        Ordering::Amd => "amd",
        _ => "other",
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
