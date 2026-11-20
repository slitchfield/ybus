//! How much parallelism one OPF solve offers, for phase O5's targets: the predictor-corrector
//! method's KKT front tree (`KktStatistics::tree`, from the solve's own LDLᵀ analysis at its
//! default ordering), and the units an AC evaluation could split by.
//!
//! Usage: `o5_parallelism <case.m>…`, or `o5_parallelism --header`. Prints one CSV line per
//! case and model (DC, AC), with the columns in `HEADER`:
//! - **The tree** (see `FrontTree`): its fronts, roots and height, the factorization's work
//!   and its critical path, and one triangular sweep's. `factor_bound` and `solve_bound` are
//!   work over critical path, the most that parallelism over independent subtrees can speed
//!   the factorization and a level-scheduled sweep.
//! - **The units:** the case's buses, in-service generators and branches, and the branches
//!   with flow limits, which an AC evaluation's per-element plans would split.
//!
//! The solve's outcome is printed too, since a failed solve's tree is still its last
//! analysis's.

use std::env;
use std::path::Path;
use std::process::ExitCode;

use ybus_powerflow::case::Case;
use ybus_powerflow::opf::ac::{ac_opf_with_options, AcOpfOptions};
use ybus_powerflow::opf::dc::{dc_opf_with_options, DcOpfOptions};
use ybus_powerflow::opf::{KktStatistics, Status};

const HEADER: &str = "case,model,outcome,fronts,roots,height,work,critical_work,factor_bound,\
                      solve_work,solve_critical_work,solve_bound,buses,generators,branches,\
                      limited_branches";

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.as_slice() == ["--header"] {
        println!("{HEADER}");
        return ExitCode::SUCCESS;
    }
    if args.is_empty() {
        eprintln!("usage: o5_parallelism <case.m>… | --header");
        return ExitCode::from(2);
    }
    for path in &args {
        if let Err(err) = report(Path::new(path)) {
            eprintln!("{path}: {err}");
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}

fn outcome(status: Status) -> &'static str {
    match status {
        Status::Converged => "converged",
        Status::MaxIterations => "max iterations",
        Status::NumericalFailure => "numerical failure",
    }
}

fn report(path: &Path) -> Result<(), String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let case = Case::from_path(path).map_err(|err| err.to_string())?;
    let buses = case.bus.len();
    let generators = case.generator.iter().filter(|g| g.status).count();
    let branches = case.branch.iter().filter(|b| b.status).count();
    let limited = case
        .branch
        .iter()
        .filter(|b| b.status && b.rate_a_mva != 0.0 && b.rate_a_mva < 1e10)
        .count();
    let line = |model: &str, status: Status, kkt: &KktStatistics| {
        let t = kkt.tree;
        let ratio = |a: u64, b: u64| if b == 0 { 0.0 } else { a as f64 / b as f64 };
        println!(
            "{name},{model},{},{},{},{},{},{},{:.3},{},{},{:.3},{buses},{generators},{branches},\
             {limited}",
            outcome(status),
            t.fronts,
            t.roots,
            t.height,
            t.work,
            t.critical_work,
            ratio(t.work, t.critical_work),
            t.solve_work,
            t.solve_critical_work,
            ratio(t.solve_work, t.solve_critical_work),
        );
    };
    let dc = dc_opf_with_options(&case, &DcOpfOptions::mehrotra()).map_err(|e| e.to_string())?;
    line("dc", dc.status, &dc.solver.kkt);
    let ac = ac_opf_with_options(&case, &AcOpfOptions::mehrotra()).map_err(|e| e.to_string())?;
    line("ac", ac.status, &ac.solver.kkt);
    Ok(())
}
