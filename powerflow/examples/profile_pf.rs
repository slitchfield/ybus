//! Runs one power flow workload on one case, repeating the whole pipeline (parse, model, solve)
//! so that a sampling profiler sees each stage in proportion to its share of one run, for
//! `scripts/profile-pf.sh`. The survey examples time each step separately, for its median, so a
//! profile of them would weight the steps by their timing budgets instead.
//!
//! Usage: `profile_pf <workload> <case.m> [<seconds>]`, with `<workload>` one of `parse`
//! (reading the case alone), `dc`, `nr`, `fdxb`, `fdbx`, `qlim2` (Newton's method with reactive
//! limits enforced one at a time), `n1` (every single-branch outage by rank-k updates) and
//! `n1-parallel` (the same in parallel, which needs the `rayon` feature), and the optimal power
//! flows: `dcopf` and `acopf` on MIPS, MATPOWER's default, and `dcopf-mehrotra` and
//! `acopf-mehrotra` on the predictor-corrector method, each at MATPOWER's OPF settings (MIPS
//! factoring the KKT matrix by LU, the predictor-corrector method by LDLᵀ), and
//! `dcopf-mehrotra-ldlt`, the DC OPF on the predictor-corrector method with the KKT matrix
//! factored by LDLᵀ explicitly, as its default now does (kept for the earlier profiles). The
//! pipeline repeats until `<seconds>` (default 5) have passed, at least once. Prints the number of
//! runs and the mean time per run to stderr.

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use powerflow::ac::{AcMethod, AcModel, AcOptions, ReactiveLimits};
use powerflow::case::Case;
use powerflow::dc::{DcModel, DcOptions};
use powerflow::opf::ac::{ac_opf_with_options, AcOpfOptions};
use powerflow::opf::dc::{dc_opf_with_options, DcOpfOptions};
use powerflow::opf::KktFactorization;

/// Outages per call to `solve_outages`, as in `n1_survey`, so a case's thousands of results are
/// never all held at once
const CHUNK: usize = 1024;

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    let (workload, case, seconds) = match args.as_slice() {
        [workload, case] => (workload.as_str(), case.as_str(), 5.0),
        [workload, case, seconds] => match seconds.parse::<f64>() {
            Ok(seconds) => (workload.as_str(), case.as_str(), seconds),
            Err(_) => return usage(),
        },
        _ => return usage(),
    };
    let budget = Duration::from_secs_f64(seconds);
    let start = Instant::now();
    let mut runs = 0;
    while runs == 0 || start.elapsed() < budget {
        if let Err(err) = run(workload, Path::new(case)) {
            eprintln!("{case}: {err}");
            return ExitCode::FAILURE;
        }
        runs += 1;
    }
    eprintln!(
        "{workload} {case}: {runs} runs, {:.3} ms per run",
        start.elapsed().as_secs_f64() * 1e3 / runs as f64
    );
    ExitCode::SUCCESS
}

fn usage() -> ExitCode {
    eprintln!(
        "usage: profile_pf parse|dc|nr|fdxb|fdbx|qlim2|n1|n1-parallel|dcopf|acopf|dcopf-mehrotra|\
         dcopf-mehrotra-ldlt|acopf-mehrotra <case.m> [<seconds>]"
    );
    ExitCode::from(2)
}

/// One run of `workload` on the case at `path`, from parsing to the solution. A solve that
/// doesn't converge still counts as a run; only a case that can't be read or modelled is an
/// error.
fn run(workload: &str, path: &Path) -> Result<(), String> {
    let case = Case::from_path(path).map_err(|err| err.to_string())?;
    let ac = |method: AcMethod, limits: ReactiveLimits| -> Result<(), String> {
        let model = AcModel::new(&case).map_err(|err| err.to_string())?;
        let mut options = AcOptions::for_method(method);
        options.reactive_limits = limits;
        // The outcome doesn't matter here, only the work
        let _ = model.solve(&case, &options);
        Ok(())
    };
    match workload {
        "parse" => Ok(()),
        "dc" => {
            let model = DcModel::new(&case).map_err(|err| err.to_string())?;
            let _ = model.solve(&case, &DcOptions::default());
            Ok(())
        }
        "nr" => ac(AcMethod::Newton, ReactiveLimits::Ignore),
        "fdxb" => ac(AcMethod::FastDecoupledXb, ReactiveLimits::Ignore),
        "fdbx" => ac(AcMethod::FastDecoupledBx, ReactiveLimits::Ignore),
        "qlim2" => ac(AcMethod::Newton, ReactiveLimits::OneAtATime),
        "n1" | "n1-parallel" => {
            let model = DcModel::new(&case).map_err(|err| err.to_string())?;
            let study = model
                .study(&DcOptions::default())
                .map_err(|err| err.to_string())?;
            for chunk in model.branch_rows().chunks(CHUNK) {
                let sets: Vec<Vec<usize>> = chunk.iter().map(|&row| vec![row]).collect();
                let _ = if workload == "n1" {
                    study.solve_outages(&case, &sets)
                } else {
                    parallel(&study, &case, &sets)?
                };
            }
            Ok(())
        }
        // The OPF's outcome doesn't matter here either, only the work
        "dcopf" | "dcopf-mehrotra" | "dcopf-mehrotra-ldlt" => {
            let mut options = if workload == "dcopf" {
                DcOpfOptions::default()
            } else {
                DcOpfOptions::mehrotra()
            };
            // LDLᵀ explicitly, as the DC OPF's default has since
            //   docs/benchmarks/ldlt-amalgamation-2026-11-08.md; the workload predates that
            if workload == "dcopf-mehrotra-ldlt" {
                options.solver = options
                    .solver
                    .with_kkt_factorization(KktFactorization::Ldlt);
            }
            dc_opf_with_options(&case, &options).map_err(|err| err.to_string())?;
            Ok(())
        }
        "acopf" | "acopf-mehrotra" => {
            let options = if workload == "acopf" {
                AcOpfOptions::default()
            } else {
                AcOpfOptions::mehrotra()
            };
            ac_opf_with_options(&case, &options).map_err(|err| err.to_string())?;
            Ok(())
        }
        _ => Err(format!("unknown workload {workload}")),
    }
}

type Outcomes = Result<
    Vec<Result<powerflow::dc::DcSolution, powerflow::powerflow_error::PowerflowError>>,
    powerflow::powerflow_error::PowerflowError,
>;

#[cfg(feature = "rayon")]
fn parallel(
    study: &powerflow::dc::DcStudy,
    case: &Case,
    sets: &[Vec<usize>],
) -> Result<Outcomes, String> {
    Ok(study.par_solve_outages(case, sets))
}

#[cfg(not(feature = "rayon"))]
fn parallel(_: &powerflow::dc::DcStudy, _: &Case, _: &[Vec<usize>]) -> Result<Outcomes, String> {
    Err("n1-parallel needs the rayon feature".to_string())
}
