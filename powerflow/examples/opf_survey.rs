//! Times the DC or AC optimal power flow on one case and compares it with MATPOWER's, for the
//! survey of every MATPOWER case with costs that `scripts/benchmark-opf.sh` runs. The script runs
//! this once per case and model, each in its own process, so that it can limit time and memory
//! and record peak memory.
//!
//! Usage: `opf_survey --model dc|ac [--solver mips|mehrotra] <case.m> [<reference-dir>]`, or
//! `opf_survey --header`. The solver is MIPS by default, as in MATPOWER, or the
//! predictor-corrector method, each at MATPOWER's OPF settings (`DcOpfOptions::mehrotra`,
//! `AcOpfOptions::mehrotra`).
//!
//! Prints one CSV line, with the columns in `HEADER`. Parsing and the OPF (model and solve
//! together, as MATPOWER's elapsed time covers its setup and solve) are each timed by their median
//! over repeated runs, for up to `BUDGET` or `MAX_RUNS` runs. Given a directory of MATPOWER's
//! results from `scripts/generate-reference.sh`, it also compares the results with MATPOWER's
//! (MP-Core) and, where they were generated, the legacy implementation's.

#[allow(dead_code)]
#[path = "../tests/common/reference.rs"]
mod reference;

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use interior_point::solution::{KktResiduals, Status};
use powerflow::case::Case;
use powerflow::opf::ac::{ac_opf_with_options, AcOpfOptions, AcOpfSolution};
use powerflow::opf::dc::{dc_opf_with_options, DcOpfOptions, DcOpfSolution};
use reference::{AcOpfReference, DcOpfReference};

/// The columns.
///
/// - `outcome`: `converged`, `max iterations`, `numerical failure`, or the error.
/// - **Against MATPOWER (MP-Core):** its outcome (`matpower_success`), iterations, whether
///   Octave's least-squares fallback was taken (`singular`) and elapsed time; then the
///   objective's relative difference, the bus prices' relative differences (`LAM_P`, and `LAM_Q`
///   for AC; relative to each price's magnitude, at least 1), and the largest absolute
///   differences in the primal solution (angles, magnitudes, dispatch and flows) and in the bound
///   and flow-limit multipliers. The local tests judge the outcome, iterations, objective and bus
///   prices; the primal and multiplier differences are recorded, since degenerate problems leave
///   them non-unique.
/// - **Against the legacy implementation,** where its results were generated: its outcome,
///   iterations and time, and the same differences of ours from it and of it from MP-Core.
/// - **The solver and its KKT residuals** (last): `mips` or `mehrotra`, and the residuals at its
///   point and multipliers on the OPF's own problem, unscaled (`DcOpfSolution::residuals`),
///   empty for an error.
const HEADER: &str = "case,buses,branches,generators,model,parse_ms,opf_ms,opf_runs,outcome,\
                      iterations,matpower_success,matpower_iterations,singular,matpower_ms,\
                      objective_rel,lam_p_rel,lam_q_rel,primal_abs,multiplier_abs,\
                      legacy_success,legacy_iterations,legacy_ms,\
                      legacy_objective_rel,legacy_lam_p_rel,legacy_lam_q_rel,legacy_primal_abs,\
                      legacy_multiplier_abs,cores_objective_rel,cores_lam_p_rel,cores_lam_q_rel,\
                      cores_primal_abs,cores_multiplier_abs,solver,kkt_stationarity,\
                      kkt_feasibility,kkt_complementarity,kkt_dual_feasibility";

/// How long to keep repeating a step for its median time
const BUDGET: Duration = Duration::from_secs(1);
/// The most runs of a step
const MAX_RUNS: usize = 101;

#[derive(Clone, Copy, PartialEq)]
enum Model {
    Dc,
    Ac,
}

#[derive(Clone, Copy, PartialEq)]
enum Solver {
    Mips,
    Mehrotra,
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.as_slice() == ["--header"] {
        println!("{HEADER}");
        return ExitCode::SUCCESS;
    }
    let (model, rest) = match args.as_slice() {
        [flag, model, rest @ ..] if flag == "--model" && model == "dc" => (Model::Dc, rest),
        [flag, model, rest @ ..] if flag == "--model" && model == "ac" => (Model::Ac, rest),
        _ => return usage(),
    };
    let (solver, rest) = match rest {
        [flag, solver, rest @ ..] if flag == "--solver" && solver == "mips" => (Solver::Mips, rest),
        [flag, solver, rest @ ..] if flag == "--solver" && solver == "mehrotra" => {
            (Solver::Mehrotra, rest)
        }
        [flag, ..] if flag == "--solver" => return usage(),
        rest => (Solver::Mips, rest),
    };
    let (case, references) = match rest {
        [case] => (case, None),
        [case, references] => (case, Some(Path::new(references))),
        _ => return usage(),
    };
    match survey(Path::new(case), references, model, solver) {
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
        "usage: opf_survey --model dc|ac [--solver mips|mehrotra] <case.m> [<reference-dir>] \
         | --header"
    );
    ExitCode::from(2)
}

/// The quantities two OPF results are compared by.
struct Fields {
    objective: f64,
    success: bool,
    /// How our solve ended, and its KKT residuals (none for a reference)
    status: Option<Status>,
    residuals: Option<KktResiduals>,
    iterations: usize,
    lam_p: Vec<f64>,
    /// `LAM_Q`, for AC only
    lam_q: Option<Vec<f64>>,
    primal: Vec<Vec<f64>>,
    multipliers: Vec<Vec<f64>>,
}

impl Fields {
    fn dc(s: &DcOpfSolution) -> Fields {
        Fields {
            objective: s.objective,
            success: s.converged(),
            status: Some(s.status),
            residuals: Some(s.residuals),
            iterations: s.iterations,
            lam_p: s.lam_p.clone(),
            lam_q: None,
            primal: vec![s.va_deg.clone(), s.pg_mw.clone(), s.pf_mw.clone()],
            multipliers: vec![
                s.mu_pmin.clone(),
                s.mu_pmax.clone(),
                s.mu_sf.clone(),
                s.mu_st.clone(),
                s.mu_angmin.clone(),
                s.mu_angmax.clone(),
            ],
        }
    }

    fn dc_reference(r: &DcOpfReference) -> Fields {
        Fields {
            objective: r.objective,
            success: r.success,
            status: None,
            residuals: None,
            iterations: r.iterations,
            lam_p: r.lam_p.clone(),
            lam_q: None,
            primal: vec![r.va_deg.clone(), r.pg_mw.clone(), r.pf_mw.clone()],
            multipliers: vec![
                r.mu_pmin.clone(),
                r.mu_pmax.clone(),
                r.mu_sf.clone(),
                r.mu_st.clone(),
                r.mu_angmin.clone(),
                r.mu_angmax.clone(),
            ],
        }
    }

    fn ac(s: &AcOpfSolution) -> Fields {
        Fields {
            objective: s.objective,
            success: s.converged(),
            status: Some(s.status),
            residuals: Some(s.residuals),
            iterations: s.iterations,
            lam_p: s.lam_p.clone(),
            lam_q: Some(s.lam_q.clone()),
            primal: vec![
                s.vm_pu.clone(),
                s.va_deg.clone(),
                s.pg_mw.clone(),
                s.qg_mvar.clone(),
                s.pf_mw.clone(),
                s.qf_mvar.clone(),
                s.pt_mw.clone(),
                s.qt_mvar.clone(),
            ],
            multipliers: vec![
                s.mu_vmin.clone(),
                s.mu_vmax.clone(),
                s.mu_pmin.clone(),
                s.mu_pmax.clone(),
                s.mu_qmin.clone(),
                s.mu_qmax.clone(),
                s.mu_sf.clone(),
                s.mu_st.clone(),
                s.mu_angmin.clone(),
                s.mu_angmax.clone(),
            ],
        }
    }

    fn ac_reference(r: &AcOpfReference) -> Fields {
        Fields {
            objective: r.objective,
            success: r.success,
            status: None,
            residuals: None,
            iterations: r.iterations,
            lam_p: r.lam_p.clone(),
            lam_q: Some(r.lam_q.clone()),
            primal: vec![
                r.vm_pu.clone(),
                r.va_deg.clone(),
                r.pg_mw.clone(),
                r.qg_mvar.clone(),
                r.pf_mw.clone(),
                r.qf_mvar.clone(),
                r.pt_mw.clone(),
                r.qt_mvar.clone(),
            ],
            multipliers: vec![
                r.mu_vmin.clone(),
                r.mu_vmax.clone(),
                r.mu_pmin.clone(),
                r.mu_pmax.clone(),
                r.mu_qmin.clone(),
                r.mu_qmax.clone(),
                r.mu_sf.clone(),
                r.mu_st.clone(),
                r.mu_angmin.clone(),
                r.mu_angmax.clone(),
            ],
        }
    }
}

/// The largest `|a − e|` over corresponding vectors, or NaN if any is NaN or a length differs.
fn largest_absolute(actual: &[Vec<f64>], expected: &[Vec<f64>]) -> f64 {
    actual
        .iter()
        .zip(expected)
        .flat_map(|(a, e)| {
            let lengths_differ = (a.len() != e.len()).then_some(f64::NAN);
            a.iter()
                .zip(e)
                .map(|(a, e)| (a - e).abs())
                .chain(lengths_differ)
                .collect::<Vec<_>>()
        })
        .fold(0.0, |largest, d| {
            if d.is_nan() || d > largest {
                d
            } else {
                largest
            }
        })
}

/// The five differences of `a` from `b`, as CSV fields: the objective's relative difference, the
/// bus prices' relative differences, and the largest primal and multiplier differences.
fn differences(a: &Fields, b: &Fields) -> String {
    let objective = (a.objective - b.objective).abs() / b.objective.abs().max(1.0);
    let lam_q = match (&a.lam_q, &b.lam_q) {
        (Some(a), Some(b)) => format!("{:e}", reference::relative_difference(a, b)),
        _ => String::new(),
    };
    format!(
        "{objective:e},{:e},{lam_q},{:e},{:e}",
        reference::relative_difference(&a.lam_p, &b.lam_p),
        largest_absolute(&a.primal, &b.primal),
        largest_absolute(&a.multipliers, &b.multipliers),
    )
}

/// The CSV line for the case at `path`. Errors only if the case can't be read; an OPF that fails
/// is reported in the line.
fn survey(
    path: &Path,
    references: Option<&Path>,
    model: Model,
    solver: Solver,
) -> Result<String, String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let (parse, _, case) = time(|| Case::from_path(path));
    let case = case.map_err(|err| err.to_string())?;
    let (opf, runs, ours) = match model {
        Model::Dc => {
            let options = match solver {
                Solver::Mips => DcOpfOptions::default(),
                Solver::Mehrotra => DcOpfOptions::mehrotra(),
            };
            let (t, n, s) = time(|| dc_opf_with_options(&case, &options));
            (t, n, s.map(|s| Fields::dc(&s)))
        }
        Model::Ac => {
            let options = match solver {
                Solver::Mips => AcOpfOptions::default(),
                Solver::Mehrotra => AcOpfOptions::mehrotra(),
            };
            let (t, n, s) = time(|| ac_opf_with_options(&case, &options));
            (t, n, s.map(|s| Fields::ac(&s)))
        }
    };
    let (outcome, iterations) = match &ours {
        Ok(fields) => (
            match fields.status.expect("our own solve has a status") {
                Status::Converged => "converged",
                Status::MaxIterations => "max iterations",
                Status::NumericalFailure => "numerical failure",
            }
            .to_string(),
            fields.iterations.to_string(),
        ),
        Err(err) => (
            format!("error: {}", err.to_string().replace(',', ";")),
            String::new(),
        ),
    };
    let (run, legacy_run) = match model {
        Model::Dc => ("dcopf", "dcopf_legacy"),
        Model::Ac => ("opf", "opf_legacy"),
    };
    let read = |dir: &Path, run: &str| -> Option<Result<(Fields, bool), String>> {
        dir.join(format!("{name}_{run}_summary.csv"))
            .is_file()
            .then(|| match model {
                Model::Dc => reference::read_dcopf(dir, name, run)
                    .map(|r| (Fields::dc_reference(&r), r.singular)),
                Model::Ac => reference::read_acopf(dir, name, run)
                    .map(|r| (Fields::ac_reference(&r), r.singular)),
            })
    };
    let milliseconds_of = |dir: &Path, run: &str| {
        reference::read_time(dir, name, run).map_or(String::new(), |s| format!("{:.3}", s * 1e3))
    };
    let empty_differences = ",,,,";
    let (core, legacy) = match references {
        Some(dir) => {
            let core = read(dir, run).transpose()?;
            let legacy = read(dir, legacy_run).transpose()?;
            let core_columns = match (&core, &ours) {
                (Some((theirs, singular)), Ok(ours)) => format!(
                    "{},{},{},{},{}",
                    theirs.success,
                    theirs.iterations,
                    singular,
                    milliseconds_of(dir, run),
                    differences(ours, theirs)
                ),
                (Some((theirs, singular)), Err(_)) => format!(
                    "{},{},{},{},{empty_differences}",
                    theirs.success,
                    theirs.iterations,
                    singular,
                    milliseconds_of(dir, run)
                ),
                (None, _) => format!(",,,,{empty_differences}"),
            };
            let legacy_columns = match (&legacy, &ours, &core) {
                (Some((legacy, _)), ours, core) => {
                    let against_ours = match ours {
                        Ok(ours) => differences(ours, legacy),
                        Err(_) => empty_differences.to_string(),
                    };
                    let cores = match core {
                        Some((core, _)) => differences(legacy, core),
                        None => empty_differences.to_string(),
                    };
                    format!(
                        "{},{},{},{against_ours},{cores}",
                        legacy.success,
                        legacy.iterations,
                        milliseconds_of(dir, legacy_run)
                    )
                }
                (None, _, _) => format!(",,,{empty_differences},{empty_differences}"),
            };
            (core_columns, legacy_columns)
        }
        None => (
            format!(",,,,{empty_differences}"),
            format!(",,,{empty_differences},{empty_differences}"),
        ),
    };
    let model_name = match model {
        Model::Dc => "dc",
        Model::Ac => "ac",
    };
    let solver_name = match solver {
        Solver::Mips => "mips",
        Solver::Mehrotra => "mehrotra",
    };
    let kkt = match ours.as_ref().ok().and_then(|fields| fields.residuals) {
        Some(r) => format!(
            "{:e},{:e},{:e},{:e}",
            r.stationarity, r.feasibility, r.complementarity, r.dual_feasibility
        ),
        None => ",,,".to_string(),
    };
    Ok(format!(
        "{name},{},{},{},{model_name},{:.3},{:.3},{runs},{outcome},{iterations},{core},{legacy},\
         {solver_name},{kkt}",
        case.bus.len(),
        case.branch.len(),
        case.generator.len(),
        milliseconds(parse),
        milliseconds(opf),
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
