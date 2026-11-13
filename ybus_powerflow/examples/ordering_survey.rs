//! Compares LU orderings on one case's power-flow matrices, for the survey of every MATPOWER
//! case that `scripts/benchmark-ordering.sh` runs: the AC Newton Jacobian at a flat start, and
//! the DC B′ (Bbus over the PV and PQ buses), each factored by LU under AMD (on `A + Aᵀ`, with
//! diagonal pivots preferred) and under COLAMD (on `A`'s columns, with plain partial pivoting).
//!
//! Usage: `ordering_survey <case.m>`, or `ordering_survey --header`.
//!
//! Prints one CSV line, with the columns in `HEADER`. Each analysis and factorization is timed
//! together by its median over repeated runs, for up to `BUDGET` or `MAX_RUNS` runs.

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use ybus_powerflow::ac::AcModel;
use ybus_powerflow::case::Case;
use ybus_powerflow::dc::DcModel;
use ybus_sparse::compressed_matrix::{CscMatrix, LuOptions, Ordering};
use ybus_sparse::Complex;

/// For the Jacobian (`jacobian_*`) and B′ (`bprime_*`): the dimension and stored entries, then
/// for each ordering the LU's stored entries (L and U, the diagonal counted once) and the
/// median time to analyze and factor (ms). `*_colamd_diagonal_fill` is the fill under
/// COLAMD's ordering with diagonal pivots preferred, as under AMD
/// (`analyze_lu_with_permutation`), which separates the ordering from the pivoting.
const HEADER: &str = "case,buses,jacobian_dim,jacobian_entries,jacobian_amd_fill,\
                      jacobian_colamd_fill,jacobian_colamd_diagonal_fill,jacobian_amd_ms,\
                      jacobian_colamd_ms,bprime_dim,bprime_entries,bprime_amd_fill,\
                      bprime_colamd_fill,bprime_colamd_diagonal_fill,bprime_amd_ms,\
                      bprime_colamd_ms";

/// How long to keep repeating a factorization for its median time
const BUDGET: Duration = Duration::from_secs(1);
/// The most runs of a factorization
const MAX_RUNS: usize = 101;

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
            eprintln!("usage: ordering_survey <case.m> | --header");
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

    let ac = AcModel::new(&case).map_err(|err| err.to_string())?;
    let flat = vec![Complex::new(1.0, 0.0); ac.bus_rows().len()];
    let jacobian = ac.jacobian(&flat);

    let dc = DcModel::new(&case).map_err(|err| err.to_string())?;
    let unknown: Vec<usize> = dc.pv_buses().iter().chain(dc.pq_buses()).copied().collect();
    let bprime = dc
        .bbus()
        .select(&unknown, &unknown)
        .map_err(|err| err.to_string())?
        .to_csc();

    Ok(format!(
        "{name},{},{},{}",
        case.bus.len(),
        compare(&jacobian)?,
        compare(&bprime)?
    ))
}

/// The matrix's dimension and entries, then each ordering's fill, then each one's time.
fn compare(a: &CscMatrix<f64>) -> Result<String, String> {
    let n = a.shape().0;
    if n == 0 {
        return Ok("0,0,,,,,".to_string());
    }
    let (amd_time, amd_fill) = factor(a, Ordering::Amd)?;
    let (colamd_time, colamd_fill) = factor(a, Ordering::Colamd)?;
    let diagonal = a
        .analyze_lu_with_permutation(a.ordering(Ordering::Colamd).map_err(|e| e.to_string())?)
        .and_then(|symbolic| symbolic.factor(a, LuOptions::default()))
        .map_err(|err| err.to_string())?;
    let diagonal_fill = diagonal.l().num_stored() + diagonal.u().num_stored() - n;
    Ok(format!(
        "{n},{},{amd_fill},{colamd_fill},{diagonal_fill},{:.3},{:.3}",
        a.num_stored(),
        milliseconds(amd_time),
        milliseconds(colamd_time),
    ))
}

/// The median time to analyze and factor `a` under `ordering`, and the factor's entries.
fn factor(a: &CscMatrix<f64>, ordering: Ordering) -> Result<(Duration, usize), String> {
    let start = Instant::now();
    let mut times = vec![];
    loop {
        let run = Instant::now();
        let lu = a
            .analyze_lu(ordering)
            .and_then(|symbolic| symbolic.factor(a, LuOptions::default()))
            .map_err(|err| err.to_string())?;
        times.push(run.elapsed());
        if start.elapsed() >= BUDGET || times.len() >= MAX_RUNS {
            times.sort();
            let fill = lu.l().num_stored() + lu.u().num_stored() - a.shape().0;
            return Ok((times[times.len() / 2], fill));
        }
    }
}

fn milliseconds(duration: Duration) -> f64 {
    duration.as_secs_f64() * 1e3
}
