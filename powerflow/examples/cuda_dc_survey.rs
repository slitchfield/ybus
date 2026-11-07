//! Times the DC power flow's GPU path against the CPU's on one case, for the CUDA spike's
//! survey (`scripts/benchmark-cuda-dc.sh` runs it once per case, each in its own process).
//! Needs the `cuda` and `rayon` features.
//!
//! Usage: `cuda_dc_survey <case.m>`, or `cuda_dc_survey --header`. Prints one CSV line per
//! workload as it finishes:
//!
//! - `batched_dc`: DC solves of copies of the case with its loads scaled (eight scalings,
//!   repeated): `DcStudy::solve` for each, serially and on rayon's pool, against
//!   `DcGpuStudy::solve_many`.
//! - `ptdf_columns`: PTDF columns for evenly spaced buses, MATPOWER's slack (skipped for a case
//!   whose islands don't all reach it):
//!   `DcSensitivity::ptdf_columns` and `par_ptdf_columns` against `DcGpuSensitivity`'s.
//! - `n1_solutions`: single-branch outages of evenly spaced in-service branches, solved whole:
//!   `DcStudy::solve_outages` and `par_solve_outages` against `DcGpuStudy::solve_outages`.
//! - `n1_screening`: the same outages' most loaded branches: the CPU's solutions reduced on the
//!   host against `DcGpuStudy::outage_loadings`.
//!
//! Each workload runs in chunks of the GPU's block width (`width`), on both paths, whose results
//! are dropped as they come: a GPU block always solves its full width, so smaller chunks would
//! pad. It's sized so the work is at most `WORK` (right-hand sides times unknowns), at most
//! `MAX_COUNT` items, and skipped if a chunk's results would take more than `CHUNK_BYTES`. Each pass is timed once, after a warm-up chunk;
//! `gpu_setup_ms` is uploading the factorization (and Bf). Then the first and a middle chunk
//! are compared: `worst_relative_difference` is the largest difference relative to each
//! result's largest magnitude (for screening, to the larger of the loading and 1), and
//! `mismatches` counts items where one path fails and the other doesn't, or their errors
//! differ.

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use powerflow::case::Case;
use powerflow::dc::{DcModel, DcSolution, DcStudy, OutageLoading, PtdfSlack};
use powerflow::powerflow_error::PowerflowError;
use powerflow::Gpu;
use rayon::prelude::*;

const HEADER: &str = "case,buses,branches,unknowns,solver,workload,count,width,\
                      cpu_serial_ms,cpu_parallel_ms,gpu_setup_ms,gpu_ms,\
                      speedup_vs_parallel,worst_relative_difference,mismatches";

/// The most bytes one chunk's results may take, on each path
const CHUNK_BYTES: usize = 2 << 30;
/// The most right-hand sides times unknowns per workload
const WORK: usize = 1 << 28;
/// The most items per workload
const MAX_COUNT: usize = 8192;

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.as_slice() {
        [flag] if flag == "--header" => {
            println!("{HEADER}");
            ExitCode::SUCCESS
        }
        [path] => match run(Path::new(path)) {
            Ok(()) => ExitCode::SUCCESS,
            Err(err) => {
                eprintln!("{path}: {err}");
                ExitCode::FAILURE
            }
        },
        _ => {
            eprintln!("usage: cuda_dc_survey <case.m> | cuda_dc_survey --header");
            ExitCode::from(2)
        }
    }
}

fn ms(seconds: f64) -> String {
    format!("{:.3}", seconds * 1e3)
}

/// The seconds `f` takes, and its result
fn time<T>(f: impl FnOnce() -> T) -> (f64, T) {
    let start = Instant::now();
    let value = f();
    (start.elapsed().as_secs_f64(), value)
}

/// `count` items out of `0..len`, evenly spaced
fn spaced(len: usize, count: usize) -> Vec<usize> {
    if len == 0 {
        return vec![];
    }
    let count = count.min(len);
    (0..count).map(|i| i * len / count).collect()
}

/// How many items a workload takes, for results of `bytes` each, solves of `unknowns` and
/// chunks of `chunk`
fn count_for(bytes: usize, unknowns: usize, chunk: usize) -> usize {
    let chunks_fit = CHUNK_BYTES / bytes.max(1) >= chunk;
    let count = (WORK / unknowns.max(1)).min(MAX_COUNT);
    if chunks_fit {
        count.max(1)
    } else {
        0
    }
}

/// The largest difference between `a` and `b`, relative to `b`'s largest magnitude
fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    let scale = b
        .iter()
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(f64::MIN_POSITIVE);
    a.iter()
        .zip(b)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()))
        / scale
}

fn solution_difference(a: &DcSolution, b: &DcSolution) -> f64 {
    relative_difference(&a.va_deg, &b.va_deg)
        .max(relative_difference(&a.pf_mw, &b.pf_mw))
        .max(relative_difference(&a.pg_mw, &b.pg_mw))
}

/// The worst difference and the mismatches between two chunks of results
fn compare<T>(
    gpu: &[Result<T, PowerflowError>],
    cpu: &[Result<T, PowerflowError>],
    difference: impl Fn(&T, &T) -> f64,
) -> (f64, usize) {
    let (mut worst, mut mismatches) = (0.0f64, 0);
    for (gpu, cpu) in gpu.iter().zip(cpu) {
        match (gpu, cpu) {
            (Ok(gpu), Ok(cpu)) => worst = worst.max(difference(gpu, cpu)),
            (Err(gpu), Err(cpu)) if gpu == cpu => {}
            _ => mismatches += 1,
        }
    }
    (worst, mismatches)
}

/// The most loaded rated branch of a CPU solution, as `OutageLoading` defines it
fn cpu_loading(case: &Case, model: &DcModel, solution: &DcSolution) -> OutageLoading {
    let mut best = OutageLoading {
        branch: None,
        loading: 0.0,
    };
    for &row in model.branch_rows() {
        let rate = case.branch[row].rate_a_mva;
        if rate > 0.0 {
            let loading = solution.pf_mw[row].abs() / rate;
            if loading > best.loading {
                best = OutageLoading {
                    branch: Some(row),
                    loading,
                };
            }
        }
    }
    best
}

/// One workload's timings and comparison
struct Row {
    workload: &'static str,
    count: usize,
    width: usize,
    serial: f64,
    parallel: f64,
    setup: f64,
    gpu: f64,
    worst: f64,
    mismatches: usize,
}

/// The case's fixed columns
struct Prefix(String);

impl Prefix {
    fn print(&self, row: &Row) {
        println!(
            "{},{},{},{},{},{},{},{},{:.2},{:.3e},{}",
            self.0,
            row.workload,
            row.count,
            row.width,
            ms(row.serial),
            ms(row.parallel),
            ms(row.setup),
            ms(row.gpu),
            row.parallel / row.gpu,
            row.worst,
            row.mismatches,
        );
    }
}

/// Times a chunked workload three ways and compares the first and a middle chunk of the GPU's
/// results with the serial CPU's
#[allow(clippy::too_many_arguments)]
fn chunked<I: Sync, T: Send>(
    items: &[I],
    chunk: usize,
    serial: impl Fn(&[I]) -> Vec<Result<T, PowerflowError>>,
    parallel: impl Fn(&[I]) -> Vec<Result<T, PowerflowError>>,
    mut gpu: impl FnMut(&[I]) -> Vec<Result<T, PowerflowError>>,
    difference: impl Fn(&T, &T) -> f64,
) -> (f64, f64, f64, f64, usize) {
    let chunks: Vec<&[I]> = items.chunks(chunk).collect();
    // Warm-up: the GPU's block and analysis, rayon's pool
    drop(gpu(&chunks[0][..chunks[0].len().min(8)]));
    drop(parallel(&chunks[0][..chunks[0].len().min(8)]));
    let (serial_time, ()) = time(|| chunks.iter().for_each(|chunk| drop(serial(chunk))));
    let (parallel_time, ()) = time(|| chunks.iter().for_each(|chunk| drop(parallel(chunk))));
    let (gpu_time, ()) = time(|| chunks.iter().for_each(|chunk| drop(gpu(chunk))));
    let (mut worst, mut mismatches) = (0.0f64, 0);
    let mut sample = vec![0, chunks.len() / 2];
    sample.dedup();
    for index in sample {
        let (w, m) = compare(&gpu(chunks[index]), &serial(chunks[index]), &difference);
        worst = worst.max(w);
        mismatches += m;
    }
    (serial_time, parallel_time, gpu_time, worst, mismatches)
}

fn run(path: &Path) -> Result<(), String> {
    let gpu = Gpu::new()
        .map_err(|err| err.to_string())?
        .ok_or("no CUDA device")?;
    let name = path
        .file_stem()
        .map_or("?".into(), |s| s.to_string_lossy().into_owned());
    let case = Case::from_path(path).map_err(|err| err.to_string())?;
    let model = DcModel::new(&case).map_err(|err| err.to_string())?;
    let study = model.study().map_err(|err| err.to_string())?;
    let unknowns = model.pv_buses().len() + model.pq_buses().len();
    let prefix = Prefix(format!(
        "{name},{},{},{unknowns},{:?}",
        case.bus.len(),
        case.branch.len(),
        study.solver()
    ));
    let solution_bytes = 8 * (case.bus.len() + case.branch.len() + case.generator.len());

    batched(&prefix, &gpu, &case, &study, unknowns, solution_bytes)?;
    ptdf(&prefix, &gpu, &case, &study, unknowns)?;
    outages(
        &prefix,
        &gpu,
        &case,
        &model,
        &study,
        unknowns,
        solution_bytes,
    )?;
    Ok(())
}

fn batched(
    prefix: &Prefix,
    gpu: &Gpu,
    case: &Case,
    study: &DcStudy,
    unknowns: usize,
    solution_bytes: usize,
) -> Result<(), String> {
    let (setup, on_gpu) = time(|| study.gpu(gpu));
    let mut on_gpu = on_gpu.map_err(|err| err.to_string())?;
    let width = on_gpu.block_width().max(1);
    let count = count_for(solution_bytes, unknowns, width);
    if count == 0 {
        return Ok(());
    }
    let variants: Vec<Case> = (0..8)
        .map(|i| {
            let mut variant = case.clone();
            for bus in &mut variant.bus {
                bus.pd_mw *= 0.8 + 0.05 * i as f64;
            }
            variant
        })
        .collect();
    let cases: Vec<&Case> = (0..count).map(|i| &variants[i % variants.len()]).collect();
    let (serial, parallel, gpu_time, worst, mismatches) = chunked(
        &cases,
        width,
        |chunk| chunk.iter().map(|case| study.solve(case)).collect(),
        |chunk| chunk.par_iter().map(|case| study.solve(case)).collect(),
        |chunk| on_gpu.solve_many(chunk).expect("the GPU works"),
        solution_difference,
    );
    prefix.print(&Row {
        workload: "batched_dc",
        count,
        width,
        serial,
        parallel,
        setup,
        gpu: gpu_time,
        worst,
        mismatches,
    });
    Ok(())
}

fn ptdf(
    prefix: &Prefix,
    gpu: &Gpu,
    case: &Case,
    study: &DcStudy,
    unknowns: usize,
) -> Result<(), String> {
    // A case with several islands, each with its own reference, has no single slack every bus
    //   reaches, so no PTDF (as on the CPU)
    let sensitivity = match study.sensitivity(PtdfSlack::Reference) {
        Ok(sensitivity) => sensitivity,
        Err(err @ PowerflowError::IslandWithoutReference { .. }) => {
            eprintln!("ptdf_columns skipped: {err}");
            return Ok(());
        }
        Err(err) => return Err(err.to_string()),
    };
    let (setup, on_gpu) = time(|| sensitivity.gpu(gpu));
    let mut on_gpu = on_gpu.map_err(|err| err.to_string())?;
    let width = on_gpu.block_width().max(1);
    let count = count_for(8 * case.branch.len(), unknowns, width);
    if count == 0 {
        return Ok(());
    }
    let buses = spaced(case.bus.len(), count);
    let wrap = |columns: Result<Vec<Vec<f64>>, PowerflowError>| -> Vec<Result<Vec<f64>, _>> {
        columns
            .expect("bus rows in range")
            .into_iter()
            .map(Ok)
            .collect()
    };
    let (serial, parallel, gpu_time, worst, mismatches) = chunked(
        &buses,
        width,
        |chunk| wrap(sensitivity.ptdf_columns(chunk)),
        |chunk| wrap(sensitivity.par_ptdf_columns(chunk)),
        |chunk| wrap(on_gpu.ptdf_columns(chunk)),
        |a, b| relative_difference(a, b),
    );
    prefix.print(&Row {
        workload: "ptdf_columns",
        count: buses.len(),
        width,
        serial,
        parallel,
        setup,
        gpu: gpu_time,
        worst,
        mismatches,
    });
    Ok(())
}

fn outages(
    prefix: &Prefix,
    gpu: &Gpu,
    case: &Case,
    model: &DcModel,
    study: &DcStudy,
    unknowns: usize,
    solution_bytes: usize,
) -> Result<(), String> {
    let (setup, on_gpu) = time(|| study.gpu(gpu));
    let mut on_gpu = on_gpu.map_err(|err| err.to_string())?;
    let width = on_gpu.block_width().max(1);
    let count = count_for(solution_bytes, unknowns, width);
    if count == 0 {
        return Ok(());
    }
    let in_service = model.branch_rows();
    let rows: Vec<usize> = spaced(in_service.len(), count)
        .into_iter()
        .map(|i| in_service[i])
        .collect();
    let sets = |chunk: &[usize]| -> Vec<Vec<usize>> { chunk.iter().map(|&r| vec![r]).collect() };
    let (serial, parallel, gpu_time, worst, mismatches) = chunked(
        &rows,
        width,
        |chunk| {
            study
                .solve_outages(case, &sets(chunk))
                .expect("the case matches")
        },
        |chunk| {
            study
                .par_solve_outages(case, &sets(chunk))
                .expect("the case matches")
        },
        |chunk| on_gpu.solve_outages(case, chunk).expect("the GPU works"),
        solution_difference,
    );
    prefix.print(&Row {
        workload: "n1_solutions",
        count: rows.len(),
        width,
        serial,
        parallel,
        setup,
        gpu: gpu_time,
        worst,
        mismatches,
    });

    let loadings = |results: Vec<Result<DcSolution, PowerflowError>>| {
        results
            .into_iter()
            .map(|result| result.map(|solution| cpu_loading(case, model, &solution)))
            .collect()
    };
    // The flow operator is made on the first screening, so it's part of the setup
    let (setup, on_gpu) = time(|| -> Result<_, PowerflowError> {
        let mut on_gpu = study.gpu(gpu)?;
        on_gpu.outage_loadings(case, &[])?;
        Ok(on_gpu)
    });
    let mut on_gpu = on_gpu.map_err(|err| err.to_string())?;
    let (serial, parallel, gpu_time, worst, mismatches) = chunked(
        &rows,
        width,
        |chunk| {
            loadings(
                study
                    .solve_outages(case, &sets(chunk))
                    .expect("the case matches"),
            )
        },
        |chunk| {
            loadings(
                study
                    .par_solve_outages(case, &sets(chunk))
                    .expect("the case matches"),
            )
        },
        |chunk| on_gpu.outage_loadings(case, chunk).expect("the GPU works"),
        |a: &OutageLoading, b: &OutageLoading| (a.loading - b.loading).abs() / b.loading.max(1.0),
    );
    prefix.print(&Row {
        workload: "n1_screening",
        count: rows.len(),
        width,
        serial,
        parallel,
        setup,
        gpu: gpu_time,
        worst,
        mismatches,
    });
    Ok(())
}
