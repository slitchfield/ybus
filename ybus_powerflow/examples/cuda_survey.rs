//! Times solving many right-hand sides against one case's B′ (Cholesky) and Newton Jacobian at
//! the case voltages (LU), on the CPU (`solve_many`, and `par_solve_many` on rayon's threads)
//! and on the GPU (`ybus_cuda`'s `GpuFactor::solve_many`) at several block widths, for the CUDA
//! spike. Needs the `cuda` and `rayon` features.
//!
//! Usage: `cuda_survey <case.m> [<count>]`, or `cuda_survey --header`. `<count>` right-hand sides
//! (default: as many as fit in 1 GB, at most 4096) are solved once each way after one warm-up
//! block; the GPU's time is split into uploading the blocks, the solves and downloading the
//! results, with the host's permuting outside all three. Prints one CSV line per matrix and
//! block width.

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use ybus_cuda::{Gpu, GpuFactor};
use ybus_powerflow::ac::AcModel;
use ybus_powerflow::case::Case;
use ybus_powerflow::dc::DcModel;
use ybus_sparse::compressed_matrix::{CsrMatrix, LuOptions, Ordering};
use ybus_sparse::Complex;
use ybus_sparse::Factorization;

const HEADER: &str = "case,matrix,n,factor_entries,rhs,cpu_serial_ms,cpu_parallel_ms,upload_factor_ms,width,gpu_ms,gpu_upload_ms,gpu_solve_ms,gpu_download_ms,worst_relative_difference";

const WIDTHS: &[usize] = &[64, 256, 1024, 4096];

/// Bytes one set of right-hand sides may take
const BUDGET: usize = 1 << 30;

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    let (path, count) = match args.as_slice() {
        [flag] if flag == "--header" => {
            println!("{HEADER}");
            return ExitCode::SUCCESS;
        }
        [path] => (path, None),
        [path, count] => match count.parse::<usize>() {
            Ok(count) if count > 0 => (path, Some(count)),
            _ => return usage(),
        },
        _ => return usage(),
    };
    match run(Path::new(path), count) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{path}: {err}");
            ExitCode::FAILURE
        }
    }
}

fn usage() -> ExitCode {
    eprintln!("usage: cuda_survey <case.m> [<count>] | cuda_survey --header");
    ExitCode::from(2)
}

/// B′: Bbus over the PV then PQ buses, as the DC power flow solves it
fn b_prime(model: &DcModel) -> CsrMatrix<f64> {
    let unknown: Vec<usize> = model
        .pv_buses()
        .iter()
        .chain(model.pq_buses())
        .copied()
        .collect();
    let mut position = vec![None; model.bus_rows().len()];
    for (k, &bus) in unknown.iter().enumerate() {
        position[bus] = Some(k);
    }
    let triplets = model
        .bbus()
        .iter()
        .filter_map(|(i, j, v)| Some((position[i]?, position[j]?, v)));
    CsrMatrix::from_triplets(unknown.len(), unknown.len(), triplets).expect("indices in range")
}

/// `count` deterministic right-hand sides of length `n`, with entries in [-1, 1]
fn right_hand_sides(n: usize, count: usize) -> Vec<Vec<f64>> {
    (0..count)
        .map(|j| {
            (0..n)
                .map(|i| (((i * 7919 + j * 104729) % 2001) as f64) / 1000.0 - 1.0)
                .collect()
        })
        .collect()
}

fn worst_relative_difference(gpu: &[Vec<f64>], cpu: &[Vec<f64>]) -> f64 {
    gpu.iter()
        .zip(cpu)
        .map(|(g, c)| {
            let scale = c
                .iter()
                .fold(0.0f64, |m, v| m.max(v.abs()))
                .max(f64::MIN_POSITIVE);
            g.iter()
                .zip(c)
                .fold(0.0f64, |m, (g, c)| m.max((g - c).abs()))
                / scale
        })
        .fold(0.0, f64::max)
}

fn ms(seconds: f64) -> String {
    format!("{:.3}", seconds * 1e3)
}

/// A CPU solve of many right-hand sides
type SolveMany<'a> = &'a dyn Fn(&[Vec<f64>]) -> Vec<Vec<f64>>;

/// The CPU's solves of one factorization against the GPU's at each block width
struct Solver<'a> {
    name: &'a str,
    n: usize,
    entries: usize,
    serial: SolveMany<'a>,
    parallel: SolveMany<'a>,
    upload: &'a dyn Fn() -> Result<GpuFactor<'a>, ybus_cuda::CudaError>,
}

fn survey(case: &str, solver: &Solver, count: Option<usize>) -> Result<(), String> {
    let n = solver.n;
    let count = count.unwrap_or((BUDGET / 8 / n.max(1)).clamp(1, 4096));
    let rhs = right_hand_sides(n, count);
    (solver.serial)(&rhs[..1]);
    let start = Instant::now();
    let cpu = (solver.serial)(&rhs);
    let serial = start.elapsed().as_secs_f64();
    let start = Instant::now();
    drop((solver.parallel)(&rhs));
    let parallel = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let mut factor = (solver.upload)().map_err(|err| err.to_string())?;
    let upload = start.elapsed().as_secs_f64();
    for &width in WIDTHS.iter().filter(|&&width| width <= count) {
        factor.set_block_width(width);
        // A warm-up block, which also makes the block and its analysis
        factor
            .solve_many(&rhs[..width])
            .map_err(|err| err.to_string())?;
        factor.reset_timings();
        let start = Instant::now();
        let solved = factor.solve_many(&rhs).map_err(|err| err.to_string())?;
        let total = start.elapsed().as_secs_f64();
        let timings = factor.timings();
        println!(
            "{case},{},{n},{},{count},{},{},{},{width},{},{},{},{},{:.3e}",
            solver.name,
            solver.entries,
            ms(serial),
            ms(parallel),
            ms(upload),
            ms(total),
            ms(timings.upload.as_secs_f64()),
            ms(timings.solve.as_secs_f64()),
            ms(timings.download.as_secs_f64()),
            worst_relative_difference(&solved, &cpu),
        );
    }
    Ok(())
}

fn run(path: &Path, count: Option<usize>) -> Result<(), String> {
    let gpu = Gpu::new()
        .map_err(|err| err.to_string())?
        .ok_or("no CUDA device")?;
    let case_name = path
        .file_stem()
        .map_or("?".into(), |s| s.to_string_lossy().into_owned());
    let case = Case::from_path(path).map_err(|err| err.to_string())?;

    // As the DC power flow factors B′: Cholesky when it's positive definite, otherwise (as with
    //   negative reactances) LU
    let b = b_prime(&DcModel::new(&case).map_err(|err| err.to_string())?);
    match b
        .analyze_cholesky(Ordering::Amd)
        .and_then(|symbolic| symbolic.factor(&b))
    {
        Ok(cholesky) => survey(
            &case_name,
            &Solver {
                name: "bprime-cholesky",
                n: b.shape().0,
                entries: cholesky.l().num_stored(),
                serial: &|rhs| cholesky.solve_many(rhs).expect("one entry per row"),
                parallel: &|rhs| cholesky.par_solve_many(rhs).expect("one entry per row"),
                upload: &|| GpuFactor::from_cholesky(&gpu, &cholesky),
            },
            count,
        )?,
        Err(_) => {
            let b = b.to_csc();
            let lu = b
                .analyze_lu(Ordering::Amd)
                .and_then(|symbolic| symbolic.factor(&b, LuOptions::default()))
                .map_err(|err| format!("B′: {err}"))?;
            survey(
                &case_name,
                &Solver {
                    name: "bprime-lu",
                    n: b.shape().0,
                    entries: lu.l().num_stored() + lu.u().num_stored(),
                    serial: &|rhs| lu.solve_many(rhs).expect("one entry per row"),
                    parallel: &|rhs| lu.par_solve_many(rhs).expect("one entry per row"),
                    upload: &|| GpuFactor::from_lu(&gpu, &lu),
                },
                count,
            )?
        }
    }

    let model = AcModel::new(&case).map_err(|err| err.to_string())?;
    let v: Vec<Complex<f64>> = model
        .bus_rows()
        .iter()
        .map(|&row| Complex::from_polar(case.bus[row].vm_pu, case.bus[row].va_deg.to_radians()))
        .collect();
    let jacobian = model.jacobian(&v);
    let lu = jacobian
        .analyze_lu(Ordering::Amd)
        .and_then(|symbolic| symbolic.factor(&jacobian, LuOptions::default()))
        .map_err(|err| format!("Jacobian: {err}"))?;
    survey(
        &case_name,
        &Solver {
            name: "jacobian",
            n: jacobian.shape().0,
            entries: lu.l().num_stored() + lu.u().num_stored(),
            serial: &|rhs| lu.solve_many(rhs).expect("one entry per row"),
            parallel: &|rhs| lu.par_solve_many(rhs).expect("one entry per row"),
            upload: &|| GpuFactor::from_lu(&gpu, &lu),
        },
        count,
    )
}
