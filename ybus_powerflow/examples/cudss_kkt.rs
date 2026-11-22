//! Times one KKT system's work, a refactorization and two refined solves (an interior-point
//! iteration's predictor and corrector), on the CPU's pivoted LDLᵀ and on the GPU's cuDSS
//! (`ybus_cuda::GpuLdlt`), for O5 target 4's second gate. Needs the `cuda` and `rayon` features.
//!
//! Usage: `cudss_kkt <label> <kkt.mtx>...`, or `cudss_kkt --header`. Each file is a KKT matrix,
//! stored in full, as an interior-point solve factored it. For each:
//!
//! - **The CPU:** as `ybus_ipm`'s KKT solver runs it: AMD, threshold 1e-4, refactored into the
//!   last factorization's storage, each solve refined normwise to 1e-12 in at most 3 steps,
//!   stopping early when a step fails to halve the residual. It runs on rayon's pool
//!   (`RAYON_NUM_THREADS`).
//! - **The GPU:** cuDSS's analysis once, then refactorizations, with the same refinement against
//!   the host's matrix, the residuals computed on the host. The transfers are counted in each
//!   stage. Separately, the corrections the GPU needs to reach 1e-12 at all, up to 10.
//! - **Reported:** the time per stage (the median of `RUNS`), the inertia of each and whether
//!   they agree, the refined residuals, and the device memory in use after factoring.
//!
//! Prints one CSV line per matrix.

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::process::ExitCode;
use std::time::{Duration, Instant};

use ybus_cuda::{Gpu, GpuLdlt};
use ybus_sparse::compressed_matrix::{CsrMatrix, LdltOptions, Ordering};
use ybus_sparse::SparseMatrix;

const HEADER: &str = "label,matrix,n,entries,threads,cpu_analysis_ms,cpu_factor_ms,cpu_solves_ms,cpu_steps,cpu_residual,cpu_positive,cpu_negative,gpu_analysis_ms,gpu_factor_ms,gpu_solves_ms,gpu_steps,gpu_residual,gpu_steps_to_tolerance,gpu_positive,gpu_negative,inertia_matches,cpu_iteration_ms,gpu_iteration_ms,speedup,device_mib";

/// The KKT solver's settings (`ybus_ipm::kkt`)
const THRESHOLD: f64 = 1e-4;
const STEPS: usize = 3;
const TOLERANCE: f64 = 1e-12;
/// The most corrections tried for the GPU to reach `TOLERANCE`
const MORE_STEPS: usize = 10;
/// Timed repetitions of each stage; the median is reported
const RUNS: usize = 5;

/// `‖b − A x‖∞ / ‖b‖∞`, and the residual
fn relative_residual(a: &CsrMatrix<f64>, x: &[f64], b: &[f64]) -> (f64, Vec<f64>) {
    let ax = a.mul_vec(x).expect("square");
    let r: Vec<f64> = b.iter().zip(&ax).map(|(b, ax)| b - ax).collect();
    // A NaN counts as infinite: `f64::max` would drop it, and a NaN solution would look exact
    let inf = |v: &[f64]| {
        v.iter().fold(0.0f64, |m, x| {
            if x.is_nan() {
                f64::INFINITY
            } else {
                m.max(x.abs())
            }
        })
    };
    (inf(&r) / inf(b), r)
}

/// The CPU's refinement rule around any solve: at most `max_steps` corrections, stopping at
///   `TOLERANCE` or when a step fails to halve the residual, keeping the best solution. Returns
///   it, the corrections applied, and its residual
fn refine(
    a: &CsrMatrix<f64>,
    b: &[f64],
    max_steps: usize,
    mut solve: impl FnMut(&[f64]) -> Vec<f64>,
) -> (Vec<f64>, usize, f64) {
    let mut x = solve(b);
    let (mut ratio, mut r) = relative_residual(a, &x, b);
    let mut steps = 0;
    for step in 1..=max_steps {
        if ratio <= TOLERANCE {
            break;
        }
        let d = solve(&r);
        let candidate: Vec<f64> = x.iter().zip(&d).map(|(x, d)| x + d).collect();
        let (next, next_r) = relative_residual(a, &candidate, b);
        if next >= ratio {
            break;
        }
        let halved = next <= ratio / 2.0;
        (x, ratio, r, steps) = (candidate, next, next_r, step);
        if !halved {
            break;
        }
    }
    (x, steps, ratio)
}

fn median(mut times: Vec<Duration>) -> f64 {
    times.sort();
    times[times.len() / 2].as_secs_f64() * 1e3
}

/// Two deterministic right-hand sides of order `n`
fn right_hand_sides(n: usize) -> [Vec<f64>; 2] {
    let mut state = 12345u64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    [
        (0..n).map(|_| next()).collect(),
        (0..n).map(|_| next()).collect(),
    ]
}

fn run(label: &str, path: &str, gpu: &Gpu) -> Result<(), Box<dyn std::error::Error>> {
    let a: CsrMatrix<f64> =
        SparseMatrix::read_matrix_market(BufReader::new(File::open(path)?))?.to_csr();
    let n = a.shape().0;
    let rhs = right_hand_sides(n);
    let options = LdltOptions {
        threshold: THRESHOLD,
    };
    // The CPU
    let start = Instant::now();
    let symbolic = a.analyze_ldlt(Ordering::Amd)?;
    let cpu_analysis = start.elapsed().as_secs_f64() * 1e3;
    let mut ldlt = symbolic.factor(&a, options)?;
    let mut factor_times = vec![];
    for _ in 0..RUNS {
        let start = Instant::now();
        ldlt = symbolic.refactor(&a, options, ldlt)?;
        factor_times.push(start.elapsed());
    }
    let mut solve_times = vec![];
    let (mut cpu_steps, mut cpu_residual) = (0, 0.0f64);
    for _ in 0..RUNS {
        let start = Instant::now();
        (cpu_steps, cpu_residual) = (0, 0.0);
        for b in &rhs {
            let refined = ldlt.solve_refined_normwise(&a, b, STEPS, TOLERANCE)?;
            cpu_steps = cpu_steps.max(refined.steps);
            cpu_residual = cpu_residual.max(refined.relative_residual);
        }
        solve_times.push(start.elapsed());
    }
    let cpu_inertia = ldlt.inertia();
    // The GPU
    let start = Instant::now();
    let mut device = GpuLdlt::new(gpu, &a)?.ok_or("no cuDSS")?;
    let gpu_analysis = start.elapsed().as_secs_f64() * 1e3;
    device.factor(&a)?;
    let mut gpu_factor_times = vec![];
    for _ in 0..RUNS {
        let start = Instant::now();
        device.factor(&a)?;
        gpu_factor_times.push(start.elapsed());
    }
    let (free, total) = gpu.memory()?;
    let mut gpu_solve_times = vec![];
    let (mut gpu_steps, mut gpu_residual) = (0, 0.0f64);
    for _ in 0..RUNS {
        let start = Instant::now();
        (gpu_steps, gpu_residual) = (0, 0.0);
        for b in &rhs {
            let (_, steps, ratio) = refine(&a, b, STEPS, |v| device.solve(v).expect("factored"));
            gpu_steps = gpu_steps.max(steps);
            gpu_residual = gpu_residual.max(ratio);
        }
        gpu_solve_times.push(start.elapsed());
    }
    // The corrections the GPU needs to reach the tolerance at all
    let mut to_tolerance = 0;
    for b in &rhs {
        let (_, steps, ratio) = refine(&a, b, MORE_STEPS, |v| device.solve(v).expect("factored"));
        to_tolerance = to_tolerance.max(if ratio <= TOLERANCE {
            steps
        } else {
            usize::MAX
        });
    }
    let gpu_inertia = device.inertia()?;
    let (cpu_factor, cpu_solves) = (median(factor_times), median(solve_times));
    let (gpu_factor, gpu_solves) = (median(gpu_factor_times), median(gpu_solve_times));
    let (cpu_iteration, gpu_iteration) = (cpu_factor + cpu_solves, gpu_factor + gpu_solves);
    let matrix = path.rsplit('/').next().unwrap_or(path);
    println!(
        "{label},{matrix},{n},{},{},{cpu_analysis:.3},{cpu_factor:.3},{cpu_solves:.3},{cpu_steps},{cpu_residual:e},{},{},{gpu_analysis:.3},{gpu_factor:.3},{gpu_solves:.3},{gpu_steps},{gpu_residual:e},{},{},{},{},{cpu_iteration:.3},{gpu_iteration:.3},{:.3},{}",
        a.indices().len(),
        rayon::current_num_threads(),
        cpu_inertia.positive,
        cpu_inertia.negative,
        if to_tolerance == usize::MAX { "never".to_string() } else { to_tolerance.to_string() },
        gpu_inertia.positive,
        gpu_inertia.negative,
        cpu_inertia == gpu_inertia,
        cpu_iteration / gpu_iteration,
        (total - free) >> 20,
    );
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("--header") {
        println!("{HEADER}");
        return ExitCode::SUCCESS;
    }
    let [label, paths @ ..] = args.as_slice() else {
        eprintln!("usage: cudss_kkt <label> <kkt.mtx>... | --header");
        return ExitCode::FAILURE;
    };
    let gpu = match Gpu::new() {
        Ok(Some(gpu)) => gpu,
        Ok(None) => {
            eprintln!("no CUDA device");
            return ExitCode::FAILURE;
        }
        Err(err) => {
            eprintln!("{err}");
            return ExitCode::FAILURE;
        }
    };
    for path in paths {
        if let Err(err) = run(label, path, &gpu) {
            eprintln!("{label} {path}: {err}");
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}
