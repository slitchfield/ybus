//! The C interface's own cost: each solve through the exports against the same solve in Rust,
//! timed side by side in one process (`docs/benchmarks/ffi-boundary-2026-11-10.md`).
//!
//! The exports are called through the rlib, which runs exactly the code a C program runs: the
//! guard, the solution's handle, a copy of every column into buffers the caller allocated, and
//! the free. A real C caller adds only a plain function call. The Rust side calls `powerflow`'s
//! own function on a case it loaded itself, and keeps the result. Before anything is timed, the
//! two results are checked to agree bit for bit.
//!
//! Usage: `cargo run --release -p ybus_ffi --example boundary_cost -- [--opf] <case.m>...`.
//! With `--opf`, the DC and AC OPFs are timed too. It prints one Markdown table row per case
//! and solve: each side's median, and the difference, absolute and as a fraction of Rust's.
//! A solve takes milliseconds, against which the boundary's microseconds are lost in the
//! machine's noise, so the last column times the boundary's work alone: reading every column
//! of one kept solution, the median of 1001 reads.

use std::ffi::CString;
use std::hint::black_box;
use std::ptr;
use std::time::{Duration, Instant};

use powerflow::ac::{ac_power_flow_with_options, AcOptions};
use powerflow::dc::{dc_power_flow_with_options, DcOptions};
use powerflow::opf::ac::{ac_opf_with_options, AcOpfOptions};
use powerflow::opf::dc::{dc_opf_with_options, DcOpfOptions};
use powerflow::Case;
use ybus::ac::*;
use ybus::case::columns::*;
use ybus::case::{ybus_case_free, ybus_case_from_path, YbusCase};
use ybus::dc::*;
use ybus::opf::ac::*;
use ybus::opf::dc::*;
use ybus::status::{ybus_status, YBUS_OK};

/// How many timed runs of each side, after one warm-up: the power flows take milliseconds, so
/// more runs steady their medians; the OPFs take far longer
const POWER_FLOW_RUNS: usize = 21;
const OPF_RUNS: usize = 5;

/// A table of the case, for the length of its columns
#[derive(Clone, Copy)]
enum Table {
    Bus,
    Gen,
    Branch,
}

/// A solution column's getter
type Get<S> = unsafe extern "C" fn(*const S, ybus_field, *mut f64, usize) -> ybus_status;

/// Reads a solution's columns, each `(table, field)` in `columns`, into `buffers`, which have
/// one entry per row of their tables
fn read<S>(
    solution: *const S,
    getters: (Get<S>, Get<S>, Get<S>),
    columns: &[(Table, ybus_field)],
    buffers: &mut [Vec<f64>],
) {
    for (&(table, field), buffer) in columns.iter().zip(buffers) {
        let get = match table {
            Table::Bus => getters.0,
            Table::Gen => getters.1,
            Table::Branch => getters.2,
        };
        // SAFETY: solution is live, and buffer has room for the table's rows
        let status = unsafe { get(solution, field, buffer.as_mut_ptr(), buffer.len()) };
        assert_eq!(status, YBUS_OK);
    }
}

/// The median of `samples`
fn median(mut samples: Vec<Duration>) -> Duration {
    samples.sort();
    samples[samples.len() / 2]
}

/// A Rust solve, returning its result's columns
type RustSolve<'a> = Box<dyn Fn() -> Vec<Vec<f64>> + 'a>;
/// A solve through the exports, reading its columns into the buffers
type CSolve<'a> = Box<dyn Fn(&mut [Vec<f64>]) + 'a>;
/// A solve through the exports, then the median time of reading all its columns, many times
type CReads<'a> = Box<dyn Fn(&mut [Vec<f64>]) -> Duration + 'a>;

/// One solve to time: its name, its columns, the Rust solve (returning those columns, which
/// are read only for the check), and the solve through the exports into the buffers
struct Solve<'a> {
    name: &'static str,
    columns: Vec<(Table, ybus_field)>,
    rust: RustSolve<'a>,
    through_c: CSolve<'a>,
    reads: CReads<'a>,
    runs: usize,
}

/// How many times the columns of one kept solution are read, for the boundary's work alone
const READ_RUNS: usize = 1001;

/// The median time of reading all of a solution's columns, `READ_RUNS` times
fn time_reads<S>(
    solution: *const S,
    getters: (Get<S>, Get<S>, Get<S>),
    columns: &[(Table, ybus_field)],
    buffers: &mut [Vec<f64>],
) -> Duration {
    median(
        (0..READ_RUNS)
            .map(|_| {
                let start = Instant::now();
                read(solution, getters, columns, black_box(&mut *buffers));
                start.elapsed()
            })
            .collect(),
    )
}

fn solves<'a>(case: &'a Case, handle: *const YbusCase, opf: bool) -> Vec<Solve<'a>> {
    let mut solves = vec![
        Solve {
            name: "DC power flow",
            columns: dc_columns(),
            rust: Box::new(move || {
                let s = black_box(dc_power_flow_with_options(case, &DcOptions::default()).unwrap());
                vec![s.va_deg, s.pg_mw, s.pf_mw]
            }),
            through_c: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(
                        ybus_dc_power_flow(handle, ptr::null(), &mut solution),
                        YBUS_OK
                    );
                    let getters = (
                        ybus_dc_solution_get_bus as Get<_>,
                        ybus_dc_solution_get_gen as Get<_>,
                        ybus_dc_solution_get_branch as Get<_>,
                    );
                    read(solution, getters, &dc_columns(), buffers);
                    ybus_dc_solution_free(solution);
                }
            }),
            reads: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(
                        ybus_dc_power_flow(handle, ptr::null(), &mut solution),
                        YBUS_OK
                    );
                    let getters = (
                        ybus_dc_solution_get_bus as Get<_>,
                        ybus_dc_solution_get_gen as Get<_>,
                        ybus_dc_solution_get_branch as Get<_>,
                    );
                    let time = time_reads(solution, getters, &dc_columns(), buffers);
                    ybus_dc_solution_free(solution);
                    time
                }
            }),
            runs: POWER_FLOW_RUNS,
        },
        Solve {
            name: "AC power flow",
            columns: ac_columns(),
            rust: Box::new(move || {
                let s = black_box(ac_power_flow_with_options(case, &AcOptions::default()).unwrap());
                vec![
                    s.vm_pu, s.va_deg, s.pg_mw, s.qg_mvar, s.pf_mw, s.qf_mvar, s.pt_mw, s.qt_mvar,
                ]
            }),
            through_c: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                let mut iterations = 0;
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(
                        ybus_ac_power_flow(handle, ptr::null(), &mut solution),
                        YBUS_OK
                    );
                    let getters = (
                        ybus_ac_solution_get_bus as Get<_>,
                        ybus_ac_solution_get_gen as Get<_>,
                        ybus_ac_solution_get_branch as Get<_>,
                    );
                    read(solution, getters, &ac_columns(), buffers);
                    assert_eq!(
                        ybus_ac_solution_iterations(solution, &mut iterations),
                        YBUS_OK
                    );
                    ybus_ac_solution_free(solution);
                }
            }),
            reads: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(
                        ybus_ac_power_flow(handle, ptr::null(), &mut solution),
                        YBUS_OK
                    );
                    let getters = (
                        ybus_ac_solution_get_bus as Get<_>,
                        ybus_ac_solution_get_gen as Get<_>,
                        ybus_ac_solution_get_branch as Get<_>,
                    );
                    let time = time_reads(solution, getters, &ac_columns(), buffers);
                    ybus_ac_solution_free(solution);
                    time
                }
            }),
            runs: POWER_FLOW_RUNS,
        },
    ];
    if opf {
        solves.push(Solve {
            name: "DC OPF",
            columns: dc_opf_columns(),
            rust: Box::new(move || {
                let s = black_box(dc_opf_with_options(case, &DcOpfOptions::default()).unwrap());
                vec![
                    s.va_deg,
                    s.lam_p,
                    s.pg_mw,
                    s.mu_pmin,
                    s.mu_pmax,
                    s.pf_mw,
                    s.mu_sf,
                    s.mu_st,
                    s.mu_angmin,
                    s.mu_angmax,
                ]
            }),
            through_c: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                let mut objective = 0.0;
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(ybus_dc_opf(handle, ptr::null(), &mut solution), YBUS_OK);
                    let getters = (
                        ybus_dc_opf_solution_get_bus as Get<_>,
                        ybus_dc_opf_solution_get_gen as Get<_>,
                        ybus_dc_opf_solution_get_branch as Get<_>,
                    );
                    read(solution, getters, &dc_opf_columns(), buffers);
                    assert_eq!(
                        ybus_dc_opf_solution_objective(solution, &mut objective),
                        YBUS_OK
                    );
                    ybus_dc_opf_solution_free(solution);
                }
            }),
            reads: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(ybus_dc_opf(handle, ptr::null(), &mut solution), YBUS_OK);
                    let getters = (
                        ybus_dc_opf_solution_get_bus as Get<_>,
                        ybus_dc_opf_solution_get_gen as Get<_>,
                        ybus_dc_opf_solution_get_branch as Get<_>,
                    );
                    let time = time_reads(solution, getters, &dc_opf_columns(), buffers);
                    ybus_dc_opf_solution_free(solution);
                    time
                }
            }),
            runs: OPF_RUNS,
        });
        solves.push(Solve {
            name: "AC OPF",
            columns: ac_opf_columns(),
            rust: Box::new(move || {
                let s = black_box(ac_opf_with_options(case, &AcOpfOptions::default()).unwrap());
                vec![
                    s.vm_pu,
                    s.va_deg,
                    s.lam_p,
                    s.lam_q,
                    s.mu_vmax,
                    s.mu_vmin,
                    s.pg_mw,
                    s.qg_mvar,
                    s.vg_pu,
                    s.mu_pmax,
                    s.mu_pmin,
                    s.mu_qmax,
                    s.mu_qmin,
                    s.pf_mw,
                    s.qf_mvar,
                    s.pt_mw,
                    s.qt_mvar,
                    s.mu_sf,
                    s.mu_st,
                    s.mu_angmin,
                    s.mu_angmax,
                ]
            }),
            through_c: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                let mut objective = 0.0;
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(ybus_ac_opf(handle, ptr::null(), &mut solution), YBUS_OK);
                    let getters = (
                        ybus_ac_opf_solution_get_bus as Get<_>,
                        ybus_ac_opf_solution_get_gen as Get<_>,
                        ybus_ac_opf_solution_get_branch as Get<_>,
                    );
                    read(solution, getters, &ac_opf_columns(), buffers);
                    assert_eq!(
                        ybus_ac_opf_solution_objective(solution, &mut objective),
                        YBUS_OK
                    );
                    ybus_ac_opf_solution_free(solution);
                }
            }),
            reads: Box::new(move |buffers| {
                let mut solution = ptr::null_mut();
                // SAFETY: handle is live, and solution writable and freed once
                unsafe {
                    assert_eq!(ybus_ac_opf(handle, ptr::null(), &mut solution), YBUS_OK);
                    let getters = (
                        ybus_ac_opf_solution_get_bus as Get<_>,
                        ybus_ac_opf_solution_get_gen as Get<_>,
                        ybus_ac_opf_solution_get_branch as Get<_>,
                    );
                    let time = time_reads(solution, getters, &ac_opf_columns(), buffers);
                    ybus_ac_opf_solution_free(solution);
                    time
                }
            }),
            runs: OPF_RUNS,
        });
    }
    solves
}

fn dc_columns() -> Vec<(Table, ybus_field)> {
    vec![
        (Table::Bus, YBUS_BUS_VA_DEG),
        (Table::Gen, YBUS_GEN_PG_MW),
        (Table::Branch, YBUS_BRANCH_PF_MW),
    ]
}

fn ac_columns() -> Vec<(Table, ybus_field)> {
    use Table::{Branch, Bus, Gen};
    vec![
        (Bus, YBUS_BUS_VM_PU),
        (Bus, YBUS_BUS_VA_DEG),
        (Gen, YBUS_GEN_PG_MW),
        (Gen, YBUS_GEN_QG_MVAR),
        (Branch, YBUS_BRANCH_PF_MW),
        (Branch, YBUS_BRANCH_QF_MVAR),
        (Branch, YBUS_BRANCH_PT_MW),
        (Branch, YBUS_BRANCH_QT_MVAR),
    ]
}

fn dc_opf_columns() -> Vec<(Table, ybus_field)> {
    use Table::{Branch, Bus, Gen};
    vec![
        (Bus, YBUS_BUS_VA_DEG),
        (Bus, YBUS_BUS_LAM_P),
        (Gen, YBUS_GEN_PG_MW),
        (Gen, YBUS_GEN_MU_PMIN),
        (Gen, YBUS_GEN_MU_PMAX),
        (Branch, YBUS_BRANCH_PF_MW),
        (Branch, YBUS_BRANCH_MU_SF),
        (Branch, YBUS_BRANCH_MU_ST),
        (Branch, YBUS_BRANCH_MU_ANGMIN),
        (Branch, YBUS_BRANCH_MU_ANGMAX),
    ]
}

fn ac_opf_columns() -> Vec<(Table, ybus_field)> {
    use Table::{Branch, Bus, Gen};
    vec![
        (Bus, YBUS_BUS_VM_PU),
        (Bus, YBUS_BUS_VA_DEG),
        (Bus, YBUS_BUS_LAM_P),
        (Bus, YBUS_BUS_LAM_Q),
        (Bus, YBUS_BUS_MU_VMAX),
        (Bus, YBUS_BUS_MU_VMIN),
        (Gen, YBUS_GEN_PG_MW),
        (Gen, YBUS_GEN_QG_MVAR),
        (Gen, YBUS_GEN_VG_PU),
        (Gen, YBUS_GEN_MU_PMAX),
        (Gen, YBUS_GEN_MU_PMIN),
        (Gen, YBUS_GEN_MU_QMAX),
        (Gen, YBUS_GEN_MU_QMIN),
        (Branch, YBUS_BRANCH_PF_MW),
        (Branch, YBUS_BRANCH_QF_MVAR),
        (Branch, YBUS_BRANCH_PT_MW),
        (Branch, YBUS_BRANCH_QT_MVAR),
        (Branch, YBUS_BRANCH_MU_SF),
        (Branch, YBUS_BRANCH_MU_ST),
        (Branch, YBUS_BRANCH_MU_ANGMIN),
        (Branch, YBUS_BRANCH_MU_ANGMAX),
    ]
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let opf = args.iter().any(|arg| arg == "--opf");
    args.retain(|arg| arg != "--opf");
    if args.is_empty() {
        eprintln!("usage: boundary_cost [--opf] <case.m>...");
        std::process::exit(2);
    }
    println!(
        "| Case | Solve | Rust (median) | Through C (median) | Difference | Relative | \
         Reading the columns |"
    );
    println!("|---|---|---:|---:|---:|---:|---:|");
    for path in &args {
        let case = Case::from_path(path).unwrap();
        let name = std::path::Path::new(path)
            .file_stem()
            .unwrap()
            .to_string_lossy();
        let c_path = CString::new(path.as_str()).unwrap();
        let mut handle = ptr::null_mut();
        // SAFETY: c_path is NUL-terminated, and handle writable
        assert_eq!(
            unsafe { ybus_case_from_path(c_path.as_ptr(), &mut handle) },
            YBUS_OK
        );
        let rows = |table: Table| match table {
            Table::Bus => case.bus.len(),
            Table::Gen => case.generator.len(),
            Table::Branch => case.branch.len(),
        };
        for solve in solves(&case, handle, opf) {
            let mut buffers: Vec<Vec<f64>> = solve
                .columns
                .iter()
                .map(|&(table, _)| vec![0.0; rows(table)])
                .collect();
            // The check, which is also each side's warm-up
            (solve.through_c)(&mut buffers);
            let expected = (solve.rust)();
            assert!(
                buffers.iter().zip(&expected).all(|(c, rust)| c
                    .iter()
                    .map(|v| v.to_bits())
                    .eq(rust.iter().map(|v| v.to_bits()))),
                "{name} {}: the results differ",
                solve.name
            );
            let (mut rust, mut through_c) = (vec![], vec![]);
            for _ in 0..solve.runs {
                let start = Instant::now();
                black_box((solve.rust)());
                rust.push(start.elapsed());
                let start = Instant::now();
                (solve.through_c)(black_box(&mut buffers));
                through_c.push(start.elapsed());
            }
            let (rust, through_c) = (median(rust), median(through_c));
            let difference = through_c.as_secs_f64() - rust.as_secs_f64();
            let reads = (solve.reads)(&mut buffers);
            println!(
                "| {name} | {} | {:.3} ms | {:.3} ms | {:+.1} µs | {:+.2}% | {:.1} µs |",
                solve.name,
                rust.as_secs_f64() * 1e3,
                through_c.as_secs_f64() * 1e3,
                difference * 1e6,
                difference / rust.as_secs_f64() * 100.0,
                reads.as_secs_f64() * 1e6
            );
        }
        // SAFETY: handle is live until this free
        unsafe { ybus_case_free(handle) };
    }
}
