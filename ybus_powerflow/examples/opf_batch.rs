//! Times batches of OPF solves serially and in parallel (`opf::batch`, feature `rayon`), for
//! phase O5's target 1 survey that `scripts/benchmark-opf-batch.sh` runs.
//!
//! Usage: `opf_batch --model dc|ac --solver mips|mehrotra --batch outages|hours
//! [--outages <n>] [--threads <t>[,<t>…]] [--line-search off|monotone|adaptive] <case.m>`, or
//! `opf_batch --header`. Each solver runs at MATPOWER's OPF settings with its own defaults
//! (`warm_shift` included). `--line-search` sets the predictor-corrector method's `LineSearch`
//! (`Off` by default; MIPS ignores it), and is its own column.
//!
//! **The batches:**
//! - `outages`: a seeded sample of `n` branch outages (not bridges) and `n` generator outages
//!   (default 14 each, 28 in all), as `opf_sequence` samples them, through `dc_opf_outages` or
//!   `ac_opf_outages`
//! - `hours`: the walk's 24 load levels after the case's own, through `solve_many` on the
//!   case's model
//!
//! **The timing.** Each batch is solved cold, then warm (from the base solution), each once
//! serially. Then the `par_` form runs on a rayon pool of each thread count (default 1, 2, 4,
//! 8, 14, 28), timed once each. Every parallel result is compared with the serial one, bit for
//! bit, through its `Debug` form. A difference prints `false` in `bit_for_bit`, and the example
//! exits nonzero, so a survey can't record a wrong speedup.
//!
//! Prints one CSV line per batch, start and thread count, with the columns in `HEADER`. The
//! serial time is on every line, and `speedup` is serial over parallel.

#[path = "shared/sequences.rs"]
mod sequences;

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use sequences::{loaded, outage_rows, sample, steps};
use ybus_powerflow::case::Case;
use ybus_powerflow::opf::ac::{AcOpfModel, AcOpfOptions};
use ybus_powerflow::opf::batch::{
    ac_opf_outages, dc_opf_outages, par_ac_opf_outages, par_dc_opf_outages, BatchStart, OpfOutage,
};
use ybus_powerflow::opf::dc::{DcOpfModel, DcOpfOptions};
use ybus_powerflow::opf::{LineSearch, OpfSolver, OpfStart};

const HEADER: &str = "case,model,solver,line_search,batch,start,threads,solves,converged,\
                      serial_ms,parallel_ms,speedup,bit_for_bit";

/// The AC OPF's options: MIPS's, or the predictor-corrector method's with `line_search`
fn ac_options(mehrotra: bool, line_search: LineSearch) -> AcOpfOptions {
    if !mehrotra {
        return AcOpfOptions::default();
    }
    let mut options = AcOpfOptions::mehrotra();
    if let OpfSolver::Mehrotra(m) = &mut options.solver {
        m.line_search = line_search;
    }
    options
}

/// The DC OPF's options, as `ac_options`
fn dc_options(mehrotra: bool, line_search: LineSearch) -> DcOpfOptions {
    if !mehrotra {
        return DcOpfOptions::default();
    }
    let mut options = DcOpfOptions::mehrotra();
    if let OpfSolver::Mehrotra(m) = &mut options.solver {
        m.line_search = line_search;
    }
    options
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.as_slice() == ["--header"] {
        println!("{HEADER}");
        return ExitCode::SUCCESS;
    }
    let (ac, rest) = match args.as_slice() {
        [flag, model, rest @ ..] if flag == "--model" && model == "dc" => (false, rest),
        [flag, model, rest @ ..] if flag == "--model" && model == "ac" => (true, rest),
        _ => return usage(),
    };
    let (mehrotra, rest) = match rest {
        [flag, solver, rest @ ..] if flag == "--solver" && solver == "mips" => (false, rest),
        [flag, solver, rest @ ..] if flag == "--solver" && solver == "mehrotra" => (true, rest),
        _ => return usage(),
    };
    let (outages, rest) = match rest {
        [flag, batch, rest @ ..] if flag == "--batch" && batch == "outages" => (true, rest),
        [flag, batch, rest @ ..] if flag == "--batch" && batch == "hours" => (false, rest),
        _ => return usage(),
    };
    let (per_kind, rest) = match rest {
        [flag, n, rest @ ..] if flag == "--outages" => match n.parse() {
            Ok(n) => (n, rest),
            Err(_) => return usage(),
        },
        rest => (14, rest),
    };
    let (threads, rest) = match rest {
        [flag, list, rest @ ..] if flag == "--threads" => {
            match list
                .split(',')
                .map(str::parse)
                .collect::<Result<Vec<usize>, _>>()
            {
                Ok(threads) => (threads, rest),
                Err(_) => return usage(),
            }
        }
        rest => (vec![1, 2, 4, 8, 14, 28], rest),
    };
    let (line_search, rest) = match rest {
        [flag, line_search, rest @ ..] if flag == "--line-search" => match line_search.as_str() {
            "off" => (LineSearch::Off, rest),
            "monotone" => (LineSearch::MonotoneFilter, rest),
            "adaptive" => (LineSearch::AdaptiveFilter, rest),
            _ => return usage(),
        },
        rest => (LineSearch::Off, rest),
    };
    let [path] = rest else { return usage() };
    match run(
        Path::new(path),
        ac,
        mehrotra,
        line_search,
        outages,
        per_kind,
        &threads,
    ) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => {
            eprintln!("{path}: a parallel batch differed from the serial one");
            ExitCode::FAILURE
        }
        Err(err) => {
            eprintln!("{path}: {err}");
            ExitCode::FAILURE
        }
    }
}

fn usage() -> ExitCode {
    eprintln!(
        "usage: opf_batch --model dc|ac --solver mips|mehrotra --batch outages|hours \
         [--outages <n>] [--threads <t>[,<t>…]] [--line-search off|monotone|adaptive] \
         <case.m> | --header"
    );
    ExitCode::from(2)
}

/// `f` timed, in milliseconds
fn timed<R>(f: impl FnOnce() -> R) -> (f64, R) {
    let begin = Instant::now();
    let result = f();
    (begin.elapsed().as_secs_f64() * 1e3, result)
}

/// `f` on a rayon pool of `threads` threads
fn on<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> Result<R, String> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|err| err.to_string())?;
    Ok(pool.install(f))
}

/// One batch, solved serially and then in parallel at each thread count: whether every
///   parallel result was the serial one's
struct Batch<'a> {
    prefix: String,
    threads: &'a [usize],
}

impl Batch<'_> {
    /// Times `serial` once and `parallel` on each pool, comparing them by their `Debug` forms;
    ///   `converged` counts a result's converged solves
    fn time<R: std::fmt::Debug + Send>(
        &self,
        start: &str,
        serial: impl Fn() -> R,
        parallel: impl Fn() -> R + Sync,
        solves: usize,
        converged: impl Fn(&R) -> usize,
    ) -> Result<bool, String> {
        let (serial_ms, expected) = timed(&serial);
        let expected_bits = format!("{expected:?}");
        let mut all = true;
        for &threads in self.threads {
            let (parallel_ms, result) = on(threads, || timed(&parallel))?;
            let same = format!("{result:?}") == expected_bits;
            all &= same;
            println!(
                "{},{start},{threads},{solves},{},{serial_ms:.3},{parallel_ms:.3},{:.3},{same}",
                self.prefix,
                converged(&expected),
                serial_ms / parallel_ms,
            );
        }
        Ok(all)
    }
}

fn run(
    path: &Path,
    ac: bool,
    mehrotra: bool,
    line_search: LineSearch,
    outages: bool,
    per_kind: usize,
    threads: &[usize],
) -> Result<bool, String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let case = Case::from_path(path).map_err(|err| err.to_string())?;
    let batch = Batch {
        prefix: format!(
            "{name},{},{},{},{}",
            if ac { "ac" } else { "dc" },
            if mehrotra { "mehrotra" } else { "mips" },
            match line_search {
                LineSearch::MonotoneFilter => "monotone",
                LineSearch::AdaptiveFilter => "adaptive",
                _ => "off",
            },
            if outages { "outages" } else { "hours" },
        ),
        threads,
    };
    let starts = [("cold", BatchStart::Cold), ("warm", BatchStart::Warm)];
    let mut all = true;
    if outages {
        let list: Vec<OpfOutage> = sample(outage_rows(&case, true)?, per_kind)
            .into_iter()
            .map(OpfOutage::branch)
            .chain(
                sample(outage_rows(&case, false)?, per_kind)
                    .into_iter()
                    .map(OpfOutage::generator),
            )
            .collect();
        for (label, start) in starts {
            all &= if ac {
                let options = ac_options(mehrotra, line_search);
                batch.time(
                    label,
                    || ac_opf_outages(&case, &list, &options, start).map_err(|e| e.to_string()),
                    || par_ac_opf_outages(&case, &list, &options, start).map_err(|e| e.to_string()),
                    list.len(),
                    |r| {
                        r.as_ref().map_or(0, |b| {
                            b.outages
                                .iter()
                                .filter(|o| o.as_ref().is_ok_and(|s| s.converged()))
                                .count()
                        })
                    },
                )?
            } else {
                let options = dc_options(mehrotra, line_search);
                batch.time(
                    label,
                    || dc_opf_outages(&case, &list, &options, start).map_err(|e| e.to_string()),
                    || par_dc_opf_outages(&case, &list, &options, start).map_err(|e| e.to_string()),
                    list.len(),
                    |r| {
                        r.as_ref().map_or(0, |b| {
                            b.outages
                                .iter()
                                .filter(|o| o.as_ref().is_ok_and(|s| s.converged()))
                                .count()
                        })
                    },
                )?
            };
        }
    } else {
        // The walk's load levels after the case's own, all from the base solution
        let cases: Vec<Case> = steps("walk", case.bus.len())
            .into_iter()
            .skip(1)
            .map(|(factors, _)| loaded(&case, &factors))
            .collect();
        if ac {
            let options = ac_options(mehrotra, line_search);
            let model = AcOpfModel::new(&case).map_err(|e| e.to_string())?;
            let base = model.solve(&case, &options).map_err(|e| e.to_string())?;
            for (label, start) in [
                ("cold", OpfStart::Cold),
                ("warm", OpfStart::Warm(&base.solver)),
            ] {
                all &= batch.time(
                    label,
                    || model.solve_many(&cases, &options, start),
                    || model.par_solve_many(&cases, &options, start),
                    cases.len(),
                    |r| {
                        r.iter()
                            .filter(|s| s.as_ref().is_ok_and(|s| s.converged()))
                            .count()
                    },
                )?;
            }
        } else {
            let options = dc_options(mehrotra, line_search);
            let model = DcOpfModel::new(&case).map_err(|e| e.to_string())?;
            let base = model.solve(&case, &options).map_err(|e| e.to_string())?;
            for (label, start) in [
                ("cold", OpfStart::Cold),
                ("warm", OpfStart::Warm(&base.solver)),
            ] {
                all &= batch.time(
                    label,
                    || model.solve_many(&cases, &options, start),
                    || model.par_solve_many(&cases, &options, start),
                    cases.len(),
                    |r| {
                        r.iter()
                            .filter(|s| s.as_ref().is_ok_and(|s| s.converged()))
                            .count()
                    },
                )?;
            }
        }
    }
    Ok(all)
}
