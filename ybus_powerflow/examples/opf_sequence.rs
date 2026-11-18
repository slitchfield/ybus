//! Solves a deterministic sequence of load levels on one case with a kept OPF model, three ways
//! per step, for the warm-start survey that `scripts/benchmark-opf-sequence.sh` runs: `cold`
//! (MATPOWER's starting point), `primal` (the previous step's point, cold multipliers) and `warm`
//! (its point and multipliers, shifted by `warm_shift`), with `solve_from`'s `OpfStart`.
//!
//! Usage: `opf_sequence --model dc|ac --solver mips|mehrotra --sequence sweep|walk
//! [--warm-shift <δ>] <case.m>`, or `opf_sequence --header`. Each solver runs at MATPOWER's OPF
//! settings (`DcOpfOptions`, `AcOpfOptions`, and their `mehrotra()`), with its own KKT default,
//! and `--warm-shift` sets its `warm_shift` (default: the solver's).
//!
//! **The sequences** scale every bus's load, active and reactive, from the case's:
//! - `sweep`: by 0.90, 0.91, …, 1.10, 21 steps
//! - `walk`: 25 steps from the case's own load, each bus's then moved by an independent uniform
//!   step in ±2% of its current value per step, from a SplitMix64 stream with a fixed seed, as a
//!   day's hours
//!
//! **The chains.** Step 0 is solved cold, and starts all three chains. From step 1 each mode
//! starts from its own chain's previous solution, so a warm chain's drift shows. A chain whose
//! previous solve didn't converge starts its step cold, and the row says so (`start`).
//!
//! **Prints** one CSV line per step and mode, with the columns in `HEADER`: the step's outcome,
//! iterations and solve time (one run: a chain can't repeat a step), its objective and KKT
//! statistics, and for `primal` and `warm` the differences from the same step's cold solve and
//! O3's bar against it (`acceptable`). The bar: the same outcome, and where the cold solve
//! converged, the objective within 1e-6 (relative, to at least 1), the prices `LAM_P` and, for
//! AC, `LAM_Q` within 1e-3 (relative, to at least 1), and each KKT residual within 100× the
//! cold solve's or below 1e-6.

use std::env;
use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use ybus_powerflow::case::Case;
use ybus_powerflow::opf::ac::{AcOpfModel, AcOpfOptions, AcOpfSolution};
use ybus_powerflow::opf::dc::{DcOpfModel, DcOpfOptions, DcOpfSolution};
use ybus_powerflow::opf::{KktResiduals, KktStatistics, OpfSolver, OpfStart, Solution, Status};

const HEADER: &str = "case,model,solver,sequence,warm_shift,step,load_factor,mode,start,outcome,\
                      iterations,solve_ms,objective,objective_rel,lam_p_rel,lam_q_rel,\
                      kkt_stationarity,kkt_feasibility,kkt_complementarity,kkt_dual_feasibility,\
                      kkt_factorizations,kkt_delays,kkt_factor_entries,acceptable";

/// O3's bar
const OBJECTIVE: f64 = 1e-6;
const PRICE: f64 = 1e-3;
const RESIDUAL_FACTOR: f64 = 100.0;
const RESIDUAL_FLOOR: f64 = 1e-6;

/// The walk's seed and step
const WALK_SEED: u64 = 0x0005_EED0_FDA7;
const WALK_STEP: f64 = 0.02;
const WALK_STEPS: usize = 25;

const MODES: [Mode; 3] = [Mode::Cold, Mode::Primal, Mode::Warm];

#[derive(Clone, Copy, PartialEq)]
enum Model {
    Dc,
    Ac,
}

#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Cold,
    Primal,
    Warm,
}

impl Mode {
    fn name(self) -> &'static str {
        match self {
            Mode::Cold => "cold",
            Mode::Primal => "primal",
            Mode::Warm => "warm",
        }
    }
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
    let (mehrotra, rest) = match rest {
        [flag, solver, rest @ ..] if flag == "--solver" && solver == "mips" => (false, rest),
        [flag, solver, rest @ ..] if flag == "--solver" && solver == "mehrotra" => (true, rest),
        _ => return usage(),
    };
    let (sequence, rest) = match rest {
        [flag, sequence, rest @ ..]
            if flag == "--sequence" && (sequence == "sweep" || sequence == "walk") =>
        {
            (sequence.as_str(), rest)
        }
        _ => return usage(),
    };
    let (shift, rest) = match rest {
        [flag, shift, rest @ ..] if flag == "--warm-shift" => match shift.parse::<f64>() {
            Ok(shift) => (Some(shift), rest),
            Err(_) => return usage(),
        },
        rest => (None, rest),
    };
    let [case] = rest else { return usage() };
    match run(Path::new(case), model, mehrotra, sequence, shift) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{case}: {err}");
            ExitCode::FAILURE
        }
    }
}

fn usage() -> ExitCode {
    eprintln!(
        "usage: opf_sequence --model dc|ac --solver mips|mehrotra --sequence sweep|walk \
         [--warm-shift <δ>] <case.m> | --header"
    );
    ExitCode::from(2)
}

/// SplitMix64 (Steele, Lea and Flood, 2014), for a walk that is the same on every machine
struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [−1, 1)
    fn symmetric(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    }
}

/// The sequence's steps: each bus's load factor, and the mean factor for the CSV
fn steps(sequence: &str, buses: usize) -> Vec<(Vec<f64>, f64)> {
    match sequence {
        "sweep" => (0..=20)
            .map(|k| {
                let factor = 0.90 + 0.01 * k as f64;
                (vec![factor; buses], factor)
            })
            .collect(),
        _ => {
            let mut random = SplitMix64(WALK_SEED);
            let mut factors = vec![1.0; buses];
            let mut steps = vec![(factors.clone(), 1.0)];
            for _ in 1..WALK_STEPS {
                for factor in &mut factors {
                    *factor *= 1.0 + WALK_STEP * random.symmetric();
                }
                let mean = factors.iter().sum::<f64>() / buses as f64;
                steps.push((factors.clone(), mean));
            }
            steps
        }
    }
}

fn loaded(case: &Case, factors: &[f64]) -> Case {
    let mut case = case.clone();
    for (bus, &factor) in case.bus.iter_mut().zip(factors) {
        bus.pd_mw *= factor;
        bus.qd_mvar *= factor;
    }
    case
}

/// One solve's results, from either model
struct Fields {
    status: Status,
    iterations: usize,
    objective: f64,
    lam_p: Vec<f64>,
    lam_q: Option<Vec<f64>>,
    residuals: KktResiduals,
    kkt: KktStatistics,
    solver: Solution,
}

impl Fields {
    fn dc(s: DcOpfSolution) -> Fields {
        Fields {
            status: s.status,
            iterations: s.iterations,
            objective: s.objective,
            lam_p: s.lam_p,
            lam_q: None,
            residuals: s.residuals,
            kkt: s.solver.kkt,
            solver: s.solver,
        }
    }

    fn ac(s: AcOpfSolution) -> Fields {
        Fields {
            status: s.status,
            iterations: s.iterations,
            objective: s.objective,
            lam_p: s.lam_p,
            lam_q: Some(s.lam_q),
            residuals: s.residuals,
            kkt: s.solver.kkt,
            solver: s.solver,
        }
    }

    fn converged(&self) -> bool {
        self.status == Status::Converged
    }
}

/// The largest `|a − e| / max(|e|, 1)`, NaN if any is or the lengths differ, as the surveys
///   compare prices
fn relative_difference(actual: &[f64], expected: &[f64]) -> f64 {
    if actual.len() != expected.len() {
        return f64::NAN;
    }
    actual
        .iter()
        .zip(expected)
        .map(|(a, e)| (a - e).abs() / e.abs().max(1.0))
        .fold(0.0, |largest, d| {
            if d.is_nan() || d > largest {
                d
            } else {
                largest
            }
        })
}

/// The differences of `step` from `cold` (objective, LAM_P, LAM_Q, as CSV fields) and O3's bar
fn compare(step: &Fields, cold: &Fields) -> (String, bool) {
    let objective = (step.objective - cold.objective).abs() / cold.objective.abs().max(1.0);
    let lam_p = relative_difference(&step.lam_p, &cold.lam_p);
    let lam_q = match (&step.lam_q, &cold.lam_q) {
        (Some(a), Some(b)) => Some(relative_difference(a, b)),
        _ => None,
    };
    let fields = format!(
        "{objective:e},{lam_p:e},{}",
        lam_q.map(|d| format!("{d:e}")).unwrap_or_default()
    );
    let within = |d: f64, bar: f64| d <= bar;
    let (a, b) = (step.residuals, cold.residuals);
    let residuals = [
        (a.stationarity, b.stationarity),
        (a.feasibility, b.feasibility),
        (a.complementarity, b.complementarity),
        (a.dual_feasibility, b.dual_feasibility),
    ]
    .iter()
    .all(|&(a, b)| a <= RESIDUAL_FACTOR * b || a <= RESIDUAL_FLOOR);
    let acceptable = step.converged() == cold.converged()
        && (!cold.converged()
            || (within(objective, OBJECTIVE)
                && within(lam_p, PRICE)
                && lam_q.is_none_or(|d| within(d, PRICE))
                && residuals));
    (fields, acceptable)
}

fn outcome(status: Status) -> &'static str {
    match status {
        Status::Converged => "converged",
        Status::MaxIterations => "max iterations",
        Status::NumericalFailure => "numerical failure",
    }
}

/// A kept model of either kind, with its options (boxed: a model is large)
enum Kept {
    Dc(Box<DcOpfModel>, DcOpfOptions),
    Ac(Box<AcOpfModel>, AcOpfOptions),
}

impl Kept {
    fn solve(&self, case: &Case, start: OpfStart<'_>) -> Result<Fields, String> {
        match self {
            Kept::Dc(model, options) => model
                .solve_from(case, options, start)
                .map(Fields::dc)
                .map_err(|err| err.to_string()),
            Kept::Ac(model, options) => model
                .solve_from(case, options, start)
                .map(Fields::ac)
                .map_err(|err| err.to_string()),
        }
    }
}

fn with_shift(solver: OpfSolver, shift: Option<f64>) -> (OpfSolver, f64) {
    match solver {
        OpfSolver::Mips(mut options) => {
            options.warm_shift = shift.unwrap_or(options.warm_shift);
            (OpfSolver::Mips(options), options.warm_shift)
        }
        OpfSolver::Mehrotra(mut options) => {
            options.warm_shift = shift.unwrap_or(options.warm_shift);
            (OpfSolver::Mehrotra(options), options.warm_shift)
        }
        _ => unreachable!("a solver this example knows"),
    }
}

fn run(
    path: &Path,
    model: Model,
    mehrotra: bool,
    sequence: &str,
    shift: Option<f64>,
) -> Result<(), String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let case = Case::from_path(path).map_err(|err| err.to_string())?;
    let (kept, shift) = match model {
        Model::Dc => {
            let mut options = if mehrotra {
                DcOpfOptions::mehrotra()
            } else {
                DcOpfOptions::default()
            };
            let (solver, shift) = with_shift(options.solver, shift);
            options.solver = solver;
            let model = DcOpfModel::new(&case).map_err(|err| err.to_string())?;
            (Kept::Dc(Box::new(model), options), shift)
        }
        Model::Ac => {
            let mut options = if mehrotra {
                AcOpfOptions::mehrotra()
            } else {
                AcOpfOptions::default()
            };
            let (solver, shift) = with_shift(options.solver, shift);
            options.solver = solver;
            let model = AcOpfModel::new(&case).map_err(|err| err.to_string())?;
            (Kept::Ac(Box::new(model), options), shift)
        }
    };
    let prefix = format!(
        "{name},{},{},{sequence},{shift:e}",
        if model == Model::Dc { "dc" } else { "ac" },
        if mehrotra { "mehrotra" } else { "mips" },
    );
    // Each chain's previous solution, if it converged
    let mut previous: [Option<Solution>; 3] = [None, None, None];
    for (step, (factors, mean)) in steps(sequence, case.bus.len()).iter().enumerate() {
        let changed = loaded(&case, factors);
        let modes: &[Mode] = if step == 0 { &MODES[..1] } else { &MODES };
        let mut results: Vec<(Mode, &str, f64, Result<Fields, String>)> = vec![];
        for &mode in modes {
            let chain = &previous[mode as usize];
            let (start, from) = match (mode, chain) {
                (Mode::Primal, Some(p)) => (OpfStart::Primal(p), "previous"),
                (Mode::Warm, Some(p)) => (OpfStart::Warm(p), "previous"),
                _ => (OpfStart::Cold, "cold"),
            };
            let begin = Instant::now();
            let result = kept.solve(&changed, start);
            let ms = begin.elapsed().as_secs_f64() * 1e3;
            results.push((mode, from, ms, result));
        }
        let cold = results[0].3.as_ref().ok();
        let mut lines = vec![];
        for (mode, from, ms, result) in &results {
            let line = match result {
                Ok(fields) => {
                    let (differences, acceptable) = match (mode, cold) {
                        (Mode::Cold, _) | (_, None) => (",,".to_string(), String::new()),
                        (_, Some(cold)) => {
                            let (d, ok) = compare(fields, cold);
                            (d, ok.to_string())
                        }
                    };
                    let r = fields.residuals;
                    format!(
                        "{prefix},{step},{mean:.6},{},{from},{},{},{ms:.3},{:e},{differences},\
                         {:e},{:e},{:e},{:e},{},{},{},{acceptable}",
                        mode.name(),
                        outcome(fields.status),
                        fields.iterations,
                        fields.objective,
                        r.stationarity,
                        r.feasibility,
                        r.complementarity,
                        r.dual_feasibility,
                        fields.kkt.factorizations,
                        fields.kkt.delays,
                        fields.kkt.factor_entries,
                    )
                }
                Err(err) => format!(
                    "{prefix},{step},{mean:.6},{},{from},error: {},,{ms:.3},,,,,,,,,,,,",
                    mode.name(),
                    err.replace(',', ";")
                ),
            };
            lines.push(line);
        }
        for line in lines {
            println!("{line}");
        }
        // The chains move on: cold's to nothing (it always starts cold), the others to their
        //   step's solution if it converged
        let solutions: Vec<(Mode, Option<Solution>)> = results
            .into_iter()
            .map(|(mode, _, _, result)| {
                let solution = result.ok().filter(Fields::converged).map(|f| f.solver);
                (mode, solution)
            })
            .collect();
        for (mode, solution) in solutions {
            if step == 0 {
                // Step 0's cold solve starts every chain
                previous = [None, solution.clone(), solution];
                break;
            }
            previous[mode as usize] = solution;
        }
    }
    Ok(())
}
