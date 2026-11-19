//! Solves a deterministic sequence of load levels on one case with a kept OPF model, several
//! ways per step, for the warm-start survey that `scripts/benchmark-opf-sequence.sh` runs: `cold`
//! (MATPOWER's starting point), `primal` (the previous step's point, cold multipliers) and `warm`
//! at each warm shift (its point and multipliers, shifted by `warm_shift`), with `solve_from`'s
//! `OpfStart`.
//!
//! Usage: `opf_sequence --model dc|ac --solver mips|mehrotra --sequence sweep|walk
//! [--warm-shifts <δ>[,<δ>…]] <case.m>`, or `opf_sequence --header`. Each solver runs at
//! MATPOWER's OPF settings (`DcOpfOptions`, `AcOpfOptions`, and their `mehrotra()`), with its
//! own KKT default. `--warm-shifts` lists the warm shifts, one warm chain each (default: the
//! solver's own `warm_shift`). The cold and primal chains don't depend on the shift, so one run
//! covers every shift with one of each.
//!
//! **The sequences** scale every bus's load, active and reactive, from the case's:
//! - `sweep`: by 0.90, 0.91, …, 1.10, 21 steps
//! - `walk`: 25 steps from the case's own load, each bus's then moved by an independent uniform
//!   step in ±2% of its current value per step, from a SplitMix64 stream with a fixed seed
//!
//! **The chains.** Step 0 is solved cold, and starts every chain. From step 1 each chain starts
//! from its own previous solution, so a warm chain's drift shows. A chain whose previous solve
//! didn't converge starts its step cold, and the row says so (`start`).
//!
//! **Prints** one CSV line per step and chain, with the columns in `HEADER`: the step's outcome,
//! iterations and solve time (one run: a chain can't repeat a step), its objective and KKT
//! statistics, and for `primal` and `warm` the differences from the same step's cold solve and
//! the warm-start survey's verdict against it (`verdict`, the rule in
//! `docs/OptimizationChecklist.md`, as corrected on 2026-11-18):
//! - `regression`: the cold solve converged and this one didn't
//! - `rescue`: this one converged and the cold solve didn't; `both failed` if neither did
//! - where both converged, `miss: …` naming what O3's bar without `LAM_Q` misses: the objective
//!   within 1e-6 and `LAM_P` within 1e-3 of the cold solve's (relative, to at least 1), and each
//!   KKT residual within 100× MIPS-LU's cold solve of the step, or below 1e-6, as O3's bar
//!   compared them. A predictor-corrector run solves that reference too (`reference_outcome`);
//!   for MIPS it's the cold solve. Residuals aren't compared when the reference didn't converge.
//! - `local optimum`: on AC, a miss of the objective alone (with the prices that follow it), the
//!   residuals within the bar
//! - `acceptable` otherwise
//!
//! `LAM_Q`'s difference is reported (`lam_q_rel`) but not barred: two cold solvers already
//! differ in it by up to 9.4e-3 (MIPS and the predictor-corrector method on case9241pegase), so
//! a difference near 1e-3 is within what the AC tolerance determines.

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
                      kkt_factorizations,kkt_delays,kkt_factor_entries,reference_outcome,\
                      verdict";

/// O3's bar, but for `LAM_Q`
const OBJECTIVE: f64 = 1e-6;
const PRICE: f64 = 1e-3;
const RESIDUAL_FACTOR: f64 = 100.0;
const RESIDUAL_FLOOR: f64 = 1e-6;

/// The walk's seed and step
const WALK_SEED: u64 = 0x0005_EED0_FDA7;
const WALK_STEP: f64 = 0.02;
const WALK_STEPS: usize = 25;

#[derive(Clone, Copy, PartialEq)]
enum Model {
    Dc,
    Ac,
}

/// A chain: how each of its steps starts
#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Cold,
    Primal,
    /// Warm, at this shift
    Warm(f64),
}

impl Mode {
    fn name(self) -> &'static str {
        match self {
            Mode::Cold => "cold",
            Mode::Primal => "primal",
            Mode::Warm(_) => "warm",
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
    let (shifts, rest) = match rest {
        [flag, shifts, rest @ ..] if flag == "--warm-shifts" => {
            match shifts.split(',').map(str::parse::<f64>).collect() {
                Ok(shifts) => (Some(shifts), rest),
                Err(_) => return usage(),
            }
        }
        rest => (None, rest),
    };
    let [case] = rest else { return usage() };
    match run(Path::new(case), model, mehrotra, sequence, shifts) {
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
         [--warm-shifts <δ>[,<δ>…]] <case.m> | --header"
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

/// The differences of `step` from `cold` (objective, LAM_P, LAM_Q, as CSV fields), and the
///   survey's verdict, with residuals against `reference`, MIPS-LU's cold solve if it converged
fn compare(
    step: &Fields,
    cold: &Fields,
    reference: Option<&KktResiduals>,
    ac: bool,
) -> (String, String) {
    let objective = (step.objective - cold.objective).abs() / cold.objective.abs().max(1.0);
    let lam_p = relative_difference(&step.lam_p, &cold.lam_p);
    let lam_q = match (&step.lam_q, &cold.lam_q) {
        (Some(a), Some(b)) => format!("{:e}", relative_difference(a, b)),
        _ => String::new(),
    };
    let fields = format!("{objective:e},{lam_p:e},{lam_q}");
    let verdict = match (cold.converged(), step.converged()) {
        (true, false) => "regression".to_string(),
        (false, true) => "rescue".to_string(),
        (false, false) => "both failed".to_string(),
        (true, true) => {
            let mut misses = vec![];
            if objective > OBJECTIVE || objective.is_nan() {
                misses.push("objective");
            }
            if lam_p > PRICE || lam_p.is_nan() {
                misses.push("LAM_P");
            }
            let mut residual_misses = vec![];
            if let Some(b) = reference {
                let a = step.residuals;
                for (name, a, b) in [
                    ("stationarity", a.stationarity, b.stationarity),
                    ("feasibility", a.feasibility, b.feasibility),
                    ("complementarity", a.complementarity, b.complementarity),
                    ("dual feasibility", a.dual_feasibility, b.dual_feasibility),
                ] {
                    if a > RESIDUAL_FACTOR * b && a > RESIDUAL_FLOOR {
                        residual_misses.push(name);
                    }
                }
            }
            if ac && residual_misses.is_empty() && misses.first() == Some(&"objective") {
                "local optimum".to_string()
            } else if misses.is_empty() && residual_misses.is_empty() {
                "acceptable".to_string()
            } else {
                misses.extend(residual_misses);
                format!("miss: {}", misses.join("; "))
            }
        }
    };
    (fields, verdict)
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

/// `solver` with its `warm_shift` set to `shift`
fn shifted(solver: OpfSolver, shift: f64) -> OpfSolver {
    match solver {
        OpfSolver::Mips(mut options) => {
            options.warm_shift = shift;
            OpfSolver::Mips(options)
        }
        OpfSolver::Mehrotra(mut options) => {
            options.warm_shift = shift;
            OpfSolver::Mehrotra(options)
        }
        _ => unreachable!("a solver this example knows"),
    }
}

/// `solver`'s own `warm_shift`
fn own_shift(solver: OpfSolver) -> f64 {
    match solver {
        OpfSolver::Mips(options) => options.warm_shift,
        OpfSolver::Mehrotra(options) => options.warm_shift,
        _ => unreachable!("a solver this example knows"),
    }
}

impl Kept {
    /// MIPS-LU's cold solve of `case`, the residuals' reference, for a predictor-corrector run
    fn reference(&self, case: &Case) -> Option<Result<Fields, String>> {
        if matches!(self.solver(), OpfSolver::Mips(_)) {
            return None;
        }
        Some(match self {
            Kept::Dc(model, _) => model
                .solve(case, &DcOpfOptions::default())
                .map(Fields::dc)
                .map_err(|err| err.to_string()),
            Kept::Ac(model, _) => model
                .solve(case, &AcOpfOptions::default())
                .map(Fields::ac)
                .map_err(|err| err.to_string()),
        })
    }

    fn solver(&self) -> OpfSolver {
        match self {
            Kept::Dc(_, options) => options.solver,
            Kept::Ac(_, options) => options.solver,
        }
    }

    /// Solves `case` from `start`, a warm start at `shift`
    fn solve(
        &self,
        case: &Case,
        start: OpfStart<'_>,
        shift: Option<f64>,
    ) -> Result<Fields, String> {
        match self {
            Kept::Dc(model, options) => {
                let mut options = *options;
                if let Some(shift) = shift {
                    options.solver = shifted(options.solver, shift);
                }
                model
                    .solve_from(case, &options, start)
                    .map(Fields::dc)
                    .map_err(|err| err.to_string())
            }
            Kept::Ac(model, options) => {
                let mut options = *options;
                if let Some(shift) = shift {
                    options.solver = shifted(options.solver, shift);
                }
                model
                    .solve_from(case, &options, start)
                    .map(Fields::ac)
                    .map_err(|err| err.to_string())
            }
        }
    }
}

fn run(
    path: &Path,
    model: Model,
    mehrotra: bool,
    sequence: &str,
    shifts: Option<Vec<f64>>,
) -> Result<(), String> {
    let name = path
        .file_stem()
        .and_then(|stem| stem.to_str())
        .ok_or("the case path has no file name")?;
    let case = Case::from_path(path).map_err(|err| err.to_string())?;
    let kept = match model {
        Model::Dc => {
            let options = if mehrotra {
                DcOpfOptions::mehrotra()
            } else {
                DcOpfOptions::default()
            };
            let model = DcOpfModel::new(&case).map_err(|err| err.to_string())?;
            Kept::Dc(Box::new(model), options)
        }
        Model::Ac => {
            let options = if mehrotra {
                AcOpfOptions::mehrotra()
            } else {
                AcOpfOptions::default()
            };
            let model = AcOpfModel::new(&case).map_err(|err| err.to_string())?;
            Kept::Ac(Box::new(model), options)
        }
    };
    let shifts = shifts.unwrap_or_else(|| vec![own_shift(kept.solver())]);
    let mut chains = vec![Mode::Cold, Mode::Primal];
    chains.extend(shifts.iter().map(|&shift| Mode::Warm(shift)));
    let prefix = format!(
        "{name},{},{},{sequence}",
        if model == Model::Dc { "dc" } else { "ac" },
        if mehrotra { "mehrotra" } else { "mips" },
    );
    // Each chain's previous solution, if it converged
    let mut previous: Vec<Option<Solution>> = vec![None; chains.len()];
    for (step, (factors, mean)) in steps(sequence, case.bus.len()).iter().enumerate() {
        let changed = loaded(&case, factors);
        let modes = if step == 0 { &chains[..1] } else { &chains[..] };
        let mut results: Vec<(&str, f64, Result<Fields, String>)> = vec![];
        for (&mode, chain) in modes.iter().zip(&previous) {
            let (start, from) = match (mode, chain) {
                (Mode::Primal, Some(p)) => (OpfStart::Primal(p), "previous"),
                (Mode::Warm(_), Some(p)) => (OpfStart::Warm(p), "previous"),
                _ => (OpfStart::Cold, "cold"),
            };
            let shift = match mode {
                Mode::Warm(shift) => Some(shift),
                _ => None,
            };
            let begin = Instant::now();
            let result = kept.solve(&changed, start, shift);
            let ms = begin.elapsed().as_secs_f64() * 1e3;
            results.push((from, ms, result));
        }
        let cold = results[0].2.as_ref().ok();
        let reference_solve = if step == 0 {
            None
        } else {
            kept.reference(&changed)
        };
        let (reference, reference_outcome) = match &reference_solve {
            None => (cold, String::new()),
            Some(Ok(fields)) => (Some(fields), outcome(fields.status).to_string()),
            Some(Err(err)) => (None, format!("error: {}", err.replace(',', ";"))),
        };
        let reference = reference
            .filter(|fields| fields.converged())
            .map(|fields| &fields.residuals);
        for (&mode, (from, ms, result)) in modes.iter().zip(&results) {
            let shift = match mode {
                Mode::Warm(shift) => format!("{shift:e}"),
                _ => String::new(),
            };
            let line = match result {
                Ok(fields) => {
                    let (differences, verdict) = match (mode, cold) {
                        (Mode::Cold, _) | (_, None) => (",,".to_string(), String::new()),
                        (_, Some(cold)) => compare(fields, cold, reference, model == Model::Ac),
                    };
                    let reference_outcome = match mode {
                        Mode::Cold => "",
                        _ => reference_outcome.as_str(),
                    };
                    let r = fields.residuals;
                    format!(
                        "{prefix},{shift},{step},{mean:.6},{},{from},{},{},{ms:.3},{:e},\
                         {differences},{:e},{:e},{:e},{:e},{},{},{},{reference_outcome},{verdict}",
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
                Err(err) => {
                    // The objective through reference_outcome left empty, then the verdict
                    let verdict = match (mode, cold) {
                        (Mode::Cold, _) => "",
                        (_, Some(cold)) if cold.converged() => "regression",
                        _ => "both failed",
                    };
                    format!(
                        "{prefix},{shift},{step},{mean:.6},{},{from},error: {},,{ms:.3}{}{verdict}",
                        mode.name(),
                        err.replace(',', ";"),
                        ",".repeat(13),
                    )
                }
            };
            println!("{line}");
        }
        // The chains move on, each to its step's solution if it converged; step 0's cold
        //   solve starts every chain
        let solutions: Vec<Option<Solution>> = results
            .into_iter()
            .map(|(_, _, result)| result.ok().filter(Fields::converged).map(|f| f.solver))
            .collect();
        if step == 0 {
            previous = vec![solutions[0].clone(); chains.len()];
        } else {
            previous = solutions;
        }
    }
    Ok(())
}
