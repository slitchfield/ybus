//! MATPOWER's DC and AC power flow results, as `scripts/generate-reference.sh` writes them: committed
//! for case14, case30 and case118 in `tests/data/matpower-reference/`, whose README describes
//! the files, and generated locally for any other case, whose results can't be committed.

use std::fs;
use std::path::{Path, PathBuf};

use powerflow::case::Case;
use powerflow::dc::DcSolution;

/// The committed reference results.
pub fn committed_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join("matpower-reference")
}

/// One case's `rundcpf` results, each table in the case's row order.
#[derive(Debug, Clone, PartialEq)]
pub struct Reference {
    /// Each bus's number and voltage angle (degrees).
    pub bus: Vec<(usize, f64)>,
    /// Each branch's from and to bus numbers, status, and real power flow at the from end (MW).
    pub branch: Vec<(usize, usize, bool, f64)>,
    /// Each generator's bus number and real power output (MW).
    pub gen: Vec<(usize, f64)>,
}

/// The largest absolute differences between a solution and a reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Deviation {
    /// Voltage angles (degrees).
    pub va_deg: f64,
    /// Branch flows at the from end (MW).
    pub pf_mw: f64,
    /// Generator outputs (MW).
    pub pg_mw: f64,
}

impl Deviation {
    /// Whether every difference is within `tolerance`'s. False if any is NaN.
    pub fn within(&self, tolerance: &Deviation) -> bool {
        self.va_deg <= tolerance.va_deg
            && self.pf_mw <= tolerance.pf_mw
            && self.pg_mw <= tolerance.pg_mw
    }
}

/// Reads `<case>_dcpf_bus.csv`, `_branch.csv` and `_gen.csv` from `dir`, checking their headers
/// and row numbers.
pub fn read(dir: &Path, case: &str) -> Result<Reference, String> {
    let bus = table(dir, case, "dcpf", "bus", "bus,va_deg")?;
    let branch = table(dir, case, "dcpf", "branch", "row,from,to,status,pf_mw")?;
    let gen = table(dir, case, "dcpf", "gen", "row,bus,pg_mw")?;
    // Bus numbers and statuses are written with %d, so they convert exactly
    Ok(Reference {
        bus: bus.iter().map(|f| (f[0] as usize, f[1])).collect(),
        branch: branch
            .iter()
            .map(|f| (f[0] as usize, f[1] as usize, f[2] != 0.0, f[3]))
            .collect(),
        gen: gen.iter().map(|f| (f[0] as usize, f[1])).collect(),
    })
}

/// The rows of `<dir>/<case>_<run>_<kind>.csv` as numbers, without the leading `row` column
/// where the header has one, after checking that it counts from 1.
fn table(
    dir: &Path,
    case: &str,
    run: &str,
    kind: &str,
    header: &str,
) -> Result<Vec<Vec<f64>>, String> {
    let path = dir.join(format!("{case}_{run}_{kind}.csv"));
    let text = fs::read_to_string(&path).map_err(|err| format!("{}: {err}", path.display()))?;
    let mut lines = text.lines();
    if lines.next() != Some(header) {
        return Err(format!("{}: the header isn't {header}", path.display()));
    }
    let width = header.split(',').count();
    let numbered = header.starts_with("row,");
    lines
        .enumerate()
        .map(|(index, line)| {
            let malformed = || format!("{}: line {} is malformed", path.display(), index + 2);
            let fields: Vec<f64> = line
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()
                .map_err(|_| malformed())?;
            if fields.len() != width || (numbered && fields[0] != (index + 1) as f64) {
                return Err(malformed());
            }
            Ok(fields[usize::from(numbered)..].to_vec())
        })
        .collect()
}

/// One case's `runpf` results (Newton's method, MATPOWER's defaults), each table in the case's
/// row order.
#[derive(Debug, Clone, PartialEq)]
pub struct AcReference {
    /// Each bus's number, voltage magnitude (p.u.) and angle (degrees).
    pub bus: Vec<(usize, f64, f64)>,
    /// Each branch's from and to bus numbers, status, and real and reactive power flows at the
    /// from end and at the to end (MW, MVAr).
    pub branch: Vec<(usize, usize, bool, [f64; 4])>,
    /// Each generator's bus number and real and reactive power output (MW, MVAr).
    pub gen: Vec<(usize, f64, f64)>,
    /// How many Newton iterations `runpf` took.
    pub iterations: usize,
    /// Whether `runpf` converged.
    pub success: bool,
}

/// Reads `<case>_acpf_bus.csv`, `_branch.csv`, `_gen.csv` and `_summary.csv` from `dir`,
/// checking their headers and row numbers.
pub fn read_ac(dir: &Path, case: &str) -> Result<AcReference, String> {
    let bus = table(dir, case, "acpf", "bus", "bus,vm_pu,va_deg")?;
    let branch = table(
        dir,
        case,
        "acpf",
        "branch",
        "row,from,to,status,pf_mw,qf_mvar,pt_mw,qt_mvar",
    )?;
    let gen = table(dir, case, "acpf", "gen", "row,bus,pg_mw,qg_mvar")?;
    let summary = table(dir, case, "acpf", "summary", "iterations,success")?;
    let [run] = summary.as_slice() else {
        return Err(format!("{case}: the AC summary should have one row"));
    };
    Ok(AcReference {
        bus: bus.iter().map(|f| (f[0] as usize, f[1], f[2])).collect(),
        branch: branch
            .iter()
            .map(|f| {
                (
                    f[0] as usize,
                    f[1] as usize,
                    f[2] != 0.0,
                    [f[3], f[4], f[5], f[6]],
                )
            })
            .collect(),
        gen: gen.iter().map(|f| (f[0] as usize, f[1], f[2])).collect(),
        iterations: run[0] as usize,
        success: run[1] != 0.0,
    })
}

/// The largest absolute differences between `solution` and `reference`, after checking that
/// `reference` describes `case`: the same buses, branches (ends and status) and generators, in
/// the same order.
pub fn deviation(
    case: &Case,
    solution: &DcSolution,
    reference: &Reference,
) -> Result<Deviation, String> {
    let same_buses = case
        .bus
        .iter()
        .map(|bus| bus.bus_label)
        .eq(reference.bus.iter().map(|bus| bus.0));
    let same_branches = case
        .branch
        .iter()
        .map(|branch| (branch.from_label, branch.to_label, branch.status))
        .eq(reference
            .branch
            .iter()
            .map(|branch| (branch.0, branch.1, branch.2)));
    let same_gens = case
        .generator
        .iter()
        .map(|generator| generator.bus_label)
        .eq(reference.gen.iter().map(|generator| generator.0));
    if !(same_buses && same_branches && same_gens) {
        return Err("the reference's tables don't match the case's".to_string());
    }
    Ok(Deviation {
        va_deg: largest_difference(&solution.va_deg, reference.bus.iter().map(|bus| bus.1)),
        pf_mw: largest_difference(&solution.pf_mw, reference.branch.iter().map(|b| b.3)),
        pg_mw: largest_difference(&solution.pg_mw, reference.gen.iter().map(|g| g.1)),
    })
}

/// The largest `|a − e|`, or NaN if any difference is NaN, which `f64::max` would skip.
fn largest_difference(actual: &[f64], expected: impl Iterator<Item = f64>) -> f64 {
    actual
        .iter()
        .zip(expected)
        .map(|(a, e)| (a - e).abs())
        .fold(0.0, |largest, difference| {
            if difference.is_nan() || difference > largest {
                difference
            } else {
                largest
            }
        })
}
