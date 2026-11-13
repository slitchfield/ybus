//! The power system case model, read from MATPOWER case files.
//!
//! A [`Case`] holds the data exactly as the case file gives it: powers in MW and MVAr, angles in
//! degrees, impedances in per unit on the system base, and buses identified by their *labels*
//! (the bus numbers in the file), not by position. Solver-ready views (consecutive indices,
//! per-unit powers, radians) are derived from it when needed, as MATPOWER's `ext2int` does,
//! so the case itself stays a faithful copy of the file.
//!
//! See `caseformat.m` in MATPOWER for the format's definition, and [`Case::parse`] for what
//! the reader supports.

use std::borrow::Cow;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use crate::powerflow_error::PowerflowError;

mod matpower;

/// A power system case: buses, generators, branches and optional cost data.
#[derive(Debug, Clone, PartialEq)]
pub struct Case {
    /// The system MVA base that per-unit values are relative to.
    pub base_mva: f64,
    /// One entry per bus, in file order.
    pub bus: Vec<Bus>,
    /// One entry per generator, in file order.
    pub generator: Vec<Gen>,
    /// One entry per branch (line or transformer), in file order.
    pub branch: Vec<Branch>,
    /// Generator costs, if the file has them: one per generator for real power, optionally
    /// followed by one per generator for reactive power.
    pub gencost: Option<Vec<GenCost>>,
    /// Bus names, if the file has them, parallel to `bus`.
    pub busnames: Option<Vec<String>>,
}

/// A bus (node) of the network.
#[derive(Clone, Debug, PartialEq)]
pub struct Bus {
    /// The bus number from the file: a label other rows refer to, not a position.
    pub bus_label: usize,
    pub bus_type: BusType,
    /// Real power demand (MW).
    pub pd_mw: f64,
    /// Reactive power demand (MVAr).
    pub qd_mvar: f64,
    /// Shunt conductance, as the MW it draws at 1.0 p.u. voltage.
    pub gs_mw: f64,
    /// Shunt susceptance, as the MVAr it injects at 1.0 p.u. voltage.
    pub bs_mvar: f64,
    pub area: usize,
    /// Voltage magnitude (p.u.).
    pub vm_pu: f64,
    /// Voltage angle (degrees).
    pub va_deg: f64,
    /// Base voltage (kV).
    pub base_kv: f64,
    /// Loss zone.
    pub zone: usize,
    /// Maximum voltage magnitude (p.u.).
    pub vmax_pu: f64,
    /// Minimum voltage magnitude (p.u.).
    pub vmin_pu: f64,
}

/// A bus's role in the power flow, from MATPOWER's codes 1 to 4.
#[derive(Clone, Debug, PartialEq, Eq, Copy)]
pub enum BusType {
    /// Code 1: real and reactive power are specified.
    PQ,
    /// Code 2: real power and voltage magnitude are specified.
    PV,
    /// Code 3: the reference (slack) bus, whose angle is the reference.
    Ref,
    /// Code 4: not connected to the network.
    Isolated,
}

impl BusType {
    /// The bus type for a MATPOWER code, or `None` if the code isn't 1 to 4.
    pub fn from_code(code: f64) -> Option<BusType> {
        match code {
            1.0 => Some(BusType::PQ),
            2.0 => Some(BusType::PV),
            3.0 => Some(BusType::Ref),
            4.0 => Some(BusType::Isolated),
            _ => None,
        }
    }
}

/// A generator, with the 21 columns of the version 2 format.
#[derive(Clone, Debug, PartialEq)]
pub struct Gen {
    /// The label of the bus it's connected to.
    pub bus_label: usize,
    /// Real power output (MW).
    pub pg_mw: f64,
    /// Reactive power output (MVAr).
    pub qg_mvar: f64,
    /// Maximum reactive power output (MVAr).
    pub qmax_mvar: f64,
    /// Minimum reactive power output (MVAr).
    pub qmin_mvar: f64,
    /// Voltage magnitude setpoint (p.u.).
    pub vg_pu: f64,
    /// The machine's MVA base.
    pub mbase_mva: f64,
    /// Whether it's in service. The file's status counts as in service when it's above 0.
    pub status: bool,
    /// Maximum real power output (MW).
    pub pmax_mw: f64,
    /// Minimum real power output (MW).
    pub pmin_mw: f64,
    /// Lower real power output of the PQ capability curve (MW).
    pub pc1_mw: f64,
    /// Upper real power output of the PQ capability curve (MW).
    pub pc2_mw: f64,
    /// Minimum reactive power output at `pc1_mw` (MVAr).
    pub qc1min_mvar: f64,
    /// Maximum reactive power output at `pc1_mw` (MVAr).
    pub qc1max_mvar: f64,
    /// Minimum reactive power output at `pc2_mw` (MVAr).
    pub qc2min_mvar: f64,
    /// Maximum reactive power output at `pc2_mw` (MVAr).
    pub qc2max_mvar: f64,
    /// Ramp rate for load following/AGC (MW/min).
    pub ramp_agc_mw_per_min: f64,
    /// Ramp rate for 10 minute reserves (MW).
    pub ramp_10_mw: f64,
    /// Ramp rate for 30 minute reserves (MW).
    pub ramp_30_mw: f64,
    /// Ramp rate for reactive power, 2 second timescale (MVAr/min).
    pub ramp_q_mvar_per_min: f64,
    /// Area participation factor.
    pub apf: f64,
}

/// A branch: a transmission line, or a transformer when `tap_ratio` or `shift_deg` is set.
#[derive(Clone, Debug, PartialEq)]
pub struct Branch {
    /// The label of the "from" bus, where a transformer's tap is.
    pub from_label: usize,
    /// The label of the "to" bus.
    pub to_label: usize,
    /// Resistance (p.u.).
    pub r: f64,
    /// Reactance (p.u.).
    pub x: f64,
    /// Total line charging susceptance (p.u.).
    pub b: f64,
    /// Long-term MVA rating; 0 means unlimited.
    pub rate_a_mva: f64,
    /// Short-term MVA rating; 0 means unlimited.
    pub rate_b_mva: f64,
    /// Emergency MVA rating; 0 means unlimited.
    pub rate_c_mva: f64,
    /// Transformer off-nominal turns ratio as in the file, where 0 means a line (equivalent to
    /// 1). See [`Branch::effective_tap`].
    pub tap_ratio: f64,
    /// Transformer phase shift angle (degrees); positive is a delay.
    pub shift_deg: f64,
    /// Whether it's in service. The file's status counts as in service when it isn't 0.
    pub status: bool,
    /// Minimum angle difference, angle(Vf) − angle(Vt) (degrees).
    pub angmin_deg: f64,
    /// Maximum angle difference, angle(Vf) − angle(Vt) (degrees).
    pub angmax_deg: f64,
}

impl Branch {
    /// The turns ratio to compute with: `tap_ratio`, except that the file's 0 (a line) is 1.
    pub fn effective_tap(&self) -> f64 {
        if self.tap_ratio == 0.0 {
            1.0
        } else {
            self.tap_ratio
        }
    }
}

/// A generator cost function, in $/hr of output in MW (or MVAr for reactive power costs).
#[derive(Clone, Debug, PartialEq)]
pub enum GenCost {
    /// Model 2: `cost(p) = coeffs[0]·p^n + … + coeffs[n]`, highest order first.
    Polynomial {
        /// Startup cost (US dollars).
        startup: f64,
        /// Shutdown cost (US dollars).
        shutdown: f64,
        coeffs: Vec<f64>,
    },
    /// Model 1: a piecewise linear function through `(p, cost)` points, in increasing `p`.
    PiecewiseLinear {
        /// Startup cost (US dollars).
        startup: f64,
        /// Shutdown cost (US dollars).
        shutdown: f64,
        points: Vec<(f64, f64)>,
    },
}

impl Case {
    /// Reads and parses a MATPOWER case file; see [`Case::from_bytes`] for how its encoding is
    /// handled and [`Case::parse`] for what's supported.
    ///
    /// Returns `PowerflowError::NoCaseFile` if the file can't be read.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Case, PowerflowError> {
        let path = path.as_ref();
        let bytes = fs::read(path).map_err(|err| PowerflowError::NoCaseFile {
            bad_path: path.to_path_buf(),
            reason: err.to_string(),
        })?;
        Case::from_bytes(&bytes)
    }

    /// Parses the raw bytes of a MATPOWER case file; see [`Case::parse`].
    ///
    /// The bytes are decoded as UTF-8, ignoring a leading byte order mark. Bytes that aren't
    /// valid UTF-8 are decoded as Latin-1 (ISO 8859-1) instead, where every byte is a character:
    /// MATLAB has long saved files in the system's native encoding, often Windows-1252 or
    /// Latin-1, and the syntax itself is ASCII either way. Only text in comments and bus names
    /// can be affected, and a Windows-1252 byte outside Latin-1 becomes the corresponding
    /// control character rather than an error.
    pub fn from_bytes(bytes: &[u8]) -> Result<Case, PowerflowError> {
        let bytes = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
        let text = match std::str::from_utf8(bytes) {
            Ok(text) => Cow::Borrowed(text),
            Err(_) => Cow::Owned(bytes.iter().map(|&byte| char::from(byte)).collect()),
        };
        Case::parse(&text)
    }

    /// Parses the text of a MATPOWER case file (version 2 format) and validates it.
    ///
    /// The fields the model reads must be assigned literals: matrices of plain numbers, cell
    /// arrays of strings, strings and numbers. Assignments to any other field are skipped,
    /// whatever they hold, and MATLAB and Octave comments, double-quoted strings, several
    /// statements on a line and a struct named by the `function` line are all accepted. Some
    /// MATPOWER cases compute or rescale their data with code after the literals; those are
    /// rejected rather than read with the wrong values.
    ///
    /// Returns `PowerflowError::Parse` (with a line number) for malformed or unsupported
    /// text, `MissingField` or `UnsupportedVersion` for a file that isn't a version 2 case,
    /// and the errors of [`Case::validate`].
    pub fn parse(text: &str) -> Result<Case, PowerflowError> {
        let case = matpower::parse(text)?;
        case.validate()?;
        Ok(case)
    }

    /// Checks that the case is internally consistent: at least one bus, a positive, finite
    /// `base_mva`, unique bus labels, generators and branches that refer to existing buses, and
    /// bus names and generator costs that match the bus and generator counts.
    ///
    /// [`Case::parse`] runs this; call it again after editing a case.
    pub fn validate(&self) -> Result<(), PowerflowError> {
        if self.bus.is_empty() {
            return Err(PowerflowError::NoBuses);
        }
        if !(self.base_mva.is_finite() && self.base_mva > 0.0) {
            return Err(PowerflowError::InvalidBaseMva {
                base_mva: self.base_mva,
            });
        }

        let mut seen = HashMap::with_capacity(self.bus.len());
        for (row, bus) in self.bus.iter().enumerate() {
            if let Some(first_row) = seen.insert(bus.bus_label, row + 1) {
                return Err(PowerflowError::DuplicateBus {
                    label: bus.bus_label,
                    first_row,
                    second_row: row + 1,
                });
            }
        }

        let check = |table: &'static str, row: usize, label: usize| {
            if seen.contains_key(&label) {
                Ok(())
            } else {
                Err(PowerflowError::UnknownBus { table, row, label })
            }
        };
        for (row, generator) in self.generator.iter().enumerate() {
            check("gen", row + 1, generator.bus_label)?;
        }
        for (row, branch) in self.branch.iter().enumerate() {
            check("branch", row + 1, branch.from_label)?;
            check("branch", row + 1, branch.to_label)?;
        }

        if let Some(names) = &self.busnames {
            if names.len() != self.bus.len() {
                return Err(PowerflowError::CountMismatch {
                    field: "bus_name",
                    expected: self.bus.len().to_string(),
                    found: names.len(),
                });
            }
        }
        if let Some(costs) = &self.gencost {
            let ng = self.generator.len();
            if costs.len() != ng && costs.len() != 2 * ng {
                return Err(PowerflowError::CountMismatch {
                    field: "gencost",
                    expected: format!("{ng} or {}", 2 * ng),
                    found: costs.len(),
                });
            }
        }
        Ok(())
    }

    /// Maps each bus label to its position in `bus`, the index solvers use.
    pub fn bus_indices(&self) -> HashMap<usize, usize> {
        self.bus
            .iter()
            .enumerate()
            .map(|(index, bus)| (bus.bus_label, index))
            .collect()
    }
}

/// Parses a MATPOWER case file's text; the same as [`Case::parse`].
impl FromStr for Case {
    type Err = PowerflowError;

    fn from_str(text: &str) -> Result<Case, PowerflowError> {
        Case::parse(text)
    }
}
