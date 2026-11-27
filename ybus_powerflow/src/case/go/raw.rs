//! Reading `case.raw`: the sections and fields Challenge 1 uses (formulation, A.2), faithfully.

use std::fs;
use std::path::Path;

use super::text::{self, classify, Kind, Record};
use super::GoFile;
use crate::case::{decode, BusType};
use crate::powerflow_error::PowerflowError;

/// The PSS/E revision Challenge 1's `case.raw` files are in, the only one read.
const REVISION: i64 = 33;

/// A `case.raw` file's network: the fields Challenge 1 uses of each section it reads, in file
/// order and in the file's units (MW, MVAr, MVA, kV, degrees and per unit on the system base).
///
/// Nothing is converted or cross-checked here: a record that refers to a bus that doesn't exist,
/// or two records with the same key, read without error. [`RawCase::parse`] says what's read
/// and what's checked.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawCase {
    /// The system MVA base (SBASE).
    pub sbase_mva: f64,
    /// The buses.
    pub buses: Vec<RawBus>,
    /// The loads; a bus may have several.
    pub loads: Vec<RawLoad>,
    /// The fixed shunts; a bus may have several.
    pub fixed_shunts: Vec<RawFixedShunt>,
    /// The generators.
    pub generators: Vec<RawGenerator>,
    /// The non-transformer branches.
    pub lines: Vec<RawLine>,
    /// The two-winding transformers.
    pub transformers: Vec<RawTransformer>,
    /// The switched shunts, at most one per bus.
    pub switched_shunts: Vec<RawSwitchedShunt>,
}

/// A bus record (A.2.2).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawBus {
    /// The line of `case.raw` the record is on.
    pub source_line: usize,
    /// I: the bus number, a label other records refer to.
    pub bus_label: usize,
    /// NAME, trimmed. Challenge 1 doesn't use it, and the data leaves it blank.
    pub name: String,
    /// BASKV: the base voltage (kV).
    pub base_kv: f64,
    /// IDE: the bus type, from PSS/E's codes, which are MATPOWER's (1 PQ, 2 PV, 3 reference,
    /// 4 isolated).
    pub bus_type: BusType,
    /// AREA: the area number.
    pub area: usize,
    /// ZONE: the zone number.
    pub zone: usize,
    /// VM: the voltage magnitude (p.u.), part of the operating point.
    pub vm_pu: f64,
    /// VA: the voltage angle (degrees), part of the operating point.
    pub va_deg: f64,
    /// NVHI: the base case maximum voltage magnitude (p.u.).
    pub nvhi_pu: f64,
    /// NVLO: the base case minimum voltage magnitude (p.u.).
    pub nvlo_pu: f64,
    /// EVHI: the contingency maximum voltage magnitude (p.u.).
    pub evhi_pu: f64,
    /// EVLO: the contingency minimum voltage magnitude (p.u.).
    pub evlo_pu: f64,
}

/// A load record (A.2.3).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawLoad {
    /// The line of `case.raw` the record is on.
    pub source_line: usize,
    /// I: the bus it's at.
    pub bus_label: usize,
    /// ID: its identifier among the bus's loads, trimmed.
    pub id: String,
    /// STATUS: whether it's in service.
    pub status: bool,
    /// PL: the constant power real load (MW).
    pub pl_mw: f64,
    /// QL: the constant power reactive load (MVAr).
    pub ql_mvar: f64,
}

/// A fixed shunt record (A.2.4).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawFixedShunt {
    /// The line of `case.raw` the record is on.
    pub source_line: usize,
    /// I: the bus it's at.
    pub bus_label: usize,
    /// ID: its identifier among the bus's fixed shunts, trimmed.
    pub id: String,
    /// STATUS: whether it's in service.
    pub status: bool,
    /// GL: the conductance, as the MW it draws at 1.0 p.u. voltage.
    pub gl_mw: f64,
    /// BL: the susceptance, as the MVAr it injects at 1.0 p.u. voltage.
    pub bl_mvar: f64,
}

/// A generator record (A.2.5).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawGenerator {
    /// The line of `case.raw` the record is on.
    pub source_line: usize,
    /// I: the bus it's at.
    pub bus_label: usize,
    /// ID: its identifier among the bus's generators, trimmed.
    pub id: String,
    /// PG: the real power output (MW), part of the operating point.
    pub pg_mw: f64,
    /// QG: the reactive power output (MVAr), part of the operating point.
    pub qg_mvar: f64,
    /// QT: the maximum reactive power output (MVAr).
    pub qt_mvar: f64,
    /// QB: the minimum reactive power output (MVAr).
    pub qb_mvar: f64,
    /// MBASE: the machine's MVA base. Challenge 1 doesn't use it, and the data has it 0 on some
    /// generators.
    pub mbase_mva: f64,
    /// STAT: whether it's in service.
    pub status: bool,
    /// PT: the maximum real power output (MW).
    pub pt_mw: f64,
    /// PB: the minimum real power output (MW).
    pub pb_mw: f64,
}

/// A non-transformer branch record (A.2.6).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawLine {
    /// The line of `case.raw` the record is on.
    pub source_line: usize,
    /// I: the origin bus.
    pub from_label: usize,
    /// J: the destination bus.
    pub to_label: usize,
    /// CKT: the circuit identifier, trimmed.
    pub circuit: String,
    /// R: the resistance (p.u.).
    pub r_pu: f64,
    /// X: the reactance (p.u.).
    pub x_pu: f64,
    /// B: the total charging susceptance (p.u.).
    pub b_pu: f64,
    /// RATEA: the base case rating (MVA).
    pub rate_a_mva: f64,
    /// RATEB: a rating Challenge 1 doesn't use (MVA).
    pub rate_b_mva: f64,
    /// RATEC: the contingency rating (MVA).
    pub rate_c_mva: f64,
    /// GI: the line-end shunt conductance at bus I (p.u.). Challenge 1 doesn't use it.
    pub gi_pu: f64,
    /// BI: the line-end shunt susceptance at bus I (p.u.). Challenge 1 doesn't use it.
    pub bi_pu: f64,
    /// GJ: the line-end shunt conductance at bus J (p.u.). Challenge 1 doesn't use it.
    pub gj_pu: f64,
    /// BJ: the line-end shunt susceptance at bus J (p.u.). Challenge 1 doesn't use it.
    pub bj_pu: f64,
    /// ST: whether it's in service.
    pub status: bool,
}

/// A two-winding transformer record, four lines long (A.2.7).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawTransformer {
    /// The line of `case.raw` the record starts on.
    pub source_line: usize,
    /// I: the winding 1 bus, where the tap and the magnetizing admittance are.
    pub from_label: usize,
    /// J: the winding 2 bus.
    pub to_label: usize,
    /// CKT: the circuit identifier, trimmed.
    pub circuit: String,
    /// CW: the winding data's units code. Challenge 1's model reads the ratios as with 1.
    pub cw: i64,
    /// CZ: the impedance data's units code. Challenge 1's model reads R12 and X12 as with 1.
    pub cz: i64,
    /// CM: the magnetizing admittance's units code. Challenge 1's model reads MAG1 and MAG2 as
    /// with 1.
    pub cm: i64,
    /// MAG1: the magnetizing conductance (p.u.).
    pub mag1_pu: f64,
    /// MAG2: the magnetizing susceptance (p.u.).
    pub mag2_pu: f64,
    /// STAT: whether it's in service.
    pub status: bool,
    /// R1-2: the resistance (p.u.).
    pub r12_pu: f64,
    /// X1-2: the reactance (p.u.).
    pub x12_pu: f64,
    /// SBASE1-2: the winding MVA base.
    pub sbase12_mva: f64,
    /// WINDV1: the winding 1 off-nominal turns ratio (p.u.).
    pub windv1: f64,
    /// NOMV1: the winding 1 nominal voltage (kV), 0 meaning the bus's base voltage.
    pub nomv1_kv: f64,
    /// ANG1: the phase shift angle (degrees).
    pub ang1_deg: f64,
    /// RATA1: the base case rating (MVA).
    pub rata1_mva: f64,
    /// RATB1: a rating Challenge 1 doesn't use (MVA).
    pub ratb1_mva: f64,
    /// RATC1: the contingency rating (MVA).
    pub ratc1_mva: f64,
    /// WINDV2: the winding 2 off-nominal turns ratio (p.u.).
    pub windv2: f64,
    /// NOMV2: the winding 2 nominal voltage (kV), 0 meaning the bus's base voltage.
    pub nomv2_kv: f64,
}

/// A switched shunt record (A.2.8).
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct RawSwitchedShunt {
    /// The line of `case.raw` the record is on.
    pub source_line: usize,
    /// I: the bus it's at.
    pub bus_label: usize,
    /// STAT: whether it's in service.
    pub status: bool,
    /// BINIT: the initial susceptance (MVAr at 1.0 p.u. voltage), part of the operating point.
    pub binit_mvar: f64,
    /// N1, B1 to N8, B8: each block's number of steps, and each step's susceptance (MVAr at
    /// 1.0 p.u. voltage).
    pub blocks: [(u32, f64); 8],
}

/// The sections of a revision 33 `case.raw` after the case identification, in order, with the
/// field counts of the ones Challenge 1 uses.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Section {
    Bus,
    Load,
    FixedShunt,
    Generator,
    Line,
    Transformer,
    SwitchedShunt,
    /// A section Challenge 1 doesn't use, named for errors
    Unused(&'static str),
}

const SECTIONS: [Section; 19] = [
    Section::Bus,
    Section::Load,
    Section::FixedShunt,
    Section::Generator,
    Section::Line,
    Section::Transformer,
    Section::Unused("area interchange"),
    Section::Unused("two-terminal DC line"),
    Section::Unused("VSC DC line"),
    Section::Unused("transformer impedance correction"),
    Section::Unused("multi-terminal DC line"),
    Section::Unused("multi-section line grouping"),
    Section::Unused("zone"),
    Section::Unused("interarea transfer"),
    Section::Unused("owner"),
    Section::Unused("FACTS device"),
    Section::SwitchedShunt,
    Section::Unused("GNE device"),
    Section::Unused("induction machine"),
];

/// The position of the switched shunt section, the last one read: once it's ended, `Q` or the
/// end of the file may stand in for the rest.
const LAST_USED: usize = 16;

fn error(line: usize, message: impl Into<String>) -> PowerflowError {
    text::error(GoFile::Raw, line, message)
}

impl RawCase {
    /// Reads and parses a `case.raw` file; see [`RawCase::parse`].
    ///
    /// Returns `PowerflowError::NoCaseFile` if the file can't be read.
    pub fn from_path(path: impl AsRef<Path>) -> Result<RawCase, PowerflowError> {
        let path = path.as_ref();
        let bytes = fs::read(path).map_err(|err| PowerflowError::NoCaseFile {
            bad_path: path.to_path_buf(),
            reason: err.to_string(),
        })?;
        RawCase::from_bytes(&bytes)
    }

    /// Parses the raw bytes of a `case.raw` file, decoded as
    /// [`Case::from_bytes`](crate::case::Case::from_bytes) decodes a MATPOWER file; see
    /// [`RawCase::parse`]. Challenge 1's files are ASCII.
    pub fn from_bytes(bytes: &[u8]) -> Result<RawCase, PowerflowError> {
        RawCase::parse(&decode(bytes))
    }

    /// Parses the text of a Challenge 1 `case.raw` file: PSS/E's raw format, revision 33, in
    /// the subset the formulation's A.1.3 and A.2 allow.
    ///
    /// - **Syntax:** fields separated by commas, strings in single quotes, `/` comments, and
    ///   every field present; see the formulation's A.1.2 and A.1.3. Lines may end with LF or
    ///   CRLF.
    /// - **Case identification:** three lines. The first has six fields, of which SBASE must be
    ///   positive and REV must be 33. The other two are free text.
    /// - **Sections:** the 19 of revision 33, in order, each ended by a line that reads `0`
    ///   once its comment is removed. The bus, load, fixed shunt, generator, non-transformer
    ///   branch, two-winding transformer and switched shunt sections are read, with exactly
    ///   their records' field counts. The others are skipped to their terminators, whatever
    ///   their records hold. After the switched shunt section's terminator, a `Q` line or the
    ///   end of the file may end the file early.
    /// - **Each field:** numbers finite, bus numbers and other labels whole numbers from 1,
    ///   statuses 0 or 1, IDE 1 to 4, step counts whole numbers from 0, and IDs and circuits
    ///   one or two upper case letters or digits once trimmed. A transformer's K must be 0: a
    ///   three-winding transformer is five lines, and Challenge 1 has none. Fields Challenge 1
    ///   doesn't use are counted but not read.
    ///
    /// Returns `PowerflowError::GoParse`, naming `case.raw` and the line, for anything else.
    pub fn parse(text: &str) -> Result<RawCase, PowerflowError> {
        let mut lines = text::lines(text);

        let (number, first) = lines.next().ok_or_else(|| error(1, "the file is empty"))?;
        let (sbase_mva, _) = case_identification(number, first)?;
        for expected in 2..=3 {
            if lines.next().is_none() {
                return Err(error(
                    expected,
                    "the file ends within its three-line case identification",
                ));
            }
        }

        let mut case = RawCase {
            sbase_mva,
            buses: Vec::new(),
            loads: Vec::new(),
            fixed_shunts: Vec::new(),
            generators: Vec::new(),
            lines: Vec::new(),
            transformers: Vec::new(),
            switched_shunts: Vec::new(),
        };
        let mut last = 3;
        for (position, &section) in SECTIONS.iter().enumerate() {
            loop {
                let Some((number, line)) = lines.next() else {
                    // The end of the file ends the rest, once the sections read have ended
                    return finish(case, position, last);
                };
                last = number;
                match classify(GoFile::Raw, number, line)? {
                    Kind::Terminator => break,
                    Kind::End => return finish(case, position, number),
                    Kind::Record(content) => {
                        read_record(&mut case, section, number, content, &mut lines)?
                    }
                }
            }
        }
        // Lines after the last terminator, such as the Q line, end the file
        Ok(case)
    }
}

/// `case` once the file has ended in section `position` at `line`: an error if the switched
/// shunt section hasn't ended.
fn finish(case: RawCase, position: usize, line: usize) -> Result<RawCase, PowerflowError> {
    if position > LAST_USED {
        Ok(case)
    } else {
        Err(error(line, "the file ends before its switched shunt data"))
    }
}

/// The case identification's first line: SBASE and REV.
fn case_identification(number: usize, line: &str) -> Result<(f64, i64), PowerflowError> {
    let content = match classify(GoFile::Raw, number, line)? {
        Kind::Record(content) => content,
        _ => {
            return Err(error(
                number,
                "the file must start with its case identification (IC, SBASE, REV, …)",
            ))
        }
    };
    let record = Record::split(GoFile::Raw, number, content)?;
    record.expect_len("the case identification", 6)?;
    let sbase_mva = record.float(1, "SBASE")?;
    if sbase_mva <= 0.0 {
        return Err(record.error(format!("SBASE must be positive, not {sbase_mva}")));
    }
    let revision = record.int(2, "REV")?;
    if revision != REVISION {
        return Err(record.error(format!(
            "REV is {revision}; only PSS/E revision {REVISION}, Challenge 1's, is read"
        )));
    }
    Ok((sbase_mva, revision))
}

/// Reads one record of `section`, starting at line `number` with `content`, into `case`. A
/// transformer takes three more lines from `lines`.
fn read_record<'a>(
    case: &mut RawCase,
    section: Section,
    number: usize,
    content: &'a str,
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
) -> Result<(), PowerflowError> {
    if let Section::Unused(_) = section {
        return Ok(());
    }
    let record = Record::split(GoFile::Raw, number, content)?;
    match section {
        Section::Bus => case.buses.push(bus(&record)?),
        Section::Load => case.loads.push(load(&record)?),
        Section::FixedShunt => case.fixed_shunts.push(fixed_shunt(&record)?),
        Section::Generator => case.generators.push(generator(&record)?),
        Section::Line => case.lines.push(line_record(&record)?),
        Section::Transformer => case.transformers.push(transformer(&record, lines)?),
        Section::SwitchedShunt => case.switched_shunts.push(switched_shunt(&record)?),
        Section::Unused(_) => unreachable!("returned above"),
    }
    Ok(())
}

fn bus(r: &Record) -> Result<RawBus, PowerflowError> {
    r.expect_len("a bus record", 13)?;
    let code = r.int(3, "IDE")?;
    let bus_type = BusType::from_code(code as f64).ok_or_else(|| {
        r.error(format!(
            "IDE (field 4) is {code}, not 1 (PQ), 2 (PV), 3 (reference) or 4 (isolated)"
        ))
    })?;
    Ok(RawBus {
        source_line: r.line(),
        bus_label: r.label(0, "I")?,
        name: r.string(1, "NAME")?,
        base_kv: r.float(2, "BASKV")?,
        bus_type,
        area: r.label(4, "AREA")?,
        zone: r.label(5, "ZONE")?,
        vm_pu: r.float(7, "VM")?,
        va_deg: r.float(8, "VA")?,
        nvhi_pu: r.float(9, "NVHI")?,
        nvlo_pu: r.float(10, "NVLO")?,
        evhi_pu: r.float(11, "EVHI")?,
        evlo_pu: r.float(12, "EVLO")?,
    })
}

fn load(r: &Record) -> Result<RawLoad, PowerflowError> {
    r.expect_len("a load record", 14)?;
    Ok(RawLoad {
        source_line: r.line(),
        bus_label: r.label(0, "I")?,
        id: r.short_id(1, "ID")?,
        status: r.status(2, "STATUS")?,
        pl_mw: r.float(5, "PL")?,
        ql_mvar: r.float(6, "QL")?,
    })
}

fn fixed_shunt(r: &Record) -> Result<RawFixedShunt, PowerflowError> {
    r.expect_len("a fixed shunt record", 5)?;
    Ok(RawFixedShunt {
        source_line: r.line(),
        bus_label: r.label(0, "I")?,
        id: r.short_id(1, "ID")?,
        status: r.status(2, "STATUS")?,
        gl_mw: r.float(3, "GL")?,
        bl_mvar: r.float(4, "BL")?,
    })
}

fn generator(r: &Record) -> Result<RawGenerator, PowerflowError> {
    r.expect_len("a generator record", 28)?;
    Ok(RawGenerator {
        source_line: r.line(),
        bus_label: r.label(0, "I")?,
        id: r.short_id(1, "ID")?,
        pg_mw: r.float(2, "PG")?,
        qg_mvar: r.float(3, "QG")?,
        qt_mvar: r.float(4, "QT")?,
        qb_mvar: r.float(5, "QB")?,
        mbase_mva: r.float(8, "MBASE")?,
        status: r.status(14, "STAT")?,
        pt_mw: r.float(16, "PT")?,
        pb_mw: r.float(17, "PB")?,
    })
}

fn line_record(r: &Record) -> Result<RawLine, PowerflowError> {
    r.expect_len("a non-transformer branch record", 24)?;
    Ok(RawLine {
        source_line: r.line(),
        from_label: r.label(0, "I")?,
        to_label: r.label(1, "J")?,
        circuit: r.short_id(2, "CKT")?,
        r_pu: r.float(3, "R")?,
        x_pu: r.float(4, "X")?,
        b_pu: r.float(5, "B")?,
        rate_a_mva: r.float(6, "RATEA")?,
        rate_b_mva: r.float(7, "RATEB")?,
        rate_c_mva: r.float(8, "RATEC")?,
        gi_pu: r.float(9, "GI")?,
        bi_pu: r.float(10, "BI")?,
        gj_pu: r.float(11, "GJ")?,
        bj_pu: r.float(12, "BJ")?,
        status: r.status(13, "ST")?,
    })
}

/// A transformer from its first line, `first`, and the three lines after it.
fn transformer<'a>(
    first: &Record,
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
) -> Result<RawTransformer, PowerflowError> {
    first.expect_len("a transformer's first line", 21)?;
    let k = first.int(2, "K")?;
    if k != 0 {
        return Err(first.error(format!(
            "K (field 3) is {k}: three-winding transformers aren't part of Challenge 1"
        )));
    }
    // Each further line is read whole before any of it is used, so its errors name it
    let mut next = |what: &str, count: usize| -> Result<Record<'a>, PowerflowError> {
        let Some((number, line)) = lines.next() else {
            return Err(error(
                first.line(),
                format!("the file ends within this transformer record, before its {what}"),
            ));
        };
        match classify(GoFile::Raw, number, line)? {
            Kind::Record(content) => {
                let record = Record::split(GoFile::Raw, number, content)?;
                record.expect_len(&format!("a transformer's {what}"), count)?;
                Ok(record)
            }
            _ => Err(error(
                number,
                format!(
                    "a transformer record is four lines, but its {what} ends the section or the file"
                ),
            )),
        }
    };
    let impedance = next("second line", 3)?;
    let winding1 = next("third line", 17)?;
    let winding2 = next("fourth line", 2)?;
    Ok(RawTransformer {
        source_line: first.line(),
        from_label: first.label(0, "I")?,
        to_label: first.label(1, "J")?,
        circuit: first.short_id(3, "CKT")?,
        cw: first.int(4, "CW")?,
        cz: first.int(5, "CZ")?,
        cm: first.int(6, "CM")?,
        mag1_pu: first.float(7, "MAG1")?,
        mag2_pu: first.float(8, "MAG2")?,
        status: first.status(11, "STAT")?,
        r12_pu: impedance.float(0, "R1-2")?,
        x12_pu: impedance.float(1, "X1-2")?,
        sbase12_mva: impedance.float(2, "SBASE1-2")?,
        windv1: winding1.float(0, "WINDV1")?,
        nomv1_kv: winding1.float(1, "NOMV1")?,
        ang1_deg: winding1.float(2, "ANG1")?,
        rata1_mva: winding1.float(3, "RATA1")?,
        ratb1_mva: winding1.float(4, "RATB1")?,
        ratc1_mva: winding1.float(5, "RATC1")?,
        windv2: winding2.float(0, "WINDV2")?,
        nomv2_kv: winding2.float(1, "NOMV2")?,
    })
}

fn switched_shunt(r: &Record) -> Result<RawSwitchedShunt, PowerflowError> {
    r.expect_len("a switched shunt record", 26)?;
    let mut blocks = [(0, 0.0); 8];
    for (k, block) in blocks.iter_mut().enumerate() {
        let n = 10 + 2 * k;
        *block = (
            r.count(n, &format!("N{}", k + 1))?,
            r.float(n + 1, &format!("B{}", k + 1))?,
        );
    }
    Ok(RawSwitchedShunt {
        source_line: r.line(),
        bus_label: r.label(0, "I")?,
        status: r.status(3, "STAT")?,
        binit_mvar: r.float(9, "BINIT")?,
        blocks,
    })
}
