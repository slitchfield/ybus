//! Reading MATPOWER case files (version 2 format).
//!
//! A case file is MATLAB code that assigns fields of a struct, `mpc` by convention. Rather than
//! interpret MATLAB, this reads the subset that data files use and rejects everything else:
//!
//! - comments: `%` and Octave's `#` to the end of the line (outside strings), and `%{ … %}` or
//!   `#{ … #}` blocks, whose markers stand alone on their lines and may nest;
//! - blank lines, `end`, and a leading `function <name> = case` line, whose output names the
//!   struct (`mpc` when there's no function line);
//! - statements separated by new lines, `;` or `,`, so several can share a line;
//! - for the fields the model reads (`version`, `baseMVA`, `bus`, `gen`, `branch`, `gencost`,
//!   `bus_name`), `<struct>.<field> = <literal>`, where the literal is a matrix of plain numbers
//!   (`[ … ]`, rows ended by `;` or a new line, values separated by spaces, tabs or commas,
//!   `...` continuing a row onto the next line), a cell array of strings, a single- or
//!   double-quoted string, or a number;
//! - for any other field, including nested ones such as `mpc.if.lims`, an assignment of
//!   anything at all, which is skipped. MATLAB's value semantics mean such an assignment can't
//!   change the fields the model reads, so its right-hand side needn't be understood, only its
//!   brackets and strings balanced to find where it ends.
//!
//! Strings follow MATLAB's rule for `'`: it starts a string unless it directly follows a name,
//! number, closing bracket, `.` or `'`, where it's the transpose operator.
//!
//! Version 1 files, which assign `baseMVA`, `bus` and so on as separate variables, are
//! reported as such. Values the model reads must not be NaN; `Inf` is allowed.
//!
//! Some MATPOWER cases rescale their data with code after the literals (for example
//! `mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;`). Ignoring that code would silently
//! give the wrong values, so any statement outside the subset is an error, as are numeric
//! expressions such as `50/3` inside a matrix. To read such a case, evaluate it with MATPOWER's
//! `loadcase` and write it out as literals first; in this repository,
//! `scripts/normalize-cases.sh` does that for every MATPOWER case.

use std::collections::HashMap;

use super::{Branch, Bus, BusType, Case, Gen, GenCost};
use crate::powerflow_error::PowerflowError;

/// The columns each table must have. Extra columns, such as the results MATPOWER appends to
/// solved cases, are ignored.
///
/// Generator columns 11-21 and branch columns 12-13 were added in version 2 of the format, and
/// MATPOWER itself accepts version 2 files without them (case533mt has 18-column generator
/// rows). Missing ones get the values MATPOWER's own version 1 to 2 conversion (`mpc_1to2` in
/// `loadcase.m`) fills in: 0 for the generator columns, and -360 and 360 (unconstrained) for
/// the branch angle limits.
const BUS_COLUMNS: usize = 13;
const GEN_COLUMNS: usize = 10;
const BRANCH_COLUMNS: usize = 11;

/// A matrix row, with the line it starts on for error messages.
struct Row {
    line: usize,
    values: Vec<f64>,
}

/// The literal assigned to an `mpc` field.
enum Value {
    Matrix(Vec<Row>),
    Cell(Vec<String>),
    Str(String),
    Number(f64),
}

impl Value {
    fn kind(&self) -> &'static str {
        match self {
            Value::Matrix(_) => "a matrix",
            Value::Cell(_) => "a cell array",
            Value::Str(_) => "a string",
            Value::Number(_) => "a number",
        }
    }
}

/// An `mpc` field's value and the line its assignment starts on.
struct Field {
    line: usize,
    value: Value,
}

fn parse_error(line: usize, message: impl Into<String>) -> PowerflowError {
    PowerflowError::Parse {
        line,
        message: message.into(),
    }
}

/// Parses case file text into a [`Case`], without the cross-table checks of
/// [`Case::validate`].
pub(super) fn parse(text: &str) -> Result<Case, PowerflowError> {
    build_case(read_fields(text)?)
}

/// The fields the case model reads. Their assignments must be literals; any other field's
/// assignment is skipped, whatever its right-hand side, because MATLAB's value semantics mean
/// `mpc.other = …` can't change these.
const USED_FIELDS: &[&str] = &[
    "version", "baseMVA", "bus", "gen", "branch", "gencost", "bus_name",
];

/// The variables a version 1 case file assigns, instead of fields of one struct.
const VERSION_1_VARIABLES: &[&str] = &["baseMVA", "bus", "gen", "branch", "areas", "gencost"];

/// One line with its comments removed, and whether it ends in a `...` continuation.
struct CodeLine<'a> {
    number: usize,
    code: &'a str,
    continued: bool,
}

/// Whether a `'` at byte `i` of `text` starts a string. MATLAB reads it as the transpose
/// operator instead when it directly follows a name, a number, a closing bracket, a `.`, or
/// another `'` (as in `x''`).
fn quote_starts_string(text: &str, i: usize) -> bool {
    match text[..i].chars().next_back() {
        Some(c) => !(c.is_alphanumeric() || matches!(c, '_' | ')' | ']' | '}' | '.' | '\'')),
        None => true,
    }
}

/// Finds the first byte outside strings in `text` (one line of code) for which `stop` returns
/// true, and returns its index. Single- and double-quoted strings are skipped, with doubled
/// quotes as escapes, and a `'` that transposes isn't taken for a string. Returns `Err(start)`
/// for a string that isn't closed on the line, since MATLAB strings can't span lines.
///
/// The scan is by bytes, not characters: the quotes and every byte a caller stops on are ASCII,
/// and in UTF-8 an ASCII byte only ever encodes that character, so the indices found are those
/// of a scan by characters. A character outside ASCII reaches `stop` as each of its bytes, all
/// at least 0x80, which the callers treat as they would the character.
fn find_code(text: &str, mut stop: impl FnMut(u8) -> bool) -> Result<Option<usize>, usize> {
    let bytes = text.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let quote = match bytes[i] {
            b'\'' if quote_starts_string(text, i) => b'\'',
            b'"' => b'"',
            b if stop(b) => return Ok(Some(i)),
            _ => {
                i += 1;
                continue;
            }
        };
        // Skip the string, where a doubled quote is an escaped one
        let start = i;
        i += 1;
        loop {
            match bytes[i..].iter().position(|&b| b == quote) {
                Some(offset) if bytes.get(i + offset + 1) == Some(&quote) => i += offset + 2,
                Some(offset) => {
                    i += offset + 1;
                    break;
                }
                None => return Err(start),
            }
        }
    }
    Ok(None)
}

/// Removes a `%` or `#` comment and a `...` continuation (and the rest of its line) from
/// `line`, outside strings. A line with an unclosed string is kept whole, for the parser to
/// report.
fn strip_line(number: usize, line: &str) -> CodeLine<'_> {
    let mut previous_dot_run = 0;
    let found = find_code(line, |b| {
        if b == b'.' {
            previous_dot_run += 1;
        } else {
            previous_dot_run = 0;
        }
        matches!(b, b'%' | b'#') || previous_dot_run == 3
    });
    match found {
        Ok(Some(i)) if line.as_bytes()[i] == b'.' => CodeLine {
            number,
            // The index is the third dot of `...`
            code: &line[..i - 2],
            continued: true,
        },
        Ok(Some(i)) => CodeLine {
            number,
            code: &line[..i],
            continued: false,
        },
        Ok(None) | Err(_) => CodeLine {
            number,
            code: line,
            continued: false,
        },
    }
}

/// Splits text into code lines, removing line comments and `%{ … %}` (or Octave's
/// `#{ … #}`) block comments, whose markers must be alone on their lines and may nest.
fn code_lines(text: &str) -> Vec<CodeLine<'_>> {
    let mut depth = 0usize;
    text.lines()
        .enumerate()
        .map(|(index, line)| {
            let number = index + 1;
            let marker = line.trim();
            let blank = CodeLine {
                number,
                code: "",
                continued: false,
            };
            if marker == "%{" || marker == "#{" {
                depth += 1;
                blank
            } else if depth > 0 {
                if marker == "%}" || marker == "#}" {
                    depth -= 1;
                }
                blank
            } else {
                strip_line(number, line)
            }
        })
        .collect()
}

/// A position in the code: a line index and a byte offset into that line's code.
#[derive(Clone, Copy)]
struct Cursor {
    index: usize,
    pos: usize,
}

/// Reads every assignment, returning the literals assigned to the fields the model uses.
/// Assignments to other fields are skipped; any other statement is an error.
fn read_fields(text: &str) -> Result<HashMap<String, Field>, PowerflowError> {
    let lines = code_lines(text);
    let mut fields: HashMap<String, Field> = HashMap::new();
    // The struct's name: the function's output, or `mpc` in a file without a function line
    let mut variable = "mpc".to_string();
    let mut seen_statement = false;
    let mut cursor = Cursor { index: 0, pos: 0 };

    while cursor.index < lines.len() {
        let line = &lines[cursor.index];
        let rest = &line.code[cursor.pos..];
        let statement = rest.trim_start();
        if statement.is_empty() {
            cursor = Cursor {
                index: cursor.index + 1,
                pos: 0,
            };
            continue;
        }
        let start = cursor.pos + (rest.len() - statement.len());
        // Empty statements between separators
        if statement.starts_with([';', ',']) {
            cursor.pos = start + 1;
            continue;
        }

        let first_word = leading_identifier(statement);
        if first_word == "function" && !seen_statement {
            variable = function_output(line.number, statement)?;
            seen_statement = true;
            cursor = Cursor {
                index: cursor.index + 1,
                pos: 0,
            };
            continue;
        }
        seen_statement = true;
        if matches!(first_word, "end" | "endfunction") {
            let after = statement[first_word.len()..].trim_start();
            if after.is_empty() || after.starts_with([';', ',']) {
                cursor.pos = start + first_word.len();
                continue;
            }
        }

        let Some((path, rhs_offset)) = assignment_target(statement, &variable) else {
            if VERSION_1_VARIABLES.contains(&first_word)
                && assignment_target(statement, first_word).is_none()
                && is_plain_assignment(&statement[first_word.len()..])
            {
                return Err(PowerflowError::UnsupportedVersion {
                    found: "1".to_string(),
                });
            }
            return Err(unsupported_statement(line.number, statement));
        };
        let rhs = Cursor {
            index: cursor.index,
            pos: start + rhs_offset,
        };

        let name = path[0];
        if !USED_FIELDS.contains(&name) {
            cursor = skip_expression(&lines, rhs)?;
            continue;
        }
        if path.len() > 1 {
            // Such as `mpc.bus.x = …`, which isn't valid for a matrix anyway
            return Err(unsupported_statement(line.number, statement));
        }

        let (value, after) = read_literal(&lines, rhs)?;
        cursor = expect_statement_end(&lines, after)?;
        if let Some(previous) = fields.insert(
            name.to_string(),
            Field {
                line: line.number,
                value,
            },
        ) {
            return Err(parse_error(
                line.number,
                format!(
                    "mpc.{name} is assigned twice (first at line {})",
                    previous.line
                ),
            ));
        }
    }
    Ok(fields)
}

/// The identifier at the start of `text`, or `""`.
fn leading_identifier(text: &str) -> &str {
    let end = text
        .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
        .unwrap_or(text.len());
    &text[..end]
}

/// Whether `text` (after a name) is `= …`, an assignment rather than a comparison.
fn is_plain_assignment(text: &str) -> bool {
    text.trim_start()
        .strip_prefix('=')
        .is_some_and(|after| !after.starts_with('='))
}

/// The output variable named by a `function` line, such as `mpc` in `function mpc = case14`.
/// A version 1 file's `function [baseMVA, bus, gen, branch, …] = name` is reported as such.
fn function_output(line: usize, statement: &str) -> Result<String, PowerflowError> {
    let unsupported = || {
        parse_error(
            line,
            format!("unsupported function line `{statement}`: it must return the case struct"),
        )
    };
    let signature = statement["function".len()..].trim_start();
    let (outputs, _) = signature.split_once('=').ok_or_else(unsupported)?;
    let outputs = outputs.trim();
    let names: Vec<&str> = outputs
        .trim_start_matches('[')
        .trim_end_matches(']')
        .split([',', ' ', '\t'])
        .filter(|name| !name.is_empty())
        .collect();
    if names.len() > 1 && names.iter().any(|name| VERSION_1_VARIABLES.contains(name)) {
        return Err(PowerflowError::UnsupportedVersion {
            found: "1".to_string(),
        });
    }
    match names.as_slice() {
        [name] if !name.is_empty() && leading_identifier(name) == *name => Ok(name.to_string()),
        _ => Err(unsupported()),
    }
}

/// If `statement` assigns to a field of `variable` (`variable.a = …` or `variable.a.b = …`),
/// returns the field path and the byte offset just past the `=`.
fn assignment_target<'a>(statement: &'a str, variable: &str) -> Option<(Vec<&'a str>, usize)> {
    let mut rest = statement.strip_prefix(variable)?;
    let mut path = Vec::new();
    while let Some(after_dot) = rest.strip_prefix('.') {
        let name = leading_identifier(after_dot);
        if name.is_empty() {
            return None;
        }
        path.push(name);
        rest = &after_dot[name.len()..];
    }
    if path.is_empty() || !is_plain_assignment(rest) {
        // Such as an indexed assignment, `mpc.bus(:, PD) = …`
        return None;
    }
    let equals = statement.len() - rest.trim_start().len();
    Some((path, equals + 1))
}

fn unsupported_statement(line: usize, statement: &str) -> PowerflowError {
    parse_error(
        line,
        format!(
            "unsupported statement `{}`: only `mpc.<field> = <literal>;` assignments are \
             supported, and this case file may compute its data with code",
            statement.trim_end()
        ),
    )
}

/// Moves past whitespace and one optional `;` or `,` ending a statement at `cursor`, which must
/// otherwise be at the end of its line.
fn expect_statement_end(lines: &[CodeLine], cursor: Cursor) -> Result<Cursor, PowerflowError> {
    let line = &lines[cursor.index];
    let rest = &line.code[cursor.pos..];
    let trimmed = rest.trim_start();
    let at = cursor.pos + (rest.len() - trimmed.len());
    if trimmed.is_empty() {
        Ok(Cursor {
            index: cursor.index + 1,
            pos: 0,
        })
    } else if trimmed.starts_with([';', ',']) {
        Ok(Cursor {
            index: cursor.index,
            pos: at + 1,
        })
    } else {
        Err(parse_error(
            line.number,
            format!("unexpected `{}` after the literal", trimmed.trim_end()),
        ))
    }
}

/// Skips an expression assigned to a field the model doesn't use, through the `;` or `,` that
/// ends it, or the end of its line. Brackets must balance, and inside them new lines don't
/// end the statement.
fn skip_expression(lines: &[CodeLine], start: Cursor) -> Result<Cursor, PowerflowError> {
    let mut depth = 0usize;
    let mut cursor = start;
    loop {
        let line = &lines[cursor.index];
        let text = &line.code[cursor.pos..];
        let mut unbalanced = false;
        let found = find_code(text, |b| match b {
            b'(' | b'[' | b'{' => {
                depth += 1;
                false
            }
            b')' | b']' | b'}' => match depth.checked_sub(1) {
                Some(d) => {
                    depth = d;
                    false
                }
                None => {
                    unbalanced = true;
                    true
                }
            },
            b';' | b',' => depth == 0,
            _ => false,
        })
        .map_err(|_| parse_error(line.number, "string is never closed"))?;
        if unbalanced {
            return Err(parse_error(
                line.number,
                "closing bracket without an opening one",
            ));
        }
        if let Some(i) = found {
            return Ok(Cursor {
                index: cursor.index,
                pos: cursor.pos + i + 1,
            });
        }
        if depth == 0 && !line.continued {
            return Ok(Cursor {
                index: cursor.index + 1,
                pos: 0,
            });
        }
        if cursor.index + 1 == lines.len() {
            return Err(parse_error(
                lines[start.index].number,
                "this statement's brackets never close",
            ));
        }
        cursor = Cursor {
            index: cursor.index + 1,
            pos: 0,
        };
    }
}

/// Reads the literal at `start` (just past the `=`), which may continue over later lines.
/// Returns the value and the position just after it.
fn read_literal(lines: &[CodeLine], start: Cursor) -> Result<(Value, Cursor), PowerflowError> {
    let line = &lines[start.index];
    let rest = &line.code[start.pos..];
    let rhs = rest.trim_start();
    let at = start.pos + (rest.len() - rhs.len());
    if rhs.starts_with('[') || rhs.starts_with('{') {
        let close = if rhs.starts_with('[') { ']' } else { '}' };
        let body = Cursor {
            index: start.index,
            pos: at + 1,
        };
        let (segments, after) = collect_until(lines, body, close)?;
        let value = if close == ']' {
            Value::Matrix(parse_matrix(&segments)?)
        } else {
            Value::Cell(parse_cell(&segments)?)
        };
        return Ok((value, after));
    }
    if rhs.starts_with(['\'', '"']) {
        let (string, after) = read_string(line.number, rhs)?;
        let pos = at + (rhs.len() - after.len());
        return Ok((
            Value::Str(string),
            Cursor {
                index: start.index,
                pos,
            },
        ));
    }

    let token_len = rhs.find([';', ',']).unwrap_or(rhs.len());
    let value = parse_number(line.number, rhs[..token_len].trim(), "a scalar")?;
    Ok((
        Value::Number(value),
        Cursor {
            index: start.index,
            pos: at + token_len,
        },
    ))
}

/// A piece of a multi-line literal: its line, its text, and whether the line was continued
/// with `...` (which doesn't end a matrix row).
struct Segment<'a> {
    line: usize,
    text: &'a str,
    continued: bool,
}

/// Gathers a literal's body from `start` (just past the opening bracket) until the closing
/// `close` (ASCII: `]` or `}`) outside strings. Returns the body's segments and the position just after `close`.
fn collect_until<'a>(
    lines: &'a [CodeLine<'a>],
    start: Cursor,
    close: char,
) -> Result<(Vec<Segment<'a>>, Cursor), PowerflowError> {
    let mut segments = Vec::new();
    let mut cursor = start;
    loop {
        let line = &lines[cursor.index];
        let text = &line.code[cursor.pos..];
        let found = find_code(text, |b| b == close as u8)
            .map_err(|_| parse_error(line.number, "string is never closed"))?;
        if let Some(end) = found {
            segments.push(Segment {
                line: line.number,
                text: &text[..end],
                continued: false,
            });
            return Ok((
                segments,
                Cursor {
                    index: cursor.index,
                    pos: cursor.pos + end + close.len_utf8(),
                },
            ));
        }
        segments.push(Segment {
            line: line.number,
            text,
            continued: line.continued,
        });
        if cursor.index + 1 == lines.len() {
            return Err(parse_error(
                lines[start.index].number,
                format!("`{close}` never closes this literal"),
            ));
        }
        cursor = Cursor {
            index: cursor.index + 1,
            pos: 0,
        };
    }
}

/// Parses a plain number, including MATLAB's `Inf` and `NaN`. Expressions are rejected.
fn parse_number(line: usize, token: &str, what: &str) -> Result<f64, PowerflowError> {
    token.parse::<f64>().map_err(|_| {
        parse_error(
            line,
            format!("`{token}` in {what} is not a plain number (expressions aren't supported)"),
        )
    })
}

/// Splits a matrix body into rows: `;` always ends a row, and so does the end of a line
/// unless it was continued with `...`. Every row must have the same number of columns.
fn parse_matrix(segments: &[Segment]) -> Result<Vec<Row>, PowerflowError> {
    let mut rows: Vec<Row> = Vec::new();
    let mut current: Option<Row> = None;
    let finish = |row: Option<Row>, rows: &mut Vec<Row>| -> Result<(), PowerflowError> {
        if let Some(row) = row {
            if let Some(first) = rows.first() {
                if row.values.len() != first.values.len() {
                    return Err(parse_error(
                        row.line,
                        format!(
                            "matrix row has {} columns, but the first row (line {}) has {}",
                            row.values.len(),
                            first.line,
                            first.values.len()
                        ),
                    ));
                }
            }
            rows.push(row);
        }
        Ok(())
    };

    for segment in segments {
        let mut parts = segment.text.split(';').peekable();
        while let Some(part) = parts.next() {
            for token in part
                .split(|c: char| c.is_whitespace() || c == ',')
                .filter(|token| !token.is_empty())
            {
                let value = parse_number(segment.line, token, "a matrix")?;
                current
                    .get_or_insert_with(|| Row {
                        line: segment.line,
                        values: Vec::new(),
                    })
                    .values
                    .push(value);
            }
            // A `;` ends the row; so does the end of a line that isn't continued
            if parts.peek().is_some() || !segment.continued {
                finish(current.take(), &mut rows)?;
            }
        }
    }
    finish(current.take(), &mut rows)?;
    Ok(rows)
}

/// Splits a cell array body into its strings, separated by spaces, commas, `;` or new lines.
fn parse_cell(segments: &[Segment]) -> Result<Vec<String>, PowerflowError> {
    let mut strings = Vec::new();
    for segment in segments {
        let mut rest = segment.text;
        loop {
            rest = rest.trim_start_matches(|c: char| c.is_whitespace() || c == ',' || c == ';');
            if rest.is_empty() {
                break;
            }
            if !rest.starts_with(['\'', '"']) {
                return Err(parse_error(
                    segment.line,
                    format!(
                        "cell arrays may only hold quoted strings, not `{}`",
                        rest.trim()
                    ),
                ));
            }
            let (string, after) = read_string(segment.line, rest)?;
            strings.push(string);
            rest = after;
        }
    }
    Ok(strings)
}

/// Reads a single- or double-quoted string at the start of `text`, where a doubled quote is
/// an escaped one. Returns the string and the text after its closing quote.
fn read_string(line: usize, text: &str) -> Result<(String, &str), PowerflowError> {
    let quote = text.chars().next().filter(|&c| c == '\'' || c == '"');
    let Some(quote) = quote else {
        return Err(parse_error(line, "expected a quoted string"));
    };
    let body = &text[1..];
    let mut string = String::new();
    let mut chars = body.char_indices().peekable();
    while let Some((i, c)) = chars.next() {
        if c == quote {
            if chars.peek().is_some_and(|&(_, next)| next == quote) {
                string.push(quote);
                chars.next();
            } else {
                return Ok((string, &body[i + 1..]));
            }
        } else {
            string.push(c);
        }
    }
    Err(parse_error(line, "string is never closed"))
}

/// Builds a [`Case`] from the assigned fields, converting each table's rows.
fn build_case(mut fields: HashMap<String, Field>) -> Result<Case, PowerflowError> {
    let version = fields
        .remove("version")
        .ok_or(PowerflowError::MissingField { field: "version" })?;
    match version.value {
        Value::Str(found) if found == "2" => {}
        Value::Str(found) => return Err(PowerflowError::UnsupportedVersion { found }),
        other => {
            return Err(parse_error(
                version.line,
                format!("mpc.version must be a string, not {}", other.kind()),
            ))
        }
    }

    let base_mva = match take(&mut fields, "baseMVA")? {
        Field {
            value: Value::Number(value),
            ..
        } => value,
        Field { line, value } => {
            return Err(parse_error(
                line,
                format!("mpc.baseMVA must be a number, not {}", value.kind()),
            ))
        }
    };

    let bus = matrix(take(&mut fields, "bus")?, "bus")?
        .iter()
        .map(bus_from_row)
        .collect::<Result<_, _>>()?;
    let generator = matrix(take(&mut fields, "gen")?, "gen")?
        .iter()
        .map(gen_from_row)
        .collect::<Result<_, _>>()?;
    let branch = matrix(take(&mut fields, "branch")?, "branch")?
        .iter()
        .map(branch_from_row)
        .collect::<Result<_, _>>()?;
    let gencost = match fields.remove("gencost") {
        Some(field) => Some(
            matrix(field, "gencost")?
                .iter()
                .map(gencost_from_row)
                .collect::<Result<_, _>>()?,
        ),
        None => None,
    };
    let busnames = match fields.remove("bus_name") {
        Some(Field {
            value: Value::Cell(names),
            ..
        }) => Some(names),
        Some(Field { line, value }) => {
            return Err(parse_error(
                line,
                format!("mpc.bus_name must be a cell array, not {}", value.kind()),
            ))
        }
        None => None,
    };
    // Any other field (gentype, genfuel, dcline, areas, user fields) isn't used

    Ok(Case {
        base_mva,
        bus,
        generator,
        branch,
        gencost,
        busnames,
    })
}

fn take(fields: &mut HashMap<String, Field>, name: &'static str) -> Result<Field, PowerflowError> {
    fields
        .remove(name)
        .ok_or(PowerflowError::MissingField { field: name })
}

fn matrix(field: Field, name: &str) -> Result<Vec<Row>, PowerflowError> {
    match field.value {
        Value::Matrix(rows) => Ok(rows),
        other => Err(parse_error(
            field.line,
            format!("mpc.{name} must be a matrix, not {}", other.kind()),
        )),
    }
}

/// Column `column` of a row, or `default` for a row too short to have it.
fn optional(values: &[f64], column: usize, default: f64) -> f64 {
    values.get(column).copied().unwrap_or(default)
}

/// Checks that the first `count` values of a row (the ones the model reads) aren't NaN. No
/// MATPOWER case uses NaN, and one in, say, a reactance would only surface later as a
/// confusing solver failure. `Inf` is allowed: limits use it.
fn reject_nan(row: &Row, table: &str, count: usize) -> Result<(), PowerflowError> {
    match row
        .values
        .iter()
        .take(count)
        .position(|value| value.is_nan())
    {
        Some(column) => Err(parse_error(
            row.line,
            format!("mpc.{table} column {} is NaN", column + 1),
        )),
        None => Ok(()),
    }
}

/// Checks a row has at least `min` columns.
fn require_columns(row: &Row, table: &str, min: usize) -> Result<(), PowerflowError> {
    if row.values.len() < min {
        Err(parse_error(
            row.line,
            format!(
                "mpc.{table} row has {} columns, but needs at least {min}",
                row.values.len()
            ),
        ))
    } else {
        Ok(())
    }
}

/// Converts a value that must be a whole number of at least `min`, such as a bus label.
///
/// The largest accepted value is 2^53 (or `usize::MAX` if smaller): past it an `f64` can't
/// represent every whole number, so the value read may not be the one the file wrote, and a
/// cast to `usize` would saturate rather than fail.
fn whole(row: &Row, column: usize, what: &str, min: usize) -> Result<usize, PowerflowError> {
    const MAX_EXACT: f64 = 9_007_199_254_740_992.0; // 2^53
    let max = MAX_EXACT.min(usize::MAX as f64);
    let value = row.values[column];
    if value.is_finite() && value.fract() == 0.0 && value >= min as f64 && value <= max {
        Ok(value as usize)
    } else {
        Err(parse_error(
            row.line,
            format!("{what} must be a whole number from {min} to {max}, not {value}"),
        ))
    }
}

fn bus_from_row(row: &Row) -> Result<Bus, PowerflowError> {
    require_columns(row, "bus", BUS_COLUMNS)?;
    reject_nan(row, "bus", BUS_COLUMNS)?;
    let v = &row.values;
    let bus_type = BusType::from_code(v[1]).ok_or(PowerflowError::UnknownBusType {
        line: row.line,
        code: v[1],
    })?;
    Ok(Bus {
        bus_label: whole(row, 0, "bus number", 1)?,
        bus_type,
        pd_mw: v[2],
        qd_mvar: v[3],
        gs_mw: v[4],
        bs_mvar: v[5],
        // The format calls area and zone positive, but real cases use 0
        area: whole(row, 6, "bus area", 0)?,
        vm_pu: v[7],
        va_deg: v[8],
        base_kv: v[9],
        zone: whole(row, 10, "bus zone", 0)?,
        vmax_pu: v[11],
        vmin_pu: v[12],
    })
}

fn gen_from_row(row: &Row) -> Result<Gen, PowerflowError> {
    require_columns(row, "gen", GEN_COLUMNS)?;
    // Up to the 21 columns of the version 2 format, the ones read
    reject_nan(row, "gen", 21)?;
    let v = &row.values;
    Ok(Gen {
        bus_label: whole(row, 0, "generator bus number", 1)?,
        pg_mw: v[1],
        qg_mvar: v[2],
        qmax_mvar: v[3],
        qmin_mvar: v[4],
        vg_pu: v[5],
        mbase_mva: v[6],
        status: v[7] > 0.0,
        pmax_mw: v[8],
        pmin_mw: v[9],
        pc1_mw: optional(v, 10, 0.0),
        pc2_mw: optional(v, 11, 0.0),
        qc1min_mvar: optional(v, 12, 0.0),
        qc1max_mvar: optional(v, 13, 0.0),
        qc2min_mvar: optional(v, 14, 0.0),
        qc2max_mvar: optional(v, 15, 0.0),
        ramp_agc_mw_per_min: optional(v, 16, 0.0),
        ramp_10_mw: optional(v, 17, 0.0),
        ramp_30_mw: optional(v, 18, 0.0),
        ramp_q_mvar_per_min: optional(v, 19, 0.0),
        apf: optional(v, 20, 0.0),
    })
}

fn branch_from_row(row: &Row) -> Result<Branch, PowerflowError> {
    require_columns(row, "branch", BRANCH_COLUMNS)?;
    reject_nan(row, "branch", 13)?;
    let v = &row.values;
    Ok(Branch {
        from_label: whole(row, 0, "branch from-bus number", 1)?,
        to_label: whole(row, 1, "branch to-bus number", 1)?,
        r: v[2],
        x: v[3],
        b: v[4],
        rate_a_mva: v[5],
        rate_b_mva: v[6],
        rate_c_mva: v[7],
        tap_ratio: v[8],
        shift_deg: v[9],
        status: v[10] != 0.0,
        angmin_deg: optional(v, 11, -360.0),
        angmax_deg: optional(v, 12, 360.0),
    })
}

/// Converts a `gencost` row: model, startup, shutdown, `N`, then `N` polynomial coefficients
/// (model 2) or `N` `(p, cost)` points (model 1). Rows can be longer, padded with zeros so
/// the matrix is rectangular.
fn gencost_from_row(row: &Row) -> Result<GenCost, PowerflowError> {
    require_columns(row, "gencost", 4)?;
    let v = &row.values;
    let n = whole(row, 3, "gencost N", 1)?;
    let (startup, shutdown) = (v[1], v[2]);
    // Values after the first four: N coefficients, or N (p, cost) pairs. Checked, since N
    //   comes from the file and 4 + 2N could otherwise overflow
    let values_per_point = match v[0] {
        2.0 => 1,
        1.0 => 2,
        model => {
            return Err(parse_error(
                row.line,
                format!("gencost model {model} is not 1 (piecewise linear) or 2 (polynomial)"),
            ))
        }
    };
    let end = n
        .checked_mul(values_per_point)
        .and_then(|count| count.checked_add(4))
        .filter(|&end| end <= v.len())
        .ok_or_else(|| {
            parse_error(
                row.line,
                format!(
                    "gencost N is {n}, but the row has only {} values after the first four",
                    v.len() - 4
                ),
            )
        })?;
    reject_nan(row, "gencost", end)?;
    let values = &v[4..end];
    Ok(if values_per_point == 1 {
        GenCost::Polynomial {
            startup,
            shutdown,
            coeffs: values.to_vec(),
        }
    } else {
        GenCost::PiecewiseLinear {
            startup,
            shutdown,
            points: values
                .chunks_exact(2)
                .map(|pair| (pair[0], pair[1]))
                .collect(),
        }
    })
}
