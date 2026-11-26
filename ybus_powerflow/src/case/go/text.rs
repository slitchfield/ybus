//! The tokenizer for the comma-separated GO files (`case.raw`, `case.rop` and `case.inl`), for
//! the subset Challenge 1 allows (formulation, A.1.2 and A.1.3):
//!
//! - **Lines:** a line ends with LF or CRLF, mixed within a file, as the data has them.
//! - **Comments:** a `/` outside single quotes starts a comment, to the end of the line.
//! - **Line kinds:** after its comment is removed and it's trimmed, a line that reads `0` ends
//!   a section (so `0 / END OF BUS DATA`, `0 /END …` and a bare `0` all do), one that reads `Q`
//!   ends the file, and any other is a record.
//! - **Fields** are separated by commas, with blanks around them ignored. A string is in single
//!   quotes, with no escapes. Every field is present: Challenge 1 leaves none to a default, so an
//!   empty field is an error.
//! - **Strictness:** what full PSS/E allows but Challenge 1 rules out (fields separated by
//!   blanks, double-quoted strings) is an error, not read with a guess.

// The reader in the next commit is its first user
#![cfg_attr(not(test), allow(dead_code))]

use super::GoFile;
use crate::powerflow_error::PowerflowError;

/// A `GoParse` error in `file` at `line`.
pub(super) fn error(file: GoFile, line: usize, message: impl Into<String>) -> PowerflowError {
    PowerflowError::GoParse {
        file,
        line,
        message: message.into(),
    }
}

/// The lines of `text`, numbered from 1, without their line endings.
pub(super) fn lines(text: &str) -> impl Iterator<Item = (usize, &str)> {
    // str::lines ends a line at LF, removing a CR before it
    text.lines().enumerate().map(|(i, line)| (i + 1, line))
}

/// What a line is, once its comment is removed.
#[derive(Debug, PartialEq)]
pub(super) enum Kind<'a> {
    /// `0`: the end of a section
    Terminator,
    /// `Q`: the end of the file
    End,
    /// Anything else: a record's text, trimmed, without its comment
    Record(&'a str),
}

/// Classifies line `number` of `file`, whose text is `line`.
pub(super) fn classify(
    file: GoFile,
    number: usize,
    line: &str,
) -> Result<Kind<'_>, PowerflowError> {
    let content = strip_comment(file, number, line)?.trim();
    Ok(match content {
        "0" => Kind::Terminator,
        "Q" => Kind::End,
        _ => Kind::Record(content),
    })
}

/// `line` without its `/` comment. A `/` inside single quotes is part of a string.
fn strip_comment(file: GoFile, number: usize, line: &str) -> Result<&str, PowerflowError> {
    let mut quoted = false;
    for (i, byte) in line.bytes().enumerate() {
        match byte {
            b'\'' => quoted = !quoted,
            b'/' if !quoted => return Ok(&line[..i]),
            _ => {}
        }
    }
    if quoted {
        Err(error(file, number, "a string's closing quote is missing"))
    } else {
        Ok(line)
    }
}

/// One field of a record, as written
#[derive(Clone, Copy, Debug, PartialEq)]
enum Field<'a> {
    /// A single-quoted string, without its quotes
    Quoted(&'a str),
    /// Anything else: a number, if the record is well formed
    Bare(&'a str),
}

/// A record's fields, with where it came from for errors
#[derive(Debug)]
pub(super) struct Record<'a> {
    file: GoFile,
    line: usize,
    fields: Vec<Field<'a>>,
}

impl<'a> Record<'a> {
    /// Splits `content`, a record's text from [`classify`], into its fields.
    pub(super) fn split(
        file: GoFile,
        line: usize,
        content: &'a str,
    ) -> Result<Record<'a>, PowerflowError> {
        let mut fields = Vec::new();
        let mut rest = content;
        loop {
            let number = fields.len() + 1;
            let trimmed = rest.trim_start();
            let (field, after) = if let Some(quoted) = trimmed.strip_prefix('\'') {
                // classify has checked that every quote is closed
                let end = quoted
                    .find('\'')
                    .ok_or_else(|| error(file, line, "a string's closing quote is missing"))?;
                let after = quoted[end + 1..].trim_start();
                if !(after.is_empty() || after.starts_with(',')) {
                    return Err(error(
                        file,
                        line,
                        format!("field {number} has text after its closing quote"),
                    ));
                }
                (Field::Quoted(&quoted[..end]), after)
            } else {
                let end = trimmed.find(',').unwrap_or(trimmed.len());
                let bare = trimmed[..end].trim_end();
                if bare.is_empty() {
                    return Err(error(
                        file,
                        line,
                        format!(
                            "field {number} is empty; Challenge 1 leaves no field to a default"
                        ),
                    ));
                }
                if bare.starts_with('"') {
                    return Err(error(
                        file,
                        line,
                        format!(
                            "field {number} is in double quotes; Challenge 1 strings use single quotes"
                        ),
                    ));
                }
                if bare.contains(char::is_whitespace) {
                    return Err(error(
                        file,
                        line,
                        format!(
                            "field {number}, {bare:?}, has a blank in it; Challenge 1 separates fields with commas"
                        ),
                    ));
                }
                (Field::Bare(bare), &trimmed[end..])
            };
            fields.push(field);
            match after.strip_prefix(',') {
                Some(next) => rest = next,
                None => break,
            }
        }
        Ok(Record { file, line, fields })
    }

    /// The line the record is on.
    pub(super) fn line(&self) -> usize {
        self.line
    }

    /// An error about this record.
    pub(super) fn error(&self, message: impl Into<String>) -> PowerflowError {
        error(self.file, self.line, message)
    }

    /// Checks that a `what` record (such as "bus") has exactly `count` fields.
    pub(super) fn expect_len(&self, what: &str, count: usize) -> Result<(), PowerflowError> {
        if self.fields.len() == count {
            Ok(())
        } else {
            Err(self.error(format!(
                "a {what} record has {count} fields, not {}",
                self.fields.len()
            )))
        }
    }

    /// Field `index` (from 0), a number, for messages naming it `name` (a PSS/E field name such
    /// as `PG`).
    fn number(&self, index: usize, name: &str) -> Result<f64, PowerflowError> {
        let at = |what: String| self.error(format!("{name} (field {}) {what}", index + 1));
        match self.fields.get(index) {
            Some(Field::Bare(text)) => match text.parse::<f64>() {
                Ok(value) if value.is_finite() => Ok(value),
                _ => Err(at(format!("must be a finite number, not {text:?}"))),
            },
            Some(Field::Quoted(_)) => Err(at("must be a number, not a string".to_string())),
            None => Err(at("is missing".to_string())),
        }
    }

    /// Field `index`, a finite number.
    pub(super) fn float(&self, index: usize, name: &str) -> Result<f64, PowerflowError> {
        self.number(index, name)
    }

    /// Field `index`, a whole number from `min` to `max`.
    fn whole(&self, index: usize, name: &str, min: i64, max: i64) -> Result<i64, PowerflowError> {
        let value = self.number(index, name)?;
        if value.fract() == 0.0 && value >= min as f64 && value <= max as f64 {
            Ok(value as i64)
        } else {
            Err(self.error(format!(
                "{name} (field {}) must be a whole number from {min} to {max}, not {value}",
                index + 1
            )))
        }
    }

    /// Field `index`, a whole number of either sign.
    pub(super) fn int(&self, index: usize, name: &str) -> Result<i64, PowerflowError> {
        // 2^53: every whole number up to it is exact in an f64
        const MAX: i64 = 1 << 53;
        self.whole(index, name, -MAX, MAX)
    }

    /// Field `index`, a bus number or other label: a whole number from 1.
    pub(super) fn label(&self, index: usize, name: &str) -> Result<usize, PowerflowError> {
        // A bus number in PSS/E is at most 999,997; this only rules out values usize and f64
        //   can't both hold exactly
        let max = (1i64 << 53).min(i64::try_from(usize::MAX).unwrap_or(i64::MAX));
        self.whole(index, name, 1, max)
            .map(|value| usize::try_from(value).expect("checked to be in range"))
    }

    /// Field `index`, a count: a whole number from 0 that fits a `u32`.
    pub(super) fn count(&self, index: usize, name: &str) -> Result<u32, PowerflowError> {
        self.whole(index, name, 0, i64::from(u32::MAX))
            .map(|value| u32::try_from(value).expect("checked to be in range"))
    }

    /// Field `index`, a status: 1 in service, 0 out of service.
    pub(super) fn status(&self, index: usize, name: &str) -> Result<bool, PowerflowError> {
        self.whole(index, name, 0, 1).map(|value| value == 1)
    }

    /// Field `index`, a single-quoted string, trimmed of the blanks inside its quotes.
    pub(super) fn string(&self, index: usize, name: &str) -> Result<String, PowerflowError> {
        match self.fields.get(index) {
            Some(Field::Quoted(text)) => Ok(text.trim().to_string()),
            Some(Field::Bare(text)) => Err(self.error(format!(
                "{name} (field {}) must be a string in single quotes, not {text:?}",
                index + 1
            ))),
            None => Err(self.error(format!("{name} (field {}) is missing", index + 1))),
        }
    }

    /// Field `index`, a short identifier (A.1.4) such as a circuit or a load's ID: a string of
    /// one or two characters, each an upper case letter or a digit, once trimmed. The data pads
    /// them to two, as in `'1 '`.
    pub(super) fn short_id(&self, index: usize, name: &str) -> Result<String, PowerflowError> {
        let id = self.string(index, name)?;
        let valid = (1..=2).contains(&id.len())
            && id
                .bytes()
                .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit());
        if valid {
            Ok(id)
        } else {
            Err(self.error(format!(
                "{name} (field {}) must be one or two upper case letters or digits, not {id:?}",
                index + 1
            )))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(content: &str) -> Result<Record<'_>, PowerflowError> {
        Record::split(GoFile::Raw, 4, content)
    }

    fn message(result: Result<impl std::fmt::Debug, PowerflowError>) -> String {
        match result {
            Err(PowerflowError::GoParse {
                file: GoFile::Raw,
                line: 4,
                message,
            }) => message,
            other => panic!("expected a case.raw line 4 error, got {other:?}"),
        }
    }

    #[test]
    fn lines_end_with_lf_or_crlf() {
        let text = "a\r\nb\nc\r\nd";
        let got: Vec<_> = lines(text).collect();
        assert_eq!(got, [(1, "a"), (2, "b"), (3, "c"), (4, "d")]);
    }

    #[test]
    fn terminators_ends_and_records() {
        for line in [
            "0 / END OF BUS DATA, BEGIN LOAD DATA",
            "0 /END OF GNE DEVICE DATA",
            "0",
            "  0  ",
            "0 / end bus section",
        ] {
            assert_eq!(
                classify(GoFile::Raw, 1, line),
                Ok(Kind::Terminator),
                "{line}"
            );
        }
        assert_eq!(classify(GoFile::Raw, 1, "Q"), Ok(Kind::End));
        // A record whose first field is 0, such as a cost point at 0 MW, isn't a terminator
        assert_eq!(
            classify(GoFile::Rop, 1, "0.0, 1500.5"),
            Ok(Kind::Record("0.0, 1500.5"))
        );
        // A / in a string isn't a comment
        assert_eq!(
            classify(GoFile::Raw, 1, "1,'A/B', 1 / note"),
            Ok(Kind::Record("1,'A/B', 1"))
        );
        assert!(matches!(
            classify(GoFile::Raw, 4, "1,'open"),
            Err(PowerflowError::GoParse { line: 4, .. })
        ));
    }

    #[test]
    fn fields_and_their_types() {
        let r = record("   1,'1 ',   1.5, -2e-3, 0, 'TL', 'A NAME  ', 1.0").unwrap();
        r.expect_len("test", 8).unwrap();
        assert_eq!(r.line(), 4);
        assert_eq!(r.label(0, "I"), Ok(1));
        assert_eq!(r.short_id(1, "ID"), Ok("1".to_string()));
        assert_eq!(r.float(2, "PL"), Ok(1.5));
        assert_eq!(r.float(3, "QL"), Ok(-0.002));
        assert_eq!(r.status(4, "STATUS"), Ok(false));
        assert_eq!(r.short_id(5, "CKT"), Ok("TL".to_string()));
        assert_eq!(r.string(6, "NAME"), Ok("A NAME".to_string()));
        // A whole number may be written with a decimal point
        assert_eq!(r.int(7, "N1"), Ok(1));
        assert_eq!(r.count(7, "N1"), Ok(1));
        assert!(message(r.expect_len("bus", 13)).contains("a bus record has 13 fields, not 8"));
    }

    #[test]
    fn strictness() {
        // Blanks as separators, double quotes, empty fields and text after a quote
        assert!(message(record("1 2, 3")).contains("has a blank in it"));
        assert!(message(record("1, \"A\"")).contains("double quotes"));
        assert!(message(record("1,,3")).contains("field 2 is empty"));
        assert!(message(record("1, 2,")).contains("field 3 is empty"));
        assert!(message(record("'A' B, 1")).contains("text after its closing quote"));

        let r = record("1.5, 'x', NaN, inf, 2, -1, 'abc', 'a', ''").unwrap();
        assert!(message(r.label(0, "I")).contains("I (field 1) must be a whole number"));
        assert!(message(r.float(1, "PL")).contains("must be a number, not a string"));
        assert!(message(r.float(2, "PL")).contains("must be a finite number"));
        assert!(message(r.float(3, "PL")).contains("must be a finite number"));
        assert!(message(r.status(4, "STATUS")).contains("from 0 to 1, not 2"));
        assert!(message(r.count(5, "N1")).contains("from 0 to"));
        assert!(message(r.string(0, "NAME")).contains("in single quotes"));
        assert!(message(r.short_id(6, "CKT")).contains("one or two upper case"));
        assert!(message(r.short_id(7, "CKT")).contains("one or two upper case"));
        assert!(message(r.short_id(8, "CKT")).contains("one or two upper case"));
        assert!(message(r.float(9, "PL")).contains("is missing"));
    }
}
