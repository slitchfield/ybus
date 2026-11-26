// What the GO case.raw reader refuses: each edit of the hand-written fixture must give a GoParse
//   error naming case.raw, the line the problem is on, and what's wrong. The reader is strict
//   about the subset Challenge 1 allows, so full PSS/E syntax outside it is refused too.

use std::fs;
use std::path::PathBuf;

use ybus_powerflow::case::go::{GoFile, RawCase};
use ybus_powerflow::powerflow_error::PowerflowError;

fn fixture_lines() -> Vec<String> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/go/small/case.raw");
    let text = fs::read_to_string(path).expect("the fixture is committed");
    text.lines().map(str::to_string).collect()
}

/// The fixture with its 1-based line `line` replaced by `text`
fn replaced(line: usize, text: &str) -> String {
    let mut lines = fixture_lines();
    lines[line - 1] = text.to_string();
    lines.join("\n") + "\n"
}

/// The fixture with `from` replaced by `to` on its 1-based line `line`, which must contain it
fn edited(line: usize, from: &str, to: &str) -> String {
    let lines = fixture_lines();
    assert!(
        lines[line - 1].contains(from),
        "line {line} has no {from:?}"
    );
    replaced(line, &lines[line - 1].replacen(from, to, 1))
}

/// The fixture's first `count` lines
fn truncated(count: usize) -> String {
    fixture_lines()[..count].join("\n") + "\n"
}

/// Checks that `text` is refused at `line` with a message containing `fragment`
fn refused(what: &str, text: &str, line: usize, fragment: &str) {
    match RawCase::parse(text) {
        Err(PowerflowError::GoParse {
            file: GoFile::Raw,
            line: got,
            message,
        }) => {
            assert_eq!(got, line, "{what}: wrong line, message {message:?}");
            assert!(
                message.contains(fragment),
                "{what}: {message:?} doesn't contain {fragment:?}"
            );
        }
        other => panic!("{what}: expected a case.raw error at line {line}, got {other:?}"),
    }
}

#[test]
fn malformed_fields_name_their_line() {
    refused(
        "a field too few",
        &edited(5, "0.95000,1.10000,0.90000", "0.95000,1.10000"),
        5,
        "a bus record has 13 fields, not 12",
    );
    refused(
        "fields separated by blanks",
        &replaced(13, "     3 '1 ' 1 0.000 15.000"),
        13,
        "has a blank in it",
    );
    refused(
        "a double-quoted string",
        &edited(9, "'1 '", "\"1 \""),
        9,
        "double quotes",
    );
    refused(
        "an unclosed quote",
        &edited(4, "'BUS ONE     '", "'BUS ONE     "),
        4,
        "closing quote is missing",
    );
    refused(
        "an empty field",
        &edited(9, "'1 ',1,", "'1 ',,"),
        9,
        "field 3 is empty",
    );
    refused(
        "text after a string",
        &edited(9, "'1 '", "'1 'X"),
        9,
        "text after its closing quote",
    );
    refused(
        "NaN",
        &edited(16, "    80.000", "       NaN"),
        16,
        "PG (field 3) must be a finite number",
    );
    refused(
        "a string for a number",
        &edited(16, "    80.000", "      '80'"),
        16,
        "PG (field 3) must be a number, not a string",
    );
    refused(
        "a bus number that isn't whole",
        &edited(20, "     1,     2,'1 '", "   1.5,     2,'1 '"),
        20,
        "I (field 1) must be a whole number",
    );
    refused(
        "a status of 2",
        &edited(9, "'1 ',1,", "'1 ',2,"),
        9,
        "STATUS (field 3) must be a whole number from 0 to 1, not 2",
    );
    refused(
        "IDE 5",
        &edited(6, "115.0000,1,", "115.0000,5,"),
        6,
        "IDE (field 4) is 5",
    );
    refused(
        "an ID that isn't short",
        &edited(10, "'2 '", "'2a'"),
        10,
        "one or two upper case letters or digits",
    );
}

#[test]
fn the_case_identification_is_checked() {
    refused(
        "REV 32",
        &edited(1, " 33,", " 32,"),
        1,
        "REV is 32; only PSS/E revision 33",
    );
    refused(
        "SBASE 0",
        &edited(1, "100.00", "0.00"),
        1,
        "SBASE must be positive",
    );
    refused("an empty file", "", 1, "the file is empty");
    refused(
        "two lines",
        &truncated(2),
        3,
        "ends within its three-line case identification",
    );
}

#[test]
fn transformers_are_four_lines_of_two_windings() {
    refused(
        "K 1",
        &edited(24, "     3,     0,'1 '", "     3,     1,'1 '"),
        24,
        "K (field 3) is 1: three-winding transformers aren't part of Challenge 1",
    );
    refused(
        "a terminator in a transformer",
        &replaced(26, "0 / END OF TRANSFORMER DATA"),
        26,
        "a transformer record is four lines, but its third line ends",
    );
    refused(
        "a short second line",
        &edited(25, ",   100.00", ""),
        25,
        "a transformer's second line has 3 fields, not 2",
    );
    refused(
        "the end of the file in a transformer",
        &truncated(25),
        24,
        "the file ends within this transformer record, before its third line",
    );
}

#[test]
fn the_file_must_reach_its_switched_shunts() {
    refused(
        "Q in the line section",
        &replaced(23, "Q"),
        23,
        "the file ends before its switched shunt data",
    );
    refused(
        "the end of the file in the zone section",
        &truncated(42),
        42,
        "the file ends before its switched shunt data",
    );
    refused(
        "the end of the file before the switched shunts' terminator",
        &truncated(48),
        48,
        "the file ends before its switched shunt data",
    );
}
