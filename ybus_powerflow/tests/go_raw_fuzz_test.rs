//! Fuzz-style tests for the GO Challenge 1 `case.raw` reader: truncated, malformed and oversized
//! input must give an `Err`, never a panic. As in `matpower_fuzz_test.rs`, inputs are generated
//! deterministically from a fixed seed, so a failure always reproduces, and the seeds are the
//! hand-written fixture, since the datasets aren't committed.

use std::fs;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::PathBuf;

use ybus_powerflow::case::go::RawCase;

/// The fixture, with LF and with CRLF line endings
fn seeds() -> Vec<String> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/go/small/case.raw");
    let lf = fs::read_to_string(path).expect("the fixture is committed");
    let crlf = lf.replace('\n', "\r\n");
    vec![lf, crlf]
}

/// Tokens that tend to break this reader: its structural characters, terminators and the file
/// end, numbers at and past the limits of the integer and float types, and awkward text.
const TOKENS: &[&str] = &[
    "0",
    "0 /",
    "Q",
    "/",
    ",",
    ",,",
    "'",
    "''",
    "\"",
    " ",
    "\t",
    "\n",
    "\r\n",
    "\r",
    "-1",
    "1.5",
    "4294967296",
    "9007199254740993",
    "18446744073709551616",
    "1e30",
    "1e400",
    "-1e400",
    "NaN",
    "inf",
    "1.5.5",
    "'AB'",
    "'ABC'",
    "'é'",
    "\u{0}",
    "\u{feff}",
];

/// The LCG used by the other test files.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next() % bound as u64) as usize
    }
}

/// Applies one to three random edits to `input`: replacing, inserting or deleting a byte,
/// inserting a token, duplicating or deleting a line, or truncating.
fn mutate(input: &[u8], rng: &mut Lcg) -> Vec<u8> {
    let mut bytes = input.to_vec();
    for _ in 0..1 + rng.below(3) {
        let at = rng.below(bytes.len() + 1);
        match rng.below(7) {
            0 if at < bytes.len() => bytes[at] = rng.next() as u8,
            1 => bytes.insert(at, rng.next() as u8),
            2 if at < bytes.len() => {
                bytes.remove(at);
            }
            3 | 4 => {
                let token = TOKENS[rng.below(TOKENS.len())].as_bytes();
                bytes.splice(at..at, token.iter().copied());
            }
            5 => {
                let text = String::from_utf8_lossy(&bytes).into_owned();
                let mut lines: Vec<&str> = text.lines().collect();
                if !lines.is_empty() {
                    let line = rng.below(lines.len());
                    if rng.below(2) == 0 {
                        lines.insert(line, lines[line]);
                    } else {
                        lines.remove(line);
                    }
                }
                bytes = lines.join("\n").into_bytes();
            }
            _ => bytes.truncate(at),
        }
    }
    bytes
}

/// Parses `bytes` through `RawCase::from_bytes`, the path `from_path` uses, and fails the test
/// if parsing panics.
fn assert_no_panic(bytes: &[u8]) {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let _ = RawCase::from_bytes(bytes);
    }));
    if result.is_err() {
        panic!(
            "parsing panicked on input {:?}",
            String::from_utf8_lossy(bytes)
        );
    }
}

#[test]
fn go_raw_seeds_parse() {
    // The mutations only probe the reader if they start from valid input
    for seed in seeds() {
        if let Err(err) = RawCase::parse(&seed) {
            panic!("seed doesn't parse: {err}");
        }
    }
}

#[test]
fn go_raw_fuzz_truncation_never_panics() {
    for seed in seeds() {
        for end in 0..=seed.len() {
            assert_no_panic(&seed.as_bytes()[..end]);
        }
    }
}

#[test]
fn go_raw_fuzz_mutation_never_panics() {
    let seeds = seeds();
    let mut rng = Lcg(0x60_ca5e);
    for round in 0..20_000 {
        let input = mutate(seeds[round % seeds.len()].as_bytes(), &mut rng);
        assert_no_panic(&input);
    }
}

#[test]
fn go_raw_fuzz_oversized_input_never_panics() {
    let head = "0, 100.00, 33, 0, 0, 60.00\nline 2\nline 3\n";
    let bus = "1,'',230.0,3,1,1,1,1.0,0.0,1.1,0.9,1.1,0.9\n";
    let cases = [
        // Bus numbers and step counts too large for usize, u32 or exactly in f64
        format!("{head}18446744073709551616,'',230.0,3,1,1,1,1.0,0.0,1.1,0.9,1.1,0.9\n"),
        format!("{head}1e300,'',230.0,3,1,1,1,1.0,0.0,1.1,0.9,1.1,0.9\n"),
        // Very long lines, fields and strings, and many records
        format!("{head}{}\n", "9".repeat(100_000)),
        format!("{head}{}\n", "1,".repeat(100_000)),
        format!("{head}'{}'\n", "A".repeat(100_000)),
        format!("{head}{}0\n", bus.repeat(10_000)),
        // A quote that never closes, and terminators only
        format!("{head}'{}", "x".repeat(10_000)),
        format!("{head}{}", "0\n".repeat(10_000)),
    ];
    for case in &cases {
        assert_no_panic(case.as_bytes());
    }
}
