//! Fuzz-style tests for the MATPOWER case reader: truncated, malformed and oversized input must
//! give an `Err`, never a panic, since users may hand it their own files. Inputs are generated
//! deterministically from a fixed seed, so a failure always reproduces, and no extra tooling
//! (cargo-fuzz, nightly) is needed. The seeds are synthetic: MATPOWER's own case files can't
//! be committed.

use std::panic::{catch_unwind, AssertUnwindSafe};

use ybus_powerflow::case::Case;

/// Valid case files covering the syntax the reader accepts: comments, cell arrays, strings,
/// `...` continuations, commas, extra columns, both cost models and ignored fields.
const SEEDS: &[&str] = &[
    "function mpc = seed1\n\
     %% a comment\n\
     mpc.version = '2';\n\
     mpc.baseMVA = 100;\n\
     mpc.bus = [\n\
     \t1\t3\t0\t0\t0\t0\t1\t1.06\t0\t135\t1\t1.06\t0.94;\n\
     \t2\t1\t21.7\t12.7\t0\t0\t1\t1\t-4.98\t135\t1\t1.06\t0.94;\n\
     ];\n\
     mpc.gen = [\n\
     \t1\t232.4\t-16.9\t10\t0\t1.06\t100\t1\t332.4\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n\
     ];\n\
     mpc.branch = [\n\
     \t1\t2\t0.01938\t0.05917\t0.0528\t0\t0\t0\t0\t0\t1\t-360\t360;\n\
     ];\n\
     mpc.gencost = [\n\
     \t2\t0\t0\t3\t0.0430292599\t20\t0;\n\
     ];\n\
     mpc.bus_name = {\n\
     \t'Bus 1 % HV';\n\
     \t'Bus ''2''';\n\
     };\n",
    "function mpc = seed2\n\
     mpc.version = '2';\n\
     mpc.baseMVA = 50;\n\
     mpc.bus = [ 1, 3, 0, 0, 0, 0, 1, 1, 0, 23, 1, 1.1, 0.9, 99, 99; 2, 2, 1, 2, 0, 0, 1, 1, ...\n\
     \t0, 23, 1, Inf, 0.9, 99, 99 ];\n\
     mpc.gen = [ 2 10 0 10 -10 1 100 1 20 0 ];\n\
     mpc.branch = [ 1 2 1e-2 0.1 0 0 0 0 0.978 -3 1 ];\n\
     mpc.gencost = [ 1 0 0 2 0 0 20 400; 2 1 1 2 30 0 0 0 ];\n\
     mpc.gentype = { 'ST' };\n\
     mpc.dcline = [ 1 2 3 ];\n\
     end\n",
    "function s = seed3\n\
     %{\n\
     block comment with mpc.bus(:, 3) = 0; inside\n\
     %}\n\
     # Octave comment\n\
     s.version = \"2\"; s.baseMVA = 100;\n\
     s.if.lims = [ 1 -100 100 ];\n\
     s.note = sprintf('made %d', 2026); s.t = [1 2]'; % it's a comment\n\
     s.bus = [ 1 3 0 0 0 0 1 1 0 135 1 1.05 0.95; 2 1 50 10 0 0 1 1 0 135 1 1.05 0.95 ];\n\
     s.gen = [ 1 0 0 1 -1 1 100 1 200 0 ];\n\
     s.branch = [ 1 2 0.01 0.1 0.02 0 0 0 0 0 1 ];\n\
     s.bus_name = { \"Bus \"\"1\"\"\"; 'two' };\n",
];

/// Tokens that tend to break parsers: numbers at and past the limits of `usize` and `f64`,
/// non-finite and malformed numbers, structural characters of the syntax, and awkward text.
const TOKENS: &[&str] = &[
    "0",
    "function s = x",
    "s.",
    "#{",
    "#}",
    "\"\"",
    "(",
    ")",
    "-1",
    "1.5",
    "4294967296",
    "18446744073709551615",
    "18446744073709551616",
    "1e30",
    "1e300",
    "1e400",
    "-1e400",
    "NaN",
    "Inf",
    "-Inf",
    "50/3",
    "1.5.5",
    "",
    " ",
    "\t",
    "\n",
    "\r\n",
    "\r",
    ";",
    ",",
    "'",
    "''",
    "\"",
    "[",
    "]",
    "{",
    "}",
    "=",
    "...",
    "%",
    "%{",
    "%}",
    "#",
    "mpc.",
    "mpc.bus",
    "mpc.gencost = [",
    "function",
    "end",
    "\u{0}",
    "\u{feff}",
    "é",
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

/// Parses `bytes` and fails the test if parsing panics. Bytes that aren't UTF-8 go through
/// `Case::from_bytes`, the same path `Case::from_path` uses.
fn assert_no_panic(bytes: &[u8]) {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let _ = Case::from_bytes(bytes);
    }));
    if result.is_err() {
        panic!(
            "parsing panicked on input {:?}",
            String::from_utf8_lossy(bytes)
        );
    }
}

#[test]
fn matpower_seeds_parse() {
    // The mutations only probe the reader if they start from valid input
    for seed in SEEDS {
        if let Err(err) = Case::parse(seed) {
            panic!("seed doesn't parse: {err}\n{seed}");
        }
    }
}

#[test]
fn matpower_fuzz_truncation_never_panics() {
    for seed in SEEDS {
        for end in 0..=seed.len() {
            assert_no_panic(&seed.as_bytes()[..end]);
        }
    }
}

#[test]
fn matpower_fuzz_mutation_never_panics() {
    let mut rng = Lcg(0x5eed_ca5e);
    for round in 0..20_000 {
        let seed = SEEDS[round % SEEDS.len()].as_bytes();
        let input = mutate(seed, &mut rng);
        assert_no_panic(&input);
    }
}

#[test]
fn matpower_fuzz_oversized_input_never_panics() {
    let header = "function mpc = big\nmpc.version = '2';\nmpc.baseMVA = 100;\n";
    let bus = "mpc.bus = [ 1 3 0 0 0 0 1 1 0 135 1 1.05 0.95 ];\n";
    let gen = "mpc.gen = [ 1 0 0 1 -1 1 100 1 2 0 ];\n";
    let branch = "mpc.branch = [ 1 1 0 0.1 0 0 0 0 0 0 1 ];\n";
    let cases = [
        // A cost N so large that 4 + N, or 4 + 2N, overflows usize
        format!("{header}{bus}{gen}{branch}mpc.gencost = [ 2 0 0 18446744073709551615 1 ];\n"),
        format!("{header}{bus}{gen}{branch}mpc.gencost = [ 1 0 0 9223372036854775808 1 ];\n"),
        format!("{header}{bus}{gen}{branch}mpc.gencost = [ 2 0 0 1e300 1 ];\n"),
        // Labels too large for usize, which a cast would saturate rather than reject
        format!("{header}mpc.bus = [ 1e30 3 0 0 0 0 1 1 0 135 1 1.05 0.95 ];\n{gen}{branch}"),
        format!(
            "{header}mpc.bus = [ 18446744073709551616 3 0 0 0 0 1 1 0 135 1 1 1 ];\n{gen}{branch}"
        ),
        // Very long lines and tokens, and many rows
        format!("{header}mpc.baseMVA = {};\n", "9".repeat(100_000)),
        format!("{header}mpc.bus = [ {} ];\n", "1 ".repeat(100_000)),
        format!(
            "{header}mpc.bus = [\n{}];\n{gen}{branch}",
            "1 1 0 0 0 0 1 1 0 135 1 1 1;\n".repeat(10_000)
        ),
        // A literal that never closes, and deep runs of brackets
        format!("{header}mpc.bus = [ {}", "[".repeat(10_000)),
        format!("{header}mpc.bus_name = {{ {}", "'a' ".repeat(10_000)),
    ];
    for case in &cases {
        assert_no_panic(case.as_bytes());
    }
}
