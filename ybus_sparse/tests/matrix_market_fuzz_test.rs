//! Fuzz-style tests for the Matrix Market reader: truncated, malformed and oversized input
//! must give an `Err`, never a panic. The inputs are generated deterministically from a fixed
//! seed, so a failure always reproduces, and no extra tooling (cargo-fuzz, nightly) is needed.

use std::fmt::Debug;
use std::panic::{catch_unwind, AssertUnwindSafe};

use ybus_sparse::sparse_matrix::SparseMatrix;
use ybus_sparse::Scalar;

/// Valid files covering every supported header, comments, blank lines and duplicates. The
/// mutations start from these, so most mutated inputs get past the header.
const SEEDS: &[&str] = &[
    "%%MatrixMarket matrix coordinate real general\n% comment\n3 4 4\n1 1 1.5\n2 3 -2\n3 4 1e3\n1 1 0.5\n",
    "%%MatrixMarket matrix coordinate integer symmetric\n\n3 3 3\n1 1 4\n2 1 -1\n3 2 7\n",
    "%%MatrixMarket matrix coordinate real skew-symmetric\n2 2 1\n2 1 3.25\n",
    "%%MatrixMarket matrix coordinate pattern general\n2 3 2\n1 3\n2 1\n",
    "%%MATRIXMARKET Matrix Coordinate Integer General\n0 0 0\n",
];

/// Tokens that tend to break parsers: zero and negative indices, values at and past the
/// limits of each integer type, non-finite and malformed floats, and structural characters.
const TOKENS: &[&str] = &[
    "0",
    "-1",
    "1",
    "2147483647",
    "-2147483648",
    "4294967295",
    "9223372036854775807",
    "-9223372036854775808",
    "18446744073709551615",
    "18446744073709551616",
    "1e400",
    "-1e400",
    "nan",
    "inf",
    "1.5.5",
    "0x10",
    "",
    " ",
    "\t",
    "\n",
    "\r\n",
    "%",
    "%%MatrixMarket",
    "pattern",
    "skew-symmetric",
    "hermitian",
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

/// Parses `bytes` as a matrix of `T` and fails the test if parsing panics. A successful parse
/// must also survive a write and re-read unchanged (unless it holds NaN).
fn assert_no_panic<T: Scalar + Debug + std::fmt::Display>(bytes: &[u8]) {
    let result = catch_unwind(AssertUnwindSafe(|| {
        let parsed = SparseMatrix::<T>::read_matrix_market(bytes);
        // Writing iterates the CSR form, which can't be allocated for dimensions near
        //   usize::MAX (a deliberate panic in compressed storage), so only round-trip shapes
        //   that compressed storage can hold
        let small = |m: &SparseMatrix<T>| m.shape.0 < 1 << 20 && m.shape.1 < 1 << 20;
        if let Some(matrix) = parsed.as_ref().ok().filter(|m| small(m)) {
            let text = matrix.to_matrix_market();
            let reparsed = SparseMatrix::<T>::from_matrix_market(&text).unwrap();
            // NaN != NaN, so a matrix holding NaN can't compare equal to its round trip
            #[allow(clippy::eq_op)]
            let holds_nan = matrix != matrix;
            if !holds_nan {
                assert!(reparsed == *matrix, "round trip changed {text:?}");
            }
        }
    }));
    if result.is_err() {
        panic!(
            "parsing as {} panicked on input {:?}",
            std::any::type_name::<T>(),
            String::from_utf8_lossy(bytes)
        );
    }
}

/// Every element type the reader is likely to be used with, including unsigned and narrow
/// integers, where negation and overflow are the risks.
fn assert_no_panic_any_type(bytes: &[u8]) {
    assert_no_panic::<f64>(bytes);
    assert_no_panic::<f32>(bytes);
    assert_no_panic::<i64>(bytes);
    assert_no_panic::<i32>(bytes);
    assert_no_panic::<u32>(bytes);
    assert_no_panic::<u8>(bytes);
}

#[test]
fn matrix_market_fuzz_truncation_never_panics() {
    for seed in SEEDS {
        for end in 0..=seed.len() {
            assert_no_panic_any_type(&seed.as_bytes()[..end]);
        }
    }
}

#[test]
fn matrix_market_fuzz_mutation_never_panics() {
    let mut rng = Lcg(0x5eed_f0f0);
    for round in 0..20_000 {
        let seed = SEEDS[round % SEEDS.len()].as_bytes();
        let input = mutate(seed, &mut rng);
        assert_no_panic_any_type(&input);
    }
}

#[test]
fn matrix_market_fuzz_oversized_input_errors() {
    let header = "%%MatrixMarket matrix coordinate real general\n";
    let cases = [
        // Huge entry count with few entries: must not pre-allocate it, and must report the
        //   shortfall
        format!("{header}2 2 {}\n1 1 1\n", usize::MAX),
        format!("{header}2 2 1000000000000\n1 1 1\n"),
        // Huge dimensions are allowed as a shape, but entries must still be in range
        format!("{header}{} {} 1\n1 1 1\n", usize::MAX, usize::MAX),
        format!("{header}{} 1 1\n{} 1 1\n", usize::MAX, usize::MAX),
        // Sizes past usize
        format!("{header}18446744073709551616 1 0\n"),
        format!("{header}1 1 18446744073709551616\n"),
        // Index past usize, and past the shape
        format!("{header}2 2 1\n18446744073709551616 1 1\n"),
        format!("{header}2 2 1\n3 1 1\n"),
        // A very long line, and a very long value token
        format!("{header}2 2 1\n1 1 {}\n", "9".repeat(100_000)),
        format!("{header}2 2 1\n{}\n", "1 ".repeat(100_000)),
        // More entries than declared
        format!("{header}2 2 1\n{}", "1 1 1\n".repeat(10_000)),
    ];
    for case in &cases {
        assert_no_panic_any_type(case.as_bytes());
    }

    // Shapes this large have to fail cleanly in the reader rather than later, since a CSR
    //   form of them could never be allocated. Only the entry-count cases are required to
    //   error; the rest only must not panic.
    for case in &cases[..2] {
        assert!(SparseMatrix::<f64>::from_matrix_market(case).is_err());
    }
}

#[test]
fn matrix_market_fuzz_invalid_utf8_errors() {
    let mut bytes = b"%%MatrixMarket matrix coordinate real general\n2 2 1\n1 1 ".to_vec();
    bytes.extend_from_slice(&[0xff, 0xfe, b'\n']);
    assert_no_panic_any_type(&bytes);
    assert!(SparseMatrix::<f64>::read_matrix_market(bytes.as_slice()).is_err());
    assert!(SparseMatrix::<f64>::read_matrix_market(&[0xc3u8][..]).is_err());
}

#[test]
fn matrix_market_skew_symmetric_values_without_a_negation_error() {
    // The mirrored entry of a skew-symmetric file is the negated value. This used to be
    //   computed as 0 - value, which overflowed (panicking in debug builds) for unsigned types
    //   and for i32::MIN.
    let skew = |value: &str| {
        format!("%%MatrixMarket matrix coordinate integer skew-symmetric\n2 2 1\n2 1 {value}\n")
    };
    assert!(SparseMatrix::<u32>::from_matrix_market(&skew("3")).is_err());
    assert!(SparseMatrix::<i32>::from_matrix_market(&skew("-2147483648")).is_err());
    let pattern = "%%MatrixMarket matrix coordinate pattern skew-symmetric\n2 2 1\n2 1\n";
    assert!(SparseMatrix::<u8>::from_matrix_market(pattern).is_err());

    // Values that do have a negation still mirror correctly, including an explicit '+'
    let local = SparseMatrix::<i32>::from_matrix_market(&skew("2147483647")).unwrap();
    assert!(local[(0, 1)] == -2147483647);
    let local = SparseMatrix::<f64>::from_matrix_market(&skew("+2.5")).unwrap();
    assert!(local[(1, 0)] == 2.5 && local[(0, 1)] == -2.5);
    let local = SparseMatrix::<i8>::from_matrix_market(&skew("-127")).unwrap();
    assert!(local[(0, 1)] == 127);
}
