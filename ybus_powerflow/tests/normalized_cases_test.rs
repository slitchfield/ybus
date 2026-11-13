//! Tests for the MATPOWER cases normalized by `scripts/normalize-cases.sh`: every case
//! MATPOWER's `loadcase` can load, including those that compute their data with code, written
//! as literals the reader accepts. Skips when the cases haven't been normalized, unless
//! `YBUS_REQUIRE_CASES=1`.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{locate_normalized, Located, NORMALIZED_MARKER};
use ybus_powerflow::case::Case;

/// The normalized cases' directory, or `None` to skip.
fn normalized_dir() -> Option<PathBuf> {
    Some(common::normalized_case("case14")?.parent()?.to_path_buf())
}

/// MATPOWER's own data directory, or `None` to skip.
fn matpower_data_dir() -> Option<PathBuf> {
    Some(common::matpower_case("case14")?.parent()?.to_path_buf())
}

/// The `.m` files in `dir`, sorted.
fn case_files(dir: &Path) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "m"))
        .collect();
    files.sort();
    files
}

/// Parses the literals at the start of a MATPOWER case that computes data with code, by
/// cutting the text where the code begins. This gives the file's raw values (in kW, ohms, ...)
/// before its own code converts them.
fn raw_literals(path: &Path, code_starts_with: &str) -> Case {
    let text = fs::read_to_string(path).unwrap();
    let cut = text
        .find(code_starts_with)
        .unwrap_or_else(|| panic!("no {code_starts_with:?} in {}", path.display()));
    Case::parse(&text[..cut]).unwrap()
}

#[test]
fn every_normalized_case_parses() {
    let Some(dir) = normalized_dir() else {
        return;
    };
    let files = case_files(&dir);
    for path in &files {
        if let Err(err) = Case::from_path(path) {
            panic!("{}: {err}", path.display());
        }
    }
    // Including cases the reader rejects in their original form: rescaled with code, and
    //   with arithmetic in their literals
    let names: Vec<String> = files
        .iter()
        .map(|path| path.file_stem().unwrap().to_string_lossy().into_owned())
        .collect();
    for name in ["case10ba", "case69", "case8387pegase", "case533mt_hi"] {
        assert!(names.iter().any(|n| n == name), "{name} wasn't normalized");
    }
}

#[test]
fn normalizing_a_literal_case_changes_nothing() {
    // For every case the reader accepts directly, the normalized copy must parse to exactly the
    //   same case. This checks the writer loses no precision (%.17g round-trips every double)
    //   and that loadcase adds nothing the reader would see differently
    let (Some(normalized), Some(original)) = (normalized_dir(), matpower_data_dir()) else {
        return;
    };
    let mut compared = 0;
    for path in case_files(&original) {
        let Ok(direct) = Case::from_path(&path) else {
            continue;
        };
        let copy = Case::from_path(normalized.join(path.file_name().unwrap())).unwrap();
        assert!(copy == direct, "{} changed when normalized", path.display());
        compared += 1;
    }
    assert!(compared >= 50, "only {compared} cases compared");
}

#[test]
fn case10ba_values_match_its_own_conversion_code() {
    // case10ba states loads in kW and impedances in ohms, then converts them with:
    //   Vbase = mpc.bus(1, BASE_KV) * 1e3;
    //   Sbase = mpc.baseMVA * 1e6;
    //   mpc.branch(:, [BR_R BR_X]) = mpc.branch(:, [BR_R BR_X]) / (Vbase^2 / Sbase);
    //   mpc.bus(:, [PD, QD]) = mpc.bus(:, [PD, QD]) / 1e3;
    // Recompute that from the raw literals and compare with what loadcase produced
    let (Some(normalized), Some(original)) = (normalized_dir(), matpower_data_dir()) else {
        return;
    };
    let raw = raw_literals(&original.join("case10ba.m"), "[PQ, PV, REF");
    let converted = Case::from_path(normalized.join("case10ba.m")).unwrap();

    let vbase = raw.bus[0].base_kv * 1e3;
    let sbase = raw.base_mva * 1e6;
    let z_base = vbase * vbase / sbase;
    for (raw_bus, bus) in raw.bus.iter().zip(&converted.bus) {
        assert_eq!(bus.pd_mw, raw_bus.pd_mw / 1e3);
        assert_eq!(bus.qd_mvar, raw_bus.qd_mvar / 1e3);
    }
    for (raw_branch, branch) in raw.branch.iter().zip(&converted.branch) {
        assert_eq!(branch.r, raw_branch.r / z_base);
        assert_eq!(branch.x, raw_branch.x / z_base);
    }
    // The values the reader refuses to load as-is differ by a factor of 1000 from the real ones
    assert_eq!(raw.bus[1].pd_mw, 1840.0);
    assert_eq!(converted.bus[1].pd_mw, 1.84);
}

#[test]
fn case8387pegase_code_is_evaluated_not_assumed() {
    // case8387pegase sets `fixed = 0;` and only pins the limits of 615 generators inside
    //   `if fixed ... end`. MATPOWER doesn't run that block by default, so those generators keep
    //   their infinite limits. Guessing at the code's effect instead of evaluating it would
    //   get this wrong
    let Some(normalized) = normalized_dir() else {
        return;
    };
    let case = Case::from_path(normalized.join("case8387pegase.m")).unwrap();
    assert_eq!(case.bus.len(), 8387);
    let unlimited = case
        .generator
        .iter()
        .filter(|g| {
            [g.qmin_mvar, g.qmax_mvar, g.pmin_mw, g.pmax_mw]
                .iter()
                .all(|limit| limit.is_infinite())
        })
        .count();
    assert_eq!(unlimited, 615);
}

#[test]
fn case533mt_arithmetic_is_evaluated() {
    // case533mt_hi writes `mpc.baseMVA = 50/3;`, which loadcase evaluates
    let Some(normalized) = normalized_dir() else {
        return;
    };
    let case = Case::from_path(normalized.join("case533mt_hi.m")).unwrap();
    assert_eq!(case.base_mva, 50.0 / 3.0);
}

/// A fake normalized directory made from `source`, with `cases`, optionally complete.
fn fake_normalized(name: &str, source: &Path, cases: &[&str], complete: bool) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("normalized_cases_test")
        .join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let source = fs::canonicalize(source).unwrap();
    fs::write(
        dir.join("SOURCE"),
        format!("source {}\nMATPOWER 8.1\n", source.display()),
    )
    .unwrap();
    for case in cases {
        fs::write(dir.join(format!("{case}.m")), "function mpc = x\n").unwrap();
    }
    if complete {
        fs::write(dir.join(NORMALIZED_MARKER), "").unwrap();
    }
    dir
}

#[test]
fn locate_normalized_finds_complete_output_from_the_same_source() {
    let source = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let dir = fake_normalized("complete", source, &["case14"], true);
    assert_eq!(
        locate_normalized("case14", source, &dir),
        Located::Found {
            path: dir.join("case14.m"),
            warning: None
        }
    );
}

#[test]
fn locate_normalized_needs_the_marker() {
    // An interrupted run is moved into place only when complete, but a stray directory without
    //   the marker mustn't count either
    let source = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let dir = fake_normalized("incomplete", source, &["case14"], false);
    match locate_normalized("case14", source, &dir) {
        Located::Missing(reason) => {
            assert!(reason.contains("scripts/normalize-cases.sh"), "{reason}");
            assert!(reason.contains("Octave"), "{reason}");
        }
        other => panic!("expected Missing, got {other:?}"),
    }
}

#[test]
fn locate_normalized_rejects_output_from_another_source() {
    // Output made from a different MATPOWER directory (say, before a version bump) is stale
    let made_from = Path::new(env!("CARGO_MANIFEST_DIR"));
    let dir = fake_normalized("stale", made_from, &["case14"], true);
    let current = Path::new(env!("CARGO_TARGET_TMPDIR"));
    match locate_normalized("case14", current, &dir) {
        Located::Missing(reason) => assert!(reason.contains("weren't made from"), "{reason}"),
        other => panic!("expected Missing, got {other:?}"),
    }
}

#[test]
fn locate_normalized_reports_a_missing_case_as_broken() {
    let source = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let dir = fake_normalized("missing-case", source, &["case14"], true);
    assert!(matches!(
        locate_normalized("case_nonexistent", source, &dir),
        Located::Broken(_)
    ));
}
