//! Tests for the MATPOWER case locator in `common`, plus a check that the cases the power flow
//! tests use are available whenever `YBUS_REQUIRE_CASES=1` (as in CI).

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{
    installed_version, locate, parse_require, pinned_version, resolve, Located, VERIFIED_MARKER,
};

/// The pinned version the locator tests pass in; any value works, as long as it's consistent.
const PINNED: &str = "8.1";

/// A fresh, empty directory under Cargo's per-target scratch space, unique to `name`.
fn scratch_dir(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("matpower_cases_test")
        .join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// A fake MATPOWER installation holding `cases`, optionally with the verified marker.
fn fake_install(name: &str, cases: &[&str], verified: bool) -> PathBuf {
    let dir = scratch_dir(name);
    fs::create_dir_all(dir.join("data")).unwrap();
    for case in cases {
        fs::write(
            dir.join("data").join(format!("{case}.m")),
            "function mpc = x\n",
        )
        .unwrap();
    }
    if verified {
        fs::write(dir.join(VERIFIED_MARKER), "").unwrap();
    }
    dir
}

/// Records `version` in a fake installation's `lib/mpver.m`, in MATPOWER's own layout,
/// including the doc comment that mentions ``Version`` and must not be mistaken for it.
fn set_version(dir: &Path, version: &str) {
    fs::create_dir_all(dir.join("lib")).unwrap();
    let mpver = format!(
        "function rv = mpver(varargin)\n\
         % ``Version``, ``Release``, and ``Date`` *(all char arrays)*.\n\
         v{{1}} = struct(  'Name',     'MATPOWER', ... \n\
         \x20               'Version',  '{version}', ...\n\
         \x20               'Release',  '', ...\n\
         \x20               'Date',     '12-Jul-2025' );\n"
    );
    fs::write(dir.join("lib").join("mpver.m"), mpver).unwrap();
}

fn found(path: PathBuf) -> Located {
    Located::Found {
        path,
        warning: None,
    }
}

#[test]
fn finds_a_case_in_a_verified_download() {
    let dir = fake_install("verified", &["case14"], true);
    assert_eq!(
        locate("case14", None, &dir, PINNED),
        found(dir.join("data").join("case14.m"))
    );
}

#[test]
fn an_unverified_download_is_missing_not_found() {
    // An interrupted extraction leaves files but no marker; it must not be used
    let dir = fake_install("unverified", &["case14"], false);
    assert!(matches!(
        locate("case14", None, &dir, PINNED),
        Located::Missing(_)
    ));
}

#[test]
fn no_download_is_missing_with_a_pointer_to_the_fetch_script() {
    let dir = scratch_dir("absent").join("matpower8.1");
    match locate("case14", None, &dir, PINNED) {
        Located::Missing(reason) => {
            assert!(reason.contains("scripts/fetch-matpower.sh"), "{reason}");
            assert!(reason.contains("MATPOWER_DIR"), "{reason}");
        }
        other => panic!("expected Missing, got {other:?}"),
    }
}

#[test]
fn a_verified_download_without_the_case_is_broken() {
    let dir = fake_install("verified-incomplete", &["case14"], true);
    assert!(matches!(
        locate("case30", None, &dir, PINNED),
        Located::Broken(_)
    ));
}

#[test]
fn matpower_dir_takes_precedence_and_needs_no_marker() {
    let env_dir = fake_install("env-dir", &["case30"], false);
    set_version(&env_dir, PINNED);
    let fetched = fake_install("env-dir-fetched", &["case30"], true);
    assert_eq!(
        locate("case30", Some(&env_dir), &fetched, PINNED),
        found(env_dir.join("data").join("case30.m"))
    );
}

#[test]
fn a_matpower_dir_of_another_version_is_used_with_a_warning() {
    let env_dir = fake_install("env-dir-old", &["case30"], false);
    set_version(&env_dir, "7.1");
    match locate(
        "case30",
        Some(&env_dir),
        &scratch_dir("env-dir-old-fetched"),
        PINNED,
    ) {
        Located::Found {
            path,
            warning: Some(warning),
        } => {
            assert_eq!(path, env_dir.join("data").join("case30.m"));
            assert!(warning.contains("is MATPOWER 7.1"), "{warning}");
            assert!(warning.contains("generated with MATPOWER 8.1"), "{warning}");
            assert!(warning.contains("your responsibility"), "{warning}");
        }
        other => panic!("expected Found with a warning, got {other:?}"),
    }
}

#[test]
fn a_matpower_dir_without_a_version_is_used_with_a_warning() {
    let env_dir = fake_install("env-dir-unversioned", &["case30"], false);
    match locate(
        "case30",
        Some(&env_dir),
        &scratch_dir("env-dir-unversioned-fetched"),
        PINNED,
    ) {
        Located::Found {
            warning: Some(warning),
            ..
        } => assert!(
            warning.contains("can't determine the MATPOWER version"),
            "{warning}"
        ),
        other => panic!("expected Found with a warning, got {other:?}"),
    }
}

#[test]
fn installed_version_parses_mpver() {
    let dir = scratch_dir("mpver-parse");
    set_version(&dir, "8.1");
    let mpver = fs::read_to_string(dir.join("lib").join("mpver.m")).unwrap();
    assert_eq!(installed_version(&mpver).as_deref(), Some("8.1"));
    // A comment mentioning 'Version' in quotes is skipped
    assert_eq!(
        installed_version("% 'Version', 'wrong'\nv = struct('Version', '7.0');").as_deref(),
        Some("7.0")
    );
    assert_eq!(installed_version("no version here"), None);
}

#[test]
fn installed_version_reads_the_real_mpver() {
    // Checks the parser against MATPOWER's actual mpver.m whenever cases are available
    let Some(case) = common::matpower_case("case14") else {
        return;
    };
    let install = case.parent().and_then(Path::parent).unwrap();
    let mpver = fs::read_to_string(install.join("lib").join("mpver.m")).unwrap();
    let parsed = installed_version(&mpver);
    assert!(
        parsed.is_some(),
        "no version parsed from {}/lib/mpver.m; has MATPOWER changed its layout?",
        install.display()
    );
    // A MATPOWER_DIR installation of another version is allowed (the locator only warns), so
    //   the exact version is only checked for the fetch script's download of the pinned release
    if std::env::var_os(common::MATPOWER_DIR_VAR).is_none() {
        let root = common::workspace_root();
        assert_eq!(parsed, Some(common::pinned_matpower_version(&root)));
    }
}

#[test]
fn a_matpower_dir_without_the_case_is_broken_even_if_fetched() {
    // An explicitly chosen directory that's wrong is a misconfiguration, so no fallback
    let env_dir = scratch_dir("env-dir-empty");
    let fetched = fake_install("env-dir-empty-fetched", &["case118"], true);
    match locate("case118", Some(&env_dir), &fetched, PINNED) {
        Located::Broken(reason) => assert!(reason.contains("MATPOWER_DIR="), "{reason}"),
        other => panic!("expected Broken, got {other:?}"),
    }
}

#[test]
fn parse_require_accepts_only_0_and_1() {
    assert!(!parse_require(None));
    assert!(!parse_require(Some("")));
    assert!(!parse_require(Some("0")));
    assert!(parse_require(Some("1")));
}

#[test]
#[should_panic(expected = "must be 0 or 1")]
fn parse_require_rejects_other_values() {
    parse_require(Some("yes"));
}

#[test]
fn resolve_skips_missing_cases_unless_required() {
    let path = PathBuf::from("data/case14.m");
    assert_eq!(resolve(found(path.clone()), true), Some(path.clone()));
    // A warning doesn't stop the case being used
    let warned = Located::Found {
        path: path.clone(),
        warning: Some("different version".into()),
    };
    assert_eq!(resolve(warned, true), Some(path));
    assert_eq!(resolve(Located::Missing("not fetched".into()), false), None);
}

#[test]
#[should_panic(expected = "YBUS_REQUIRE_CASES=1")]
fn resolve_fails_on_missing_cases_when_required() {
    resolve(Located::Missing("not fetched".into()), true);
}

#[test]
#[should_panic(expected = "wrong install")]
fn resolve_always_fails_on_a_broken_install() {
    resolve(Located::Broken("wrong install".into()), false);
}

#[test]
fn pinned_version_reads_the_fetch_script() {
    assert_eq!(
        pinned_version("set -e\nVERSION=\"8.1\"\nSHA256=\"abc\"\n").as_deref(),
        Some("8.1")
    );
    assert_eq!(pinned_version("no version here"), None);
}

#[test]
fn reference_results_match_the_pinned_matpower_version() {
    // The committed reference results must come from the release the tests fetch, or the
    // power flow tests would compare against results from different case data
    let root = common::workspace_root();
    let script = fs::read_to_string(root.join("scripts/fetch-matpower.sh")).unwrap();
    let pinned = pinned_version(&script).expect("the fetch script pins a version");
    let versions =
        fs::read_to_string(root.join("ybus_powerflow/tests/data/matpower-reference/VERSIONS"))
            .unwrap();
    assert!(
        versions
            .lines()
            .any(|line| line.trim() == format!("MATPOWER {pinned}")),
        "fetch-matpower.sh pins MATPOWER {pinned}, but the reference results' VERSIONS file \
         says:\n{versions}\nRegenerate them with scripts/generate-reference.sh"
    );
}

#[test]
fn power_flow_cases_are_available() {
    // Skips locally when MATPOWER isn't fetched; with YBUS_REQUIRE_CASES=1 it checks the
    // fetch worked, so the data-driven tests can't all skip silently
    for case in ["case14", "case30", "case118", "case2383wp"] {
        let Some(path) = common::matpower_case(case) else {
            return;
        };
        let text = fs::read_to_string(&path).unwrap();
        assert!(
            text.contains("mpc.bus") && text.contains("mpc.branch"),
            "{} doesn't look like a MATPOWER case file",
            path.display()
        );
    }
}
