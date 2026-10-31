//! Shared test support: locating MATPOWER case files.
//!
//! MATPOWER's case files are not committed (see `docs/PowerFlowChecklist.md`, Phase 1.5).
//! `scripts/fetch-matpower.sh` downloads them, so tests that need a case look it up here:
//!
//! ```ignore
//! mod common;
//!
//! #[test]
//! fn reads_case14() {
//!     let Some(path) = common::matpower_case("case14") else {
//!         return; // not fetched; the skip was reported on stderr
//!     };
//!     // ... use path ...
//! }
//! ```
//!
//! Lookup order:
//! 1. `MATPOWER_DIR`, if set, names an existing MATPOWER installation, as for the fetch script.
//!    A case missing from it is an error, not a skip, since the directory was chosen on purpose.
//!    It needs no `.ybus-verified` marker, since the fetch script didn't produce it. Its version
//!    (from `lib/mpver.m`) is compared with the pinned one, and a mismatch is reported as a
//!    warning, not an error: the committed reference results come from the pinned release, so
//!    comparisons may fail, but keeping a supplied installation compatible is the user's
//!    responsibility.
//! 2. Otherwise `testdata/matpower/matpower<VERSION>/` from the fetch script, at the version the
//!    script pins. It only counts once the script's `.ybus-verified` marker shows a complete,
//!    checksum-verified extraction.
//!
//! When the cases aren't there, tests skip by returning early. Set `YBUS_REQUIRE_CASES=1` to
//! make a missing case fail instead. CI does this, so a broken fetch can't turn into silently
//! passing tests.
//!
//! [`normalized_case`] finds the same cases as rewritten by `scripts/normalize-cases.sh`, which
//! evaluates the ones that compute their data with code (with MATPOWER's `loadcase`, under
//! Octave) and writes every case as plain literals. It follows the same rules, and its output
//! only counts once complete (`.ybus-normalized`) and made from the MATPOWER installation
//! currently in use (its `SOURCE` file).

// Each integration test file compiles this module separately and uses only part of it
#![allow(dead_code)]

/// The backward error of a DC power flow's angles.
pub mod accuracy;
/// Hand-built cases for the DC power flow tests.
pub mod dc_fixture;
/// MATPOWER's DC power flow results, to compare against.
pub mod reference;

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

/// Names an existing MATPOWER installation to use instead of `testdata/`.
pub const MATPOWER_DIR_VAR: &str = "MATPOWER_DIR";
/// When `1`, a missing case fails the test instead of skipping it.
pub const REQUIRE_VAR: &str = "YBUS_REQUIRE_CASES";
/// Written by `scripts/fetch-matpower.sh` after a complete, verified extraction.
pub const VERIFIED_MARKER: &str = ".ybus-verified";

/// Written by `scripts/normalize-cases.sh` once its output is complete.
pub const NORMALIZED_MARKER: &str = ".ybus-normalized";

/// The outcome of looking for a case file.
#[derive(Debug, PartialEq)]
pub enum Located {
    /// The case file exists at this path. `warning` says why results from it may not match the
    /// committed reference results, such as a different MATPOWER version.
    Found {
        path: PathBuf,
        warning: Option<String>,
    },
    /// MATPOWER hasn't been fetched (or not completely). Tests skip, unless cases are required.
    Missing(String),
    /// The installation exists but is wrong. Tests always fail, so a misconfiguration is fixed
    /// rather than skipped.
    Broken(String),
}

/// The workspace root, the parent of this crate's directory.
pub fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("the powerflow crate is inside the workspace")
        .to_path_buf()
}

/// The MATPOWER version pinned by `scripts/fetch-matpower.sh`, read from its `VERSION="…"`
/// line so the version is written down in one place.
pub fn pinned_version(script: &str) -> Option<String> {
    script.lines().find_map(|line| {
        line.trim()
            .strip_prefix("VERSION=\"")?
            .strip_suffix('"')
            .map(str::to_string)
    })
}

/// The MATPOWER version `scripts/fetch-matpower.sh` pins, which the reference results were
/// generated with.
pub fn pinned_matpower_version(root: &Path) -> String {
    let script = root.join("scripts").join("fetch-matpower.sh");
    let text = fs::read_to_string(&script)
        .unwrap_or_else(|err| panic!("can't read {}: {err}", script.display()));
    pinned_version(&text)
        .unwrap_or_else(|| panic!("no VERSION=\"...\" line in {}", script.display()))
}

/// The directory the fetch script extracts release `version` into.
pub fn fetched_dir(root: &Path, version: &str) -> PathBuf {
    root.join("testdata")
        .join("matpower")
        .join(format!("matpower{version}"))
}

/// The version recorded in a MATPOWER installation's `lib/mpver.m`, which holds it as
/// `'Version',  '8.1', ...` inside a struct literal. Comment lines are skipped.
pub fn installed_version(mpver: &str) -> Option<String> {
    mpver
        .lines()
        .filter(|line| !line.trim_start().starts_with('%'))
        .find_map(|line| {
            let after_key = line.split_once("'Version'")?.1;
            let value = after_key.split_once('\'')?.1;
            value
                .split_once('\'')
                .map(|(version, _)| version.to_string())
        })
}

/// A warning if the MATPOWER installation at `dir` is not version `pinned`, or if its version
/// can't be determined.
fn version_warning(dir: &Path, pinned: &str) -> Option<String> {
    let mpver = dir.join("lib").join("mpver.m");
    let advice = format!(
        "the reference results were generated with MATPOWER {pinned}, so tests comparing against \
         them may fail. Keeping a {MATPOWER_DIR_VAR} installation compatible is your \
         responsibility; unset {MATPOWER_DIR_VAR} to use the pinned download from \
         scripts/fetch-matpower.sh"
    );
    match fs::read_to_string(&mpver)
        .ok()
        .as_deref()
        .and_then(installed_version)
    {
        Some(version) if version == pinned => None,
        Some(version) => Some(format!(
            "{MATPOWER_DIR_VAR}={} is MATPOWER {version}, but {advice}",
            dir.display()
        )),
        None => Some(format!(
            "can't determine the MATPOWER version of {MATPOWER_DIR_VAR}={} (no readable version \
             in {}); {advice}",
            dir.display(),
            mpver.display()
        )),
    }
}

/// Looks for `data/<case>.m`, in `env_dir` if given (the value of `MATPOWER_DIR`) and
/// otherwise in `fetched_dir`, the download of release `pinned`. A case found in `env_dir`
/// carries a warning if that installation isn't release `pinned`. Reads no environment
/// variables, so it can be tested directly.
pub fn locate(case: &str, env_dir: Option<&Path>, fetched_dir: &Path, pinned: &str) -> Located {
    let file = format!("{case}.m");
    if let Some(dir) = env_dir {
        let path = dir.join("data").join(&file);
        return if path.is_file() {
            Located::Found {
                path,
                warning: version_warning(dir, pinned),
            }
        } else {
            Located::Broken(format!(
                "{MATPOWER_DIR_VAR}={} has no data/{file}; is it a MATPOWER directory?",
                dir.display()
            ))
        };
    }

    if !fetched_dir.join(VERIFIED_MARKER).is_file() {
        return Located::Missing(format!(
            "MATPOWER is not fetched (no complete download at {}); run scripts/fetch-matpower.sh \
             or set {MATPOWER_DIR_VAR}",
            fetched_dir.display()
        ));
    }
    let path = fetched_dir.join("data").join(&file);
    if path.is_file() {
        // The fetch script downloaded exactly the pinned release, so no version check
        Located::Found {
            path,
            warning: None,
        }
    } else {
        Located::Broken(format!(
            "the verified MATPOWER download at {} has no data/{file}",
            fetched_dir.display()
        ))
    }
}

/// Parses `YBUS_REQUIRE_CASES`: unset, empty or `0` means skip, `1` means require. Anything
/// else is rejected, so a typo such as `yes` can't silently mean "skip".
pub fn parse_require(value: Option<&str>) -> bool {
    match value {
        None | Some("") | Some("0") => false,
        Some("1") => true,
        Some(other) => panic!("{REQUIRE_VAR} must be 0 or 1, not {other:?}"),
    }
}

/// Turns a lookup into the path to use, `None` to skip, or a panic. A warning goes to stderr,
/// which the test harness shows only for failing tests: where a version mismatch matters.
pub fn resolve(located: Located, require: bool) -> Option<PathBuf> {
    match located {
        Located::Found { path, warning } => {
            if let Some(warning) = warning {
                eprintln!("warning: {warning}");
            }
            Some(path)
        }
        Located::Missing(reason) if require => {
            panic!("{reason} ({REQUIRE_VAR}=1, so missing cases are an error)")
        }
        Located::Missing(reason) => {
            eprintln!("skipping: {reason}");
            None
        }
        Located::Broken(reason) => panic!("{reason}"),
    }
}

/// The directory `scripts/normalize-cases.sh` writes the cases normalized from the MATPOWER
/// installation at `source` into: `testdata/matpower/<source's name>-normalized/`.
pub fn normalized_dir(root: &Path, source: &Path) -> PathBuf {
    let name = source
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| "matpower".to_string());
    root.join("testdata")
        .join("matpower")
        .join(format!("{name}-normalized"))
}

/// Looks for the normalized `<case>.m` in `normalized`, which must be complete and made from
/// the MATPOWER installation at `source`. Reads no environment variables.
pub fn locate_normalized(case: &str, source: &Path, normalized: &Path) -> Located {
    if !normalized.join(NORMALIZED_MARKER).is_file() {
        return Located::Missing(format!(
            "MATPOWER cases aren't normalized (no complete output at {}); run \
             scripts/normalize-cases.sh, which needs Octave",
            normalized.display()
        ));
    }

    // The script records the physical path of the MATPOWER directory it read
    let expected = fs::canonicalize(source).unwrap_or_else(|_| source.to_path_buf());
    let recorded = fs::read_to_string(normalized.join("SOURCE"))
        .ok()
        .and_then(|text| text.lines().next().map(str::to_string));
    if recorded.as_deref() != Some(format!("source {}", expected.display()).as_str()) {
        return Located::Missing(format!(
            "the normalized cases at {} weren't made from {}; rerun scripts/normalize-cases.sh",
            normalized.display(),
            expected.display()
        ));
    }

    let path = normalized.join(format!("{case}.m"));
    if path.is_file() {
        Located::Found {
            path,
            warning: None,
        }
    } else {
        Located::Broken(format!(
            "the normalized cases at {} have no {case}.m; is it a MATPOWER case that loadcase \
             can load?",
            normalized.display()
        ))
    }
}

/// The path of the normalized `<case>.m`, written as literals by `scripts/normalize-cases.sh`,
/// or `None` (after reporting why on stderr) when the cases haven't been normalized and cases
/// aren't required. Panics if cases are required but missing. Normalized from `MATPOWER_DIR`
/// when that's set, with the same version warning as [`matpower_case`].
pub fn normalized_case(case: &str) -> Option<PathBuf> {
    let env_dir = env::var_os(MATPOWER_DIR_VAR)
        .filter(|dir| !dir.is_empty())
        .map(PathBuf::from);
    let require = parse_require(env::var(REQUIRE_VAR).ok().as_deref());
    let root = workspace_root();
    let pinned = pinned_matpower_version(&root);
    let source = env_dir
        .clone()
        .unwrap_or_else(|| fetched_dir(&root, &pinned));
    let located = match locate_normalized(case, &source, &normalized_dir(&root, &source)) {
        Located::Found { path, .. } => Located::Found {
            path,
            warning: env_dir.and_then(|dir| version_warning(&dir, &pinned)),
        },
        other => other,
    };
    resolve(located, require)
}

/// The path of MATPOWER's `data/<case>.m`, or `None` (after reporting why on stderr) when
/// MATPOWER hasn't been fetched and cases aren't required. Panics if cases are required but
/// missing, or if the installation is broken. A `MATPOWER_DIR` installation of a different
/// version than the pinned one is used, with a warning on stderr.
pub fn matpower_case(case: &str) -> Option<PathBuf> {
    let env_dir = env::var_os(MATPOWER_DIR_VAR)
        .filter(|dir| !dir.is_empty())
        .map(PathBuf::from);
    let require = parse_require(env::var(REQUIRE_VAR).ok().as_deref());
    let root = workspace_root();
    let pinned = pinned_matpower_version(&root);
    resolve(
        locate(
            case,
            env_dir.as_deref(),
            &fetched_dir(&root, &pinned),
            &pinned,
        ),
        require,
    )
}
