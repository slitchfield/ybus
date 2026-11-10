//! For the tests: MATPOWER's case files, which aren't committed but fetched by
//! `scripts/fetch-matpower.sh`. The rules are `powerflow/tests/common`'s, in miniature: the
//! `MATPOWER_DIR` installation if set, otherwise the pinned download; a missing case is skipped,
//! saying so, unless `YBUS_REQUIRE_CASES=1`, as CI sets.

use std::path::{Path, PathBuf};

/// The path of MATPOWER's `<name>.m`, or `None` (reported on stderr) if it isn't there
pub(crate) fn matpower_case(name: &str) -> Option<PathBuf> {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let dir = match std::env::var_os("MATPOWER_DIR").filter(|dir| !dir.is_empty()) {
        Some(dir) => PathBuf::from(dir),
        None => {
            let script = std::fs::read_to_string(workspace.join("scripts/fetch-matpower.sh"))
                .expect("scripts/fetch-matpower.sh is readable");
            let version = script
                .lines()
                .find_map(|line| line.strip_prefix("VERSION=\""))
                .and_then(|rest| rest.strip_suffix('"'))
                .expect("fetch-matpower.sh has a VERSION=\"...\" line");
            workspace.join(format!("testdata/matpower/matpower{version}"))
        }
    };
    let path = dir.join("data").join(format!("{name}.m"));
    if path.is_file() {
        return Some(path);
    }
    let missing = format!(
        "{} is missing; run scripts/fetch-matpower.sh",
        path.display()
    );
    if std::env::var("YBUS_REQUIRE_CASES").as_deref() == Ok("1") {
        panic!("{missing}");
    }
    eprintln!("skipped: {missing}");
    None
}

/// The cases with committed MATPOWER references
pub(crate) const REFERENCE_CASES: [&str; 3] = ["case14", "case30", "case118"];
