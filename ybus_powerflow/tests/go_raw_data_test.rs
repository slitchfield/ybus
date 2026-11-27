//! The GO Challenge 1 `case.raw` reader on the fetched data (`scripts/fetch-go-c1.sh`): every
//! scenario of both bundles reads, and the totals and sizes match what G2's throwaway reader
//! counted (`docs/GoParserChecklist.md`, "What the data holds"), so each checks the other. It
//! skips when the data isn't fetched, unless `YBUS_REQUIRE_CASES=1`, as in CI. Also tests the
//! data's locator in `common`.

mod common;

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use common::{go_c1_bundles, go_c1_scenarios, locate_go_c1, Located, VERIFIED_MARKER};
use ybus_powerflow::case::go::RawCase;

/// A fresh, empty directory under Cargo's per-target scratch space, unique to `name`.
fn scratch_dir(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("go_raw_data_test")
        .join(name);
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// A fake bundle `name` under `root`, with one scenario, optionally with the verified marker.
fn fake_bundle(root: &Path, name: &str, verified: bool) {
    let scenario = root
        .join(name)
        .join("Part")
        .join("Network_01R-10")
        .join("scenario_1");
    fs::create_dir_all(&scenario).unwrap();
    fs::write(scenario.join("case.raw"), "").unwrap();
    if verified {
        fs::write(root.join(name).join(VERIFIED_MARKER), "").unwrap();
    }
}

fn bundles() -> Vec<String> {
    vec!["one".to_string(), "two".to_string()]
}

#[test]
fn the_bundles_come_from_the_fetch_script() {
    let script = fs::read_to_string(common::workspace_root().join("scripts/fetch-go-c1.sh"))
        .expect("the fetch script is committed");
    assert_eq!(
        go_c1_bundles(&script),
        [
            "Challenge_1_Original_Dataset_1_Scenarios",
            "Challenge_1_Original_Dataset_2_Scenarios"
        ]
    );
}

#[test]
fn the_fetched_data_counts_once_every_bundle_is_verified() {
    let root = scratch_dir("verified");
    fake_bundle(&root, "one", true);
    fake_bundle(&root, "two", true);
    assert_eq!(
        locate_go_c1(None, &root, &bundles()),
        Located::Found {
            path: root.clone(),
            warning: None
        }
    );
    assert_eq!(go_c1_scenarios(&root).len(), 2);

    let root = scratch_dir("one unverified");
    fake_bundle(&root, "one", true);
    fake_bundle(&root, "two", false);
    match locate_go_c1(None, &root, &bundles()) {
        Located::Missing(reason) => {
            assert!(reason.contains("no complete download of two"), "{reason}");
            assert!(reason.contains("scripts/fetch-go-c1.sh"), "{reason}");
        }
        other => panic!("expected Missing, got {other:?}"),
    }
}

#[test]
fn go_c1_dir_takes_precedence_and_must_hold_scenarios() {
    let fetched = scratch_dir("unfetched");
    let copy = scratch_dir("copy");
    fake_bundle(&copy, "anything", false);
    assert_eq!(
        locate_go_c1(Some(&copy), &fetched, &bundles()),
        Located::Found {
            path: copy.clone(),
            warning: None
        }
    );
    let empty = scratch_dir("empty");
    assert!(matches!(
        locate_go_c1(Some(&empty), &fetched, &bundles()),
        Located::Broken(reason) if reason.contains("GO_C1_DIR=")
    ));
}

/// The smallest and largest of each of the seven counts below
type Ranges = [(usize, usize); 7];

/// The sizes G2 recorded for each network, over all its scenarios: how many scenarios, then the
/// smallest and largest count of buses, loads, fixed shunts, generators, lines, transformers and
/// switched shunts. A network's scenarios aren't all the same network.
const SIZES: [(&str, usize, Ranges); 5] = [
    (
        "Network_01",
        60,
        [
            (500, 500),
            (200, 200),
            (0, 0),
            (90, 90),
            (468, 468),
            (131, 131),
            (17, 19),
        ],
    ),
    (
        "Network_03",
        60,
        [
            (793, 793),
            (568, 568),
            (49, 49),
            (210, 210),
            (769, 769),
            (143, 143),
            (50, 50),
        ],
    ),
    (
        "Network_05",
        40,
        [
            (2000, 2000),
            (1125, 1125),
            (0, 3),
            (544, 544),
            (2359, 2359),
            (847, 847),
            (153, 160),
        ],
    ),
    (
        "Network_07",
        60,
        [
            (2312, 2312),
            (1529, 1529),
            (121, 121),
            (617, 617),
            (2156, 2156),
            (857, 857),
            (201, 201),
        ],
    ),
    (
        "Network_10",
        20,
        [
            (7977, 7977),
            (6671, 6671),
            (0, 4587),
            (560, 570),
            (8896, 8900),
            (2548, 2693),
            (6, 5082),
        ],
    ),
];

#[test]
fn every_fetched_case_raw_reads_and_matches_g2s_counts() {
    let Some(dir) = common::go_c1_dir() else {
        return; // not fetched; the skip was reported on stderr
    };
    let scenarios = go_c1_scenarios(&dir);

    // Totals over every scenario: buses, generators, lines, transformers, nonzero MAG2, CZ 2,
    //   and generators out of service
    let mut totals = [0usize; 7];
    // Each network's scenarios, and the range of each of its counts
    let mut sizes: BTreeMap<String, (usize, Ranges)> = BTreeMap::new();
    for scenario in &scenarios {
        let path = scenario.join("case.raw");
        let case = RawCase::from_path(&path)
            .unwrap_or_else(|err| panic!("{} doesn't read: {err}", path.display()));
        let transformers = &case.transformers;
        for (total, count) in totals.iter_mut().zip([
            case.buses.len(),
            case.generators.len(),
            case.lines.len(),
            transformers.len(),
            transformers.iter().filter(|t| t.mag2_pu != 0.0).count(),
            transformers.iter().filter(|t| t.cz == 2).count(),
            case.generators.iter().filter(|g| !g.status).count(),
        ]) {
            *total += count;
        }

        let network = scenario
            .parent()
            .and_then(|network| network.file_name())
            .and_then(|name| name.to_str())
            .expect("a scenario is in a network folder");
        let counts = [
            case.buses.len(),
            case.loads.len(),
            case.fixed_shunts.len(),
            case.generators.len(),
            case.lines.len(),
            case.transformers.len(),
            case.switched_shunts.len(),
        ];
        let (seen, ranges) = sizes
            .entry(network[..network.len().min(10)].to_string())
            .or_insert((0, counts.map(|count| (count, count))));
        *seen += 1;
        for ((low, high), count) in ranges.iter_mut().zip(counts) {
            *low = (*low).min(count);
            *high = (*high).max(count);
        }
    }
    assert_eq!(scenarios.len(), 240, "scenarios under {}", dir.display());
    assert_eq!(
        totals,
        [455_840, 88_160, 475_878, 153_006, 15_180, 5_700, 35_720],
        "buses, generators, lines, transformers, nonzero MAG2, CZ 2, generators out of service"
    );
    let want: BTreeMap<String, (usize, Ranges)> = SIZES
        .into_iter()
        .map(|(network, count, ranges)| (network.to_string(), (count, ranges)))
        .collect();
    assert_eq!(sizes, want);
}
