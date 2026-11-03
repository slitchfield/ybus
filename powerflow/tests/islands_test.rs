// Tests for `find_islands` and `find_bridges`: hand-built fixtures, MATPOWER's results for
//   case14, case30 and case118 (committed), and every local case whose results
//   `scripts/generate-reference.sh` has generated. Bridges are also checked by brute force:
//   a branch is a bridge exactly when taking it out adds a component.

mod common;

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use common::dc_fixture::{fixture, B, G, L};
use powerflow::case::Case;
use powerflow::islands::{find_bridges, find_islands, Islands};

/// A line in service
fn line(from: usize, to: usize) -> L {
    L(from, to, 0.1, 0.0, 0.0, true)
}

/// The in-service branches as (case branch row, from bus row, to bus row)
fn edges(case: &Case) -> Vec<(usize, usize, usize)> {
    let row_of: HashMap<usize, usize> = case
        .bus
        .iter()
        .enumerate()
        .map(|(row, bus)| (bus.bus_label, row))
        .collect();
    case.branch
        .iter()
        .enumerate()
        .filter(|(_, branch)| branch.status)
        .map(|(row, branch)| (row, row_of[&branch.from_label], row_of[&branch.to_label]))
        .collect()
}

/// The number of connected components of the graph on `n` buses with `edges`, leaving out
/// the branch at case row `skip`, by union–find: independent of the code under test.
fn components(n: usize, edges: &[(usize, usize, usize)], skip: Option<usize>) -> usize {
    fn root(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }
    let mut parent: Vec<usize> = (0..n).collect();
    let mut count = n;
    for &(row, from, to) in edges {
        if Some(row) == skip {
            continue;
        }
        let (a, b) = (root(&mut parent, from), root(&mut parent, to));
        if a != b {
            parent[a] = b;
            count -= 1;
        }
    }
    count
}

/// The bridges by brute force over `rows` (case branch rows): those in service whose outage
/// adds a component.
fn brute_force_bridges(case: &Case, rows: impl IntoIterator<Item = usize>) -> Vec<usize> {
    let edges = edges(case);
    let n = case.bus.len();
    let base = components(n, &edges, None);
    rows.into_iter()
        .filter(|&row| case.branch[row].status)
        .filter(|&row| components(n, &edges, Some(row)) > base)
        .collect()
}

#[test]
fn islands_are_listed_largest_first() {
    // Buses 1–3 a triangle, 4–5 a pair, 6–7 a pair joined by an out-of-service branch only
    //   to the rest, 8 with no branch, and 9 isolated (type 4) but still joined to 8
    let case = fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 0.0, 0.0, 0.0),
            B(3, 1, 0.0, 0.0, 0.0),
            B(4, 2, 0.0, 0.0, 0.0),
            B(5, 1, 0.0, 0.0, 0.0),
            B(6, 1, 0.0, 0.0, 0.0),
            B(7, 1, 0.0, 0.0, 0.0),
            B(8, 1, 0.0, 0.0, 0.0),
            B(9, 4, 0.0, 0.0, 0.0),
            B(10, 1, 0.0, 0.0, 0.0),
        ],
        &[G(1, 0.0, true), G(4, 0.0, true)],
        &[
            line(1, 2),
            line(2, 3),
            line(3, 1),
            line(6, 7),
            line(4, 5),
            L(5, 6, 0.1, 0.0, 0.0, false),
            line(9, 8),
            // A branch from a bus to itself joins nothing, so bus 10 has no branch
            line(10, 10),
        ],
    );
    assert_eq!(
        find_islands(&case).unwrap(),
        Islands {
            // Equal sizes keep the order of their first bus: 4–5, then 6–7, then 8–9
            groups: vec![vec![0, 1, 2], vec![3, 4], vec![5, 6], vec![7, 8]],
            isolated: vec![9],
        }
    );
}

#[test]
fn bridges_are_the_branches_whose_outage_splits_an_island() {
    // A triangle 1–2–3, a pendant bus 4 off bus 3 (a bridge), bus 5 joined to 4 by two
    //   parallel branches (neither a bridge), and a path 5–6–7 (both bridges)
    let case = fixture(
        &[
            B(1, 3, 0.0, 0.0, 0.0),
            B(2, 1, 0.0, 0.0, 0.0),
            B(3, 1, 0.0, 0.0, 0.0),
            B(4, 1, 0.0, 0.0, 0.0),
            B(5, 1, 0.0, 0.0, 0.0),
            B(6, 1, 0.0, 0.0, 0.0),
            B(7, 1, 0.0, 0.0, 0.0),
        ],
        &[G(1, 0.0, true)],
        &[
            line(1, 2),
            line(2, 3),
            line(3, 1),
            line(3, 4),
            line(4, 5),
            line(5, 4),
            line(5, 6),
            line(6, 7),
            L(7, 1, 0.1, 0.0, 0.0, false),
        ],
    );
    assert_eq!(find_bridges(&case).unwrap(), vec![3, 6, 7]);
    assert_eq!(
        brute_force_bridges(&case, 0..case.branch.len()),
        vec![3, 6, 7]
    );

    // Putting 7–1 in service closes the loop: only the parallel pair's neighbours remain
    let mut closed = case.clone();
    closed.branch[8].status = true;
    assert_eq!(find_bridges(&closed).unwrap(), Vec::<usize>::new());
}

#[test]
fn a_long_path_doesnt_overflow_the_stack() {
    // 20,000 buses in a line: every branch is a bridge, found without recursion
    let n = 20_000;
    let buses: Vec<B> = (1..=n)
        .map(|label| B(label, if label == 1 { 3 } else { 1 }, 0.0, 0.0, 0.0))
        .collect();
    let branches: Vec<L> = (1..n).map(|label| line(label, label + 1)).collect();
    let case = fixture(&buses, &[G(1, 0.0, true)], &branches);
    assert_eq!(find_bridges(&case).unwrap().len(), n - 1);
    assert_eq!(find_islands(&case).unwrap().groups.len(), 1);
}

/// MATPOWER's islands for `name` in `dir`, as `generate_reference.m` writes them (each bus
/// row's island, 1-based in MATPOWER's order, or 0), or `None` without the file.
fn reference_islands(dir: &Path, name: &str) -> Option<Islands> {
    let text = fs::read_to_string(dir.join(format!("{name}_islands.csv"))).ok()?;
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("bus,island"));
    let mut groups: Vec<Vec<usize>> = vec![];
    let mut isolated = vec![];
    for (row, line) in lines.enumerate() {
        let island: usize = line.split(',').nth(1).unwrap().parse().unwrap();
        if island == 0 {
            isolated.push(row);
        } else {
            if groups.len() < island {
                groups.resize(island, vec![]);
            }
            groups[island - 1].push(row);
        }
    }
    Some(Islands { groups, isolated })
}

/// MATPOWER's bridges for `name` in `dir` (0-based rows), or `None` where `find_bridges`
/// failed.
fn reference_bridges(dir: &Path, name: &str) -> Option<Vec<usize>> {
    let text = fs::read_to_string(dir.join(format!("{name}_bridges.csv"))).ok()?;
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("row"));
    Some(
        lines
            .map(|line| line.parse::<usize>().unwrap() - 1)
            .collect(),
    )
}

/// Checks one case against MATPOWER's results in `dir`, and its bridges by brute force: every
/// branch on cases of up to 3,000 buses, and every m-th of about 300 otherwise. Returns
/// whether MATPOWER's bridges were compared.
fn check_case(case: &Case, dir: &Path, name: &str) -> bool {
    let islands = find_islands(case).unwrap();
    if let Some(reference) = reference_islands(dir, name) {
        assert_eq!(islands, reference, "{name}: islands differ from MATPOWER's");
    }
    let bridges = find_bridges(case).unwrap();
    let nl = case.branch.len();
    let step = if case.bus.len() <= 3000 {
        1
    } else {
        nl.div_ceil(300)
    };
    let sample: Vec<usize> = (0..nl).step_by(step).collect();
    let expected = brute_force_bridges(case, sample.iter().copied());
    let found: Vec<usize> = bridges
        .iter()
        .copied()
        .filter(|row| row % step == 0)
        .collect();
    assert_eq!(
        found, expected,
        "{name}: bridges differ from the brute force"
    );
    // Where MATPOWER's find_bridges ran: it fails on many cases (see find_bridges's docs)
    match reference_bridges(dir, name) {
        Some(reference) => {
            assert_eq!(bridges, reference, "{name}: bridges differ from MATPOWER's");
            true
        }
        None => false,
    }
}

#[test]
fn islands_and_bridges_match_matpower_on_the_committed_cases() {
    let dir = common::reference::committed_dir();
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            return; // not fetched; the skip was reported on stderr
        };
        let case = Case::from_path(&path).unwrap();
        assert!(
            reference_islands(&dir, name).is_some(),
            "{name}: no islands file"
        );
        let compared = check_case(&case, &dir, name);
        // case118's parallel branches make MATPOWER's find_bridges fail
        assert_eq!(compared, name != "case118", "{name}");
    }
}

/// The local reference results' directory, beside the normalized cases, or `None` to skip.
fn local_reference_dir() -> Option<PathBuf> {
    let normalized = common::normalized_case("case14")?.parent()?.to_path_buf();
    let name = normalized
        .file_name()?
        .to_str()?
        .strip_suffix("-normalized")?;
    let dir = normalized.with_file_name(format!("{name}-reference"));
    dir.is_dir().then_some(dir)
}

#[test]
fn islands_and_bridges_match_matpower_on_every_local_case() {
    let Some(dir) = local_reference_dir() else {
        eprintln!("skipping: no local reference results (scripts/generate-reference.sh)");
        return;
    };
    let normalized = common::normalized_case("case14").unwrap();
    let normalized = normalized.parent().unwrap();
    let (mut cases, mut compared) = (0, 0);
    let mut paths: Vec<PathBuf> = fs::read_dir(normalized)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "m"))
        .collect();
    paths.sort();
    for path in paths {
        let name = path.file_stem().unwrap().to_str().unwrap().to_string();
        if reference_islands(&dir, &name).is_none() {
            continue;
        }
        let case = Case::from_path(&path).unwrap();
        cases += 1;
        compared += usize::from(check_case(&case, &dir, &name));
    }
    eprintln!("{cases} cases' islands match MATPOWER's; bridges compared on {compared}");
}
