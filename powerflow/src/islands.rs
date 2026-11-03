//! A network's islands and bridges, as MATPOWER's `find_islands` and `find_bridges` define them.
//!
//! Both work on the graph whose nodes are every bus in the case, of any type, and whose edges
//! are the in-service branches. A branch from a bus to itself joins nothing, so it is left out.
//! [`find_islands`] lists the connected components, and [`find_bridges`] the branches whose
//! outage would split one: those an N-1 study can't take out without islanding part of the
//! network.

use std::collections::{HashMap, VecDeque};

use crate::case::Case;
use crate::powerflow_error::PowerflowError;

/// A network's islands, from [`find_islands`]: groups of case bus rows connected by in-service
/// branches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Islands {
    /// Each island's case bus rows, in increasing order. Islands are listed largest first, and
    /// islands of the same size in the order of their first bus, as MATPOWER lists them.
    pub groups: Vec<Vec<usize>>,
    /// The case bus rows with no in-service branch, in increasing order. They belong to no
    /// group.
    pub isolated: Vec<usize>,
}

/// The islands of `case`, as MATPOWER's `find_islands` (with `connected_components`) finds
/// them: every bus with an in-service branch belongs to exactly one group, whatever its type,
/// so an isolated bus (type 4) that still has an in-service branch is grouped too.
///
/// Validates `case` first with [`Case::validate`], returning its errors.
pub fn find_islands(case: &Case) -> Result<Islands, PowerflowError> {
    let ends: Vec<(usize, usize)> = branch_edges(case)?
        .iter()
        .map(|&(_, from, to)| (from, to))
        .collect();
    let neighbours = adjacency(case.bus.len(), &ends);
    let mut island = vec![None; case.bus.len()];
    let mut groups: Vec<Vec<usize>> = vec![];
    let mut isolated = vec![];
    let mut queue = VecDeque::new();
    for start in 0..case.bus.len() {
        if neighbours[start].is_empty() {
            isolated.push(start);
            continue;
        }
        if island[start].is_some() {
            continue;
        }
        let id = groups.len();
        island[start] = Some(id);
        let mut group = vec![start];
        queue.push_back(start);
        while let Some(bus) = queue.pop_front() {
            for &(next, _) in &neighbours[bus] {
                if island[next].is_none() {
                    island[next] = Some(id);
                    group.push(next);
                    queue.push_back(next);
                }
            }
        }
        group.sort_unstable();
        groups.push(group);
    }
    // Largest first; a stable sort keeps islands of equal size in discovery order
    groups.sort_by_key(|group| std::cmp::Reverse(group.len()));
    Ok(Islands { groups, isolated })
}

/// The bridges of `case`: the in-service branches whose outage would split an island, as case
/// branch rows in increasing order. Found by Tarjan's algorithm (Tarjan 1974) over the graph of
/// [`find_islands`].
///
/// Changed from MATPOWER's `find_bridges`, which identifies each branch by its pair of buses:
/// here parallel branches are distinct edges, so neither of a pair of parallel branches is a
/// bridge, since the other still joins its buses. MATPOWER's, on a case with parallel branches,
/// adds their row numbers together and so marks the wrong branch or fails, and it fails on any
/// case with more buses than branches; on cases it handles, the results agree.
///
/// Validates `case` first with [`Case::validate`], returning its errors.
pub fn find_bridges(case: &Case) -> Result<Vec<usize>, PowerflowError> {
    let edges = branch_edges(case)?;
    let rows: Vec<usize> = edges.iter().map(|&(row, _, _)| row).collect();
    let ends: Vec<(usize, usize)> = edges.iter().map(|&(_, f, t)| (f, t)).collect();
    let mut bridges: Vec<usize> = bridges(case.bus.len(), &ends)
        .into_iter()
        .map(|edge| rows[edge])
        .collect();
    bridges.sort_unstable();
    Ok(bridges)
}

/// The in-service branches of `case` as (case branch row, from bus row, to bus row), leaving
/// out a branch from a bus to itself, after validating `case`.
fn branch_edges(case: &Case) -> Result<Vec<(usize, usize, usize)>, PowerflowError> {
    case.validate()?;
    let row_of: HashMap<usize, usize> = case
        .bus
        .iter()
        .enumerate()
        .map(|(row, bus)| (bus.bus_label, row))
        .collect();
    Ok(case
        .branch
        .iter()
        .enumerate()
        .filter(|(_, branch)| branch.status)
        .map(|(row, branch)| (row, row_of[&branch.from_label], row_of[&branch.to_label]))
        .filter(|&(_, from, to)| from != to)
        .collect())
}

/// Each node's (neighbour, edge index) pairs, for edges given as (from, to) pairs. A loop joins
/// nothing, so it is left out.
fn adjacency(n: usize, ends: &[(usize, usize)]) -> Vec<Vec<(usize, usize)>> {
    let mut neighbours = vec![vec![]; n];
    for (edge, &(from, to)) in ends.iter().enumerate() {
        if from != to {
            neighbours[from].push((to, edge));
            neighbours[to].push((from, edge));
        }
    }
    neighbours
}

/// The bridges of the graph on nodes `0..n` with the edges `ends`, as indices into `ends`, in
/// the order found. Tarjan's algorithm, iteratively, so a long path can't overflow the stack:
/// a depth-first search numbers the nodes in discovery order, and an edge to a child is a
/// bridge when nothing in the child's subtree has an edge back to the node or above it. Only
/// the edge to the parent is skipped, not every edge to the parent node, so of two parallel
/// edges each is the other's way back. A loop (`from == to`) is never a bridge.
pub(crate) fn bridges(n: usize, ends: &[(usize, usize)]) -> Vec<usize> {
    let neighbours = adjacency(n, ends);
    let mut discovered: Vec<Option<usize>> = vec![None; n];
    let mut low = vec![0; n];
    let mut time = 0;
    let mut found = vec![];
    // The search's stack of (node, the edge it was reached by, next neighbour slot)
    let mut stack: Vec<(usize, Option<usize>, usize)> = vec![];
    for root in 0..n {
        if discovered[root].is_some() {
            continue;
        }
        discovered[root] = Some(time);
        low[root] = time;
        time += 1;
        stack.push((root, None, 0));
        while let Some(&mut (node, via, ref mut next)) = stack.last_mut() {
            if let Some(&(child, edge)) = neighbours[node].get(*next) {
                *next += 1;
                if via == Some(edge) {
                    continue;
                }
                match discovered[child] {
                    Some(order) => low[node] = low[node].min(order),
                    None => {
                        discovered[child] = Some(time);
                        low[child] = time;
                        time += 1;
                        stack.push((child, Some(edge), 0));
                    }
                }
            } else {
                stack.pop();
                if let (Some(&(parent, _, _)), Some(edge)) = (stack.last(), via) {
                    low[parent] = low[parent].min(low[node]);
                    let parent_order = discovered[parent].expect("a parent is discovered");
                    if low[node] > parent_order {
                        found.push(edge);
                    }
                }
            }
        }
    }
    found
}
