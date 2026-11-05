//! Whether taking out branches strands a bus from every reference bus, without searching the
//! network on every outage.
//!
//! An outage can only strand buses if it cuts the network, and a single branch cuts it only if
//! it's a bridge (`islands::bridges`, which counts parallel branches separately, as the search
//! does). The bridges join the network's 2-edge-connected components into a forest; each tree
//! holds a reference bus, since the model checked that every bus reaches one. Rooting each tree
//! at a component with a reference bus, a bridge's removal strands exactly the subtree below it
//! when that subtree holds no reference bus, and the first stranded bus, in internal order, is
//! the subtree's smallest. So one pass over the forest answers every single-branch outage, as
//! the search would, error included. Outages of several branches search breadth-first from the
//! reference buses, over an adjacency built once.

use std::collections::VecDeque;

use crate::islands::bridges;
use crate::topology::Topology;

/// The answers for single-branch outages, and the adjacency for the rest.
#[derive(Debug, Clone)]
pub(crate) struct Reachability {
    /// For each internal branch, the smallest internal bus its outage strands, if any
    single: Vec<Option<usize>>,
    /// Each bus's (neighbour, internal branch), in CSR form: bus `i`'s are
    /// `adjacency[start[i]..start[i + 1]]`
    start: Vec<usize>,
    adjacency: Vec<(usize, usize)>,
}

/// Scratch for the breadth-first search, reused across outages.
#[derive(Debug, Default)]
pub(crate) struct ReachScratch {
    reached: Vec<bool>,
    queue: VecDeque<usize>,
}

impl Reachability {
    /// The answers for `topology`'s network, every bus of which reaches a reference bus.
    pub(crate) fn new(topology: &Topology) -> Reachability {
        let nb = topology.bus_rows.len();
        let ends = &topology.ends;

        // Each bus's neighbours, with the branch, in CSR form
        let mut degree = vec![0; nb + 1];
        for &(from, to) in ends {
            degree[from + 1] += 1;
            degree[to + 1] += 1;
        }
        let mut start = degree;
        for i in 0..nb {
            start[i + 1] += start[i];
        }
        let mut next = start.clone();
        let mut adjacency = vec![(0, 0); start[nb]];
        for (k, &(from, to)) in ends.iter().enumerate() {
            adjacency[next[from]] = (to, k);
            next[from] += 1;
            adjacency[next[to]] = (from, k);
            next[to] += 1;
        }

        // The 2-edge-connected components: buses joined by branches that aren't bridges
        let mut is_bridge = vec![false; ends.len()];
        for k in bridges(nb, ends) {
            is_bridge[k] = true;
        }
        let mut component = vec![usize::MAX; nb];
        let mut count = 0;
        let mut queue = VecDeque::new();
        for root in 0..nb {
            if component[root] != usize::MAX {
                continue;
            }
            component[root] = count;
            queue.push_back(root);
            while let Some(bus) = queue.pop_front() {
                for &(next, k) in &adjacency[start[bus]..start[bus + 1]] {
                    if !is_bridge[k] && component[next] == usize::MAX {
                        component[next] = count;
                        queue.push_back(next);
                    }
                }
            }
            count += 1;
        }

        // Each component's smallest bus and whether it holds a reference bus
        let mut smallest = vec![usize::MAX; count];
        let mut has_ref = vec![false; count];
        for (bus, &c) in component.iter().enumerate() {
            smallest[c] = smallest[c].min(bus);
        }
        for &bus in &topology.ref_buses {
            has_ref[component[bus]] = true;
        }

        // The forest of components joined by bridges, each tree searched from a component with
        //   a reference bus, recording each component's parent bridge in discovery order
        let mut tree: Vec<Vec<(usize, usize)>> = vec![vec![]; count];
        for (k, &(from, to)) in ends.iter().enumerate() {
            if is_bridge[k] {
                let (a, b) = (component[from], component[to]);
                tree[a].push((b, k));
                tree[b].push((a, k));
            }
        }
        let mut parent = vec![usize::MAX; count];
        let mut parent_bridge = vec![usize::MAX; count];
        let mut visited = vec![false; count];
        let mut order = Vec::with_capacity(count);
        let roots = topology.ref_buses.iter().map(|&bus| component[bus]);
        for root in roots {
            if visited[root] {
                continue;
            }
            visited[root] = true;
            let mut stack = vec![root];
            while let Some(c) = stack.pop() {
                order.push(c);
                for &(next, k) in &tree[c] {
                    if !visited[next] {
                        visited[next] = true;
                        parent[next] = c;
                        parent_bridge[next] = k;
                        stack.push(next);
                    }
                }
            }
        }

        // Each subtree's smallest bus and whether it holds a reference bus, children first
        let mut subtree_smallest = smallest;
        let mut subtree_has_ref = has_ref;
        for &c in order.iter().rev() {
            let p = parent[c];
            if p != usize::MAX {
                subtree_smallest[p] = subtree_smallest[p].min(subtree_smallest[c]);
                subtree_has_ref[p] = subtree_has_ref[p] || subtree_has_ref[c];
            }
        }

        // A bridge strands the subtree below it when that subtree has no reference bus
        let mut single = vec![None; ends.len()];
        for c in 0..count {
            let k = parent_bridge[c];
            if k != usize::MAX && !subtree_has_ref[c] {
                single[k] = Some(subtree_smallest[c]);
            }
        }
        Reachability {
            single,
            start,
            adjacency,
        }
    }

    /// The smallest internal bus left unable to reach a reference bus without the internal
    /// branches `outaged` (sorted, distinct), if any: what a new model of the case with those
    /// branches out of service would report.
    pub(crate) fn stranded(
        &self,
        topology: &Topology,
        outaged: &[usize],
        scratch: &mut ReachScratch,
    ) -> Option<usize> {
        match outaged {
            [] => None,
            [k] => self.single[*k],
            _ => {
                let nb = topology.bus_rows.len();
                let ReachScratch { reached, queue } = scratch;
                reached.clear();
                reached.resize(nb, false);
                queue.clear();
                for &bus in &topology.ref_buses {
                    reached[bus] = true;
                    queue.push_back(bus);
                }
                while let Some(bus) = queue.pop_front() {
                    for &(next, k) in &self.adjacency[self.start[bus]..self.start[bus + 1]] {
                        if !reached[next] && outaged.binary_search(&k).is_err() {
                            reached[next] = true;
                            queue.push_back(next);
                        }
                    }
                }
                reached.iter().position(|&reached| !reached)
            }
        }
    }
}
