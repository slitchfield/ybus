//! What a power flow solves over, shared by the DC and AC models: the connected buses, the
//! in-service generators and branches between them, each bus's role, and the check that every
//! connected bus reaches a reference bus. It follows MATPOWER's `ext2int` and `bustypes`.

use std::collections::{HashMap, VecDeque};

use crate::case::{Branch, BusType, Case};
use crate::powerflow_error::PowerflowError;

/// A case's topology for a power flow, in internal numbering: internal bus `i` is case bus row
/// `bus_rows[i]`, and internal branch `k` is case branch row `branch_rows[k]`; both keep the
/// case's order.
#[derive(Debug, Clone)]
pub(crate) struct Topology {
    pub(crate) bus_rows: Vec<usize>,
    pub(crate) branch_rows: Vec<usize>,
    /// The internal (from, to) buses of each internal branch, parallel to `branch_rows`
    pub(crate) ends: Vec<(usize, usize)>,
    pub(crate) gen_rows: Vec<usize>,
    /// The internal bus of each in-service generator, parallel to `gen_rows`
    pub(crate) gen_buses: Vec<usize>,
    pub(crate) ref_buses: Vec<usize>,
    pub(crate) pv_buses: Vec<usize>,
    pub(crate) pq_buses: Vec<usize>,
    /// The table lengths of the case the topology was built from: buses, generators, branches
    pub(crate) case_shape: (usize, usize, usize),
}

impl Topology {
    /// The topology of `case`, after validating it with [`Case::validate`].
    ///
    /// - **`ext2int`:** buses that aren't isolated (type 4) are connected, in case order.
    ///   Generators and branches count only when in service and attached to connected buses.
    ///   `check_branch` is called on each such branch, with its 1-based case row, as it's found,
    ///   so a model can reject a branch it can't represent before any later check runs.
    /// - **`bustypes`:** a REF or PV bus counts as one only if it has an in-service generator,
    ///   and is otherwise treated as PQ. With no REF bus, the first PV bus becomes the
    ///   reference, and with neither, it's `PowerflowError::NoReferenceBus`.
    /// - **Islands:** every connected bus must reach a reference bus through in-service
    ///   branches, or it's `PowerflowError::IslandWithoutReference`, naming a stranded bus.
    ///   That's checked here rather than inferred from a failed factorization, which can't name
    ///   the island and may not fail at all, since rounding can leave a tiny positive pivot.
    pub(crate) fn new(
        case: &Case,
        mut check_branch: impl FnMut(usize, &Branch) -> Result<(), PowerflowError>,
    ) -> Result<Topology, PowerflowError> {
        // Only Case::parse validates, so a case edited since needs checking here
        case.validate()?;

        // ext2int: connected buses, in case order, and the label → internal index map
        let bus_rows: Vec<usize> = (0..case.bus.len())
            .filter(|&row| case.bus[row].bus_type != BusType::Isolated)
            .collect();
        let mut internal = HashMap::with_capacity(bus_rows.len());
        for (index, &row) in bus_rows.iter().enumerate() {
            internal.insert(case.bus[row].bus_label, index);
        }
        let nb = bus_rows.len();

        // Generators and branches count only when in service and fully connected
        let (mut gen_rows, mut gen_buses) = (vec![], vec![]);
        for (row, generator) in case.generator.iter().enumerate() {
            if let (true, Some(&bus)) = (generator.status, internal.get(&generator.bus_label)) {
                gen_rows.push(row);
                gen_buses.push(bus);
            }
        }
        let mut branch_rows = vec![];
        let mut ends = vec![];
        for (row, branch) in case.branch.iter().enumerate() {
            let from = internal.get(&branch.from_label);
            let to = internal.get(&branch.to_label);
            if let (true, Some(&from), Some(&to)) = (branch.status, from, to) {
                check_branch(row + 1, branch)?;
                branch_rows.push(row);
                ends.push((from, to));
            }
        }

        // bustypes
        let mut has_gen = vec![false; nb];
        for &bus in &gen_buses {
            has_gen[bus] = true;
        }
        let (mut ref_buses, mut pv_buses, mut pq_buses) = (vec![], vec![], vec![]);
        for (index, &row) in bus_rows.iter().enumerate() {
            match case.bus[row].bus_type {
                BusType::Ref if has_gen[index] => ref_buses.push(index),
                BusType::PV if has_gen[index] => pv_buses.push(index),
                _ => pq_buses.push(index),
            }
        }
        if ref_buses.is_empty() {
            if pv_buses.is_empty() {
                return Err(PowerflowError::NoReferenceBus);
            }
            ref_buses.push(pv_buses.remove(0));
        }

        // Every connected bus must reach a reference bus
        let mut neighbours = vec![vec![]; nb];
        for &(from, to) in &ends {
            neighbours[from].push(to);
            neighbours[to].push(from);
        }
        let mut reached = vec![false; nb];
        let mut queue: VecDeque<usize> = ref_buses.iter().copied().collect();
        for &bus in &ref_buses {
            reached[bus] = true;
        }
        while let Some(bus) = queue.pop_front() {
            for &next in &neighbours[bus] {
                if !reached[next] {
                    reached[next] = true;
                    queue.push_back(next);
                }
            }
        }
        if let Some(stranded) = reached.iter().position(|&reached| !reached) {
            return Err(PowerflowError::IslandWithoutReference {
                bus_label: case.bus[bus_rows[stranded]].bus_label,
            });
        }

        Ok(Topology {
            bus_rows,
            branch_rows,
            ends,
            gen_rows,
            gen_buses,
            ref_buses,
            pv_buses,
            pq_buses,
            case_shape: (case.bus.len(), case.generator.len(), case.branch.len()),
        })
    }

    /// The PV then PQ buses, the unknowns of a power flow, in MATPOWER's order.
    pub(crate) fn pv_pq(&self) -> Vec<usize> {
        self.pv_buses
            .iter()
            .chain(&self.pq_buses)
            .copied()
            .collect()
    }

    /// `PowerflowError::CaseMismatch` unless `case` has the table lengths of the case this
    /// topology came from.
    pub(crate) fn check_case(&self, case: &Case) -> Result<(), PowerflowError> {
        let found = (case.bus.len(), case.generator.len(), case.branch.len());
        if found == self.case_shape {
            Ok(())
        } else {
            Err(PowerflowError::CaseMismatch {
                expected: self.case_shape,
                found,
            })
        }
    }
}
