//! Sensitivities of the DC power flow: power transfer distribution factors (PTDF) and line outage
//! distribution factors (LODF), as MATPOWER's `makePTDF` and `makeLODF` compute them.
//!
//! A PTDF is the change in a branch's from-end flow per unit of power injected at a bus and
//! withdrawn at the slack, so `PTDF[k][j] = b_k (θ_f − θ_t)` where `θ` solves `B θ = e_j` with the
//! slack's row and column removed and its angle fixed at zero. Only angle differences matter,
//! so the angle reference doesn't. MATPOWER returns the whole dense matrix (branches × buses),
//! which large cases can't hold, so this module also gives selected rows and columns, each from
//! one solve with a sparse right-hand side.
//!
//! An LODF is the change in each branch's flow, as a fraction of a branch's flow before its
//! outage, when that branch is taken out: `LODF[l][k] = H[l][k] / (1 − H[k][k])`, where `H[·][k]`
//! is the PTDF of a transfer from branch `k`'s from bus to its to bus, and `LODF[k][k] = −1`.
//! A transfer doesn't depend on the slack, so neither does the LODF. Taking out a bridge
//! splits the network, and `1 − H[k][k]` is zero in exact arithmetic, so rounding decides what
//! the division gives; here a bridge's column is NaN, as with `makeLODF`'s `mask_bridge`, with
//! bridges found from the graph rather than the arithmetic.

use std::borrow::Cow;
use std::collections::VecDeque;

use sparse_mat::compressed_matrix::Ordering;
use sparse_mat::sparse_matrix::SparseVector;

use super::{DcModel, DcStudy, Factor};
use crate::islands::bridges;
use crate::powerflow_error::PowerflowError;

/// How a PTDF's injections are balanced.
#[derive(Debug, Clone, PartialEq, Default)]
#[non_exhaustive]
pub enum PtdfSlack {
    /// At MATPOWER's default slack, as `makePTDF` chooses it: the first connected bus of type
    /// REF in the case (with or without a generator), or if there is none, the model's first
    /// reference bus.
    #[default]
    Reference,
    /// At one bus: a case bus row, which must be connected.
    Bus(usize),
    /// Spread over buses in proportion to these weights, one per case bus row, as `makePTDF`'s
    /// distributed slack: `PTDF − (PTDF w) 1ᵀ` with the weights normalized to sum to 1. A
    /// bus that isn't connected must have zero weight.
    Weights(Vec<f64>),
}

/// A DC network's sensitivities for one slack, from [`DcStudy::sensitivity`]: PTDF rows,
/// columns, transfers and the full matrix.
///
/// Everything is indexed like the case's tables: rows by case branch row, columns by case bus
/// row. A branch that's out of service, or a bus that isn't connected, gets zeros.
#[derive(Debug, Clone)]
pub struct DcSensitivity<'a> {
    model: &'a DcModel,
    /// The single slack bus (internal) the factorization leaves out
    slack: usize,
    /// For a distributed slack, each internal bus's normalized weight
    weights: Option<Vec<f64>>,
    /// Each internal bus's position among the factorization's unknowns, `None` for the slack
    position: Vec<Option<usize>>,
    factor: Option<Cow<'a, Factor>>,
    /// Whether each internal branch is a bridge of the model's graph
    bridge: Vec<bool>,
}

impl DcStudy<'_> {
    /// The network's sensitivities with the given slack.
    ///
    /// When the slack is the model's only reference bus, B′'s factorization serves. Otherwise
    /// (several references, or another slack bus) Bbus without the slack's row and column is
    /// factored here, as `makePTDF` does: it treats any other reference bus as an ordinary bus.
    /// A distributed slack uses the default slack's factorization, since its correction makes
    /// the result independent of that choice.
    ///
    /// Returns `PowerflowError::RowOutOfRange` for a slack bus row past the bus table,
    /// `InvalidSlack` for a slack bus that isn't connected or unusable weights,
    /// `IslandWithoutReference` naming a bus that can't reach the slack (in a model with
    /// several islands), and `SingularNetwork` if the reduced Bbus is singular.
    pub fn sensitivity(&self, slack: PtdfSlack) -> Result<DcSensitivity<'_>, PowerflowError> {
        let model = self.model;
        let topology = &model.topology;
        let nb = topology.bus_rows.len();
        let mut internal = vec![None; topology.case_shape.0];
        for (index, &row) in topology.bus_rows.iter().enumerate() {
            internal[row] = Some(index);
        }
        let (slack_bus, weights) = match slack {
            PtdfSlack::Reference => (model.matpower_slack, None),
            PtdfSlack::Bus(row) => {
                let Some(&index) = internal.get(row) else {
                    return Err(PowerflowError::RowOutOfRange {
                        table: "bus",
                        row,
                        len: internal.len(),
                    });
                };
                let index = index.ok_or(PowerflowError::InvalidSlack {
                    reason: "the slack bus isn't connected",
                })?;
                (index, None)
            }
            PtdfSlack::Weights(weights) => {
                if weights.len() != internal.len() {
                    return Err(PowerflowError::InvalidSlack {
                        reason: "the weights need one per bus",
                    });
                }
                let mut normalized = vec![0.0; nb];
                for (row, &weight) in weights.iter().enumerate() {
                    match internal[row] {
                        Some(index) => normalized[index] = weight,
                        None if weight != 0.0 => {
                            return Err(PowerflowError::InvalidSlack {
                                reason: "a bus that isn't connected has a nonzero weight",
                            })
                        }
                        None => {}
                    }
                }
                let total: f64 = normalized.iter().sum();
                if total == 0.0 || !total.is_finite() {
                    return Err(PowerflowError::InvalidSlack {
                        reason: "the weights don't have a finite, nonzero sum",
                    });
                }
                for weight in &mut normalized {
                    *weight /= total;
                }
                (model.matpower_slack, Some(normalized))
            }
        };

        // Every bus must reach the slack, or the reduced Bbus is singular
        let mut neighbours = vec![vec![]; nb];
        for &(from, to) in &topology.ends {
            neighbours[from].push(to);
            neighbours[to].push(from);
        }
        let mut reached = vec![false; nb];
        reached[slack_bus] = true;
        let mut queue = VecDeque::from([slack_bus]);
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
                bus_label: model.bus_labels[stranded],
            });
        }

        let (unknown, factor) = if topology.ref_buses == [slack_bus] {
            (
                self.unknown.clone(),
                self.factor.as_ref().map(Cow::Borrowed),
            )
        } else {
            let unknown: Vec<usize> = (0..nb).filter(|&bus| bus != slack_bus).collect();
            let factor = if unknown.is_empty() {
                None
            } else {
                Some(Cow::Owned(model.factor_over(&unknown, Ordering::Amd)?))
            };
            (unknown, factor)
        };
        let mut position = vec![None; nb];
        for (index, &bus) in unknown.iter().enumerate() {
            position[bus] = Some(index);
        }
        let mut bridge = vec![false; topology.ends.len()];
        for k in bridges(nb, &topology.ends) {
            bridge[k] = true;
        }
        Ok(DcSensitivity {
            model,
            slack: slack_bus,
            weights,
            position,
            factor,
            bridge,
        })
    }
}

impl DcSensitivity<'_> {
    /// The slack bus the factorization leaves out, as a case bus row. With a distributed
    /// slack, it's the default slack, which doesn't affect the results.
    pub fn slack_bus(&self) -> usize {
        self.model.topology.bus_rows[self.slack]
    }

    /// PTDF rows for the case branch rows `branches`: each is the branch's sensitivity to an
    /// injection at every bus (one entry per case bus row), from one solve each. An
    /// out-of-service branch's row is zeros. Returns `PowerflowError::RowOutOfRange` for a
    /// row past the branch table.
    pub fn ptdf_rows(&self, branches: &[usize]) -> Result<Vec<Vec<f64>>, PowerflowError> {
        let internal = self.internal_branches(branches)?;
        let topology = &self.model.topology;
        let nb = topology.bus_rows.len();
        Ok(internal
            .into_iter()
            .map(|branch| {
                let mut row = vec![0.0; topology.case_shape.0];
                let Some(k) = branch else {
                    return row;
                };
                // B is symmetric, so row k is b_k (B⁻¹ (e_f − e_t))ᵀ
                let (from, to) = topology.ends[k];
                let theta = self.solve(&[(from, 1.0), (to, -1.0)]);
                let b = self.susceptance(k);
                let mut values: Vec<f64> = (0..nb).map(|bus| b * theta[bus]).collect();
                if let Some(weights) = &self.weights {
                    let shift: f64 = values.iter().zip(weights).map(|(h, w)| h * w).sum();
                    for value in &mut values {
                        *value -= shift;
                    }
                }
                for (&case_row, value) in topology.bus_rows.iter().zip(values) {
                    row[case_row] = value;
                }
                row
            })
            .collect())
    }

    /// PTDF columns for the case bus rows `buses`: each is every branch's sensitivity to an
    /// injection at that bus (one entry per case branch row), from one solve each. A bus that
    /// isn't connected gets zeros. Returns `PowerflowError::RowOutOfRange` for a row past the
    /// bus table.
    pub fn ptdf_columns(&self, buses: &[usize]) -> Result<Vec<Vec<f64>>, PowerflowError> {
        let internal = self.internal_buses(buses)?;
        let shift = self.distributed_shift();
        Ok(internal
            .into_iter()
            .map(|bus| self.ptdf_column(bus, shift.as_deref()))
            .collect())
    }

    /// [`DcSensitivity::ptdf_columns`] with the columns computed in parallel (feature `rayon`).
    /// Each column is computed by one thread with exactly the sequential operations, so the
    /// result is bit for bit `ptdf_columns`'s.
    #[cfg(feature = "rayon")]
    pub fn par_ptdf_columns(&self, buses: &[usize]) -> Result<Vec<Vec<f64>>, PowerflowError> {
        use rayon::prelude::*;
        let internal = self.internal_buses(buses)?;
        let shift = self.distributed_shift();
        Ok(internal
            .into_par_iter()
            .map(|bus| self.ptdf_column(bus, shift.as_deref()))
            .collect())
    }

    /// With a distributed slack, the flows the weighted withdrawal causes, `H w`, which every
    /// column loses.
    fn distributed_shift(&self) -> Option<Vec<f64>> {
        self.weights.as_ref().map(|weights| {
            let injections: Vec<(usize, f64)> = weights.iter().copied().enumerate().collect();
            self.flows(&self.solve(&injections))
        })
    }

    /// One PTDF column, for an internal bus (`None`: not connected, zeros).
    fn ptdf_column(&self, bus: Option<usize>, shift: Option<&[f64]>) -> Vec<f64> {
        let Some(bus) = bus else {
            return vec![0.0; self.model.topology.case_shape.2];
        };
        let mut column = self.flows(&self.solve(&[(bus, 1.0)]));
        if let Some(shift) = shift {
            for (value, shift) in column.iter_mut().zip(shift) {
                *value -= shift;
            }
        }
        column
    }

    /// PTDF columns for transfers from one case bus row to another: each is every branch's
    /// sensitivity to an injection at the first bus withdrawn at the second (one entry per case
    /// branch row), as `makePTDF`'s `txfr`. A transfer doesn't depend on the slack. A bus
    /// that isn't connected injects nothing. Returns `PowerflowError::RowOutOfRange` for a row
    /// past the bus table.
    pub fn ptdf_transfers(
        &self,
        transfers: &[(usize, usize)],
    ) -> Result<Vec<Vec<f64>>, PowerflowError> {
        let from: Vec<usize> = transfers.iter().map(|&(from, _)| from).collect();
        let to: Vec<usize> = transfers.iter().map(|&(_, to)| to).collect();
        let (from, to) = (self.internal_buses(&from)?, self.internal_buses(&to)?);
        Ok(from
            .into_iter()
            .zip(to)
            .map(|(from, to)| {
                let mut injections = vec![];
                injections.extend(from.map(|bus| (bus, 1.0)));
                injections.extend(to.map(|bus| (bus, -1.0)));
                self.flows(&self.solve(&injections))
            })
            .collect())
    }

    /// The whole PTDF matrix, as `makePTDF` returns it: one row per case branch row, each with
    /// one entry per case bus row. It's dense, so it suits small cases: one solve per bus.
    pub fn ptdf(&self) -> Vec<Vec<f64>> {
        let (nbus, nl) = (
            self.model.topology.case_shape.0,
            self.model.topology.case_shape.2,
        );
        let all: Vec<usize> = (0..nbus).collect();
        let columns = self.ptdf_columns(&all).expect("every bus row is in range");
        (0..nl)
            .map(|k| columns.iter().map(|column| column[k]).collect())
            .collect()
    }

    /// LODF columns for taking out the case branch rows `outaged`: each is every branch's change
    /// in flow as a fraction of the outaged branch's flow before (one entry per case branch
    /// row), with −1 for the outaged branch itself, from one solve each. A bridge's column is
    /// NaN. An out-of-service branch's entry is zero. Returns `PowerflowError::RowOutOfRange`
    /// for a row past the branch table, and `BranchOutOfService` for a branch not in service.
    pub fn lodf_columns(&self, outaged: &[usize]) -> Result<Vec<Vec<f64>>, PowerflowError> {
        let internal = self.internal_branches(outaged)?;
        outaged
            .iter()
            .zip(internal)
            .map(|(&row, branch)| self.lodf_column(row, branch))
            .collect()
    }

    /// [`DcSensitivity::lodf_columns`] with the columns computed in parallel (feature `rayon`),
    /// bit for bit `lodf_columns`'s. An error is the first in `outaged`'s order, as there.
    #[cfg(feature = "rayon")]
    pub fn par_lodf_columns(&self, outaged: &[usize]) -> Result<Vec<Vec<f64>>, PowerflowError> {
        use rayon::prelude::*;
        let internal = self.internal_branches(outaged)?;
        // Collecting a parallel iterator straight into a Result stops at some error, not
        //   necessarily the first, so collect every result and take the first error in order
        let results: Vec<Result<Vec<f64>, PowerflowError>> = outaged
            .par_iter()
            .zip(internal)
            .map(|(&row, branch)| self.lodf_column(row, branch))
            .collect();
        results.into_iter().collect()
    }

    /// One LODF column, for the outage of case branch row `row` (internal branch `branch`,
    /// `None` if it isn't in service).
    fn lodf_column(&self, row: usize, branch: Option<usize>) -> Result<Vec<f64>, PowerflowError> {
        let k = branch.ok_or(PowerflowError::BranchOutOfService { row })?;
        if self.bridge[k] {
            return Ok(vec![f64::NAN; self.model.topology.case_shape.2]);
        }
        let (from, to) = self.model.topology.ends[k];
        let mut column = self.flows(&self.solve(&[(from, 1.0), (to, -1.0)]));
        let remaining = 1.0 - column[row];
        for value in &mut column {
            *value /= remaining;
        }
        column[row] = -1.0;
        Ok(column)
    }

    /// The whole LODF matrix, as `makeLODF` returns it (with bridges masked): one row per case
    /// branch row, each with one entry per case branch row, so `lodf()[l][k]` is branch `l`'s
    /// change for the outage of branch `k`. The column of a branch that's out of service is
    /// zeros, since taking it out changes nothing. It's dense, so it suits small cases: one
    /// solve per in-service branch.
    pub fn lodf(&self) -> Vec<Vec<f64>> {
        let topology = &self.model.topology;
        let nl = topology.case_shape.2;
        let columns = self
            .lodf_columns(&topology.branch_rows)
            .expect("in-service branch rows are in range and in service");
        let mut matrix = vec![vec![0.0; nl]; nl];
        for (&k, column) in topology.branch_rows.iter().zip(&columns) {
            for (row, &value) in matrix.iter_mut().zip(column) {
                row[k] = value;
            }
        }
        matrix
    }

    /// The internal branch of each case branch row, `None` if it isn't in service, checking
    /// the rows are in range.
    pub(super) fn internal_branches(
        &self,
        rows: &[usize],
    ) -> Result<Vec<Option<usize>>, PowerflowError> {
        let topology = &self.model.topology;
        let len = topology.case_shape.2;
        rows.iter()
            .map(|&row| {
                if row >= len {
                    Err(PowerflowError::RowOutOfRange {
                        table: "branch",
                        row,
                        len,
                    })
                } else {
                    Ok(topology.branch_rows.binary_search(&row).ok())
                }
            })
            .collect()
    }

    /// The internal bus of each case bus row, `None` if it isn't connected, checking the rows
    /// are in range.
    fn internal_buses(&self, rows: &[usize]) -> Result<Vec<Option<usize>>, PowerflowError> {
        let topology = &self.model.topology;
        let len = topology.case_shape.0;
        rows.iter()
            .map(|&row| {
                if row >= len {
                    Err(PowerflowError::RowOutOfRange {
                        table: "bus",
                        row,
                        len,
                    })
                } else {
                    Ok(topology.bus_rows.binary_search(&row).ok())
                }
            })
            .collect()
    }

    /// The angles (one per internal bus, the slack's zero) for injections at internal buses,
    /// by one solve with a sparse right-hand side. The slack's injection is dropped: the
    /// slack balances it.
    pub(super) fn solve(&self, injections: &[(usize, f64)]) -> Vec<f64> {
        let nb = self.position.len();
        let mut theta = vec![0.0; nb];
        let Some(factor) = &self.factor else {
            return theta;
        };
        let mut entries: Vec<(usize, f64)> = injections
            .iter()
            .filter_map(|&(bus, value)| self.position[bus].map(|position| (position, value)))
            .collect();
        entries.sort_unstable_by_key(|&(position, _)| position);
        entries.dedup_by(|later, earlier| {
            let same = later.0 == earlier.0;
            if same {
                earlier.1 += later.1;
            }
            same
        });
        let rhs = SparseVector::new(
            nb - 1,
            entries.iter().map(|&(position, _)| position).collect(),
            entries.iter().map(|&(_, value)| value).collect(),
        )
        .expect("positions are below the number of unknowns, increasing and distinct");
        let solved = factor.solve_sparse(&rhs);
        for (bus, position) in self.position.iter().enumerate() {
            if let Some(position) = position {
                theta[bus] = solved[*position];
            }
        }
        theta
    }

    /// Each case branch's from-end flow, `b_k (θ_f − θ_t)`, for angles at the internal buses:
    /// zero for a branch out of service. No phase-shifter term, since a PTDF is a sensitivity.
    pub(super) fn flows(&self, theta: &[f64]) -> Vec<f64> {
        let topology = &self.model.topology;
        let mut flows = vec![0.0; topology.case_shape.2];
        for (k, (&row, &(from, to))) in topology.branch_rows.iter().zip(&topology.ends).enumerate()
        {
            flows[row] = self.susceptance(k) * (theta[from] - theta[to]);
        }
        flows
    }

    /// Internal branch `k`'s susceptance, as the model stores it.
    pub(super) fn susceptance(&self, k: usize) -> f64 {
        self.model.susceptance(k)
    }
}
