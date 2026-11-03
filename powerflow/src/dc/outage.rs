//! DC power flow with branches taken out, by rank-k updates of B′'s factorization instead of
//! refactoring.
//!
//! Taking out branches `K` subtracts `b_k c_k c_kᵀ` from Bbus for each, where `c_k = e_f − e_t`.
//! Over the unknowns, with `A` = B′, `C`'s columns the `c_k` there and `D = diag(b_k)`, the
//! Sherman–Morrison–Woodbury identity gives
//!
//! `(A − C D Cᵀ)⁻¹ r = x + Y W⁻¹ Cᵀ x`, with `x = A⁻¹ r`, `Y = A⁻¹ C` and `W = D⁻¹ − Cᵀ Y`,
//!
//! so the outage needs one sparse-RHS solve per outaged branch (the columns of `Y`) and a k×k
//! solve, against the kept factorization. The right-hand side changes too: an outaged phase
//! shifter's injection goes, and so does an outaged branch's coupling to a reference bus. Both
//! are combinations of the `c_k`, so `x` is the base solution plus a combination of `Y`'s
//! columns, with no further solve.

use std::collections::VecDeque;

use sparse_mat::compressed_matrix::CscMatrix;
use sparse_mat::sparse_matrix::SparseVector;

use super::{DcSolution, DcStudy};
use crate::case::Case;
use crate::powerflow_error::PowerflowError;

impl DcStudy<'_> {
    /// The DC power flow for `case`'s injections with the case branch rows `outaged` taken out,
    /// without refactoring B′: the same result, to within rounding, as building a model of the
    /// case with those branches out of service and solving it.
    ///
    /// Returns the errors of [`DcStudy::solve`], `PowerflowError::RowOutOfRange` for a row past
    /// the branch table, `BranchOutOfService` for a branch that isn't in the model (out of
    /// service, or touching an isolated bus), and `IslandWithoutReference` if the outage leaves
    /// a bus unable to reach a reference bus, as a new model would report. An outage that
    /// splits off an island with its own reference bus is solved.
    pub fn solve_outage(
        &self,
        case: &Case,
        outaged: &[usize],
    ) -> Result<DcSolution, PowerflowError> {
        self.solve_outages(case, &[outaged.to_vec()])?
            .pop()
            .expect("one outage, one result")
    }

    /// [`DcStudy::solve_outage`] for each set of case branch rows in `outages`, sharing the
    /// base solve: for a single-branch outage, each costs one sparse-RHS solve. The outer
    /// result holds the errors that apply to every outage (those of [`DcStudy::solve`]), and
    /// each inner one an outage's own.
    pub fn solve_outages(
        &self,
        case: &Case,
        outages: &[Vec<usize>],
    ) -> Result<Vec<Result<DcSolution, PowerflowError>>, PowerflowError> {
        let (pbus, theta, base, position) = self.outage_base(case)?;
        Ok(outages
            .iter()
            .map(|rows| self.outage(case, rows, &pbus, &theta, &base, &position))
            .collect())
    }

    /// [`DcStudy::solve_outages`] with the outages solved in parallel (feature `rayon`), after
    /// one base solve. Each outage is solved by one thread with exactly the sequential
    /// operations, so the results are bit for bit `solve_outages`'s, in the same order.
    #[cfg(feature = "rayon")]
    pub fn par_solve_outages(
        &self,
        case: &Case,
        outages: &[Vec<usize>],
    ) -> Result<Vec<Result<DcSolution, PowerflowError>>, PowerflowError> {
        use rayon::prelude::*;
        let (pbus, theta, base, position) = self.outage_base(case)?;
        Ok(outages
            .par_iter()
            .map(|rows| self.outage(case, rows, &pbus, &theta, &base, &position))
            .collect())
    }

    /// The base solve the outages share: the scheduled injections, the base angles (radians),
    /// the base unknowns, and each internal bus's position among the unknowns.
    #[allow(clippy::type_complexity)]
    fn outage_base(
        &self,
        case: &Case,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<Option<usize>>), PowerflowError> {
        let model = self.model;
        model.topology.check_case(case)?;
        let (pbus, mut theta) = self.scheduled(case);
        let base = match &self.factor {
            Some(factor) => {
                let base = factor.solve(&self.rhs(&pbus, &theta));
                // Checked here, as solve checks it, so a singular B′ fails every outage alike
                self.place(case, &mut theta, &base)?;
                base
            }
            None => vec![],
        };
        let mut position = vec![None; model.topology.bus_rows.len()];
        for (index, &bus) in self.unknown.iter().enumerate() {
            position[bus] = Some(index);
        }
        Ok((pbus, theta, base, position))
    }

    /// One outage, given the base solve: the scheduled injections, the base angles (radians)
    /// and the base unknowns, and each internal bus's position among the unknowns.
    fn outage(
        &self,
        case: &Case,
        rows: &[usize],
        pbus: &[f64],
        theta: &[f64],
        base: &[f64],
        position: &[Option<usize>],
    ) -> Result<DcSolution, PowerflowError> {
        let model = self.model;
        let topology = &model.topology;
        let mut outaged = Vec::with_capacity(rows.len());
        for &row in rows {
            if row >= topology.case_shape.2 {
                return Err(PowerflowError::RowOutOfRange {
                    table: "branch",
                    row,
                    len: topology.case_shape.2,
                });
            }
            let k = topology
                .branch_rows
                .binary_search(&row)
                .map_err(|_| PowerflowError::BranchOutOfService { row })?;
            outaged.push(k);
        }
        outaged.sort_unstable();
        outaged.dedup();
        self.check_reach(&outaged)?;

        // The outaged phase shifters' injections leave the scheduled injections
        let mut pbus = pbus.to_vec();
        for &k in &outaged {
            let (from, to) = topology.ends[k];
            pbus[from] += model.pfinj_pu[k];
            pbus[to] -= model.pfinj_pu[k];
        }
        let mut theta = theta.to_vec();
        if let Some(factor) = &self.factor {
            // The outaged branches that touch an unknown, with c_k over the unknowns as
            //   (position, ±1) and α_k, its share of the change in the right-hand side
            let mut columns = vec![];
            for &k in &outaged {
                let (from, to) = topology.ends[k];
                let b = model.susceptance(k);
                let mut entries = vec![];
                let mut reference_part = 0.0;
                match position[from] {
                    Some(p) => entries.push((p, 1.0)),
                    None => reference_part += theta[from],
                }
                match position[to] {
                    Some(p) => entries.push((p, -1.0)),
                    None => reference_part -= theta[to],
                }
                if from == to || entries.is_empty() {
                    continue;
                }
                entries.sort_unstable_by_key(|&(p, _)| p);
                let alpha = model.pfinj_pu[k] + b * reference_part;
                columns.push((entries, b, alpha));
            }
            let n = self.unknown.len();
            // Y's columns, A⁻¹ c_k, and x = A⁻¹ r′ = base + Σ α_k Y_k
            let y: Vec<Vec<f64>> = columns
                .iter()
                .map(|(entries, _, _)| {
                    let c = SparseVector::new(
                        n,
                        entries.iter().map(|&(p, _)| p).collect(),
                        entries.iter().map(|&(_, v)| v).collect(),
                    )
                    .expect("distinct positions below the number of unknowns, sorted");
                    factor.solve_sparse(&c)
                })
                .collect();
            let mut x = base.to_vec();
            for ((_, _, alpha), y_k) in columns.iter().zip(&y) {
                for (x_i, y_i) in x.iter_mut().zip(y_k) {
                    *x_i += alpha * y_i;
                }
            }
            // cᵀ v over the unknowns
            let dot = |entries: &[(usize, f64)], v: &[f64]| -> f64 {
                entries.iter().map(|&(p, sign)| sign * v[p]).sum()
            };
            let m = columns.len();
            if m > 0 {
                // W = D⁻¹ − Cᵀ Y, and θ′ = x + Y W⁻¹ Cᵀ x
                let mut triplets = Vec::with_capacity(m * m);
                for (i, (entries, b, _)) in columns.iter().enumerate() {
                    for (j, y_j) in y.iter().enumerate() {
                        let mut w = -dot(entries, y_j);
                        if i == j {
                            w += 1.0 / b;
                        }
                        triplets.push((i, j, w));
                    }
                }
                let w = CscMatrix::from_triplets(m, m, triplets).expect("indices below m");
                let ctx: Vec<f64> = columns
                    .iter()
                    .map(|(entries, _, _)| dot(entries, &x))
                    .collect();
                let singular = || PowerflowError::SingularNetwork {
                    bus_label: case.bus[topology.bus_rows[topology.ends[outaged[0]].0]].bus_label,
                };
                // Nonsingular in exact arithmetic once every bus still reaches a reference;
                //   only reactances that cancel out can make it singular
                let z = w
                    .lu()
                    .and_then(|lu| lu.solve(&ctx))
                    .map_err(|_| singular())?;
                for (z_j, y_j) in z.iter().zip(&y) {
                    for (x_i, y_i) in x.iter_mut().zip(y_j) {
                        *x_i += z_j * y_i;
                    }
                }
            }
            self.place(case, &mut theta, &x)?;
        }
        Ok(self.finish(case, &theta, &pbus, &outaged))
    }

    /// `PowerflowError::IslandWithoutReference` unless every connected bus still reaches a
    /// reference bus without the internal branches `outaged` (sorted), as a new model would
    /// check.
    fn check_reach(&self, outaged: &[usize]) -> Result<(), PowerflowError> {
        let topology = &self.model.topology;
        let nb = topology.bus_rows.len();
        let mut neighbours = vec![vec![]; nb];
        for (k, &(from, to)) in topology.ends.iter().enumerate() {
            if outaged.binary_search(&k).is_err() {
                neighbours[from].push(to);
                neighbours[to].push(from);
            }
        }
        let mut reached = vec![false; nb];
        let mut queue: VecDeque<usize> = topology.ref_buses.iter().copied().collect();
        for &bus in &topology.ref_buses {
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
        match reached.iter().position(|&reached| !reached) {
            Some(stranded) => Err(PowerflowError::IslandWithoutReference {
                bus_label: self.model.bus_labels[stranded],
            }),
            None => Ok(()),
        }
    }
}
