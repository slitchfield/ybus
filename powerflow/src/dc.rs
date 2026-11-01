//! DC power flow, computed exactly as MATPOWER's `rundcpf` does.
//!
//! The DC approximation keeps only real power and voltage angles: every voltage magnitude is
//! 1 p.u., resistances and line charging are ignored, and angle differences are small enough
//! that sin θ ≈ θ. Branch flows are then linear in the bus angles, `P = B θ`, and one sparse
//! solve gives the angles.
//!
//! [`DcModel::new`] builds the matrices once, following MATPOWER's `ext2int`, `bustypes` and
//! `makeBdc`; [`DcModel::solve`] computes the injections and solves, following `runpf` and
//! `dcpf`. [`dc_power_flow`] does both. Results are indexed like the case's own tables, so
//! callers never see the internal numbering.
//!
//! Elements are dropped from the solve as MATPOWER's `ext2int` drops them:
//! - **isolated buses** (type 4), whose load is then not served, and whose angle is reported
//!   unchanged from the case
//! - **out-of-service generators and branches**, and those connected to an isolated bus. Their
//!   results are reported as 0.0, as `runpf` reports them (lines 474–480 of MATPOWER 8.1's
//!   `runpf.m`).

use std::collections::{HashMap, VecDeque};
use std::f64::consts::PI;

use sparse_mat::compressed_matrix::{CsrMatrix, Ordering};
use sparse_mat::sparse_error::SparseError;

use crate::case::{BusType, Case};
use crate::powerflow_error::PowerflowError;

/// The DC model of a case: the B matrices and phase-shifter injections in per unit, over the
/// connected buses and in-service branches, with each bus's role in the solve.
///
/// Internal bus `i` is case bus row `bus_rows()[i]`, and internal branch `k` is case branch row
/// `branch_rows()[k]`; both keep the case's order.
#[derive(Debug, Clone)]
pub struct DcModel {
    base_mva: f64,
    bbus: CsrMatrix<f64>,
    bf: CsrMatrix<f64>,
    pbusinj_pu: Vec<f64>,
    pfinj_pu: Vec<f64>,
    bus_rows: Vec<usize>,
    branch_rows: Vec<usize>,
    gen_rows: Vec<usize>,
    /// The internal bus of each in-service generator, parallel to `gen_rows`
    gen_buses: Vec<usize>,
    ref_buses: Vec<usize>,
    pv_buses: Vec<usize>,
    pq_buses: Vec<usize>,
    /// The table lengths of the case the model was built from: buses, generators, branches
    case_shape: (usize, usize, usize),
}

impl DcModel {
    /// Builds the DC model of `case`, as MATPOWER's `makeBdc` does after `ext2int`.
    ///
    /// Each in-service branch has series susceptance `b = 1 / (x · tap)`, using
    /// [`Branch::effective_tap`](crate::case::Branch::effective_tap), and stamps
    /// `Bbus = Cftᵀ · diag(b) · Cft` and `Bf = diag(b) · Cft`, where `Cft` has +1 at the
    /// branch's from bus and −1 at its to bus. Parallel branches sum. A phase shift of `φ`
    /// degrees adds the injection `−b · φ · π/180` at the from end (`pfinj_pu`), and its
    /// effect on each bus (`pbusinj_pu = Cftᵀ · pfinj_pu`).
    ///
    /// Bus roles follow MATPOWER's `bustypes`: a REF or PV bus counts as one only if it has an
    /// in-service generator, and is otherwise treated as PQ. With no REF bus, the first PV bus
    /// becomes the reference. A case may have several reference buses.
    ///
    /// Validates `case` first with [`Case::validate`], returning its errors: the model looks
    /// buses up by label, so an unknown label would silently drop a generator or branch, and a
    /// duplicate one would attach it to the wrong bus. Then returns
    /// `PowerflowError::ZeroReactance` for an in-service branch with `x == 0`,
    /// `NoReferenceBus` if no bus can be the reference, and `IslandWithoutReference` if some
    /// connected bus can't reach a reference bus through in-service branches, which would make
    /// the solve singular.
    pub fn new(case: &Case) -> Result<DcModel, PowerflowError> {
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
                if branch.x == 0.0 {
                    return Err(PowerflowError::ZeroReactance { row: row + 1 });
                }
                branch_rows.push(row);
                ends.push((from, to));
            }
        }

        // makeBdc. Each branch's (from, to) and (to, from) stamps are pushed together, so both
        //   positions sum the same values in the same order and Bbus is exactly symmetric,
        //   which cholesky requires. Every stamped position stays stored even if its stamps
        //   cancel, so the pattern depends only on the topology.
        let nl = branch_rows.len();
        let mut bbus_triplets = Vec::with_capacity(4 * nl);
        let mut bf_triplets = Vec::with_capacity(2 * nl);
        let mut pfinj_pu = Vec::with_capacity(nl);
        let mut pbusinj_pu = vec![0.0; nb];
        for (k, (&row, &(from, to))) in branch_rows.iter().zip(&ends).enumerate() {
            let branch = &case.branch[row];
            let b = 1.0 / branch.x / branch.effective_tap();
            bbus_triplets.extend([(from, from, b), (from, to, -b), (to, from, -b), (to, to, b)]);
            bf_triplets.extend([(k, from, b), (k, to, -b)]);
            let pfinj = b * (-branch.shift_deg * PI / 180.0);
            pfinj_pu.push(pfinj);
            pbusinj_pu[from] += pfinj;
            pbusinj_pu[to] -= pfinj;
        }
        let bbus = CsrMatrix::from_triplets(nb, nb, bbus_triplets)
            .expect("internal bus indices are below nb");
        let bf = CsrMatrix::from_triplets(nl, nb, bf_triplets)
            .expect("internal branch and bus indices are in range");

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

        // Every connected bus must reach a reference bus, or B′ is singular: an island with
        //   no reference has no fixed angle. Checked here rather than inferred from a failed
        //   factorization, which can't name the island and may not fail at all, since
        //   rounding can leave a tiny positive pivot.
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

        Ok(DcModel {
            base_mva: case.base_mva,
            bbus,
            bf,
            pbusinj_pu,
            pfinj_pu,
            bus_rows,
            branch_rows,
            gen_rows,
            gen_buses,
            ref_buses,
            pv_buses,
            pq_buses,
            case_shape: (case.bus.len(), case.generator.len(), case.branch.len()),
        })
    }

    /// The system MVA base the per-unit values are relative to.
    pub fn base_mva(&self) -> f64 {
        self.base_mva
    }

    /// The bus susceptance matrix, connected buses × connected buses (p.u.), in CSR form with
    /// one stored entry for every position a branch stamps.
    pub fn bbus(&self) -> &CsrMatrix<f64> {
        &self.bbus
    }

    /// The branch susceptance matrix, in-service branches × connected buses (p.u.), which maps
    /// bus angles to from-end branch flows, in CSR form.
    pub fn bf(&self) -> &CsrMatrix<f64> {
        &self.bf
    }

    /// The phase shifters' equivalent injection at each connected bus (p.u.).
    pub fn pbusinj_pu(&self) -> &[f64] {
        &self.pbusinj_pu
    }

    /// The phase shifters' contribution to each in-service branch's from-end flow (p.u.).
    pub fn pfinj_pu(&self) -> &[f64] {
        &self.pfinj_pu
    }

    /// The case bus row of each internal bus.
    pub fn bus_rows(&self) -> &[usize] {
        &self.bus_rows
    }

    /// The case branch row of each internal branch.
    pub fn branch_rows(&self) -> &[usize] {
        &self.branch_rows
    }

    /// The case generator rows of the in-service generators at connected buses.
    pub fn gen_rows(&self) -> &[usize] {
        &self.gen_rows
    }

    /// The reference buses (internal indices), in case order.
    pub fn ref_buses(&self) -> &[usize] {
        &self.ref_buses
    }

    /// The PV buses (internal indices), in case order.
    pub fn pv_buses(&self) -> &[usize] {
        &self.pv_buses
    }

    /// The PQ buses (internal indices), in case order. Includes REF and PV buses that have no
    /// in-service generator.
    pub fn pq_buses(&self) -> &[usize] {
        &self.pq_buses
    }
}

/// The result of a DC power flow, indexed like the case's tables.
#[derive(Debug, Clone, PartialEq)]
pub struct DcSolution {
    /// Each bus's voltage angle (degrees), one per case bus row. An isolated bus keeps its
    /// angle from the case.
    pub va_deg: Vec<f64>,
    /// Each branch's real power flow at the from end (MW), one per case branch row. The flow at
    /// the to end is its negation, since DC power flow is lossless. 0.0 for a branch that's out
    /// of service or touches an isolated bus.
    pub pf_mw: Vec<f64>,
    /// Each generator's real power output (MW), one per case generator row. The first
    /// generator at each reference bus also picks up that bus's mismatch (the slack). 0.0 for a
    /// generator that's out of service or at an isolated bus.
    pub pg_mw: Vec<f64>,
    /// Which factorization solved the system.
    pub solver: DcSolver,
}

/// The factorization a DC power flow used for B′ (Bbus without the reference rows and
/// columns).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DcSolver {
    /// Cholesky, which applies when B′ is symmetric positive definite: always, unless some
    /// branch has a negative reactance (a series capacitor).
    Cholesky,
    /// LU with partial pivoting, the fallback when Cholesky finds B′ isn't positive definite.
    /// (It's also the fallback if B′ isn't exactly symmetric, which the model's construction
    /// rules out.)
    Lu,
    /// No solve was needed: every connected bus is a reference bus.
    NoSolve,
}

impl DcModel {
    /// Solves the DC power flow for `case`'s injections, as MATPOWER's `runpf` and `dcpf` do.
    ///
    /// `case` must be the case the model was built from, or a copy of it with different
    /// injections or reference angles: `pd_mw`, `gs_mw`, `pg_mw` and `va_deg` are read from
    /// `case`, and everything else from the model. Changing anything else, such as a status or
    /// a reactance, needs a new model.
    ///
    /// The net injection at each connected bus is `(Pg − Pd − Gs) / baseMVA − pbusinj`, with
    /// `Pg` summed over its in-service generators. Reference bus angles are fixed at their
    /// `va_deg`, and the remaining angles solve `B′ θ = P′`, by Cholesky ordered by approximate
    /// minimum degree, or by LU if B′ isn't positive definite. Branch flows are
    /// `(Bf θ + pfinj) · baseMVA`.
    ///
    /// Returns `PowerflowError::CaseMismatch` if `case`'s tables have different lengths from
    /// the model's case, and `SingularNetwork` if B′ is singular even though every bus reaches
    /// a reference, which takes reactances that cancel out.
    pub fn solve(&self, case: &Case) -> Result<DcSolution, PowerflowError> {
        self.solve_with_ordering(case, Ordering::Amd)
    }

    /// [`DcModel::solve`] with B′'s Cholesky factorization ordered by `ordering` instead of the
    /// default, approximate minimum degree, for comparing orderings. The LU fallback isn't ordered
    /// yet. Returns the errors of `solve`.
    pub fn solve_with_ordering(
        &self,
        case: &Case,
        ordering: Ordering,
    ) -> Result<DcSolution, PowerflowError> {
        let found = (case.bus.len(), case.generator.len(), case.branch.len());
        if found != self.case_shape {
            return Err(PowerflowError::CaseMismatch {
                expected: self.case_shape,
                found,
            });
        }
        let base_mva = self.base_mva;
        let to_radians = PI / 180.0;

        // makeSbus's real part, less the phase-shifter injections and shunt conductance
        let mut pbus: Vec<f64> = self
            .bus_rows
            .iter()
            .zip(&self.pbusinj_pu)
            .map(|(&row, pbusinj)| {
                let bus = &case.bus[row];
                -(bus.pd_mw + bus.gs_mw) / base_mva - pbusinj
            })
            .collect();
        for (&row, &bus) in self.gen_rows.iter().zip(&self.gen_buses) {
            pbus[bus] += case.generator[row].pg_mw / base_mva;
        }
        let mut theta: Vec<f64> = self
            .bus_rows
            .iter()
            .map(|&row| case.bus[row].va_deg * to_radians)
            .collect();

        // dcpf: θ[pv, pq] = B[pv, pq; pv, pq] \ (P[pv, pq] − B[pv, pq; ref] · θ[ref]), with
        //   the unknowns in MATPOWER's order
        let unknown: Vec<usize> = self
            .pv_buses
            .iter()
            .chain(&self.pq_buses)
            .copied()
            .collect();
        let solver = if unknown.is_empty() {
            DcSolver::NoSolve
        } else {
            let b_unknown = self
                .bbus
                .select(&unknown, &unknown)
                .expect("indices in range");
            let b_ref = self
                .bbus
                .select(&unknown, &self.ref_buses)
                .expect("indices in range");
            let theta_ref: Vec<f64> = self.ref_buses.iter().map(|&bus| theta[bus]).collect();
            let coupling = b_ref
                .mul_vec(&theta_ref)
                .expect("one angle per reference bus");
            let rhs: Vec<f64> = unknown
                .iter()
                .zip(&coupling)
                .map(|(&bus, coupling)| pbus[bus] - coupling)
                .collect();
            let singular = |position: usize| PowerflowError::SingularNetwork {
                bus_label: case.bus[self.bus_rows[unknown[position]]].bus_label,
            };
            // Bbus is built exactly symmetric, but LU is the fallback for NotSymmetric too, in
            //   case that ever changes: a slower solve beats a panic
            let cholesky = b_unknown
                .analyze_cholesky(ordering)
                .and_then(|symbolic| symbolic.factor(&b_unknown))
                .and_then(|factor| factor.solve(&rhs));
            let (solved, solver) = match cholesky {
                Ok(solved) => (solved, DcSolver::Cholesky),
                Err(SparseError::NotPositiveDefinite { .. } | SparseError::NotSymmetric) => {
                    match b_unknown.solve(&rhs) {
                        Ok(solved) => (solved, DcSolver::Lu),
                        Err(SparseError::SingularMatrix { pivot }) => return Err(singular(pivot)),
                        Err(err) => unreachable!("LU of a square system: {err}"),
                    }
                }
                Err(err) => unreachable!("Cholesky of a square system: {err}"),
            };
            // A nearly singular B′ can factor with a tiny pivot and give huge or non-finite
            //   angles instead of an error; MATPOWER's dcpf treats |θ| > 1e5 as failure.
            if let Some(position) = solved.iter().position(|t| !t.is_finite() || t.abs() > 1e5) {
                return Err(singular(position));
            }
            for (&bus, value) in unknown.iter().zip(solved) {
                theta[bus] = value;
            }
            solver
        };

        let mut va_deg: Vec<f64> = case.bus.iter().map(|bus| bus.va_deg).collect();
        for (&row, angle) in self.bus_rows.iter().zip(&theta) {
            va_deg[row] = angle / to_radians;
        }

        let mut pf_mw = vec![0.0; case.branch.len()];
        let bf_theta = self.bf.mul_vec(&theta).expect("one angle per bus");
        for ((&row, flow), pfinj) in self.branch_rows.iter().zip(bf_theta).zip(&self.pfinj_pu) {
            pf_mw[row] = (flow + pfinj) * base_mva;
        }

        // The slack: the first in-service generator at each reference bus takes up the
        //   difference between the bus's solved injection and its scheduled one
        let mut pg_mw = vec![0.0; case.generator.len()];
        for &row in &self.gen_rows {
            pg_mw[row] = case.generator[row].pg_mw;
        }
        let injection = self.bbus.mul_vec(&theta).expect("one angle per bus");
        for &bus in &self.ref_buses {
            let first = self
                .gen_buses
                .iter()
                .position(|&gen_bus| gen_bus == bus)
                .expect("a reference bus has an in-service generator");
            pg_mw[self.gen_rows[first]] += (injection[bus] - pbus[bus]) * base_mva;
        }

        Ok(DcSolution {
            va_deg,
            pf_mw,
            pg_mw,
            solver,
        })
    }
}

/// Runs a DC power flow on `case`: [`DcModel::new`] followed by [`DcModel::solve`], with
/// their errors.
pub fn dc_power_flow(case: &Case) -> Result<DcSolution, PowerflowError> {
    DcModel::new(case)?.solve(case)
}
