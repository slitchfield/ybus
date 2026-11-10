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

use std::f64::consts::PI;
use std::sync::OnceLock;

use sparse_mat::compressed_matrix::{
    Cholesky, CsrMatrix, Lu, LuOptions, Ordering, SparseSolveWorkspace,
};
use sparse_mat::sparse_error::SparseError;
use sparse_mat::sparse_matrix::SparseVector;
use sparse_mat::Factorization;

use crate::case::{BusType, Case};
use crate::powerflow_error::PowerflowError;
use crate::topology::Topology;
use reach::Reachability;

#[cfg(feature = "cuda")]
mod gpu;
mod outage;
mod reach;
mod sensitivity;
#[cfg(feature = "cuda")]
pub use gpu::{DcGpuSensitivity, DcGpuStudy, OutageLoading};
pub use sensitivity::{DcSensitivity, PtdfSlack};

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
    topology: Topology,
    /// The bus label of each internal bus, for errors that name one
    bus_labels: Vec<usize>,
    /// MATPOWER's default PTDF slack (internal): the first connected bus of type REF, or if
    /// there is none, the first reference bus `bustypes` chose
    matpower_slack: usize,
    /// Which outages strand buses, built on the first outage that asks
    reachability: OnceLock<Reachability>,
}

/// MATPOWER's `makeBdc` over a network in internal numbering: the B matrices and phase-shifter
/// injections (p.u.), shared by the DC power flow and the DC optimal power flow.
pub(crate) struct Bdc {
    pub(crate) bbus: CsrMatrix<f64>,
    pub(crate) bf: CsrMatrix<f64>,
    pub(crate) pbusinj_pu: Vec<f64>,
    pub(crate) pfinj_pu: Vec<f64>,
}

impl Bdc {
    /// The matrices over `nb` connected buses and the in-service branches at case rows
    /// `branch_rows`, whose internal (from, to) buses are `ends`. Each branch's reactance must
    /// be nonzero.
    pub(crate) fn new(
        case: &Case,
        nb: usize,
        branch_rows: &[usize],
        ends: &[(usize, usize)],
    ) -> Bdc {
        // Each branch's (from, to) and (to, from) stamps are pushed together, so both positions
        //   sum the same values in the same order and Bbus is exactly symmetric, which cholesky
        //   requires. Every stamped position stays stored even if its stamps cancel, so the
        //   pattern depends only on the topology.
        let nl = branch_rows.len();
        let mut bbus_triplets = Vec::with_capacity(4 * nl);
        let mut bf_triplets = Vec::with_capacity(2 * nl);
        let mut pfinj_pu = Vec::with_capacity(nl);
        let mut pbusinj_pu = vec![0.0; nb];
        for (k, (&row, &(from, to))) in branch_rows.iter().zip(ends).enumerate() {
            let branch = &case.branch[row];
            let b = 1.0 / branch.x / branch.effective_tap();
            bbus_triplets.extend([(from, from, b), (from, to, -b), (to, from, -b), (to, to, b)]);
            bf_triplets.extend([(k, from, b), (k, to, -b)]);
            let pfinj = b * (-branch.shift_deg * PI / 180.0);
            pfinj_pu.push(pfinj);
            pbusinj_pu[from] += pfinj;
            pbusinj_pu[to] -= pfinj;
        }
        Bdc {
            bbus: CsrMatrix::from_triplets(nb, nb, bbus_triplets)
                .expect("internal bus indices are below nb"),
            bf: CsrMatrix::from_triplets(nl, nb, bf_triplets)
                .expect("internal branch and bus indices are in range"),
            pbusinj_pu,
            pfinj_pu,
        }
    }
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
        let topology = Topology::new(case, |row, branch| {
            if branch.x == 0.0 {
                Err(PowerflowError::ZeroReactance { row })
            } else {
                Ok(())
            }
        })?;
        let Bdc {
            bbus,
            bf,
            pbusinj_pu,
            pfinj_pu,
        } = Bdc::new(
            case,
            topology.bus_rows.len(),
            &topology.branch_rows,
            &topology.ends,
        );

        let bus_labels = topology
            .bus_rows
            .iter()
            .map(|&row| case.bus[row].bus_label)
            .collect();
        let matpower_slack = topology
            .bus_rows
            .iter()
            .position(|&row| case.bus[row].bus_type == BusType::Ref)
            .unwrap_or(topology.ref_buses[0]);
        Ok(DcModel {
            base_mva: case.base_mva,
            bbus,
            bf,
            pbusinj_pu,
            pfinj_pu,
            topology,
            bus_labels,
            matpower_slack,
            reachability: OnceLock::new(),
        })
    }

    /// Which outages strand buses (see `reach`), built once, on first use.
    fn reachability(&self) -> &Reachability {
        self.reachability
            .get_or_init(|| Reachability::new(&self.topology))
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
        &self.topology.bus_rows
    }

    /// The case branch row of each internal branch.
    pub fn branch_rows(&self) -> &[usize] {
        &self.topology.branch_rows
    }

    /// The case generator rows of the in-service generators at connected buses.
    pub fn gen_rows(&self) -> &[usize] {
        &self.topology.gen_rows
    }

    /// The reference buses (internal indices), in case order.
    pub fn ref_buses(&self) -> &[usize] {
        &self.topology.ref_buses
    }

    /// The PV buses (internal indices), in case order.
    pub fn pv_buses(&self) -> &[usize] {
        &self.topology.pv_buses
    }

    /// The PQ buses (internal indices), in case order. Includes REF and PV buses that have no
    /// in-service generator.
    pub fn pq_buses(&self) -> &[usize] {
        &self.topology.pq_buses
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

/// Options for the DC power flow. `Default` gives MATPOWER's method, with B′ ordered by
/// approximate minimum degree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct DcOptions {
    /// The fill-reducing ordering B′ is factored under (by Cholesky, or by the LU fallback):
    /// [`Ordering::Amd`] by default. The others are for comparing orderings; every ordering
    /// gives the same solution, to within rounding.
    pub ordering: Ordering,
}

impl Default for DcOptions {
    fn default() -> DcOptions {
        DcOptions {
            ordering: Ordering::Amd,
        }
    }
}

/// The factorization a DC power flow used for B′ (Bbus without the reference rows and
/// columns).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DcSolver {
    /// Cholesky, which applies when B′ is symmetric positive definite: always, unless some
    /// branch has a negative reactance (a series capacitor).
    Cholesky,
    /// LU with threshold partial pivoting, the fallback when Cholesky finds B′ isn't positive
    /// definite.
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
    /// `case`, and everything else from the model. A changed topology (a bus's label or type,
    /// a generator's bus or status, or a branch's ends or status) is refused; any other change,
    /// such as a reactance, is ignored. Either needs a new model.
    ///
    /// The net injection at each connected bus is `(Pg − Pd − Gs) / baseMVA − pbusinj`, with
    /// `Pg` summed over its in-service generators. Reference bus angles are fixed at their
    /// `va_deg`, and the remaining angles solve `B′ θ = P′`, by Cholesky ordered by approximate
    /// minimum degree, or, if B′ isn't positive definite, by LU with the same ordering and
    /// threshold partial pivoting. Branch flows are
    /// `(Bf θ + pfinj) · baseMVA`. B′ is ordered as `options` choose.
    ///
    /// Returns `PowerflowError::CaseMismatch` if `case`'s tables have different lengths from
    /// the model's case, `TopologyMismatch` if its topology differs, and `SingularNetwork` if B′
    /// is singular even though every bus reaches a reference, which takes reactances that cancel
    /// out.
    pub fn solve(&self, case: &Case, options: &DcOptions) -> Result<DcSolution, PowerflowError> {
        // The case is checked before B′ is factored, as it always was
        self.topology.check_case(case)?;
        self.study(options)?.solve(case)
    }

    /// Factors B′ once, for any number of solves and studies against it: Bbus over the PV and
    /// PQ buses (MATPOWER's order), by Cholesky under the ordering `options` choose
    /// (approximate minimum degree by default), or, if B′ isn't positive definite, by LU with
    /// the same ordering and threshold partial pivoting.
    ///
    /// Returns `PowerflowError::SingularNetwork` if B′ is singular even though every bus
    /// reaches a reference, which takes reactances that cancel out.
    pub fn study(&self, options: &DcOptions) -> Result<DcStudy<'_>, PowerflowError> {
        let ordering = options.ordering;
        // dcpf's unknowns: θ[pv, pq], in MATPOWER's order
        let unknown = self.topology.pv_pq();
        if unknown.is_empty() {
            return Ok(DcStudy {
                model: self,
                unknown,
                factor: None,
            });
        }
        let factor = self.factor_over(&unknown, ordering)?;
        Ok(DcStudy {
            model: self,
            unknown,
            factor: Some(factor),
        })
    }

    /// Internal branch `k`'s susceptance `b = 1 / (x · tap)`, as Bf stores it at the from end.
    /// A branch from a bus to itself stamps b and −b in one place, so it reads 0 there, which
    /// is harmless: its angle difference is zero.
    pub(crate) fn susceptance(&self, k: usize) -> f64 {
        let (from, _) = self.topology.ends[k];
        let position = self
            .bf
            .position(k, from)
            .expect("in range")
            .expect("Bf stores each branch's from end");
        self.bf.values()[position]
    }

    /// Factors Bbus over the internal buses `unknown` (rows and columns in that order), by
    /// Cholesky under `ordering`, or by LU if that isn't positive definite. Returns
    /// `PowerflowError::SingularNetwork`, naming a bus, if it's singular.
    fn factor_over(&self, unknown: &[usize], ordering: Ordering) -> Result<Factor, PowerflowError> {
        let b_unknown = self
            .bbus
            .select(unknown, unknown)
            .expect("indices in range");
        let singular = |position: usize| PowerflowError::SingularNetwork {
            bus_label: self.bus_labels[unknown[position]],
        };
        // Bbus is built exactly symmetric, but LU is the fallback for NotSymmetric too, in
        //   case that ever changes: a slower solve beats a panic
        let cholesky = b_unknown
            .analyze_cholesky(ordering)
            .and_then(|symbolic| symbolic.factor(&b_unknown));
        let factor = match cholesky {
            Ok(factor) => Factor::Cholesky(factor),
            Err(SparseError::NotPositiveDefinite { .. } | SparseError::NotSymmetric) => {
                // If B′ is singular, LU's error names a column of the matrix it factored that
                //   has no usable pivot, numbered as in B′ itself (not in the AMD order it
                //   factored in). B′ was selected from Bbus with `unknown` as its columns, so
                //   column j of B′ is bus `unknown[j]`, which `singular` turns into that bus's
                //   label for the error
                let b_columns = b_unknown.to_csc();
                let lu = b_columns
                    .analyze_lu(ordering)
                    .and_then(|symbolic| symbolic.factor(&b_columns, LuOptions::default()));
                match lu {
                    Ok(factor) => Factor::Lu(factor),
                    Err(SparseError::SingularMatrix { pivot }) => return Err(singular(pivot)),
                    Err(err) => unreachable!("LU of a square system: {err}"),
                }
            }
            Err(err) => unreachable!("Cholesky of a square system: {err}"),
        };
        Ok(factor)
    }
}

/// A factorization of B′ (Bbus without the reference rows and columns), from
/// [`DcModel::study`], kept for any number of solves and studies of the same network.
#[derive(Debug, Clone)]
pub struct DcStudy<'a> {
    model: &'a DcModel,
    /// The PV then PQ buses (internal indices): B′'s rows and columns, in MATPOWER's order
    unknown: Vec<usize>,
    /// `None` when every connected bus is a reference bus, so there is nothing to solve
    factor: Option<Factor>,
}

/// B′'s factorization: Cholesky when it's positive definite, LU otherwise
#[derive(Debug, Clone)]
enum Factor {
    Cholesky(Cholesky<f64>),
    Lu(Lu<f64>),
}

impl Factor {
    fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        match self {
            Factor::Cholesky(factor) => factor.solve(rhs),
            Factor::Lu(factor) => factor.solve(rhs),
        }
        .expect("one entry per unknown")
    }

    fn solve_sparse_into(
        &self,
        rhs: &SparseVector<f64>,
        workspace: &mut SparseSolveWorkspace<f64>,
        x: &mut [f64],
    ) {
        match self {
            Factor::Cholesky(factor) => factor.solve_sparse_into(rhs, workspace, x),
            Factor::Lu(factor) => factor.solve_sparse_into(rhs, workspace, x),
        }
        .expect("one entry per unknown")
    }

    fn solve_sparse(&self, rhs: &SparseVector<f64>) -> Vec<f64> {
        match self {
            Factor::Cholesky(factor) => factor.solve_sparse(rhs),
            Factor::Lu(factor) => factor.solve_sparse(rhs),
        }
        .expect("one entry per unknown")
    }
}

impl DcStudy<'_> {
    /// The model this study factored.
    pub fn model(&self) -> &DcModel {
        self.model
    }

    /// Which factorization B′ has.
    pub fn solver(&self) -> DcSolver {
        match self.factor {
            Some(Factor::Cholesky(_)) => DcSolver::Cholesky,
            Some(Factor::Lu(_)) => DcSolver::Lu,
            None => DcSolver::NoSolve,
        }
    }

    /// Solves the DC power flow for `case`'s injections, exactly as [`DcModel::solve`] does
    /// (which factors B′ and calls this), with the factorization kept.
    ///
    /// Returns `PowerflowError::CaseMismatch` if `case`'s tables have different lengths from
    /// the model's case, `TopologyMismatch` if its topology differs, and `SingularNetwork` if B′
    /// is so nearly singular that some angle is non-finite or above 1e5 radians, which
    /// MATPOWER's `dcpf` treats as failure.
    pub fn solve(&self, case: &Case) -> Result<DcSolution, PowerflowError> {
        self.model.topology.check_case(case)?;
        let (pbus, mut theta) = self.scheduled(case);
        if let Some(factor) = &self.factor {
            let solved = factor.solve(&self.rhs(&pbus, &theta));
            self.place(case, &mut theta, &solved)?;
        }
        Ok(self.finish(case, &theta, &pbus, &[]))
    }

    /// makeSbus's real part at each connected bus, less the phase-shifter injections and shunt
    /// conductance (p.u.), and each connected bus's angle from the case (radians).
    fn scheduled(&self, case: &Case) -> (Vec<f64>, Vec<f64>) {
        let model = self.model;
        let base_mva = model.base_mva;
        let mut pbus: Vec<f64> = model
            .topology
            .bus_rows
            .iter()
            .zip(&model.pbusinj_pu)
            .map(|(&row, pbusinj)| {
                let bus = &case.bus[row];
                -(bus.pd_mw + bus.gs_mw) / base_mva - pbusinj
            })
            .collect();
        for (&row, &bus) in model
            .topology
            .gen_rows
            .iter()
            .zip(&model.topology.gen_buses)
        {
            pbus[bus] += case.generator[row].pg_mw / base_mva;
        }
        let theta = model
            .topology
            .bus_rows
            .iter()
            .map(|&row| case.bus[row].va_deg * (PI / 180.0))
            .collect();
        (pbus, theta)
    }

    /// dcpf's right-hand side over the unknowns: `P[pv, pq] − B[pv, pq; ref] · θ[ref]`.
    fn rhs(&self, pbus: &[f64], theta: &[f64]) -> Vec<f64> {
        self.rhs_with(&self.reference_coupling(), pbus, theta)
    }

    /// `B[pv, pq; ref]`, which couples the unknowns to the reference buses' fixed angles.
    fn reference_coupling(&self) -> CsrMatrix<f64> {
        self.model
            .bbus
            .select(&self.unknown, &self.model.topology.ref_buses)
            .expect("indices in range")
    }

    /// [`DcStudy::rhs`] with [`DcStudy::reference_coupling`] already made, for many solves.
    fn rhs_with(&self, b_ref: &CsrMatrix<f64>, pbus: &[f64], theta: &[f64]) -> Vec<f64> {
        let model = self.model;
        let theta_ref: Vec<f64> = model
            .topology
            .ref_buses
            .iter()
            .map(|&bus| theta[bus])
            .collect();
        let coupling = b_ref
            .mul_vec(&theta_ref)
            .expect("one angle per reference bus");
        self.unknown
            .iter()
            .zip(&coupling)
            .map(|(&bus, coupling)| pbus[bus] - coupling)
            .collect()
    }

    /// Writes the solved unknowns into `theta`. A nearly singular B′ can factor with a tiny
    /// pivot and give huge or non-finite angles instead of an error; MATPOWER's dcpf treats
    /// |θ| > 1e5 as failure, `PowerflowError::SingularNetwork` here.
    fn place(&self, case: &Case, theta: &mut [f64], solved: &[f64]) -> Result<(), PowerflowError> {
        if let Some(position) = solved.iter().position(|t| !t.is_finite() || t.abs() > 1e5) {
            let row = self.model.topology.bus_rows[self.unknown[position]];
            return Err(PowerflowError::SingularNetwork {
                bus_label: case.bus[row].bus_label,
            });
        }
        for (&bus, &value) in self.unknown.iter().zip(solved) {
            theta[bus] = value;
        }
        Ok(())
    }

    /// The solution from the connected buses' angles (radians) and scheduled injections, with
    /// the internal branches `outaged` taken out: their flows are zero, and they carry nothing
    /// into the reference buses' injections.
    fn finish(&self, case: &Case, theta: &[f64], pbus: &[f64], outaged: &[usize]) -> DcSolution {
        let model = self.model;
        let base_mva = model.base_mva;
        let mut va_deg: Vec<f64> = case.bus.iter().map(|bus| bus.va_deg).collect();
        for (&row, angle) in model.topology.bus_rows.iter().zip(theta) {
            va_deg[row] = angle / (PI / 180.0);
        }

        // Each in-service branch's flow, row by row as `Bf θ` computes it
        let mut pf_mw = vec![0.0; case.branch.len()];
        for (k, (&row, pfinj)) in model
            .topology
            .branch_rows
            .iter()
            .zip(&model.pfinj_pu)
            .enumerate()
        {
            pf_mw[row] = (row_dot(&model.bf, k, theta) + pfinj) * base_mva;
        }
        // `Bbus θ`, only at the reference buses, whose injections the slack needs, less what
        //   the outaged branches carried there
        let ref_buses = &model.topology.ref_buses;
        let mut injection: Vec<f64> = ref_buses
            .iter()
            .map(|&bus| row_dot(&model.bbus, bus, theta))
            .collect();
        for &k in outaged {
            pf_mw[model.topology.branch_rows[k]] = 0.0;
            let (from, to) = model.topology.ends[k];
            let flow = model.susceptance(k) * (theta[from] - theta[to]);
            for (r, &bus) in ref_buses.iter().enumerate() {
                if bus == from {
                    injection[r] -= flow;
                }
                if bus == to {
                    injection[r] += flow;
                }
            }
        }

        // The slack: the first in-service generator at each reference bus takes up the
        //   difference between the bus's solved injection and its scheduled one
        let mut pg_mw = vec![0.0; case.generator.len()];
        for &row in &model.topology.gen_rows {
            pg_mw[row] = case.generator[row].pg_mw;
        }
        for (r, &bus) in ref_buses.iter().enumerate() {
            let first = model
                .topology
                .gen_buses
                .iter()
                .position(|&gen_bus| gen_bus == bus)
                .expect("a reference bus has an in-service generator");
            pg_mw[model.topology.gen_rows[first]] += (injection[r] - pbus[bus]) * base_mva;
        }

        DcSolution {
            va_deg,
            pf_mw,
            pg_mw,
            solver: self.solver(),
        }
    }
}

/// Runs a DC power flow on `case` with the default options; see [`dc_power_flow_with_options`].
pub fn dc_power_flow(case: &Case) -> Result<DcSolution, PowerflowError> {
    dc_power_flow_with_options(case, &DcOptions::default())
}

/// Runs a DC power flow on `case` with `options`: [`DcModel::new`] followed by
/// [`DcModel::solve`], with their errors.
pub fn dc_power_flow_with_options(
    case: &Case,
    options: &DcOptions,
) -> Result<DcSolution, PowerflowError> {
    DcModel::new(case)?.solve(case, options)
}

/// Row `row` of `m` times `x`, summed from zero in the row's order, exactly as `mul_vec`
/// computes each entry of `m x`: so a single entry is bit for bit the product's.
fn row_dot(m: &CsrMatrix<f64>, row: usize, x: &[f64]) -> f64 {
    let (cols, values) = m.outer(row).expect("a row of the matrix");
    cols.iter()
        .zip(values)
        .fold(0.0, |sum, (&col, &value)| sum + value * x[col])
}
