//! AC power flow, computed as MATPOWER's `runpf` does with its default options.
//!
//! [`AcModel::new`] builds the bus admittance matrix Ybus and the branch matrices Yf and Yt once,
//! following MATPOWER's `ext2int`, `bustypes` and `makeYbus`, through assembly maps so branch
//! parameters can be re-stamped without rebuilding them ([`AcModel::restamp`]). Elements are
//! dropped as in the DC power flow: isolated buses, and out-of-service generators and branches
//! or those attached to an isolated bus.

use std::f64::consts::PI;

use sparse_mat::compressed_matrix::{AssemblyMap, CsrMatrix};
use sparse_mat::Complex;

use crate::case::{Branch, Case};
use crate::powerflow_error::PowerflowError;
use crate::topology::Topology;

/// The AC model of a case: the admittance matrices in per unit over the connected buses and
/// in-service branches, with each bus's role in the solve.
///
/// Internal bus `i` is case bus row `bus_rows()[i]`, and internal branch `k` is case branch row
/// `branch_rows()[k]`; both keep the case's order.
#[derive(Debug, Clone)]
pub struct AcModel {
    base_mva: f64,
    ybus: CsrMatrix<Complex<f64>>,
    yf: CsrMatrix<Complex<f64>>,
    yt: CsrMatrix<Complex<f64>>,
    ybus_map: AssemblyMap,
    yf_map: AssemblyMap,
    yt_map: AssemblyMap,
    topology: Topology,
}

/// The values a case stamps into Ybus, Yf and Yt, in the order `AcModel::new` lists their
/// positions.
struct Stamps {
    ybus: Vec<Complex<f64>>,
    yf: Vec<Complex<f64>>,
    yt: Vec<Complex<f64>>,
}

/// A branch's π-model stamps, as in `makeYbus`: series admittance `Ys = 1 / (r + jx)`, line
/// charging `Bc`, and the complex tap `τ = tap · e^{jφ}` (a tap of 0 meaning 1). Returns
/// `(Yff, Yft, Ytf, Ytt)`.
fn branch_stamps(branch: &Branch) -> [Complex<f64>; 4] {
    let ys = Complex::new(1.0, 0.0) / Complex::new(branch.r, branch.x);
    let tap = Complex::from_polar(branch.effective_tap(), branch.shift_deg * PI / 180.0);
    let ytt = ys + Complex::new(0.0, branch.b / 2.0);
    let yff = ytt / (tap * tap.conj());
    let yft = -ys / tap.conj();
    let ytf = -ys / tap;
    [yff, yft, ytf, ytt]
}

/// `PowerflowError::ZeroImpedance` for an in-service branch with `r = x = 0`, whose series
/// admittance is infinite.
fn check_impedance(row: usize, branch: &Branch) -> Result<(), PowerflowError> {
    if branch.r == 0.0 && branch.x == 0.0 {
        Err(PowerflowError::ZeroImpedance { row })
    } else {
        Ok(())
    }
}

impl AcModel {
    /// Builds the AC model of `case`, as MATPOWER's `makeYbus` does after `ext2int`.
    ///
    /// Each in-service branch `k` from bus `f` to bus `t` stamps its π-model: `Yff` at (f, f),
    /// `Yft` at (f, t), `Ytf` at (t, f) and `Ytt` at (t, t) of Ybus, with
    /// `Ytt = Ys + jBc/2`, `Yff = Ytt / |τ|²`, `Yft = −Ys / conj(τ)` and `Ytf = −Ys / τ`, where
    /// `Ys = 1 / (r + jx)` and `τ = tap · e^{jφ}` (a tap of 0 meaning 1). Row `k` of Yf holds
    /// `Yff` and `Yft`, and row `k` of Yt holds `Ytf` and `Ytt`, so `Yf V` and `Yt V` are the
    /// branch currents at each end. Each connected bus adds its shunt `(Gs + jBs) / baseMVA` on
    /// the diagonal. Parallel branches sum, and every stamped position stays stored, so the
    /// patterns depend only on the topology.
    ///
    /// Bus roles and the checks follow [`DcModel::new`](crate::dc::DcModel::new), except that
    /// the branch check is `PowerflowError::ZeroImpedance` for an in-service branch with
    /// `r = x = 0`. Returns the errors of `Case::validate`, then `ZeroImpedance`,
    /// `NoReferenceBus` or `IslandWithoutReference`.
    pub fn new(case: &Case) -> Result<AcModel, PowerflowError> {
        let topology = Topology::new(case, check_impedance)?;
        let nb = topology.bus_rows.len();
        let nl = topology.branch_rows.len();
        let stamps = Self::stamps(case, &topology);

        let mut ybus_positions = Vec::with_capacity(4 * nl + nb);
        let mut yf_positions = Vec::with_capacity(2 * nl);
        let mut yt_positions = Vec::with_capacity(2 * nl);
        for (k, &(f, t)) in topology.ends.iter().enumerate() {
            ybus_positions.extend([(f, f), (f, t), (t, f), (t, t)]);
            yf_positions.extend([(k, f), (k, t)]);
            yt_positions.extend([(k, f), (k, t)]);
        }
        ybus_positions.extend((0..nb).map(|i| (i, i)));

        let build = |rows, positions: &[(usize, usize)], values: &[Complex<f64>]| {
            CsrMatrix::from_triplets_mapped(
                rows,
                nb,
                positions
                    .iter()
                    .zip(values)
                    .map(|(&(i, j), &value)| (i, j, value)),
            )
            .expect("internal indices are in range")
        };
        let (ybus, ybus_map) = build(nb, &ybus_positions, &stamps.ybus);
        let (yf, yf_map) = build(nl, &yf_positions, &stamps.yf);
        let (yt, yt_map) = build(nl, &yt_positions, &stamps.yt);
        Ok(AcModel {
            base_mva: case.base_mva,
            ybus,
            yf,
            yt,
            ybus_map,
            yf_map,
            yt_map,
            topology,
        })
    }

    /// The stamp values of `case`'s branches and bus shunts, in the order `new` lists their
    /// positions: four Ybus stamps per branch then one shunt per bus, and two Yf and two Yt
    /// stamps per branch.
    fn stamps(case: &Case, topology: &Topology) -> Stamps {
        let nl = topology.branch_rows.len();
        let mut stamps = Stamps {
            ybus: Vec::with_capacity(4 * nl + topology.bus_rows.len()),
            yf: Vec::with_capacity(2 * nl),
            yt: Vec::with_capacity(2 * nl),
        };
        for &row in &topology.branch_rows {
            let [yff, yft, ytf, ytt] = branch_stamps(&case.branch[row]);
            stamps.ybus.extend([yff, yft, ytf, ytt]);
            stamps.yf.extend([yff, yft]);
            stamps.yt.extend([ytf, ytt]);
        }
        stamps.ybus.extend(topology.bus_rows.iter().map(|&row| {
            let bus = &case.bus[row];
            Complex::new(bus.gs_mw, bus.bs_mvar) / case.base_mva
        }));
        stamps
    }

    /// Re-stamps Ybus, Yf and Yt from `case`'s branch impedances, charging, taps and phase
    /// shifts and its bus shunts, through the assembly maps, without rebuilding the matrices.
    ///
    /// `case` must have the topology the model was built from: the same tables, and the same
    /// statuses and bus types, since those decide which branches and buses the matrices cover.
    /// Returns `PowerflowError::CaseMismatch` if its tables' lengths differ, or
    /// `ZeroImpedance` for an in-service branch with `r = x = 0`; on an error the model is
    /// unchanged.
    pub fn restamp(&mut self, case: &Case) -> Result<(), PowerflowError> {
        self.topology.check_case(case)?;
        for &row in &self.topology.branch_rows {
            check_impedance(row + 1, &case.branch[row])?;
        }
        let stamps = Self::stamps(case, &self.topology);
        let refilled = "the maps were built for these matrices from the same stamp list";
        self.ybus
            .refill(&self.ybus_map, &stamps.ybus)
            .expect(refilled);
        self.yf.refill(&self.yf_map, &stamps.yf).expect(refilled);
        self.yt.refill(&self.yt_map, &stamps.yt).expect(refilled);
        self.base_mva = case.base_mva;
        Ok(())
    }

    /// The complex power each connected bus injects (p.u.), as MATPOWER's `makeSbus` computes
    /// it for constant-power loads: `(Pg + jQg)` of its in-service generators less its load
    /// `(Pd + jQd)`, over `baseMVA`. Read from `case`, which must have the model's tables.
    pub fn sbus_pu(&self, case: &Case) -> Vec<Complex<f64>> {
        let mut sbus: Vec<Complex<f64>> = self
            .topology
            .bus_rows
            .iter()
            .map(|&row| {
                let bus = &case.bus[row];
                -Complex::new(bus.pd_mw, bus.qd_mvar) / self.base_mva
            })
            .collect();
        for (&row, &bus) in self.topology.gen_rows.iter().zip(&self.topology.gen_buses) {
            let generator = &case.generator[row];
            sbus[bus] += Complex::new(generator.pg_mw, generator.qg_mvar) / self.base_mva;
        }
        sbus
    }

    /// The system MVA base the per-unit values are relative to.
    pub fn base_mva(&self) -> f64 {
        self.base_mva
    }

    /// The bus admittance matrix, connected buses × connected buses (p.u.), in CSR form with
    /// one stored entry for every position a branch or shunt stamps.
    pub fn ybus(&self) -> &CsrMatrix<Complex<f64>> {
        &self.ybus
    }

    /// The from-end branch admittance matrix, in-service branches × connected buses (p.u.):
    /// `Yf V` is each branch's current at its from end.
    pub fn yf(&self) -> &CsrMatrix<Complex<f64>> {
        &self.yf
    }

    /// The to-end branch admittance matrix, in-service branches × connected buses (p.u.):
    /// `Yt V` is each branch's current at its to end.
    pub fn yt(&self) -> &CsrMatrix<Complex<f64>> {
        &self.yt
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
