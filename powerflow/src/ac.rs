//! AC power flow, computed as MATPOWER's `runpf` does with its default options.
//!
//! [`AcModel::new`] builds the bus admittance matrix Ybus and the branch matrices Yf and Yt once,
//! following MATPOWER's `ext2int`, `bustypes` and `makeYbus`, through assembly maps so branch
//! parameters can be re-stamped without rebuilding them ([`AcModel::restamp`]). Elements are
//! dropped as in the DC power flow: isolated buses, and out-of-service generators and branches
//! or those attached to an isolated bus.

use std::f64::consts::PI;

use sparse_mat::compressed_matrix::{
    AssemblyMap, CscMatrix, CsrMatrix, Lu, LuOptions, Ordering, RefactorOutcome,
};
use sparse_mat::sparse_error::SparseError;
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
    pi_stamps(
        branch.r,
        branch.x,
        branch.b,
        branch.effective_tap(),
        branch.shift_deg,
    )
}

/// The π-model stamps `(Yff, Yft, Ytf, Ytt)` of resistance `r`, reactance `x`, charging `b`, tap
/// ratio `tap` and phase shift `shift_deg`.
fn pi_stamps(r: f64, x: f64, b: f64, tap: f64, shift_deg: f64) -> [Complex<f64>; 4] {
    let ys = Complex::new(1.0, 0.0) / Complex::new(r, x);
    let tap = Complex::from_polar(tap, shift_deg * PI / 180.0);
    let ytt = ys + Complex::new(0.0, b / 2.0);
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
        let same_stamps_message = "the maps were built for these matrices from the same stamp list";
        self.ybus
            .refill(&self.ybus_map, &stamps.ybus)
            .expect(same_stamps_message);
        self.yf
            .refill(&self.yf_map, &stamps.yf)
            .expect(same_stamps_message);
        self.yt
            .refill(&self.yt_map, &stamps.yt)
            .expect(same_stamps_message);
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

/// How to solve the AC power flow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AcMethod {
    /// Newton's method in polar form, as MATPOWER's `newtonpf` (`pf.alg = 'NR'`, the default).
    Newton,
    /// Fast-decoupled, XB version, as MATPOWER's `fdpf` with `makeB` (`pf.alg = 'FDXB'`): B′
    /// ignores resistance.
    FastDecoupledXb,
    /// Fast-decoupled, BX version (`pf.alg = 'FDBX'`): B″ ignores resistance.
    FastDecoupledBx,
}

/// Options for the AC power flow. `Default` gives MATPOWER's defaults, for Newton's method.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct AcOptions {
    /// Convergence when the largest real or reactive power mismatch is below this (p.u.).
    /// MATPOWER's `pf.tol`, 1e-8 by default.
    pub tolerance: f64,
    /// The most iterations before giving up: Newton iterations, or for fast-decoupled, P
    /// half-iterations. MATPOWER's `pf.nr.max_it` (10) or `pf.fd.max_it` (30).
    pub max_iterations: usize,
    /// The solution method.
    pub method: AcMethod,
}

impl AcOptions {
    /// MATPOWER's defaults for `method`: tolerance 1e-8, and at most 10 Newton iterations or
    /// 30 fast-decoupled ones.
    pub fn for_method(method: AcMethod) -> AcOptions {
        AcOptions {
            tolerance: 1e-8,
            max_iterations: match method {
                AcMethod::Newton => 10,
                AcMethod::FastDecoupledXb | AcMethod::FastDecoupledBx => 30,
            },
            method,
        }
    }
}

impl Default for AcOptions {
    fn default() -> AcOptions {
        AcOptions::for_method(AcMethod::Newton)
    }
}

/// The result of an AC power flow, indexed like the case's tables.
#[derive(Debug, Clone, PartialEq)]
pub struct AcSolution {
    /// Each bus's voltage magnitude (p.u.), one per case bus row. An isolated bus keeps its
    /// value from the case.
    pub vm_pu: Vec<f64>,
    /// Each bus's voltage angle (degrees), one per case bus row. An isolated bus keeps its
    /// angle from the case.
    pub va_deg: Vec<f64>,
    /// Each generator's real power output (MW), one per case generator row: its scheduled
    /// output, except that the first generator at each reference bus takes up the slack. 0.0
    /// for a generator that's out of service or at an isolated bus.
    pub pg_mw: Vec<f64>,
    /// Each generator's reactive power output (MVAr): at a REF or PV bus, what holds its
    /// voltage, shared among the bus's generators as MATPOWER's `pfsoln` shares it; at a PQ
    /// bus, its scheduled output. 0.0 for a generator that's out of service or at an isolated
    /// bus.
    pub qg_mvar: Vec<f64>,
    /// Each branch's real power flowing into it at the from end (MW), one per case branch
    /// row. 0.0 for a branch that's out of service or touches an isolated bus.
    pub pf_mw: Vec<f64>,
    /// Each branch's reactive power flowing into it at the from end (MVAr).
    pub qf_mvar: Vec<f64>,
    /// Each branch's real power flowing into it at the to end (MW). `pf_mw + pt_mw` is the
    /// branch's loss.
    pub pt_mw: Vec<f64>,
    /// Each branch's reactive power flowing into it at the to end (MVAr).
    pub qt_mvar: Vec<f64>,
    /// How many iterations were taken: Newton iterations, or for fast-decoupled, P
    /// half-iterations, as MATPOWER counts them. 0 if the starting point already met the
    /// tolerance.
    pub iterations: usize,
    /// The largest real or reactive power mismatch at the solution (p.u.), below the tolerance.
    /// Fast-decoupled measures it divided by each bus's voltage magnitude, as `fdpf` does.
    pub mismatch_pu: f64,
    /// How many matrices were factored afresh. Newton's method: 1 when every later iteration
    /// refactored the Jacobian with the first factorization's pivots. Fast-decoupled: 2, B′ and
    /// B″. 0 when no iteration was needed.
    pub full_factorizations: usize,
}

/// The unknowns of a power flow: the PV buses, the PQ buses, and both, in MATPOWER's order.
struct Unknowns {
    pv: Vec<usize>,
    pq: Vec<usize>,
    pv_pq: Vec<usize>,
}

/// A converged iteration: the voltages, and what the solution reports about getting there.
struct Solved {
    v: Vec<Complex<f64>>,
    iterations: usize,
    mismatch_pu: f64,
    full_factorizations: usize,
}

impl AcModel {
    /// Solves the AC power flow for `case`'s injections and voltage setpoints, as MATPOWER's
    /// `runpf` does: by Newton's method (`newtonpf`) or fast-decoupled (`fdpf`), as
    /// `options.method` chooses, then `pfsoln`.
    ///
    /// `case` must be the case the model was built from, or a copy with different injections,
    /// setpoints or starting voltages: `pd_mw`, `qd_mvar`, `pg_mw`, `qg_mvar`, `vg_pu`,
    /// `qmax_mvar`, `qmin_mvar`, `vm_pu` and `va_deg` are read from `case`, and everything else
    /// from the model (fast-decoupled also reads the branch data, to build B′ and B″).
    ///
    /// - **Start:** each bus's case voltage, with each in-service generator's `vg_pu` setting
    ///   the magnitude at REF and PV buses (the last generator listed at a bus wins, as in
    ///   `runpf`).
    /// - **Newton:** the mismatch is `V ∘ conj(Ybus V) − Sbus` (`sbus_pu`); its real part at PV
    ///   and PQ buses and imaginary part at PQ buses must all fall below `options.tolerance`.
    ///   The Jacobian `[[Re dS/dVa, Re dS/dVm], [Im dS/dVa, Im dS/dVm]]` over those unknowns is
    ///   built from Ybus each iteration, keeps one pattern, and is factored by AMD-ordered LU
    ///   once, then refactored with the same pivots.
    /// - **Fast-decoupled:** `makeB`'s B′ (over PV and PQ buses) and B″ (over PQ buses) are each
    ///   factored once by AMD-ordered LU; each iteration corrects the angles with B′, then the
    ///   PQ magnitudes with B″, measuring mismatches divided by the voltage magnitudes.
    /// - **Results:** `pfsoln`'s: the slack at each reference bus, reactive output at REF and PV
    ///   buses, and every in-service branch's flows at both ends.
    ///
    /// Reactive power limits aren't enforced, as in `runpf`'s default. Returns
    /// `PowerflowError::CaseMismatch` if `case`'s tables differ in length from the model's,
    /// `NoConvergence` if the tolerance isn't met within `options.max_iterations`, and
    /// `SingularNetwork` if the Jacobian, B′ or B″ is singular. Fast-decoupled also returns
    /// `ZeroReactance` for an in-service branch with `x = 0`, whose B′ or B″ stamps would be
    /// infinite once resistance is ignored.
    pub fn solve(&self, case: &Case, options: AcOptions) -> Result<AcSolution, PowerflowError> {
        self.topology.check_case(case)?;
        let topology = &self.topology;
        let unknowns = Unknowns {
            pv: topology.pv_buses.clone(),
            pq: topology.pq_buses.clone(),
            pv_pq: topology.pv_pq(),
        };

        // runpf's start: case voltages, with generator setpoints at voltage-controlled buses.
        //   Each setpoint scales the case voltage, so its angle is kept
        let start: Vec<Complex<f64>> = topology
            .bus_rows
            .iter()
            .map(|&row| {
                let bus = &case.bus[row];
                Complex::from_polar(bus.vm_pu, bus.va_deg * PI / 180.0)
            })
            .collect();
        let mut v = start.clone();
        let mut voltage_controlled = vec![true; v.len()];
        for &bus in &unknowns.pq {
            voltage_controlled[bus] = false;
        }
        for (&row, &bus) in topology.gen_rows.iter().zip(&topology.gen_buses) {
            if voltage_controlled[bus] {
                v[bus] = start[bus] * (case.generator[row].vg_pu / start[bus].norm());
            }
        }
        let sbus = self.sbus_pu(case);

        let solved = match options.method {
            AcMethod::Newton => self.newton(case, v, &sbus, &unknowns, options)?,
            AcMethod::FastDecoupledXb => {
                self.fast_decoupled(case, v, &sbus, &unknowns, options, true)?
            }
            AcMethod::FastDecoupledBx => {
                self.fast_decoupled(case, v, &sbus, &unknowns, options, false)?
            }
        };
        Ok(self.results(
            case,
            &solved.v,
            solved.iterations,
            solved.mismatch_pu,
            solved.full_factorizations,
        ))
    }

    /// `newtonpf`, from the starting voltages `v`.
    fn newton(
        &self,
        case: &Case,
        mut v: Vec<Complex<f64>>,
        sbus: &[Complex<f64>],
        unknowns: &Unknowns,
        options: AcOptions,
    ) -> Result<Solved, PowerflowError> {
        let topology = &self.topology;
        let (npv, npq) = (unknowns.pv.len(), unknowns.pq.len());
        let mut vm: Vec<f64> = v.iter().map(|value| value.norm()).collect();
        let mut va: Vec<f64> = v.iter().map(|value| value.arg()).collect();
        let (mut f, mut norm_f) = self.mismatch(&v, sbus, unknowns);
        let mut iterations = 0;
        let mut full_factorizations = 0;
        let mut lu: Option<Lu<f64>> = None;
        let mut converged = norm_f < options.tolerance;
        while !converged && iterations < options.max_iterations {
            iterations += 1;
            let jacobian = self.jacobian_over(&v, unknowns);
            let singular = |column: usize| PowerflowError::SingularNetwork {
                bus_label: case.bus[topology.bus_rows[unknowns.column_bus(column)]].bus_label,
            };
            // Factor once, then reuse the pivots while they stay acceptable
            let factor = match lu.as_mut() {
                Some(lu) => match lu.refactor(&jacobian) {
                    Ok(RefactorOutcome::Reused) => Ok(()),
                    Ok(RefactorOutcome::Repivoted) => {
                        full_factorizations += 1;
                        Ok(())
                    }
                    Err(err) => Err(err),
                },
                None => jacobian
                    .analyze_lu(Ordering::Amd)
                    .and_then(|symbolic| symbolic.factor(&jacobian, LuOptions::default()))
                    .map(|factored| {
                        full_factorizations += 1;
                        lu = Some(factored);
                    }),
            };
            match factor {
                Ok(()) => {}
                Err(SparseError::SingularMatrix { pivot }) => return Err(singular(pivot)),
                Err(err) => unreachable!("LU of a square Jacobian with a fixed pattern: {err}"),
            }
            let rhs: Vec<f64> = f.iter().map(|value| -value).collect();
            let dx = lu
                .as_ref()
                .expect("factored above")
                .solve(&rhs)
                .expect("one mismatch per unknown");

            // newtonpf's update, then the magnitudes and angles again from V, in case a
            //   magnitude went negative
            for (k, &bus) in unknowns.pv.iter().enumerate() {
                va[bus] += dx[k];
            }
            for (k, &bus) in unknowns.pq.iter().enumerate() {
                va[bus] += dx[npv + k];
                vm[bus] += dx[npv + npq + k];
            }
            for ((value, &magnitude), &angle) in v.iter_mut().zip(&vm).zip(&va) {
                *value = Complex::from_polar(magnitude, angle);
            }
            for ((value, magnitude), angle) in v.iter().zip(&mut vm).zip(&mut va) {
                *magnitude = value.norm();
                *angle = value.arg();
            }
            (f, norm_f) = self.mismatch(&v, sbus, unknowns);
            converged = norm_f < options.tolerance;
        }
        if !converged {
            return Err(PowerflowError::NoConvergence {
                iterations,
                mismatch_pu: norm_f,
            });
        }
        Ok(Solved {
            v,
            iterations,
            mismatch_pu: norm_f,
            full_factorizations,
        })
    }

    /// `fdpf` with `makeB`'s matrices, XB if `xb` and BX otherwise, from the starting voltages
    /// `v`.
    fn fast_decoupled(
        &self,
        case: &Case,
        mut v: Vec<Complex<f64>>,
        sbus: &[Complex<f64>],
        unknowns: &Unknowns,
        options: AcOptions,
        xb: bool,
    ) -> Result<Solved, PowerflowError> {
        let topology = &self.topology;
        let mut vm: Vec<f64> = v.iter().map(|value| value.norm()).collect();
        let mut va: Vec<f64> = v.iter().map(|value| value.arg()).collect();
        // A P half-iteration reads only the real mismatches, and a Q half-iteration only the
        //   reactive ones computed just before it
        let (mut p, _, mut norm) = self.decoupled_mismatch(&v, &vm, sbus, unknowns);
        if norm < options.tolerance {
            return Ok(Solved {
                v,
                iterations: 0,
                mismatch_pu: norm,
                full_factorizations: 0,
            });
        }

        // B′ over the PV and PQ buses, B″ over the PQ buses, each factored once
        let (b_p, b_pp) = self.make_b(case, xb)?;
        let factor = |matrix: &CsrMatrix<f64>, buses: &[usize]| {
            let reduced = matrix
                .select(buses, buses)
                .expect("bus indices are in range")
                .to_csc();
            reduced.lu().map_err(|err| match err {
                SparseError::SingularMatrix { pivot } => PowerflowError::SingularNetwork {
                    bus_label: case.bus[topology.bus_rows[buses[pivot]]].bus_label,
                },
                err => unreachable!("LU of a square matrix: {err}"),
            })
        };
        let lu_p = factor(&b_p, &unknowns.pv_pq)?;
        let lu_pp = factor(&b_pp, &unknowns.pq)?;

        let mut iterations = 0;
        while iterations < options.max_iterations {
            iterations += 1;
            // The P half-iteration: angles at PV and PQ buses
            let dva = lu_p.solve(&p).expect("one mismatch per unknown");
            for (&bus, d) in unknowns.pv_pq.iter().zip(dva) {
                va[bus] -= d;
            }
            for ((value, &magnitude), &angle) in v.iter_mut().zip(&vm).zip(&va) {
                *value = Complex::from_polar(magnitude, angle);
            }
            let (_, q, after_p) = self.decoupled_mismatch(&v, &vm, sbus, unknowns);
            norm = after_p;
            if norm < options.tolerance {
                break;
            }
            // The Q half-iteration: magnitudes at PQ buses
            let dvm = lu_pp.solve(&q).expect("one mismatch per unknown");
            for (&bus, d) in unknowns.pq.iter().zip(dvm) {
                vm[bus] -= d;
            }
            for ((value, &magnitude), &angle) in v.iter_mut().zip(&vm).zip(&va) {
                *value = Complex::from_polar(magnitude, angle);
            }
            let (next_p, _, after_q) = self.decoupled_mismatch(&v, &vm, sbus, unknowns);
            norm = after_q;
            if norm < options.tolerance {
                break;
            }
            p = next_p;
        }
        if norm >= options.tolerance {
            return Err(PowerflowError::NoConvergence {
                iterations,
                mismatch_pu: norm,
            });
        }
        Ok(Solved {
            v,
            iterations,
            mismatch_pu: norm,
            full_factorizations: 2,
        })
    }

    /// `fdpf`'s mismatches, `mis = (V ∘ conj(Ybus V) − Sbus) / Vm`: the real parts at PV and
    /// PQ buses, the imaginary parts at PQ buses, and the largest magnitude among them (both
    /// must fall below the tolerance, so their maximum must).
    fn decoupled_mismatch(
        &self,
        v: &[Complex<f64>],
        vm: &[f64],
        sbus: &[Complex<f64>],
        unknowns: &Unknowns,
    ) -> (Vec<f64>, Vec<f64>, f64) {
        let current = self.ybus.mul_vec(v).expect("one voltage per bus");
        let mis: Vec<Complex<f64>> = v
            .iter()
            .zip(&current)
            .zip(sbus)
            .zip(vm)
            .map(|(((&v, &i), &s), &vm)| (v * i.conj() - s) / vm)
            .collect();
        let p: Vec<f64> = unknowns.pv_pq.iter().map(|&bus| mis[bus].re).collect();
        let q: Vec<f64> = unknowns.pq.iter().map(|&bus| mis[bus].im).collect();
        let norm = p
            .iter()
            .chain(&q)
            .fold(0.0, |largest: f64, value| largest.max(value.abs()));
        (p, q, norm)
    }

    /// `makeB`'s B′ and B″ over all connected buses, as `−Im` of admittance matrices built with
    /// changed branch data. B′: line charging, bus shunts and taps removed (phase shifts kept),
    /// and for XB resistance too. B″: phase shifts removed, and for BX resistance too, with
    /// charging, taps and bus shunts kept. Returns `PowerflowError::ZeroReactance` for an
    /// in-service branch with `x = 0`.
    fn make_b(
        &self,
        case: &Case,
        xb: bool,
    ) -> Result<(CsrMatrix<f64>, CsrMatrix<f64>), PowerflowError> {
        let topology = &self.topology;
        let nb = topology.bus_rows.len();
        let nl = topology.branch_rows.len();
        let mut prime = Vec::with_capacity(4 * nl);
        let mut double_prime = Vec::with_capacity(4 * nl + nb);
        for (&row, &(f, t)) in topology.branch_rows.iter().zip(&topology.ends) {
            let branch = &case.branch[row];
            if branch.x == 0.0 {
                return Err(PowerflowError::ZeroReactance { row: row + 1 });
            }
            let r_prime = if xb { 0.0 } else { branch.r };
            let r_double_prime = if xb { branch.r } else { 0.0 };
            let stamps_prime = pi_stamps(r_prime, branch.x, 0.0, 1.0, branch.shift_deg);
            let stamps_double_prime = pi_stamps(
                r_double_prime,
                branch.x,
                branch.b,
                branch.effective_tap(),
                0.0,
            );
            for (matrix, stamps) in [
                (&mut prime, stamps_prime),
                (&mut double_prime, stamps_double_prime),
            ] {
                let [yff, yft, ytf, ytt] = stamps;
                matrix.extend([
                    (f, f, -yff.im),
                    (f, t, -yft.im),
                    (t, f, -ytf.im),
                    (t, t, -ytt.im),
                ]);
            }
        }
        double_prime.extend(
            topology
                .bus_rows
                .iter()
                .enumerate()
                .map(|(i, &row)| (i, i, -case.bus[row].bs_mvar / case.base_mva)),
        );
        let build =
            |triplets| CsrMatrix::from_triplets(nb, nb, triplets).expect("indices in range");
        Ok((build(prime), build(double_prime)))
    }

    /// `newtonpf`'s mismatch vector `[Re mis[pv]; Re mis[pq]; Im mis[pq]]`, with
    /// `mis = V ∘ conj(Ybus V) − Sbus`, and its infinity norm.
    fn mismatch(
        &self,
        v: &[Complex<f64>],
        sbus: &[Complex<f64>],
        unknowns: &Unknowns,
    ) -> (Vec<f64>, f64) {
        let current = self.ybus.mul_vec(v).expect("one voltage per bus");
        let mis: Vec<Complex<f64>> = v
            .iter()
            .zip(&current)
            .zip(sbus)
            .map(|((&v, &i), &s)| v * i.conj() - s)
            .collect();
        let f: Vec<f64> = unknowns
            .pv_pq
            .iter()
            .map(|&bus| mis[bus].re)
            .chain(unknowns.pq.iter().map(|&bus| mis[bus].im))
            .collect();
        let norm = f
            .iter()
            .fold(0.0, |largest: f64, value| largest.max(value.abs()));
        (f, norm)
    }

    /// The Newton–Raphson Jacobian at the voltages `v` (complex, per unit, one per connected bus
    /// in internal order, as `bus_rows` lists them), in CSC form: the derivatives of the real
    /// power mismatch at PV and PQ buses and the reactive mismatch at PQ buses, with respect to
    /// the angles at PV and PQ buses and the magnitudes at PQ buses, in that order. Its pattern
    /// depends only on the topology. This is the matrix `solve` factors at each Newton iteration.
    ///
    /// Panics if `v` doesn't have one entry per connected bus, or if any has zero magnitude.
    pub fn jacobian(&self, v: &[Complex<f64>]) -> CscMatrix<f64> {
        assert_eq!(
            v.len(),
            self.topology.bus_rows.len(),
            "one voltage per connected bus"
        );
        let unknowns = Unknowns {
            pv: self.topology.pv_buses.clone(),
            pq: self.topology.pq_buses.clone(),
            pv_pq: self.topology.pv_pq(),
        };
        self.jacobian_over(v, &unknowns)
    }

    /// `newtonpf`'s Jacobian in CSC form, at the voltages `v`, from `dSbus_dV`:
    /// `dS/dVa = j·diag(V)·conj(diag(I) − Ybus·diag(V))` and
    /// `dS/dVm = diag(V)·conj(Ybus·diag(V/|V|)) + conj(diag(I))·diag(V/|V|)`, with `I = Ybus V`,
    /// stacked as `[[Re dS/dVa[pvpq, pvpq], Re dS/dVm[pvpq, pq]], [Im dS/dVa[pq, pvpq],
    /// Im dS/dVm[pq, pq]]]`. Ybus stores its whole diagonal, so every iteration's Jacobian has
    /// the same pattern.
    ///
    /// `v` is the current Newton iterate: the complex voltage phasor `Vm·e^{jVa}` of every
    /// connected bus, in per unit and in the model's internal bus order (entry `i` is case bus
    /// row `bus_rows()[i]`), so it has one entry per row of Ybus. Every bus is included, the
    /// reference buses too, since their voltages enter `I = Ybus V`; `unknowns` picks the rows
    /// and columns that form the Jacobian. No magnitude may be zero, since `V/|V|` divides by
    /// it; a power flow iterate never reaches zero in practice.
    fn jacobian_over(&self, v: &[Complex<f64>], unknowns: &Unknowns) -> CscMatrix<f64> {
        let shapes_match_message = "v and the factors have one entry per bus, matching Ybus";
        let current = self.ybus.mul_vec(v).expect(shapes_match_message);
        let v_unit: Vec<Complex<f64>> = v.iter().map(|value| value / value.norm()).collect();
        let j_v: Vec<Complex<f64>> = v.iter().map(|value| value * Complex::i()).collect();

        let ybus_v = self.ybus.scale_cols(v).expect(shapes_match_message);
        let ds_dva = (&CsrMatrix::from_diagonal(current.iter().copied()) - &ybus_v)
            .conj()
            .scale_rows(&j_v)
            .expect(shapes_match_message);
        let ds_dvm = &self
            .ybus
            .scale_cols(&v_unit)
            .expect(shapes_match_message)
            .conj()
            .scale_rows(v)
            .expect(shapes_match_message)
            + &CsrMatrix::from_diagonal(current.iter().zip(&v_unit).map(|(i, u)| i.conj() * u));

        let block = |matrix: &CsrMatrix<f64>, rows: &[usize], cols: &[usize]| {
            matrix.select(rows, cols).expect("bus indices are in range")
        };
        let (re_va, im_va) = (ds_dva.real(), ds_dva.imag());
        let (re_vm, im_vm) = (ds_dvm.real(), ds_dvm.imag());
        let (pvpq, pq) = (&unknowns.pv_pq, &unknowns.pq);
        let top = CsrMatrix::hstack([block(&re_va, pvpq, pvpq), block(&re_vm, pvpq, pq)]);
        let bottom = CsrMatrix::hstack([block(&im_va, pq, pvpq), block(&im_vm, pq, pq)]);
        CsrMatrix::vstack([
            top.expect(shapes_match_message),
            bottom.expect(shapes_match_message),
        ])
        .expect(shapes_match_message)
        .to_csc()
    }

    /// `pfsoln`'s results, indexed like the case's tables.
    fn results(
        &self,
        case: &Case,
        v: &[Complex<f64>],
        iterations: usize,
        mismatch_pu: f64,
        full_factorizations: usize,
    ) -> AcSolution {
        let topology = &self.topology;
        let base_mva = self.base_mva;
        let mut vm_pu: Vec<f64> = case.bus.iter().map(|bus| bus.vm_pu).collect();
        let mut va_deg: Vec<f64> = case.bus.iter().map(|bus| bus.va_deg).collect();
        for (&row, value) in topology.bus_rows.iter().zip(v) {
            vm_pu[row] = value.norm();
            va_deg[row] = value.arg() * 180.0 / PI;
        }

        // Generators: scheduled outputs for those in service at connected buses, zero for the
        //   rest, as runpf reports them
        let mut pg_mw = vec![0.0; case.generator.len()];
        let mut qg_mvar = vec![0.0; case.generator.len()];
        for &row in &topology.gen_rows {
            pg_mw[row] = case.generator[row].pg_mw;
            qg_mvar[row] = case.generator[row].qg_mvar;
        }
        // pfsoln's "on": in service at a bus whose case type isn't PQ
        let on: Vec<(usize, usize)> = topology
            .gen_rows
            .iter()
            .zip(&topology.gen_buses)
            .filter(|&(_, &bus)| {
                case.bus[topology.bus_rows[bus]].bus_type != crate::case::BusType::PQ
            })
            .map(|(&row, &bus)| (row, bus))
            .collect();
        let injection = |bus: usize| -> Complex<f64> {
            let row: Complex<f64> = self
                .ybus
                .outer(bus)
                .map(|(cols, values)| {
                    cols.iter()
                        .zip(values)
                        .fold(Complex::new(0.0, 0.0), |sum, (&j, &y)| sum + y * v[j])
                })
                .expect("bus indices are in range");
            v[bus] * row.conj()
        };
        let s_on: Vec<Complex<f64>> = on.iter().map(|&(_, bus)| injection(bus)).collect();
        for (&(row, bus), s) in on.iter().zip(&s_on) {
            let qd = case.bus[topology.bus_rows[bus]].qd_mvar;
            qg_mvar[row] = s.im * base_mva + qd;
        }
        if on.len() > 1 {
            share_reactive_power(case, &on, &mut qg_mvar);
        }
        // The slack: the first "on" generator at each reference bus takes up the difference
        for &reference in &topology.ref_buses {
            let at_bus: Vec<usize> = (0..on.len()).filter(|&k| on[k].1 == reference).collect();
            let first = *at_bus
                .first()
                .expect("a reference bus has an in-service generator");
            let pd = case.bus[topology.bus_rows[reference]].pd_mw;
            let others: f64 = at_bus[1..].iter().map(|&k| pg_mw[on[k].0]).sum();
            pg_mw[on[first].0] = s_on[first].re * base_mva + pd - others;
        }

        // Branch flows at both ends, zero for branches out of service or touching an isolated
        //   bus
        let nl = case.branch.len();
        let (mut pf_mw, mut qf_mvar, mut pt_mw, mut qt_mvar) =
            (vec![0.0; nl], vec![0.0; nl], vec![0.0; nl], vec![0.0; nl]);
        let i_from = self.yf.mul_vec(v).expect("one voltage per bus");
        let i_to = self.yt.mul_vec(v).expect("one voltage per bus");
        for (k, (&row, &(f, t))) in topology.branch_rows.iter().zip(&topology.ends).enumerate() {
            let sf = v[f] * i_from[k].conj() * base_mva;
            let st = v[t] * i_to[k].conj() * base_mva;
            (pf_mw[row], qf_mvar[row]) = (sf.re, sf.im);
            (pt_mw[row], qt_mvar[row]) = (st.re, st.im);
        }
        AcSolution {
            vm_pu,
            va_deg,
            pg_mw,
            qg_mvar,
            pf_mw,
            qf_mvar,
            pt_mw,
            qt_mvar,
            iterations,
            mismatch_pu,
            full_factorizations,
        }
    }
}

impl Unknowns {
    /// The bus (internal index) whose angle (the first `pv_pq.len()` columns) or magnitude (the
    /// rest) Jacobian column `column` stands for.
    fn column_bus(&self, column: usize) -> usize {
        match self.pv_pq.get(column) {
            Some(&bus) => bus,
            None => self.pq[column - self.pv_pq.len()],
        }
    }
}

/// `pfsoln`'s sharing of each bus's reactive output among its "on" generators (in service at a
/// REF or PV bus), when there is more than one such generator in the system. `qg_mvar` holds
/// each one's bus total on entry.
///
/// Each generator first gets an equal share. Infinite limits are replaced by a bound `M`: the
/// generator's share plus its finite limits' magnitudes, summed over the generators at its bus.
/// Then each generator gets its minimum plus the same fraction of its range, the fraction that
/// puts the bus's total where it should be. At a bus whose total range is zero (within
/// `10 ε`), the excess over the minimum is split equally.
fn share_reactive_power(case: &Case, on: &[(usize, usize)], qg_mvar: &mut [f64]) {
    let mut count_at: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
    for &(_, bus) in on {
        *count_at.entry(bus).or_default() += 1.0;
    }
    let share: Vec<f64> = on
        .iter()
        .map(|&(row, bus)| qg_mvar[row] / count_at[&bus])
        .collect();
    let finite_magnitude = |limit: f64| {
        if limit.is_infinite() {
            0.0
        } else {
            limit.abs()
        }
    };
    let mut bound_at: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
    for (k, &(row, bus)) in on.iter().enumerate() {
        let generator = &case.generator[row];
        *bound_at.entry(bus).or_default() += share[k].abs()
            + finite_magnitude(generator.qmax_mvar)
            + finite_magnitude(generator.qmin_mvar);
    }
    let limit = |value: f64, bus: usize| {
        if value == f64::INFINITY {
            bound_at[&bus]
        } else if value == f64::NEG_INFINITY {
            -bound_at[&bus]
        } else {
            value
        }
    };
    let (qmin, qmax): (Vec<f64>, Vec<f64>) = on
        .iter()
        .map(|&(row, bus)| {
            let generator = &case.generator[row];
            (
                limit(generator.qmin_mvar, bus),
                limit(generator.qmax_mvar, bus),
            )
        })
        .unzip();
    let mut total: std::collections::HashMap<usize, (f64, f64, f64)> =
        std::collections::HashMap::new();
    for (k, &(_, bus)) in on.iter().enumerate() {
        let entry = total.entry(bus).or_default();
        entry.0 += share[k];
        entry.1 += qmin[k];
        entry.2 += qmax[k];
    }
    for (k, &(row, bus)) in on.iter().enumerate() {
        let (q_total, q_min, q_max) = total[&bus];
        qg_mvar[row] = if (q_min - q_max).abs() < 10.0 * f64::EPSILON {
            qmin[k] + (q_total - q_min) / count_at[&bus]
        } else {
            qmin[k] + (q_total - q_min) / (q_max - q_min + f64::EPSILON) * (qmax[k] - qmin[k])
        };
    }
}

/// Runs an AC power flow on `case` with MATPOWER's default options (Newton's method):
/// [`AcModel::new`] followed by [`AcModel::solve`], with their errors.
pub fn ac_power_flow(case: &Case) -> Result<AcSolution, PowerflowError> {
    AcModel::new(case)?.solve(case, AcOptions::default())
}
