//! The AC optimal power flow, computed as MATPOWER's `runopf` computes it with MIPS, its
//! default solver.
//!
//! The model follows MATPOWER 8.1's `opf_setup` for the AC case with MIPS, after `ext2int`
//! (isolated buses, and generators and branches that are out of service or touch one, are left
//! out):
//!
//! - **Variables** `x = [Va, Vm, Pg, Qg, y]`: the connected buses' angles (rad) and voltage
//!   magnitudes (p.u.), the in-service generators' real and reactive outputs (p.u.), and one
//!   helper variable per piecewise linear cost. Every bus of type REF has its angle fixed at the
//!   case's by equal bounds (no `bustypes`), each magnitude is bounded by `Vmin` and `Vmax`, and
//!   each output by its limits.
//! - **Nonlinear equalities** `Pmis` and `Qmis`: the real and imaginary parts of the power
//!   balance `V ∘ conj(Ybus V) − Cg (Pg + jQg) + (Pd + jQd)/baseMVA = 0`, constant-power loads
//!   (MATPOWER's default, with no voltage-dependent loads).
//! - **Nonlinear inequalities** `Sf` and `St`: the squared apparent power at each end of the
//!   branches with `0 < |rate_a| < 1e10` (MATPOWER tests `rate_a ≠ 0`), at most
//!   `(rate_a/baseMVA)²` (`opf.flow_lim = 'S'`).
//! - **Linear constraints,** in `opf_setup`'s order: `vl`, `makeAvl`'s constant power factor of
//!   each dispatchable load (a generator with `Pmin < 0 = Pmax` and a nonzero reactive limit),
//!   `ang`, `makeAang`'s angle difference limits, and `ycon`, `makeAy`'s segment constraints.
//!   (`opf_setup`'s first linear constraints, for PQ capability curves, are refused instead.)
//! - **Objective:** the polynomial costs of `Pg` and `Qg` from [`OpfCosts`], constant terms
//!   included as `opt_model`'s `nlp_costfcn` includes them, plus `Σ y`.
//!
//! It is solved by [`mips`] from `nlpopf_solver`'s starting point, with MATPOWER's MIPS options
//! for the AC OPF (`feastol` = `opf.violation` = 5e-6, `cost_mult` = 1e-4, the rest MIPS's
//! defaults), or, if the options choose it, by [`mehrotra`] from the same point with the same
//! `feastol` and `cost_mult`, and with the constraints' derivatives from Technical Note 2's
//! formulas (see the
//! crate-private `opf::derivatives`). The results and multipliers are mapped back to the case
//! as `nlpopf_solver` and `opf` map them.

use interior_point::kkt::KktFactorization;
use interior_point::mehrotra::{mehrotra, MehrotraOptions};
use interior_point::mips::{mips, MipsOptions};
use interior_point::problem::{Constraints, NonlinearConstraints, NonlinearMultipliers, Problem};
use interior_point::solution::{kkt_residuals, KktResiduals, Status};
use sparse_mat::compressed_matrix::CsrMatrix;
use sparse_mat::Complex;

use crate::ac::{check_impedance, Admittances};
use crate::case::Case;
use crate::opf::cost::{OpfCosts, QuadraticCosts};
use crate::opf::derivatives::{
    flows, squared_magnitude_hessian, squared_magnitudes, weighted_hessian, Flows,
};
use crate::opf::{
    helper_start, per_degree, reference_buses, starting_point, LinearRows, OpfSolver, SolverResult,
};
use crate::powerflow_error::PowerflowError;
use crate::topology::Network;

type C64 = Complex<f64>;

/// The AC optimal power flow's options.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct AcOpfOptions {
    /// The solver and its options. The default is MATPOWER's for the AC OPF: MIPS, with its
    /// defaults but `feastol` set to `opf.violation`, 5e-6, and `cost_mult` 1e-4.
    pub solver: OpfSolver,
}

impl Default for AcOpfOptions {
    fn default() -> Self {
        let mut mips = MipsOptions::default();
        mips.feastol = 5e-6;
        mips.cost_mult = 1e-4;
        AcOpfOptions {
            solver: OpfSolver::Mips(mips),
        }
    }
}

impl AcOpfOptions {
    /// The predictor-corrector method with MATPOWER's AC OPF settings: its defaults, but
    /// `feastol` 5e-6 and `cost_mult` 1e-4, as MIPS runs by default, and the KKT matrix factored
    /// by the pivoted LDLᵀ (`KktFactorization::Ldlt`, without the inertia correction). On the
    /// survey of every MATPOWER case with costs it met phase O3's bar wherever LU did, with the
    /// same outcomes and iterations, 1.12× faster in the median and 2.70× in total
    /// (`docs/benchmarks/ldlt-opf-2026-11-08.md`), and more since the LDLᵀ's bookkeeping and
    /// amalgamation (`docs/benchmarks/ldlt-amalgamation-2026-11-08.md`).
    pub fn mehrotra() -> Self {
        let mut mehrotra = MehrotraOptions::default();
        mehrotra.feastol = 5e-6;
        mehrotra.cost_mult = 1e-4;
        mehrotra.kkt_factorization = KktFactorization::Ldlt;
        AcOpfOptions {
            solver: OpfSolver::Mehrotra(mehrotra),
        }
    }
}

/// An AC optimal power flow's results, indexed like the case's own tables, as `runopf` reports
/// them.
///
/// Out-of-service generators and branches, and those at an isolated bus, report 0.0 for every
/// quantity except `vg_pu`, which keeps the case's, as MATPOWER's `opf` sets them. An isolated
/// bus keeps its voltage from the case and has no prices. Results are reported whether or not
/// the solver converged: see [`AcOpfSolution::status`].
#[derive(Clone, Debug, PartialEq)]
pub struct AcOpfSolution {
    /// Each bus's voltage magnitude (p.u.)
    pub vm_pu: Vec<f64>,
    /// Each bus's voltage angle (degrees)
    pub va_deg: Vec<f64>,
    /// Each bus's real power price, `LAM_P` ($/MWh)
    pub lam_p: Vec<f64>,
    /// Each bus's reactive power price, `LAM_Q` ($/MVAr-hr)
    pub lam_q: Vec<f64>,
    /// The multiplier of each bus's lower voltage magnitude limit, `MU_VMIN` ($/p.u.-hr)
    pub mu_vmin: Vec<f64>,
    /// The multiplier of each bus's upper voltage magnitude limit, `MU_VMAX` ($/p.u.-hr)
    pub mu_vmax: Vec<f64>,
    /// Each generator's real power output (MW)
    pub pg_mw: Vec<f64>,
    /// Each generator's reactive power output (MVAr)
    pub qg_mvar: Vec<f64>,
    /// Each generator's voltage setpoint: its bus's voltage magnitude (p.u.)
    pub vg_pu: Vec<f64>,
    /// The multipliers of each generator's lower and upper real power limits ($/MWh)
    pub mu_pmin: Vec<f64>,
    pub mu_pmax: Vec<f64>,
    /// The multipliers of each generator's lower and upper reactive power limits ($/MVAr-hr)
    pub mu_qmin: Vec<f64>,
    pub mu_qmax: Vec<f64>,
    /// Each branch's real and reactive power flowing in at the from end (MW, MVAr)
    pub pf_mw: Vec<f64>,
    pub qf_mvar: Vec<f64>,
    /// Each branch's real and reactive power flowing in at the to end (MW, MVAr)
    pub pt_mw: Vec<f64>,
    pub qt_mvar: Vec<f64>,
    /// The multipliers of each branch's apparent power limit at the from and to ends, `MU_SF`
    /// and `MU_ST` ($/MVA-hr)
    pub mu_sf: Vec<f64>,
    pub mu_st: Vec<f64>,
    /// The multipliers of each branch's angle difference limits ($/degree)
    pub mu_angmin: Vec<f64>,
    pub mu_angmax: Vec<f64>,
    /// The total cost ($/hr)
    pub objective: f64,
    /// How the solver ended
    pub status: Status,
    /// The solver's iterations
    pub iterations: usize,
    /// The KKT residuals at the solver's point and multipliers, on the model's own problem
    /// (per unit, and $/hr for the objective), by
    /// [`kkt_residuals`]: one measure for either
    /// solver
    pub residuals: KktResiduals,
    /// The solver's own result, on the model's variables `[Va, Vm, Pg, Qg, y]` and
    /// constraints, with each iteration's record
    pub solver: SolverResult,
}

impl AcOpfSolution {
    /// Whether the solver converged, MATPOWER's `success`.
    pub fn converged(&self) -> bool {
        self.status == Status::Converged
    }
}

/// The AC OPF's objective and nonlinear constraints over `x = [Va, Vm, Pg, Qg, y]`, for MIPS.
struct AcOpfProblem {
    nb: usize,
    ng: usize,
    ny: usize,
    ybus: CsrMatrix<C64>,
    identity: CsrMatrix<C64>,
    /// The limited branches' from-end and to-end incidence (`Cf`, `Ct`) and admittance rows
    cf: CsrMatrix<C64>,
    yf: CsrMatrix<C64>,
    ct: CsrMatrix<C64>,
    yt: CsrMatrix<C64>,
    /// Each limited branch's `(rate_a/baseMVA)²`
    flow_max: Vec<f64>,
    /// Each connected bus's constant-power load (p.u.)
    load: Vec<C64>,
    /// Each in-service generator's bus
    gen_buses: Vec<usize>,
    active: QuadraticCosts,
    reactive: Option<QuadraticCosts>,
}

impl AcOpfProblem {
    fn nx(&self) -> usize {
        2 * self.nb + 2 * self.ng + self.ny
    }

    fn voltages(&self, x: &[f64]) -> Vec<C64> {
        let nb = self.nb;
        (0..nb).map(|i| C64::from_polar(x[nb + i], x[i])).collect()
    }

    /// The limited branches' flows at the from and to ends, if any branch is limited.
    fn branch_flows(&self, v: &[C64]) -> Option<(Flows, Flows)> {
        (!self.flow_max.is_empty())
            .then(|| (flows(&self.cf, &self.yf, v), flows(&self.ct, &self.yt, v)))
    }

    /// Places `block` (rows × 2nb, over `[Va, Vm]`) at row offset `row` of a matrix with `rows`
    /// rows and `nx` columns, adding `extra` entries.
    fn widen(
        &self,
        rows: usize,
        blocks: &[(&CsrMatrix<f64>, usize, usize)],
        extra: impl IntoIterator<Item = (usize, usize, f64)>,
    ) -> CsrMatrix<f64> {
        let mut triplets: Vec<(usize, usize, f64)> = extra.into_iter().collect();
        for &(block, row, col) in blocks {
            triplets.extend(block.iter().map(|(i, j, v)| (row + i, col + j, v)));
        }
        CsrMatrix::from_triplets(rows, self.nx(), triplets).expect("indices in range")
    }

    /// `Σ ½ Q p² + c p + k` over `costs` at `p`, its gradient into `gradient` from `offset`, as
    /// `opt_model` evaluates a quadratic cost.
    fn cost(costs: &QuadraticCosts, p: &[f64], gradient: &mut [f64], offset: usize) -> f64 {
        let mut f = 0.0;
        for (k, &p) in p.iter().enumerate() {
            let (q, c, constant) = (costs.quadratic[k], costs.linear[k], costs.constant[k]);
            f += 0.5 * q * p * p + c * p + constant;
            gradient[offset + k] = q * p + c;
        }
        f
    }

    /// The costs' Hessian: `Q` on the outputs' diagonal.
    fn cost_hessian(&self) -> Vec<(usize, usize, f64)> {
        let (pg, qg) = (2 * self.nb, 2 * self.nb + self.ng);
        let mut entries: Vec<(usize, usize, f64)> = self
            .active
            .quadratic
            .iter()
            .enumerate()
            .filter(|&(_, &q)| q != 0.0)
            .map(|(k, &q)| (pg + k, pg + k, q))
            .collect();
        if let Some(reactive) = &self.reactive {
            entries.extend(
                reactive
                    .quadratic
                    .iter()
                    .enumerate()
                    .filter(|&(_, &q)| q != 0.0)
                    .map(|(k, &q)| (qg + k, qg + k, q)),
            );
        }
        entries
    }
}

impl Problem for AcOpfProblem {
    fn dimension(&self) -> usize {
        self.nx()
    }

    fn objective(&self, x: &[f64]) -> (f64, Vec<f64>) {
        let (nb, ng) = (self.nb, self.ng);
        let (pg, qg, y) = (2 * nb, 2 * nb + ng, 2 * nb + 2 * ng);
        let mut gradient = vec![0.0; self.nx()];
        let mut f = Self::cost(&self.active, &x[pg..qg], &mut gradient, pg);
        if let Some(reactive) = &self.reactive {
            f += Self::cost(reactive, &x[qg..y], &mut gradient, qg);
        }
        f += x[y..].iter().sum::<f64>();
        gradient[y..].fill(1.0);
        (f, gradient)
    }

    fn objective_hessian(&self, _: &[f64]) -> CsrMatrix<f64> {
        CsrMatrix::from_triplets(self.nx(), self.nx(), self.cost_hessian())
            .expect("indices in range")
    }

    fn has_nonlinear_constraints(&self) -> bool {
        true
    }

    fn constraints(&self, x: &[f64]) -> NonlinearConstraints {
        let (nb, ng) = (self.nb, self.ng);
        let (pg, qg) = (2 * nb, 2 * nb + ng);
        let v = self.voltages(x);

        // The power balance and its Jacobian
        let bus = flows(&self.identity, &self.ybus, &v);
        let mut mismatch: Vec<C64> = bus.s.iter().zip(&self.load).map(|(s, d)| s + d).collect();
        for (k, &b) in self.gen_buses.iter().enumerate() {
            mismatch[b] -= C64::new(x[pg + k], x[qg + k]);
        }
        let g = mismatch
            .iter()
            .map(|m| m.re)
            .chain(mismatch.iter().map(|m| m.im))
            .collect();
        let (dp_dva, dq_dva) = (bus.ds_dva.real(), bus.ds_dva.imag());
        let (dp_dvm, dq_dvm) = (bus.ds_dvm.real(), bus.ds_dvm.imag());
        let outputs = self
            .gen_buses
            .iter()
            .enumerate()
            .flat_map(|(k, &b)| [(b, pg + k, -1.0), (nb + b, qg + k, -1.0)]);
        let dg = self.widen(
            2 * nb,
            &[
                (&dp_dva, 0, 0),
                (&dp_dvm, 0, nb),
                (&dq_dva, nb, 0),
                (&dq_dvm, nb, nb),
            ],
            outputs,
        );

        // The squared apparent power limits and their Jacobian
        let nl2 = self.flow_max.len();
        let (h, dh) = match self.branch_flows(&v) {
            Some((from, to)) => {
                let (af, daf) = squared_magnitudes(&from);
                let (at, dat) = squared_magnitudes(&to);
                let h = af
                    .iter()
                    .chain(&at)
                    .zip(self.flow_max.iter().chain(&self.flow_max))
                    .map(|(a, max)| a - max)
                    .collect();
                let dh = self.widen(2 * nl2, &[(&daf, 0, 0), (&dat, nl2, 0)], []);
                (h, dh)
            }
            None => (vec![], self.widen(0, &[], [])),
        };
        NonlinearConstraints { h, g, dh, dg }
    }

    fn lagrangian_hessian(
        &self,
        x: &[f64],
        multipliers: NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        let nb = self.nb;
        let v = self.voltages(x);
        let (lam_p, lam_q) = multipliers.eq.split_at(nb);
        let mut voltage = weighted_hessian(&self.identity, &self.ybus, &v, lam_p, lam_q);
        if let Some((from, to)) = self.branch_flows(&v) {
            let (mu_f, mu_t) = multipliers.ineq.split_at(self.flow_max.len());
            let hf = squared_magnitude_hessian(&self.cf, &self.yf, &v, &from, mu_f);
            let ht = squared_magnitude_hessian(&self.ct, &self.yt, &v, &to, mu_t);
            voltage = voltage
                .try_add(&hf)
                .and_then(|h| h.try_add(&ht))
                .expect("all 2nb × 2nb");
        }
        let costs = self
            .cost_hessian()
            .into_iter()
            .map(|(i, j, q)| (i, j, cost_mult * q));
        self.widen(self.nx(), &[(&voltage, 0, 0)], costs)
    }
}

/// Solves `case`'s AC optimal power flow with MATPOWER's defaults; see [`ac_opf_with_options`].
pub fn ac_opf(case: &Case) -> Result<AcOpfSolution, PowerflowError> {
    ac_opf_with_options(case, &AcOpfOptions::default())
}

/// Solves `case`'s AC optimal power flow as MATPOWER's `runopf` does with MIPS, or with the
/// solver `options` choose (see the [module documentation](self)). A solve that the solver
/// doesn't finish is still returned, with its [`AcOpfSolution::status`].
///
/// # Errors
///
/// - the errors of [`Case::validate`]
/// - [`PowerflowError::ZeroImpedance`] for an in-service branch with `r = x = 0`
/// - [`PowerflowError::NoReferenceBus`] if no connected bus has type REF, where MATPOWER fails
///   too
/// - [`PowerflowError::MissingField`] if the case has no costs, and
///   [`PowerflowError::UnsupportedCost`] for costs [`OpfCosts`] refuses
/// - [`PowerflowError::UnsupportedGenerator`] for an in-service generator with a PQ capability
///   curve, or a dispatchable load with both reactive limits nonzero or a reactive output off
///   its power factor by more than 1e-6 p.u., which `makeAvl` rejects
/// - [`PowerflowError::Solver`] if the solver refuses `options`, such as MIPS's `xi` outside
///   `[0.5, 1)`
pub fn ac_opf_with_options(
    case: &Case,
    options: &AcOpfOptions,
) -> Result<AcOpfSolution, PowerflowError> {
    let network = Network::new(case, check_impedance)?;
    let nb = network.bus_rows.len();
    let ng = network.gen_rows.len();
    let base = case.base_mva;
    let refs = reference_buses(case, &network)?;
    let costs = OpfCosts::new(case, &network.gen_rows, true)?;
    let ny = costs.piecewise.len();
    for &row in &network.gen_rows {
        check_capability_curve(case, row)?;
    }
    let admittances =
        Admittances::new(case, &network.bus_rows, &network.branch_rows, &network.ends);
    let (va, vm, pg, qg, y) = (0, nb, 2 * nb, 2 * nb + ng, 2 * nb + 2 * ng);
    let nx = 2 * nb + 2 * ng + ny;

    // The variables' bounds
    let mut xmin = vec![f64::NEG_INFINITY; nx];
    let mut xmax = vec![f64::INFINITY; nx];
    for &i in &refs {
        let angle = case.bus[network.bus_rows[i]].va_deg.to_radians();
        (xmin[va + i], xmax[va + i]) = (angle, angle);
    }
    for (i, &row) in network.bus_rows.iter().enumerate() {
        xmin[vm + i] = case.bus[row].vmin_pu;
        xmax[vm + i] = case.bus[row].vmax_pu;
    }
    for (k, &row) in network.gen_rows.iter().enumerate() {
        let generator = &case.generator[row];
        (xmin[pg + k], xmax[pg + k]) = (generator.pmin_mw / base, generator.pmax_mw / base);
        (xmin[qg + k], xmax[qg + k]) = (generator.qmin_mvar / base, generator.qmax_mvar / base);
    }

    // The linear constraints, block by block: vl, ang, ycon
    let mut rows = LinearRows::default();
    for (k, &row) in network.gen_rows.iter().enumerate() {
        if let Some((p, q)) = power_factor(case, row)? {
            rows.push([(pg + k, p), (qg + k, q)], 0.0, 0.0);
        }
    }
    let first_angle_row = rows.len();
    let angle_limited = rows.angle_limits(case, &network, va);
    rows.piecewise(&costs, pg, Some(qg), y);
    let (a, l, u) = rows.finish(nx);

    // The limited branches' rows of the admittances, with their incidence
    let limited: Vec<usize> = (0..network.branch_rows.len())
        .filter(|&k| {
            let rate = case.branch[network.branch_rows[k]].rate_a_mva;
            rate != 0.0 && rate < 1e10
        })
        .collect();
    let incidence = |end: fn(&(usize, usize)) -> usize| {
        CsrMatrix::from_triplets(
            limited.len(),
            nb,
            limited
                .iter()
                .enumerate()
                .map(|(n, &k)| (n, end(&network.ends[k]), C64::new(1.0, 0.0))),
        )
        .expect("internal indices in range")
    };
    let all_buses: Vec<usize> = (0..nb).collect();
    let select = |m: &CsrMatrix<C64>| m.select(&limited, &all_buses).expect("rows in range");
    let problem = AcOpfProblem {
        nb,
        ng,
        ny,
        identity: CsrMatrix::identity(nb),
        cf: incidence(|&(f, _)| f),
        yf: select(&admittances.yf),
        ct: incidence(|&(_, t)| t),
        yt: select(&admittances.yt),
        flow_max: limited
            .iter()
            .map(|&k| (case.branch[network.branch_rows[k]].rate_a_mva / base).powi(2))
            .collect(),
        load: network
            .bus_rows
            .iter()
            .map(|&row| C64::new(case.bus[row].pd_mw, case.bus[row].qd_mvar) / base)
            .collect(),
        gen_buses: network.gen_buses.clone(),
        active: costs.active.clone(),
        reactive: costs.reactive.clone(),
        ybus: admittances.ybus.clone(),
    };

    // nlpopf_solver's starting point: every angle at the first reference bus's, and each
    //   magnitude midway between its limits, clipped to [0.5, 1.5]
    let mut x0 = starting_point(&xmin, &xmax);
    x0[va..vm].fill(xmin[va + refs[0]]);
    for (i, &row) in network.bus_rows.iter().enumerate() {
        let bus = &case.bus[row];
        x0[vm + i] = (bus.vmax_pu.min(1.5) + bus.vmin_pu.max(0.5)) / 2.0;
    }
    if let Some(start) = helper_start(&costs) {
        x0[y..].fill(start);
    }

    let constraints = Constraints {
        a: Some(a),
        l: Some(l),
        u: Some(u),
        xmin: Some(xmin),
        xmax: Some(xmax),
    };
    let refused = |error| PowerflowError::Solver { error };
    let solution = match &options.solver {
        OpfSolver::Mips(options) => {
            SolverResult::Mips(mips(&problem, &x0, &constraints, options).map_err(refused)?)
        }
        OpfSolver::Mehrotra(options) => {
            SolverResult::Mehrotra(mehrotra(&problem, &x0, &constraints, options).map_err(refused)?)
        }
    };
    let residuals = kkt_residuals(&problem, &constraints, solution.x(), solution.multipliers())
        .expect("the solver's multipliers fit its problem");

    // The results, back on the case's tables
    let x = solution.x();
    let multipliers = solution.multipliers();
    let v = problem.voltages(x);
    let nbus = case.bus.len();
    let mut vm_pu: Vec<f64> = case.bus.iter().map(|bus| bus.vm_pu).collect();
    let mut va_deg: Vec<f64> = case.bus.iter().map(|bus| bus.va_deg).collect();
    let (mut lam_p, mut lam_q) = (vec![0.0; nbus], vec![0.0; nbus]);
    let (mut mu_vmin, mut mu_vmax) = (vec![0.0; nbus], vec![0.0; nbus]);
    for (i, &row) in network.bus_rows.iter().enumerate() {
        vm_pu[row] = x[vm + i];
        va_deg[row] = x[va + i].to_degrees();
        lam_p[row] = multipliers.eq_nonlin[i] / base;
        lam_q[row] = multipliers.eq_nonlin[nb + i] / base;
        mu_vmin[row] = multipliers.lower[vm + i];
        mu_vmax[row] = multipliers.upper[vm + i];
    }
    let ngen = case.generator.len();
    let zeros = || vec![0.0; ngen];
    let (mut pg_mw, mut qg_mvar) = (zeros(), zeros());
    let (mut mu_pmin, mut mu_pmax, mut mu_qmin, mut mu_qmax) = (zeros(), zeros(), zeros(), zeros());
    let mut vg_pu: Vec<f64> = case.generator.iter().map(|g| g.vg_pu).collect();
    for (k, (&row, &bus)) in network.gen_rows.iter().zip(&network.gen_buses).enumerate() {
        pg_mw[row] = x[pg + k] * base;
        qg_mvar[row] = x[qg + k] * base;
        vg_pu[row] = x[vm + bus];
        mu_pmin[row] = multipliers.lower[pg + k] / base;
        mu_pmax[row] = multipliers.upper[pg + k] / base;
        mu_qmin[row] = multipliers.lower[qg + k] / base;
        mu_qmax[row] = multipliers.upper[qg + k] / base;
    }
    let nl = case.branch.len();
    let zeros = || vec![0.0; nl];
    let (mut pf_mw, mut qf_mvar, mut pt_mw, mut qt_mvar) = (zeros(), zeros(), zeros(), zeros());
    let (mut mu_sf, mut mu_st, mut mu_angmin, mut mu_angmax) = (zeros(), zeros(), zeros(), zeros());
    let i_from = admittances.yf.mul_vec(&v).expect("Yf has a column per bus");
    let i_to = admittances.yt.mul_vec(&v).expect("Yt has a column per bus");
    for (k, (&row, &(f, t))) in network.branch_rows.iter().zip(&network.ends).enumerate() {
        let sf = v[f] * i_from[k].conj() * base;
        let st = v[t] * i_to[k].conj() * base;
        (pf_mw[row], qf_mvar[row]) = (sf.re, sf.im);
        (pt_mw[row], qt_mvar[row]) = (st.re, st.im);
    }
    // MATPOWER's MU_SF for squared apparent power limits: 2 μ rate_a/baseMVA, per MVA
    let nl2 = limited.len();
    for (n, &k) in limited.iter().enumerate() {
        let row = network.branch_rows[k];
        let rate = case.branch[row].rate_a_mva / base;
        mu_sf[row] = 2.0 * multipliers.ineq_nonlin[n] * rate / base;
        mu_st[row] = 2.0 * multipliers.ineq_nonlin[nl2 + n] * rate / base;
    }
    for (n, &k) in angle_limited.iter().enumerate() {
        let row = network.branch_rows[k];
        mu_angmin[row] = per_degree(multipliers.mu_l[first_angle_row + n]);
        mu_angmax[row] = per_degree(multipliers.mu_u[first_angle_row + n]);
    }

    Ok(AcOpfSolution {
        vm_pu,
        va_deg,
        lam_p,
        lam_q,
        mu_vmin,
        mu_vmax,
        pg_mw,
        qg_mvar,
        vg_pu,
        mu_pmin,
        mu_pmax,
        mu_qmin,
        mu_qmax,
        pf_mw,
        qf_mvar,
        pt_mw,
        qt_mvar,
        mu_sf,
        mu_st,
        mu_angmin,
        mu_angmax,
        objective: solution.f(),
        status: solution.status(),
        iterations: solution.iterations(),
        residuals,
        solver: solution,
    })
}

/// `PowerflowError::UnsupportedGenerator` if the generator at case row `row` has a PQ
/// capability curve, as MATPOWER's `hasPQcap` finds one: `Pc1` or `Pc2` nonzero, unless its
/// reactive limits and the curve's are all equal (a fixed reactive output).
fn check_capability_curve(case: &Case, row: usize) -> Result<(), PowerflowError> {
    let g = &case.generator[row];
    let fixed = [
        g.qmax_mvar,
        g.qc1max_mvar,
        g.qc1min_mvar,
        g.qc2max_mvar,
        g.qc2min_mvar,
    ]
    .iter()
    .all(|&q| q == g.qmin_mvar);
    if (g.pc1_mw != 0.0 || g.pc2_mw != 0.0) && !fixed {
        return Err(PowerflowError::UnsupportedGenerator {
            row,
            reason: "PQ capability curves aren't supported",
        });
    }
    Ok(())
}

/// `makeAvl`'s constant power factor for the generator at case row `row`, if it's a dispatchable
/// load (`Pmin < 0 = Pmax`) with a nonzero reactive limit: the coefficients `(sin θ, −cos θ)` of
/// `Pg` and `Qg` in the row `sin θ · Pg − cos θ · Qg = 0`, with `θ = atan2(Qlim, Pmin)` and
/// `Qlim` its nonzero reactive limit.
fn power_factor(case: &Case, row: usize) -> Result<Option<(f64, f64)>, PowerflowError> {
    let g = &case.generator[row];
    let base = case.base_mva;
    let is_load = g.pmin_mw < 0.0 && g.pmax_mw == 0.0;
    if !is_load || (g.qmin_mvar == 0.0 && g.qmax_mvar == 0.0) {
        return Ok(None);
    }
    let unsupported = |reason| PowerflowError::UnsupportedGenerator { row, reason };
    if g.qmin_mvar != 0.0 && g.qmax_mvar != 0.0 {
        return Err(unsupported(
            "a dispatchable load needs Qmin or Qmax equal to zero",
        ));
    }
    let (pmin, qg, pg) = (g.pmin_mw / base, g.qg_mvar / base, g.pg_mw / base);
    let qlim = if g.qmin_mvar == 0.0 {
        g.qmax_mvar
    } else {
        g.qmin_mvar
    } / base;
    if (qg - pg * qlim / pmin).abs() > 1e-6 {
        return Err(unsupported(
            "a dispatchable load's reactive output must match its power factor",
        ));
    }
    let theta = qlim.atan2(pmin);
    Ok(Some((theta.sin(), -theta.cos())))
}
