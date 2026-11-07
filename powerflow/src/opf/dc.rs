//! The DC optimal power flow, computed as MATPOWER's `rundcopf` computes it with MIPS.
//!
//! The model follows MATPOWER 8.1's `opf_setup` for the DC case, after `ext2int` (isolated
//! buses, and generators and branches that are out of service or touch one, are left out):
//!
//! - **Variables** `x = [Va, Pg, y]`: the connected buses' angles (rad), the in-service
//!   generators' outputs (p.u.), and one helper variable per piecewise linear cost. Every bus of
//!   type REF has its angle fixed at the case's by equal bounds (with no `bustypes`: a REF bus
//!   needs no generator here), and each output is bounded by `Pmin` and `Pmax`.
//! - **Linear constraints,** in `opf_setup`'s order:
//!   - `Pmis`: the power balance `B Va − Cg Pg = −(Pd + Gs)/baseMVA − Pbusinj`, one row per bus
//!   - `Pf`: the flow limits `−(rate_a/baseMVA + Pfinj) ≤ Bf Va ≤ rate_a/baseMVA − Pfinj`, on the
//!     branches with `0 < |rate_a| < 1e10` (MATPOWER tests `rate_a ≠ 0`)
//!   - `ang`: `makeAang`'s angle difference limits
//!   - `ycon`: `makeAy`'s segment constraints, `mⱼ Pg − y ≤ bⱼ`
//! - **Objective:** `½ Pgᵀ diag(Q) Pg + cᵀ Pg + Σ y`, from [`OpfCosts`]; the constant terms are
//!   added to the objective afterwards, as MATPOWER's `opt_model.solve` does, so MIPS never sees
//!   them.
//!
//! It is solved by [`qps_mips`] from `dcopf_solver`'s starting point, with MATPOWER's MIPS
//! options for the OPF (`feastol` = `opf.violation` = 5e-6, the rest MIPS's defaults), or, if
//! the options choose it, by [`qps_mehrotra`] from the same point with the same `feastol`. The
//! results and multipliers are mapped back to the case as `dcopf_solver` and `opf` map them.

use std::f64::consts::PI;

use interior_point::mehrotra::{qps_mehrotra, MehrotraOptions};
use interior_point::mips::{qps_mips, MipsOptions};
use interior_point::problem::{Constraints, QuadraticProblem};
use interior_point::solution::{kkt_residuals, KktResiduals, Status};
use sparse_mat::compressed_matrix::CsrMatrix;

use crate::case::Case;
use crate::dc::Bdc;
use crate::opf::cost::OpfCosts;
use crate::opf::{
    helper_start, per_degree, reference_buses, starting_point, LinearRows, OpfSolver, SolverResult,
};
use crate::powerflow_error::PowerflowError;
use crate::topology::Network;

/// The DC optimal power flow's options.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct DcOpfOptions {
    /// The solver and its options. The default is MATPOWER's for the OPF: MIPS, with its
    /// defaults but `feastol` set to `opf.violation`, 5e-6.
    pub solver: OpfSolver,
}

impl Default for DcOpfOptions {
    fn default() -> Self {
        let mut mips = MipsOptions::default();
        mips.feastol = 5e-6;
        DcOpfOptions {
            solver: OpfSolver::Mips(mips),
        }
    }
}

impl DcOpfOptions {
    /// The predictor-corrector method with MATPOWER's OPF settings: its defaults, but `feastol`
    /// 5e-6, as MIPS runs by default.
    pub fn mehrotra() -> Self {
        let mut mehrotra = MehrotraOptions::default();
        mehrotra.feastol = 5e-6;
        DcOpfOptions {
            solver: OpfSolver::Mehrotra(mehrotra),
        }
    }
}

/// A DC optimal power flow's results, indexed like the case's own tables, as `rundcopf`
/// reports them.
///
/// Out-of-service generators and branches, and those at an isolated bus, report 0.0 for every
/// quantity, as MATPOWER's `opf` sets them. An isolated bus keeps its angle from the case and
/// has no prices. Results are reported whether or not the solver converged: see
/// [`DcOpfSolution::status`].
#[derive(Clone, Debug, PartialEq)]
pub struct DcOpfSolution {
    /// Each bus's voltage angle (degrees)
    pub va_deg: Vec<f64>,
    /// Each bus's real power price, `LAM_P` ($/MWh): the power balance's multiplier
    pub lam_p: Vec<f64>,
    /// Each generator's real power output (MW)
    pub pg_mw: Vec<f64>,
    /// The multiplier of each generator's lower real power limit, `MU_PMIN` ($/MWh)
    pub mu_pmin: Vec<f64>,
    /// The multiplier of each generator's upper real power limit, `MU_PMAX` ($/MWh)
    pub mu_pmax: Vec<f64>,
    /// Each branch's real power flow at the from end (MW); the to end's is its negation
    pub pf_mw: Vec<f64>,
    /// The multiplier of each branch's flow limit in the from-to direction, `MU_SF` ($/MWh)
    pub mu_sf: Vec<f64>,
    /// The multiplier of each branch's flow limit in the to-from direction, `MU_ST` ($/MWh)
    pub mu_st: Vec<f64>,
    /// The multiplier of each branch's lower angle difference limit, `MU_ANGMIN` ($/degree)
    pub mu_angmin: Vec<f64>,
    /// The multiplier of each branch's upper angle difference limit, `MU_ANGMAX` ($/degree)
    pub mu_angmax: Vec<f64>,
    /// The total cost ($/hr), including the costs' constant terms
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
    /// The solver's own result, on the model's variables `[Va, Pg, y]` and constraints, with
    /// each iteration's record
    pub solver: SolverResult,
}

impl DcOpfSolution {
    /// Whether the solver converged, MATPOWER's `success`.
    pub fn converged(&self) -> bool {
        self.status == Status::Converged
    }
}

/// Solves `case`'s DC optimal power flow with MATPOWER's defaults; see [`dc_opf_with_options`].
pub fn dc_opf(case: &Case) -> Result<DcOpfSolution, PowerflowError> {
    dc_opf_with_options(case, &DcOpfOptions::default())
}

/// Solves `case`'s DC optimal power flow as MATPOWER's `rundcopf` does with MIPS, or with the
/// solver `options` choose (see the [module documentation](self)). A solve that the solver
/// doesn't finish is still returned, with its [`DcOpfSolution::status`].
///
/// # Errors
///
/// - the errors of [`Case::validate`]
/// - [`PowerflowError::ZeroReactance`] for an in-service branch with `x == 0`
/// - [`PowerflowError::NoReferenceBus`] if no connected bus has type REF, where MATPOWER fails
///   too (`dcopf_solver` starts every angle at the first reference bus's)
/// - [`PowerflowError::MissingField`] if the case has no costs, and
///   [`PowerflowError::UnsupportedCost`] for costs [`OpfCosts`] refuses
pub fn dc_opf_with_options(
    case: &Case,
    options: &DcOpfOptions,
) -> Result<DcOpfSolution, PowerflowError> {
    let network = Network::new(case, |row, branch| {
        if branch.x == 0.0 {
            Err(PowerflowError::ZeroReactance { row })
        } else {
            Ok(())
        }
    })?;
    let nb = network.bus_rows.len();
    let ng = network.gen_rows.len();
    let base = case.base_mva;
    let refs = reference_buses(case, &network)?;
    let costs = OpfCosts::new(case, &network.gen_rows, false)?;
    let ny = costs.piecewise.len();
    let Bdc {
        bbus,
        bf,
        pbusinj_pu,
        pfinj_pu,
    } = Bdc::new(case, nb, &network.branch_rows, &network.ends);
    let (va, pg, y) = (0, nb, nb + ng);
    let nx = nb + ng + ny;

    // The variables' bounds
    let mut xmin = vec![f64::NEG_INFINITY; nx];
    let mut xmax = vec![f64::INFINITY; nx];
    for &i in &refs {
        let angle = case.bus[network.bus_rows[i]].va_deg * (PI / 180.0);
        (xmin[va + i], xmax[va + i]) = (angle, angle);
    }
    for (k, &row) in network.gen_rows.iter().enumerate() {
        let generator = &case.generator[row];
        xmin[pg + k] = generator.pmin_mw / base;
        xmax[pg + k] = generator.pmax_mw / base;
    }

    // The linear constraints, block by block: Pmis, Pf, ang, ycon
    let mut rows = LinearRows::default();
    let mut generators_at: Vec<Vec<usize>> = vec![vec![]; nb];
    for (k, &bus) in network.gen_buses.iter().enumerate() {
        generators_at[bus].push(k);
    }
    for (i, &row) in network.bus_rows.iter().enumerate() {
        let bus = &case.bus[row];
        let mismatch = -(bus.pd_mw + bus.gs_mw) / base - pbusinj_pu[i];
        let (cols, values) = bbus.outer(i).expect("an internal bus");
        let angles = cols.iter().zip(values).map(|(&j, &b)| (va + j, b));
        let outputs = generators_at[i].iter().map(|&k| (pg + k, -1.0));
        rows.push(angles.chain(outputs), mismatch, mismatch);
    }
    let limited: Vec<usize> = (0..network.branch_rows.len())
        .filter(|&k| {
            let rate = case.branch[network.branch_rows[k]].rate_a_mva;
            rate != 0.0 && rate < 1e10
        })
        .collect();
    for &k in &limited {
        let (cols, values) = bf.outer(k).expect("an internal branch");
        let rate = case.branch[network.branch_rows[k]].rate_a_mva / base;
        rows.push(
            cols.iter().zip(values).map(|(&j, &b)| (va + j, b)),
            -(rate + pfinj_pu[k]),
            rate - pfinj_pu[k],
        );
    }
    let first_angle_row = rows.len();
    let angle_limited = rows.angle_limits(case, &network, va);
    rows.piecewise(&costs, pg, None, y);
    let (a, l, u) = rows.finish(nx);

    // The objective, without the constant terms
    let quadratic: Vec<(usize, usize, f64)> = costs
        .active
        .quadratic
        .iter()
        .enumerate()
        .filter(|&(_, &q)| q != 0.0)
        .map(|(k, &q)| (pg + k, pg + k, q))
        .collect();
    let h = (!quadratic.is_empty())
        .then(|| CsrMatrix::from_triplets(nx, nx, quadratic).expect("indices in range"));
    let mut c = vec![0.0; nx];
    c[pg..pg + ng].copy_from_slice(&costs.active.linear);
    c[y..].fill(1.0);

    // dcopf_solver's starting point, with every angle at the first reference bus's
    let mut x0 = starting_point(&xmin, &xmax);
    x0[va..pg].fill(xmin[va + refs[0]]);
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
    let qp = QuadraticProblem { h, c };
    let consistent = "the DC OPF model's dimensions are consistent, and it has bounds";
    let solution = match &options.solver {
        OpfSolver::Mips(mips) => {
            SolverResult::Mips(qps_mips(&qp, &constraints, Some(&x0), mips).expect(consistent))
        }
        OpfSolver::Mehrotra(mehrotra) => SolverResult::Mehrotra(
            qps_mehrotra(&qp, &constraints, Some(&x0), mehrotra).expect(consistent),
        ),
    };
    let residuals = kkt_residuals(&qp, &constraints, solution.x(), solution.multipliers())
        .expect("the solver's multipliers fit its problem");

    // The results, back on the case's tables
    let x = solution.x();
    let multipliers = solution.multipliers();
    let mut va_deg: Vec<f64> = case.bus.iter().map(|bus| bus.va_deg).collect();
    let mut lam_p = vec![0.0; case.bus.len()];
    for (i, &row) in network.bus_rows.iter().enumerate() {
        va_deg[row] = x[va + i] * 180.0 / PI;
        lam_p[row] = (multipliers.mu_u[i] - multipliers.mu_l[i]) / base;
    }
    let ngen = case.generator.len();
    let (mut pg_mw, mut mu_pmin, mut mu_pmax) = (vec![0.0; ngen], vec![0.0; ngen], vec![0.0; ngen]);
    for (k, &row) in network.gen_rows.iter().enumerate() {
        pg_mw[row] = x[pg + k] * base;
        mu_pmin[row] = multipliers.lower[pg + k] / base;
        mu_pmax[row] = multipliers.upper[pg + k] / base;
    }
    let nl = case.branch.len();
    let mut pf_mw = vec![0.0; nl];
    let (mut mu_sf, mut mu_st) = (vec![0.0; nl], vec![0.0; nl]);
    let (mut mu_angmin, mut mu_angmax) = (vec![0.0; nl], vec![0.0; nl]);
    let flows = bf.mul_vec(&x[va..pg]).expect("Bf has a column per bus");
    for (k, &row) in network.branch_rows.iter().enumerate() {
        pf_mw[row] = (flows[k] + pfinj_pu[k]) * base;
    }
    for (n, &k) in limited.iter().enumerate() {
        let row = network.branch_rows[k];
        mu_sf[row] = multipliers.mu_u[nb + n] / base;
        mu_st[row] = multipliers.mu_l[nb + n] / base;
    }
    for (n, &k) in angle_limited.iter().enumerate() {
        let row = network.branch_rows[k];
        mu_angmin[row] = per_degree(multipliers.mu_l[first_angle_row + n]);
        mu_angmax[row] = per_degree(multipliers.mu_u[first_angle_row + n]);
    }
    let constant: f64 = costs.active.constant.iter().sum();

    Ok(DcOpfSolution {
        va_deg,
        lam_p,
        pg_mw,
        mu_pmin,
        mu_pmax,
        pf_mw,
        mu_sf,
        mu_st,
        mu_angmin,
        mu_angmax,
        objective: solution.f() + constant,
        status: solution.status(),
        iterations: solution.iterations(),
        residuals,
        solver: solution,
    })
}

/// `makeAang`'s angle difference limits over the in-service branches at case rows
/// `branch_rows`: each limited branch's internal index with its lower and upper limits (rad).
///
/// A branch is limited when either limit is nonzero and inside ±360°, or exactly one of them
/// is zero; a limit at or past ±360° is then unbounded on that side, and a zero limit is
/// literally 0°, as in MATPOWER.
pub(crate) fn angle_limits(case: &Case, branch_rows: &[usize]) -> Vec<(usize, f64, f64)> {
    branch_rows
        .iter()
        .enumerate()
        .filter_map(|(k, &row)| {
            let (min, max) = (case.branch[row].angmin_deg, case.branch[row].angmax_deg);
            let limited = (min != 0.0 && min > -360.0)
                || (max != 0.0 && max < 360.0)
                || (min != 0.0 && max == 0.0)
                || (min == 0.0 && max != 0.0);
            limited.then(|| {
                let lower = if min < -360.0 { f64::NEG_INFINITY } else { min };
                let upper = if max > 360.0 { f64::INFINITY } else { max };
                (k, lower * PI / 180.0, upper * PI / 180.0)
            })
        })
        .collect()
}
