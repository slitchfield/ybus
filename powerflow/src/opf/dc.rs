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
//! options for the OPF (`feastol` = `opf.violation` = 5e-6, the rest MIPS's defaults), and the
//! results and multipliers are mapped back to the case as `dcopf_solver` and `opf` map them.

use std::f64::consts::PI;

use interior_point::mips::{qps_mips, MipsOptions, MipsSolution, MipsStatus};
use interior_point::problem::{Constraints, QuadraticProblem};
use sparse_mat::compressed_matrix::CsrMatrix;

use crate::case::{BusType, Case};
use crate::dc::Bdc;
use crate::opf::cost::{CostVariable, OpfCosts};
use crate::powerflow_error::PowerflowError;
use crate::topology::Network;

/// The DC optimal power flow's options.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct DcOpfOptions {
    /// MIPS's options. The default is MATPOWER's for the OPF: MIPS's defaults with `feastol`
    /// set to `opf.violation`, 5e-6.
    pub mips: MipsOptions,
}

impl Default for DcOpfOptions {
    fn default() -> Self {
        let mut mips = MipsOptions::default();
        mips.feastol = 5e-6;
        DcOpfOptions { mips }
    }
}

/// A DC optimal power flow's results, indexed like the case's own tables, as `rundcopf`
/// reports them.
///
/// Out-of-service generators and branches, and those at an isolated bus, report 0.0 for every
/// quantity, as MATPOWER's `opf` sets them. An isolated bus keeps its angle from the case and
/// has no prices. Results are reported whether or not MIPS converged: see
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
    /// How MIPS ended
    pub status: MipsStatus,
    /// MIPS's iterations
    pub iterations: usize,
    /// MIPS's own result, on the model's variables `[Va, Pg, y]` and constraints, with each
    /// iteration's record
    pub solver: MipsSolution,
}

impl DcOpfSolution {
    /// Whether MIPS converged, MATPOWER's `success`.
    pub fn converged(&self) -> bool {
        self.status == MipsStatus::Converged
    }
}

/// Solves `case`'s DC optimal power flow with MATPOWER's defaults; see [`dc_opf_with_options`].
pub fn dc_opf(case: &Case) -> Result<DcOpfSolution, PowerflowError> {
    dc_opf_with_options(case, &DcOpfOptions::default())
}

/// Solves `case`'s DC optimal power flow as MATPOWER's `rundcopf` does with MIPS (see the
/// [module documentation](self)). A solve that MIPS doesn't finish is still returned, with its
/// [`DcOpfSolution::status`].
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
    let refs: Vec<usize> = (0..nb)
        .filter(|&i| case.bus[network.bus_rows[i]].bus_type == BusType::Ref)
        .collect();
    if refs.is_empty() {
        return Err(PowerflowError::NoReferenceBus);
    }
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
    let mut triplets: Vec<(usize, usize, f64)> = vec![];
    let (mut l, mut u) = (vec![], vec![]);
    triplets.extend(bbus.iter().map(|(i, j, b)| (i, va + j, b)));
    for (k, &bus) in network.gen_buses.iter().enumerate() {
        triplets.push((bus, pg + k, -1.0));
    }
    for (i, &row) in network.bus_rows.iter().enumerate() {
        let bus = &case.bus[row];
        let mismatch = -(bus.pd_mw + bus.gs_mw) / base - pbusinj_pu[i];
        l.push(mismatch);
        u.push(mismatch);
    }
    let limited: Vec<usize> = (0..network.branch_rows.len())
        .filter(|&k| {
            let rate = case.branch[network.branch_rows[k]].rate_a_mva;
            rate != 0.0 && rate < 1e10
        })
        .collect();
    for &k in &limited {
        let row = l.len();
        let (cols, values) = bf.outer(k).expect("an internal branch");
        triplets.extend(cols.iter().zip(values).map(|(&j, &b)| (row, va + j, b)));
        let rate = case.branch[network.branch_rows[k]].rate_a_mva / base;
        l.push(-(rate + pfinj_pu[k]));
        u.push(rate - pfinj_pu[k]);
    }
    let angle_limited = angle_limits(case, &network.branch_rows);
    for &(k, lower, upper) in &angle_limited {
        let row = l.len();
        let (from, to) = network.ends[k];
        triplets.extend([(row, va + from, 1.0), (row, va + to, -1.0)]);
        l.push(lower);
        u.push(upper);
    }
    for (j, cost) in costs.piecewise.iter().enumerate() {
        let CostVariable::Active(k) = cost.variable else {
            unreachable!("a DC OPF has no reactive costs")
        };
        for (&m, &b) in cost.slopes.iter().zip(&cost.bounds) {
            let row = l.len();
            triplets.extend([(row, pg + k, m), (row, y + j, -1.0)]);
            l.push(f64::NEG_INFINITY);
            u.push(b);
        }
    }
    let a = CsrMatrix::from_triplets(l.len(), nx, triplets).expect("indices in range");

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

    // dcopf_solver's starting point: midway between the bounds (infinite ones replaced by
    //   ±1e10), or 1 inside a one-sided bound; every angle at the first reference bus's; and
    //   each y at the largest piecewise linear cost's last value, plus 10%
    let mut x0: Vec<f64> = xmin
        .iter()
        .zip(&xmax)
        .map(
            |(&lo, &hi)| match (lo == f64::NEG_INFINITY, hi == f64::INFINITY) {
                (true, false) => hi - 1.0,
                (false, true) => lo + 1.0,
                (true, true) => 0.0,
                (false, false) => (lo + hi) / 2.0,
            },
        )
        .collect();
    x0[va..pg].fill(xmin[va + refs[0]]);
    if let Some(largest) = costs
        .piecewise
        .iter()
        .map(|cost| cost.last_cost)
        .reduce(f64::max)
    {
        x0[y..].fill(largest + 0.1 * largest.abs());
    }

    let constraints = Constraints {
        a: Some(a),
        l: Some(l),
        u: Some(u),
        xmin: Some(xmin),
        xmax: Some(xmax),
    };
    let solution = qps_mips(
        &QuadraticProblem { h, c },
        &constraints,
        Some(&x0),
        &options.mips,
    )
    .expect("the DC OPF model's dimensions are consistent, and it has bounds");

    // The results, back on the case's tables
    let x = &solution.x;
    let multipliers = &solution.multipliers;
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
    let first_angle_row = nb + limited.len();
    for (n, &(k, _, _)) in angle_limited.iter().enumerate() {
        let row = network.branch_rows[k];
        mu_angmin[row] = multipliers.mu_l[first_angle_row + n] * PI / 180.0;
        mu_angmax[row] = multipliers.mu_u[first_angle_row + n] * PI / 180.0;
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
        objective: solution.f + constant,
        status: solution.status,
        iterations: solution.iterations,
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
