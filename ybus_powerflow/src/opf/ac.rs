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
//! formulas (in `opf::derivatives`, public only for its tests). The results and multipliers are mapped back to the case
//! as `nlpopf_solver` and `opf` map them.
//!
//! [`mips`]: ybus_ipm::mips::mips
//! [`mehrotra`]: ybus_ipm::mehrotra::mehrotra

use ybus_ipm::kkt::KktFactorization;
use ybus_ipm::mehrotra::{mehrotra_from, MehrotraOptions};
use ybus_ipm::mips::{mips_from, MipsOptions};
use ybus_ipm::problem::{
    same_pattern, Constraints, NonlinearConstraints, NonlinearMultipliers, Problem,
};
use ybus_ipm::solution::{kkt_residuals, KktResiduals, Solution, Status};
use ybus_sparse::compressed_matrix::{concurrently, AssemblyMap, CsrMatrix, SumPlan};
use ybus_sparse::Complex;

use crate::ac::{check_impedance, Admittances};
use crate::case::Case;
use crate::opf::cost::{CostVariable, OpfCosts, QuadraticCosts};
use crate::opf::derivatives::{
    squared_magnitude_hessian, weighted_hessian, Flows, FlowsPlan, SquaredMagnitudeHessianPlan,
    SquaredMagnitudesPlan, WeightedHessianPlan,
};
use crate::opf::{
    helper_start, per_degree, reference_buses, starting_point, Key, Layout, LinearRows, OpfSolver,
    OpfStart, WarmStart,
};
use crate::powerflow_error::PowerflowError;
use crate::topology::{Network, Signature};

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
    pub solver: Solution,
}

impl AcOpfSolution {
    /// Whether the solver converged, MATPOWER's `success`.
    pub fn converged(&self) -> bool {
        self.status == Status::Converged
    }
}

/// The AC OPF's objective and nonlinear constraints over `x = [Va, Vm, Pg, Qg, y]`, for MIPS.
#[derive(Clone)]
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
    /// The constraints' plans and the Lagrangian Hessian's, made once the rest of the problem
    /// is (`ConstraintPlans::new`, `HessianPlans::new`)
    constraint_plans: Option<ConstraintPlans>,
    hessian: Option<HessianPlans>,
}

/// The plans of `AcOpfProblem::constraints`, which compute the constraints and their Jacobians
/// bit for bit as the functions do
#[derive(Clone)]
struct ConstraintPlans {
    bus: FlowsPlan,
    /// The limited branches' from and to ends: their flows, and their squared magnitudes
    branches: Option<[(FlowsPlan, SquaredMagnitudesPlan); 2]>,
    /// `widen`'s assembly of the power balance's Jacobian: the outputs' entries, then the
    /// voltages' four blocks
    dg: CsrMatrix<f64>,
    dg_map: AssemblyMap,
    /// `widen`'s assembly of the flow limits' Jacobian: the from ends' rows, then the to ends'
    dh: CsrMatrix<f64>,
    dh_map: AssemblyMap,
    /// Whether the independent evaluations run concurrently (feature `rayon`): see
    /// `EVALUATION_BUSES`
    parallel: bool,
}

/// The bus count from which an AC OPF's evaluation runs its independent parts concurrently
/// (feature `rayon`): the bus and branch-end flows, Hessians and their parts. Fixed, not a thread
/// count, and the results are the same either way. From a calibration at 14 threads
/// (`docs/benchmarks/o5-parallel-evaluation-2026-11-21.md`): case300 was slower with it,
/// ×0.75–0.87, and case1354pegase and every larger case faster.
const EVALUATION_BUSES: usize = 1000;

impl ConstraintPlans {
    fn new(problem: &AcOpfProblem) -> Self {
        let nb = problem.nb;
        let bus = FlowsPlan::new(&problem.identity, &problem.ybus);
        // widen's triplets, in its order: the outputs', then each block's
        let block = bus.evaluate(&vec![C64::new(1.0, 0.0); nb]).ds_dva.real();
        let mut triplets: Vec<(usize, usize, f64)> = problem.outputs().collect();
        for (row, col) in [(0, 0), (0, nb), (nb, 0), (nb, nb)] {
            triplets.extend(block.iter().map(|(i, j, v)| (row + i, col + j, v)));
        }
        let (dg, dg_map) = CsrMatrix::from_triplets_mapped(2 * nb, problem.nx(), triplets)
            .expect("indices in range");
        let nl2 = problem.flow_max.len();
        let branches = (nl2 > 0).then(|| {
            [(&problem.cf, &problem.yf), (&problem.ct, &problem.yt)]
                .map(|(c, y)| (FlowsPlan::new(c, y), SquaredMagnitudesPlan::new(c, y)))
        });
        let mut triplets = Vec::new();
        if let Some(ends) = &branches {
            for (row, (_, magnitudes)) in [0, nl2].into_iter().zip(ends) {
                triplets.extend(magnitudes.pattern().iter().map(|(i, j, v)| (row + i, j, v)));
            }
        }
        let (dh, dh_map) = CsrMatrix::from_triplets_mapped(2 * nl2, problem.nx(), triplets)
            .expect("indices in range");
        ConstraintPlans {
            bus,
            branches,
            dg,
            dg_map,
            dh,
            dh_map,
            parallel: problem.nb >= EVALUATION_BUSES,
        }
    }
}

/// The plans of `AcOpfProblem::lagrangian_hessian`, which compute it bit for bit as the
/// functions do, from patterns that depend only on the network and the costs
#[derive(Clone)]
struct HessianPlans {
    bus: WeightedHessianPlan,
    /// The limited branches' from and to ends, and the sums `bus + from` and `+ to`
    branches: Option<BranchHessianPlans>,
    /// `widen`'s assembly of the costs' entries and the voltage block, as a kept map
    widened: CsrMatrix<f64>,
    widen: AssemblyMap,
    /// Whether the bus Hessian's pattern stores every position a branch end's Hessian can lack
    /// (each bus's `(i, i)`, `(i, nb + i)` and `(nb + i, i)`), so that summing it in restores
    /// the position as `try_add`'s union does. It does whenever Ybus stores its diagonal.
    bus_covers_diagonals: bool,
    /// Whether the independent evaluations run concurrently (feature `rayon`): see
    /// `EVALUATION_BUSES`
    parallel: bool,
}

#[derive(Clone)]
struct BranchHessianPlans {
    from: SquaredMagnitudeHessianPlan,
    to: SquaredMagnitudeHessianPlan,
    bus_from: SumPlan,
    bus_from_to: SumPlan,
}

impl HessianPlans {
    fn new(problem: &AcOpfProblem) -> Self {
        let parallel = problem.nb >= EVALUATION_BUSES;
        let bus = WeightedHessianPlan::new(&problem.identity, &problem.ybus).parallel(parallel);
        let nb = problem.nb;
        let stores = |row, col| {
            let pattern = bus.pattern();
            let start = pattern.pointers()[row];
            pattern.indices()[start..pattern.pointers()[row + 1]]
                .binary_search(&col)
                .is_ok()
        };
        let bus_covers_diagonals =
            (0..nb).all(|i| stores(i, i) && stores(i, nb + i) && stores(nb + i, i));
        let mut voltage = bus.pattern().clone();
        let branches = (!problem.flow_max.is_empty()).then(|| {
            let from =
                SquaredMagnitudeHessianPlan::new(&problem.cf, &problem.yf).parallel(parallel);
            let to = SquaredMagnitudeHessianPlan::new(&problem.ct, &problem.yt).parallel(parallel);
            let bus_from = voltage.plan_sum(from.pattern()).expect("both 2nb × 2nb");
            voltage = bus_from.add(&voltage, from.pattern()).expect("planned");
            let bus_from_to = voltage.plan_sum(to.pattern()).expect("both 2nb × 2nb");
            voltage = bus_from_to.add(&voltage, to.pattern()).expect("planned");
            BranchHessianPlans {
                from,
                to,
                bus_from,
                bus_from_to,
            }
        });
        // widen's triplets, in its order: the costs' entries, then the voltage block's
        let mut triplets = problem.cost_hessian();
        triplets.extend(voltage.iter());
        let (widened, widen) =
            CsrMatrix::from_triplets_mapped(problem.nx(), problem.nx(), triplets)
                .expect("indices in range");
        HessianPlans {
            bus,
            branches,
            widened,
            widen,
            bus_covers_diagonals,
            parallel,
        }
    }
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
        let plans = self
            .constraint_plans
            .as_ref()
            .expect("made with the problem");
        plans.branches.as_ref().map(|[(from, _), (to, _)]| {
            concurrently(plans.parallel, || from.evaluate(v), || to.evaluate(v))
        })
    }

    /// The power balance's entries for the generators' outputs, `-1` at each generator's bus
    /// row and output column, in `widen`'s order
    fn outputs(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        let (nb, ng) = (self.nb, self.ng);
        let (pg, qg) = (2 * nb, 2 * nb + ng);
        self.gen_buses
            .iter()
            .enumerate()
            .flat_map(move |(k, &b)| [(b, pg + k, -1.0), (nb + b, qg + k, -1.0)])
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
        let mut out = NonlinearConstraints::empty(self.nx());
        self.constraints_into(x, &mut out);
        out
    }

    /// The constraints written into `out`: the Jacobians' values refilled in place when `out`
    /// already holds the planned patterns, as it does from the second call of a solve on
    fn constraints_into(&self, x: &[f64], out: &mut NonlinearConstraints) {
        let (nb, ng) = (self.nb, self.ng);
        let (pg, qg) = (2 * nb, 2 * nb + ng);
        let v = self.voltages(x);

        let plans = self
            .constraint_plans
            .as_ref()
            .expect("made with the problem");

        // The bus flows and each branch end's flows and squared magnitudes are independent:
        //   concurrently if the plans are parallel, each the same computation either way
        let ends =
            |plan: &(FlowsPlan, SquaredMagnitudesPlan)| plan.1.evaluate(&plan.0.evaluate(&v));
        let (bus, ends) = concurrently(
            plans.parallel,
            || plans.bus.evaluate(&v),
            || {
                plans
                    .branches
                    .as_ref()
                    .map(|[from, to]| concurrently(plans.parallel, || ends(from), || ends(to)))
            },
        );

        // The power balance and its Jacobian
        let mut mismatch: Vec<C64> = bus.s.iter().zip(&self.load).map(|(s, d)| s + d).collect();
        for (k, &b) in self.gen_buses.iter().enumerate() {
            mismatch[b] -= C64::new(x[pg + k], x[qg + k]);
        }
        out.g.clear();
        out.g.extend(
            mismatch
                .iter()
                .map(|m| m.re)
                .chain(mismatch.iter().map(|m| m.im)),
        );
        // widen's assembly of the outputs' entries and the blocks ∂P/∂Va, ∂P/∂Vm, ∂Q/∂Va and
        //   ∂Q/∂Vm, each the real or imaginary part of a derivative
        let (dva, dvm) = (bus.ds_dva.values(), bus.ds_dvm.values());
        let values: Vec<f64> = (self.outputs().map(|(_, _, v)| v))
            .chain(dva.iter().map(|s| s.re))
            .chain(dvm.iter().map(|s| s.re))
            .chain(dva.iter().map(|s| s.im))
            .chain(dvm.iter().map(|s| s.im))
            .collect();
        refill(&mut out.dg, &plans.dg, &plans.dg_map, &values);

        // The squared apparent power limits and their Jacobian
        out.h.clear();
        match ends {
            Some(((af, daf), (at, dat))) => {
                out.h.extend(
                    af.iter()
                        .chain(&at)
                        .zip(self.flow_max.iter().chain(&self.flow_max))
                        .map(|(a, max)| a - max),
                );
                // widen's assembly: the from ends' rows, then the to ends'
                let values: Vec<f64> = daf.values().iter().chain(dat.values()).copied().collect();
                refill(&mut out.dh, &plans.dh, &plans.dh_map, &values);
            }
            None => refill(&mut out.dh, &plans.dh, &plans.dh_map, &[]),
        }
    }

    fn lagrangian_hessian(
        &self,
        x: &[f64],
        multipliers: NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        let mut out = CsrMatrix::from_triplets(0, 0, vec![]).expect("an empty matrix");
        self.lagrangian_hessian_into(x, multipliers, cost_mult, &mut out);
        out
    }

    /// The Hessian written into `out`, its values refilled in place when `out` already holds
    /// the planned pattern
    fn lagrangian_hessian_into(
        &self,
        x: &[f64],
        multipliers: NonlinearMultipliers<'_>,
        cost_mult: f64,
        out: &mut CsrMatrix<f64>,
    ) {
        let v = self.voltages(x);
        if !self.planned_lagrangian_hessian(&v, &multipliers, cost_mult, out) {
            *out = self.direct_lagrangian_hessian(&v, &multipliers, cost_mult);
        }
    }
}

/// `matrix`'s values reassembled by `map` from `values`, after it takes `template`'s pattern
/// if it doesn't have it yet: a kept `widen`'s assembly into a caller's buffer
fn refill(
    matrix: &mut CsrMatrix<f64>,
    template: &CsrMatrix<f64>,
    map: &AssemblyMap,
    values: &[f64],
) {
    if !same_pattern(matrix, template) {
        *matrix = template.clone();
    }
    matrix.refill(map, values).expect("planned");
}

impl AcOpfProblem {
    /// The Lagrangian Hessian from the kept plans, bit for bit `direct_lagrangian_hessian`'s,
    /// written into `out`; false, with `out` untouched, where the plans decline (see
    /// `HessianPlans`)
    fn planned_lagrangian_hessian(
        &self,
        v: &[C64],
        multipliers: &NonlinearMultipliers<'_>,
        cost_mult: f64,
        out: &mut CsrMatrix<f64>,
    ) -> bool {
        self.planned_lagrangian_hessian_values(v, multipliers, cost_mult)
            .map(|values| {
                let plans = self.hessian.as_ref().expect("made with the problem");
                refill(out, &plans.widened, &plans.widen, &values);
            })
            .is_some()
    }

    /// The planned Hessian's values in `widen`'s triplet order, or `None` where the plans
    /// decline
    fn planned_lagrangian_hessian_values(
        &self,
        v: &[C64],
        multipliers: &NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> Option<Vec<f64>> {
        let plans = self.hessian.as_ref().expect("made with the problem");
        let (lam_p, lam_q) = multipliers.eq.split_at(self.nb);
        // The bus Hessian must match its plan's pattern. A branch end's may lack positions,
        //   which hold -0.0: the bus Hessian stores them, so the sums give try_add's union,
        //   with the bus Hessian's values copied there
        // The bus Hessian and each branch end's are independent: concurrently if the plans are
        //   parallel, each the same computation either way
        let (voltage, ends) = concurrently(
            plans.parallel,
            || plans.bus.evaluate(v, lam_p, lam_q),
            || match (&plans.branches, self.branch_flows(v)) {
                (Some(branches), Some((from, to))) => {
                    let (mu_f, mu_t) = multipliers.ineq.split_at(self.flow_max.len());
                    Some(concurrently(
                        plans.parallel,
                        || branches.from.evaluate_with_absent(v, &from, mu_f),
                        || branches.to.evaluate_with_absent(v, &to, mu_t),
                    ))
                }
                _ => None,
            },
        );
        let mut voltage = voltage?;
        if let (Some(branches), Some((hf, ht))) = (&plans.branches, ends) {
            let ((hf, absent_f), (ht, absent_t)) = (hf?, ht?);
            if !plans.bus_covers_diagonals && !(absent_f.is_empty() && absent_t.is_empty()) {
                return None;
            }
            voltage = branches.bus_from.add(&voltage, &hf).expect("planned");
            voltage = branches.bus_from_to.add(&voltage, &ht).expect("planned");
        }
        let values: Vec<f64> = self
            .cost_hessian()
            .into_iter()
            .map(|(_, _, q)| cost_mult * q)
            .chain(voltage.values().iter().copied())
            .collect();
        Some(values)
    }

    /// The Lagrangian Hessian from the derivative functions, the reference the plans replay
    fn direct_lagrangian_hessian(
        &self,
        v: &[C64],
        multipliers: &NonlinearMultipliers<'_>,
        cost_mult: f64,
    ) -> CsrMatrix<f64> {
        let nb = self.nb;
        let v = v.to_vec();
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
    AcOpfModel::new(case)?.solve(case, options)
}

/// What an AC OPF reads from its case's values, in the order the checks run: the costs, the
/// variables' bounds, the linear rows, and the limited branches' limits, with the loads.
struct Inputs {
    costs: OpfCosts,
    xmin: Vec<f64>,
    xmax: Vec<f64>,
    a: CsrMatrix<f64>,
    l: Vec<f64>,
    u: Vec<f64>,
    /// The generators (internal) with a power factor row, the dispatchable loads
    power_factor: Vec<usize>,
    first_angle_row: usize,
    angle_limited: Vec<usize>,
    /// The branches (internal) with a flow limit
    limited: Vec<usize>,
    flow_max: Vec<f64>,
    load: Vec<C64>,
}

impl Inputs {
    fn new(case: &Case, network: &Network, refs: &[usize]) -> Result<Inputs, PowerflowError> {
        let nb = network.bus_rows.len();
        let ng = network.gen_rows.len();
        let base = case.base_mva;
        let costs = OpfCosts::new(case, &network.gen_rows, true)?;
        let ny = costs.piecewise.len();
        for &row in &network.gen_rows {
            check_capability_curve(case, row)?;
        }
        let (va, vm, pg, qg, y) = (0, nb, 2 * nb, 2 * nb + ng, 2 * nb + 2 * ng);
        let nx = 2 * nb + 2 * ng + ny;

        // The variables' bounds
        let mut xmin = vec![f64::NEG_INFINITY; nx];
        let mut xmax = vec![f64::INFINITY; nx];
        for &i in refs {
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
        let mut power_factor_rows = Vec::new();
        for (k, &row) in network.gen_rows.iter().enumerate() {
            if let Some((p, q)) = power_factor(case, row)? {
                rows.push([(pg + k, p), (qg + k, q)], 0.0, 0.0);
                power_factor_rows.push(k);
            }
        }
        let first_angle_row = rows.len();
        let angle_limited = rows.angle_limits(case, network, va);
        rows.piecewise(&costs, pg, Some(qg), y);
        let (a, l, u) = rows.finish(nx);

        let limited: Vec<usize> = (0..network.branch_rows.len())
            .filter(|&k| {
                let rate = case.branch[network.branch_rows[k]].rate_a_mva;
                rate != 0.0 && rate < 1e10
            })
            .collect();
        let flow_max = limited
            .iter()
            .map(|&k| (case.branch[network.branch_rows[k]].rate_a_mva / base).powi(2))
            .collect();
        let load = network
            .bus_rows
            .iter()
            .map(|&row| C64::new(case.bus[row].pd_mw, case.bus[row].qd_mvar) / base)
            .collect();
        Ok(Inputs {
            costs,
            xmin,
            xmax,
            a,
            l,
            u,
            power_factor: power_factor_rows,
            first_angle_row,
            angle_limited,
            limited,
            flow_max,
            load,
        })
    }
}

/// The AC optimal power flow of one network, built once and solved for many cases on it, as
/// [`AcModel`](crate::ac::AcModel) is for the power flow: each solve keeps the network, its
/// admittances and the problem's kept plans of its derivatives and Newton system, and reads only
/// the values that change.
///
/// [`AcOpfModel::solve`] takes the case the model was built from, or a copy with different
/// values, for example a load sweep or another hour of a day:
///
/// - **Read from the case:** the loads (`pd_mw`, `qd_mvar`); the generators' limits and data
///   (`pmin_mw`, `pmax_mw`, `qmin_mvar`, `qmax_mvar`, `pg_mw`, `qg_mvar`, `vg_pu`); the costs'
///   coefficients and breakpoints; the buses' voltage limits (`vmin_pu`, `vmax_pu`) and the
///   reference angles; and the branches' flow limits (`rate_a_mva`) and angle limits.
/// - **Kept from the model:** the topology and the branches' electrical data (`r`, `x`, `b`,
///   taps and shifts), which make the admittances.
/// - **Refused** with [`PowerflowError::TopologyMismatch`]: a changed topology (as
///   [`AcModel::solve`](crate::ac::AcModel::solve) refuses), or a change to which branches are
///   limited (`rate_a_mva` becoming 0 or not), which branches have angle limits, which
///   generators are dispatchable loads, a piecewise cost's number of points, or which
///   generators' quadratic cost coefficients are 0. Each of those fixes a matrix's pattern; a
///   case that changes one needs a new model.
///
/// [`ac_opf_with_options`] is `AcOpfModel::new(case)?.solve(case, options)`, so a model's
/// solve of its own case gives the same results, bit for bit.
pub struct AcOpfModel {
    signature: Signature,
    network: Network,
    refs: Vec<usize>,
    admittances: Admittances,
    /// The problem with its plans, from the model's case; each solve copies it with its case's
    /// values
    problem: AcOpfProblem,
    /// What fixes the patterns beyond the topology, from the model's case
    power_factor: Vec<usize>,
    angle_limited: Vec<usize>,
    limited: Vec<usize>,
    /// Each piecewise cost's variable and number of segments
    piecewise: Vec<(CostVariable, usize)>,
    cost_hessian: Vec<(usize, usize)>,
}

impl AcOpfModel {
    /// `nlpopf_solver`'s starting point for `case`, with the variables' bounds and the costs
    ///   from its `Inputs`: every angle at the first reference bus's, and each magnitude midway
    ///   between its limits, clipped to [0.5, 1.5]
    fn cold_point(&self, case: &Case, xmin: &[f64], xmax: &[f64], costs: &OpfCosts) -> Vec<f64> {
        let nb = self.network.bus_rows.len();
        let ng = self.network.gen_rows.len();
        let (va, vm, y) = (0, nb, 2 * nb + 2 * ng);
        let mut x0 = starting_point(xmin, xmax);
        x0[va..vm].fill(xmin[va + self.refs[0]]);
        for (i, &row) in self.network.bus_rows.iter().enumerate() {
            let bus = &case.bus[row];
            x0[vm + i] = (bus.vmax_pu.min(1.5) + bus.vmin_pu.max(0.5)) / 2.0;
        }
        if let Some(start) = helper_start(costs) {
            x0[y..].fill(start);
        }
        x0
    }

    /// The problem's layout for `case` (see `opf::Layout`), and the `Inputs` it was read from:
    ///   `[Va, Vm, Pg, Qg, y]`; the power balances `P`, then `Q`; the flow limits `Sf`, then
    ///   `St`; and the rows of `A` as `Inputs::new` pushes them, power factors, angle limits,
    ///   then the piecewise costs' segments
    fn layout(&self, case: &Case) -> Result<(Layout, Inputs), PowerflowError> {
        Network::new(case, check_impedance)?;
        self.signature.check(case)?;
        let inputs = Inputs::new(case, &self.network, &self.refs)?;
        let network = &self.network;
        let buses = &network.bus_rows;
        let gens = &network.gen_rows;
        let branch = |k: &usize| network.branch_rows[*k];
        let mut layout = Layout::default();
        layout.columns.extend(buses.iter().map(|&row| Key::Va(row)));
        layout.columns.extend(buses.iter().map(|&row| Key::Vm(row)));
        layout.columns.extend(gens.iter().map(|&row| Key::Pg(row)));
        layout.columns.extend(gens.iter().map(|&row| Key::Qg(row)));
        layout.eq.extend(buses.iter().map(|&row| Key::P(row)));
        layout.eq.extend(buses.iter().map(|&row| Key::Q(row)));
        layout
            .ineq
            .extend(inputs.limited.iter().map(|k| Key::Sf(branch(k))));
        layout
            .ineq
            .extend(inputs.limited.iter().map(|k| Key::St(branch(k))));
        layout.rows.extend(
            inputs
                .power_factor
                .iter()
                .map(|&k| Key::PowerFactor(gens[k])),
        );
        layout
            .rows
            .extend(inputs.angle_limited.iter().map(|k| Key::Angle(branch(k))));
        layout.piecewise(&inputs.costs, gens);
        Ok((layout, inputs))
    }

    /// A warm start for this model's solve of `case` ([`OpfStart::Mapped`]) from `previous`,
    /// `from`'s solution of `from_case`: another model of the same case's tables, for example
    /// with a branch or a generator out of service.
    ///
    /// Each variable and multiplier is matched by what it belongs to, by case row: a bus's
    /// voltage and power balances, a generator's outputs, a branch's flow and angle limits, a
    /// piecewise cost's helper and segments, a dispatchable load's power factor. A match takes
    /// `previous`'s value, from its final iterate's multipliers
    /// ([`Solution::iterate_multipliers`]). A variable without one, such as a generator back in
    /// service, takes its cold start's value, and a multiplier 0, which the solver's shift
    /// raises. Both models' `cost_mult` must be the same, as for [`OpfStart::Warm`].
    ///
    /// Returns `PowerflowError::CaseMismatch` if the cases' tables differ in length, `case`'s
    /// or `from_case`'s errors as [`AcOpfModel::solve`] would return them, and
    /// `PowerflowError::Solver` (`Dimension`) if `previous` doesn't fit `from`'s problem.
    pub fn map_start(
        &self,
        case: &Case,
        from: &AcOpfModel,
        from_case: &Case,
        previous: &Solution,
    ) -> Result<WarmStart, PowerflowError> {
        let shape = |case: &Case| (case.bus.len(), case.generator.len(), case.branch.len());
        if shape(case) != shape(from_case) {
            return Err(PowerflowError::CaseMismatch {
                expected: shape(from_case),
                found: shape(case),
            });
        }
        let (ours, inputs) = self.layout(case)?;
        let (theirs, _) = from.layout(from_case)?;
        let cold = self.cold_point(case, &inputs.xmin, &inputs.xmax, &inputs.costs);
        ours.map(&theirs, previous, &cold)
    }

    /// Builds the model of `case`'s network, checking the case as [`ac_opf`] does.
    ///
    /// Returns the errors [`ac_opf`] returns before solving.
    pub fn new(case: &Case) -> Result<AcOpfModel, PowerflowError> {
        let network = Network::new(case, check_impedance)?;
        let nb = network.bus_rows.len();
        let ng = network.gen_rows.len();
        let refs = reference_buses(case, &network)?;
        let inputs = Inputs::new(case, &network, &refs)?;
        let admittances =
            Admittances::new(case, &network.bus_rows, &network.branch_rows, &network.ends);

        // The limited branches' rows of the admittances, with their incidence
        let limited = &inputs.limited;
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
        let select = |m: &CsrMatrix<C64>| m.select(limited, &all_buses).expect("rows in range");
        let mut problem = AcOpfProblem {
            nb,
            ng,
            ny: inputs.costs.piecewise.len(),
            identity: CsrMatrix::identity(nb),
            cf: incidence(|&(f, _)| f),
            yf: select(&admittances.yf),
            ct: incidence(|&(_, t)| t),
            yt: select(&admittances.yt),
            flow_max: inputs.flow_max.clone(),
            load: inputs.load.clone(),
            gen_buses: network.gen_buses.clone(),
            active: inputs.costs.active.clone(),
            reactive: inputs.costs.reactive.clone(),
            ybus: admittances.ybus.clone(),
            constraint_plans: None,
            hessian: None,
        };
        problem.constraint_plans = Some(ConstraintPlans::new(&problem));
        problem.hessian = Some(HessianPlans::new(&problem));
        let cost_hessian = cost_hessian_positions(&problem);
        Ok(AcOpfModel {
            signature: Signature::new(case),
            refs,
            admittances,
            problem,
            power_factor: inputs.power_factor,
            angle_limited: inputs.angle_limited,
            limited: inputs.limited,
            piecewise: piecewise_shape(&inputs.costs),
            cost_hessian,
            network,
        })
    }

    /// Solves the AC optimal power flow for `case`, which must be the model's case or a copy
    /// with different values (see the type's documentation), as [`ac_opf_with_options`] does.
    ///
    /// Returns [`ac_opf_with_options`]'s errors, and `PowerflowError::CaseMismatch` or
    /// `TopologyMismatch` for a case the model can't solve.
    pub fn solve(
        &self,
        case: &Case,
        options: &AcOpfOptions,
    ) -> Result<AcOpfSolution, PowerflowError> {
        self.solve_from(case, options, OpfStart::Cold)
    }

    /// [`AcOpfModel::solve_from`] for each of `cases`, all from `start`, serially, in order: cases
    /// on this model's network with new values, such as a day's hours (see
    /// [`crate::opf::batch`]).
    pub fn solve_many(
        &self,
        cases: &[Case],
        options: &AcOpfOptions,
        start: OpfStart<'_>,
    ) -> Vec<Result<AcOpfSolution, PowerflowError>> {
        cases
            .iter()
            .map(|case| self.solve_from(case, options, start))
            .collect()
    }

    /// [`AcOpfModel::solve_many`] in parallel (feature `rayon`), on rayon's current pool, the model
    /// shared: bit for bit its results, in the same order, at any thread count.
    #[cfg(feature = "rayon")]
    pub fn par_solve_many(
        &self,
        cases: &[Case],
        options: &AcOpfOptions,
        start: OpfStart<'_>,
    ) -> Vec<Result<AcOpfSolution, PowerflowError>> {
        crate::opf::batch::par_map(cases, |case| self.solve_from(case, options, start))
    }

    /// [`AcOpfModel::solve`] from `start`: MATPOWER's starting point, exactly
    /// [`AcOpfModel::solve`]'s, or a previous solution of this model's (see [`OpfStart`]).
    ///
    /// Returns [`AcOpfModel::solve`]'s errors, and `PowerflowError::Solver` for a previous
    /// solution that doesn't fit the model's problem.
    pub fn solve_from(
        &self,
        case: &Case,
        options: &AcOpfOptions,
        start: OpfStart<'_>,
    ) -> Result<AcOpfSolution, PowerflowError> {
        Network::new(case, check_impedance)?;
        self.signature.check(case)?;
        let network = &self.network;
        let refs = &self.refs;
        let admittances = &self.admittances;
        let nb = network.bus_rows.len();
        let ng = network.gen_rows.len();
        let base = case.base_mva;
        let inputs = Inputs::new(case, network, refs)?;
        self.check_patterns(case, &inputs)?;
        let Inputs {
            costs,
            xmin,
            xmax,
            a,
            l,
            u,
            first_angle_row,
            angle_limited,
            limited,
            flow_max,
            load,
            ..
        } = inputs;
        let mut problem = self.problem.clone();
        problem.flow_max = flow_max;
        problem.load = load;
        problem.active = costs.active.clone();
        problem.reactive = costs.reactive.clone();
        if cost_hessian_positions(&problem) != self.cost_hessian {
            return Err(self.cost_mismatch(&problem));
        }
        let (va, vm, pg, qg) = (0, nb, 2 * nb, 2 * nb + ng);

        let x0 = self.cold_point(case, &xmin, &xmax, &costs);

        let constraints = Constraints {
            a: Some(a),
            l: Some(l),
            u: Some(u),
            xmin: Some(xmin),
            xmax: Some(xmax),
        };
        let refused = |error| PowerflowError::Solver { error };
        let warm = start.warm_start();
        let start = start.solver_start(&x0, &warm);
        let solution = match &options.solver {
            OpfSolver::Mips(options) => {
                mips_from(&problem, &constraints, start, options).map_err(refused)?
            }
            OpfSolver::Mehrotra(options) => {
                mehrotra_from(&problem, &constraints, start, options).map_err(refused)?
            }
        };
        let residuals = kkt_residuals(&problem, &constraints, &solution.x, &solution.multipliers)
            .expect("the solver's multipliers fit its problem");

        // The results, back on the case's tables
        let x = &solution.x;
        let multipliers = &solution.multipliers;
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
        let (mut mu_pmin, mut mu_pmax, mut mu_qmin, mut mu_qmax) =
            (zeros(), zeros(), zeros(), zeros());
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
        let (mut mu_sf, mut mu_st, mut mu_angmin, mut mu_angmax) =
            (zeros(), zeros(), zeros(), zeros());
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
            objective: solution.f,
            status: solution.status,
            iterations: solution.iterations,
            residuals,
            solver: solution,
        })
    }

    /// Refuses `case` if its values change which rows or entries a matrix has, naming the
    /// first such element
    fn check_patterns(&self, case: &Case, inputs: &Inputs) -> Result<(), PowerflowError> {
        let network = &self.network;
        let first_difference = |model: &[usize], found: &[usize]| -> Option<usize> {
            let (mut i, mut j) = (0, 0);
            loop {
                match (model.get(i), found.get(j)) {
                    (Some(a), Some(b)) if a == b => (i, j) = (i + 1, j + 1),
                    (Some(&a), Some(&b)) => return Some(a.min(b)),
                    (Some(&a), None) | (None, Some(&a)) => return Some(a),
                    (None, None) => return None,
                }
            }
        };
        let mismatch = |table, row: usize, field| PowerflowError::TopologyMismatch {
            table,
            row: row + 1,
            field,
        };
        if let Some(k) = first_difference(&self.limited, &inputs.limited) {
            return Err(mismatch("branch", network.branch_rows[k], "rate_a_mva"));
        }
        if let Some(k) = first_difference(&self.angle_limited, &inputs.angle_limited) {
            return Err(mismatch("branch", network.branch_rows[k], "angmin_deg"));
        }
        if let Some(k) = first_difference(&self.power_factor, &inputs.power_factor) {
            return Err(mismatch("gen", network.gen_rows[k], "pmin_mw"));
        }
        let piecewise = piecewise_shape(&inputs.costs);
        if piecewise != self.piecewise {
            let differing = (0..self.piecewise.len().max(piecewise.len()))
                .find(|&n| self.piecewise.get(n) != piecewise.get(n))
                .expect("the shapes differ");
            let variable = self
                .piecewise
                .get(differing)
                .or_else(|| piecewise.get(differing))
                .expect("one of them has it")
                .0;
            return Err(mismatch(
                "gencost",
                self.gencost_row(case, variable),
                "ncost",
            ));
        }
        Ok(())
    }

    /// The refusal for costs whose quadratic coefficients are 0 at other generators
    fn cost_mismatch(&self, problem: &AcOpfProblem) -> PowerflowError {
        let found = cost_hessian_positions(problem);
        let (nb, ng) = (self.network.bus_rows.len(), self.network.gen_rows.len());
        let column = self
            .cost_hessian
            .iter()
            .zip(&found)
            .find(|(a, b)| a != b)
            .map(|(a, b)| a.0.min(b.0))
            .or_else(|| self.cost_hessian.get(found.len()).map(|p| p.0))
            .or_else(|| found.get(self.cost_hessian.len()).map(|p| p.0))
            .expect("the positions differ");
        let k = column - 2 * nb;
        let (row, reactive) = if k < ng { (k, false) } else { (k - ng, true) };
        let ngen = self.signature.shape().1;
        let gen_row = self.network.gen_rows[row];
        PowerflowError::TopologyMismatch {
            table: "gencost",
            row: gen_row + if reactive { ngen } else { 0 } + 1,
            field: "quadratic",
        }
    }

    /// The `gencost` row of a cost: the generator's row, then the reactive costs after all the
    /// active ones
    fn gencost_row(&self, case: &Case, variable: CostVariable) -> usize {
        match variable {
            CostVariable::Active(k) => self.network.gen_rows[k],
            CostVariable::Reactive(k) => case.generator.len() + self.network.gen_rows[k],
        }
    }
}

/// The positions of the costs' Hessian, which `cost_hessian` filters by nonzero coefficients
fn cost_hessian_positions(problem: &AcOpfProblem) -> Vec<(usize, usize)> {
    problem
        .cost_hessian()
        .into_iter()
        .map(|(i, j, _)| (i, j))
        .collect()
}

/// Each piecewise cost's variable and its number of segments, which fix its rows
fn piecewise_shape(costs: &OpfCosts) -> Vec<(CostVariable, usize)> {
    costs
        .piecewise
        .iter()
        .map(|cost| (cost.variable, cost.slopes.len()))
        .collect()
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
