//! Batches of optimal power flows (phase O5, target 1): many solves sharing a base, serially
//! or in parallel (feature `rayon`).
//!
//! - **Outages** ([`ac_opf_outages`], [`dc_opf_outages`], and their `par_` forms): a base case
//!   and a list of outages ([`OpfOutage`], each a set of branches and generators taken out of
//!   service, so N-1 or N-k). The base is solved once, cold. Each outage then gets its own model
//!   and solve, cold or warm from the base solution mapped onto it
//!   ([`AcOpfModel::map_start`](super::ac::AcOpfModel::map_start)), as an N-1 study runs.
//! - **Cases on one network** ([`AcOpfModel::solve_many`](super::ac::AcOpfModel::solve_many),
//!   [`DcOpfModel::solve_many`](super::dc::DcOpfModel::solve_many), and their `par_` forms):
//!   one kept model solving many cases with new values, such as a day's hours, all from one
//!   start.
//!
//! **The parallel forms are bit for bit the serial ones, at any thread count.** Each solve runs
//! on one thread, with its own data and exactly the serial form's operations. The base
//! solution and model are only read, and the results come back in the input's order.

use super::ac::{AcOpfModel, AcOpfOptions, AcOpfSolution};
use super::dc::{DcOpfModel, DcOpfOptions, DcOpfSolution};
use super::{OpfStart, Solution};
use crate::case::Case;
use crate::powerflow_error::PowerflowError;

/// One outage of a batch: the case rows of the branches and generators taken out of service,
/// together, in a copy of the base case. A row already out of service stays out.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct OpfOutage {
    /// The branches' case rows
    pub branches: Vec<usize>,
    /// The generators' case rows
    pub generators: Vec<usize>,
}

impl OpfOutage {
    /// The outage of the branch at case row `row`.
    pub fn branch(row: usize) -> OpfOutage {
        OpfOutage::new(vec![row], vec![])
    }

    /// The outage of the generator at case row `row`.
    pub fn generator(row: usize) -> OpfOutage {
        OpfOutage::new(vec![], vec![row])
    }

    /// The outage of every branch in `branches` and generator in `generators`, together.
    pub fn new(branches: Vec<usize>, generators: Vec<usize>) -> OpfOutage {
        OpfOutage {
            branches,
            generators,
        }
    }

    /// `case` with this outage's branches and generators out of service.
    ///
    /// Returns `PowerflowError::RowOutOfRange` for a row beyond its table.
    pub fn apply(&self, case: &Case) -> Result<Case, PowerflowError> {
        let check = |table: &'static str, rows: &[usize], len: usize| match rows
            .iter()
            .find(|&&row| row >= len)
        {
            Some(&row) => Err(PowerflowError::RowOutOfRange { table, row, len }),
            None => Ok(()),
        };
        check("branch", &self.branches, case.branch.len())?;
        check("generator", &self.generators, case.generator.len())?;
        let mut outage = case.clone();
        for &row in &self.branches {
            outage.branch[row].status = false;
        }
        for &row in &self.generators {
            outage.generator[row].status = false;
        }
        Ok(outage)
    }
}

/// Where a batch's solves start.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum BatchStart {
    /// Each solve from MATPOWER's starting point
    Cold,
    /// Each outage from the base solution, mapped onto its model and shifted by the solver's
    /// `warm_shift` ([`OpfStart::Mapped`])
    Warm,
}

/// A batch of outages' results: the base solution, whether the outages started warm from it,
/// and each outage's result, in the order given.
#[derive(Debug, PartialEq)]
#[non_exhaustive]
pub struct OpfOutages<S> {
    /// The base case's solution, solved cold
    pub base: S,
    /// Whether the outages started from the base solution: a warm batch whose base didn't
    /// converge solves its outages cold
    pub warm: bool,
    /// Each outage's result, with its own errors, such as an outage that leaves a case the OPF
    /// can't solve
    pub outages: Vec<Result<S, PowerflowError>>,
}

/// What a batch needs of a model
trait BatchModel: Sized + Sync {
    type Options: Sync;
    type Solution: Send;

    fn build(case: &Case) -> Result<Self, PowerflowError>;
    fn solve(
        &self,
        case: &Case,
        options: &Self::Options,
        start: OpfStart<'_>,
    ) -> Result<Self::Solution, PowerflowError>;
    fn map(
        &self,
        case: &Case,
        from: &Self,
        from_case: &Case,
        previous: &Solution,
    ) -> Result<super::WarmStart, PowerflowError>;
    fn solver(solution: &Self::Solution) -> &Solution;
}

impl BatchModel for AcOpfModel {
    type Options = AcOpfOptions;
    type Solution = AcOpfSolution;

    fn build(case: &Case) -> Result<Self, PowerflowError> {
        AcOpfModel::new(case)
    }
    fn solve(
        &self,
        case: &Case,
        options: &AcOpfOptions,
        start: OpfStart<'_>,
    ) -> Result<AcOpfSolution, PowerflowError> {
        self.solve_from(case, options, start)
    }
    fn map(
        &self,
        case: &Case,
        from: &Self,
        from_case: &Case,
        previous: &Solution,
    ) -> Result<super::WarmStart, PowerflowError> {
        self.map_start(case, from, from_case, previous)
    }
    fn solver(solution: &AcOpfSolution) -> &Solution {
        &solution.solver
    }
}

impl BatchModel for DcOpfModel {
    type Options = DcOpfOptions;
    type Solution = DcOpfSolution;

    fn build(case: &Case) -> Result<Self, PowerflowError> {
        DcOpfModel::new(case)
    }
    fn solve(
        &self,
        case: &Case,
        options: &DcOpfOptions,
        start: OpfStart<'_>,
    ) -> Result<DcOpfSolution, PowerflowError> {
        self.solve_from(case, options, start)
    }
    fn map(
        &self,
        case: &Case,
        from: &Self,
        from_case: &Case,
        previous: &Solution,
    ) -> Result<super::WarmStart, PowerflowError> {
        self.map_start(case, from, from_case, previous)
    }
    fn solver(solution: &DcOpfSolution) -> &Solution {
        &solution.solver
    }
}

/// The base case's model and cold solution, and the solution the outages start from, if warm
fn base<M: BatchModel>(
    case: &Case,
    options: &M::Options,
) -> Result<(M, M::Solution), PowerflowError> {
    let model = M::build(case)?;
    let solution = model.solve(case, options, OpfStart::Cold)?;
    Ok((model, solution))
}

/// One outage's solve: its own model, from `previous` mapped onto it if warm
fn outage<M: BatchModel>(
    case: &Case,
    model: &M,
    previous: Option<&Solution>,
    outage: &OpfOutage,
    options: &M::Options,
) -> Result<M::Solution, PowerflowError> {
    let changed = outage.apply(case)?;
    let own = M::build(&changed)?;
    match previous {
        Some(previous) => {
            let mapped = own.map(&changed, model, case, previous)?;
            own.solve(&changed, options, OpfStart::Mapped(&mapped))
        }
        None => own.solve(&changed, options, OpfStart::Cold),
    }
}

fn outages<M: BatchModel>(
    case: &Case,
    list: &[OpfOutage],
    options: &M::Options,
    start: BatchStart,
    parallel: bool,
) -> Result<OpfOutages<M::Solution>, PowerflowError> {
    let (model, base) = base::<M>(case, options)?;
    let previous =
        (start == BatchStart::Warm && M::solver(&base).converged()).then(|| M::solver(&base));
    let one = |o: &OpfOutage| outage(case, &model, previous, o, options);
    let results = if parallel {
        par_map(list, one)
    } else {
        list.iter().map(one).collect()
    };
    Ok(OpfOutages {
        warm: previous.is_some(),
        base,
        outages: results,
    })
}

/// `f` over `items` in parallel, in their order
#[cfg(feature = "rayon")]
pub(crate) fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync + Send) -> Vec<R> {
    use rayon::prelude::*;
    // Each solve's own factorizations serially: the batch is the parallelism, and theirs would
    //   only compete with it (docs/benchmarks/o5-parallel-ldlt-2026-11-20.md)
    items
        .par_iter()
        .map(|item| ybus_sparse::compressed_matrix::serially(|| f(item)))
        .collect()
}

#[cfg(not(feature = "rayon"))]
fn par_map<T, R>(_: &[T], _: impl Fn(&T) -> R) -> Vec<R> {
    unreachable!("the parallel forms need the rayon feature")
}

/// The AC optimal power flow of `case` and of each outage in `outages`, serially (see the
/// [module documentation](self)): the base solved once, cold, and each outage on its own
/// model, from `start`.
///
/// Returns the base's errors, as
/// [`ac_opf_with_options`](super::ac::ac_opf_with_options) returns them, and otherwise each
/// outage's own result.
pub fn ac_opf_outages(
    case: &Case,
    outages: &[OpfOutage],
    options: &AcOpfOptions,
    start: BatchStart,
) -> Result<OpfOutages<AcOpfSolution>, PowerflowError> {
    self::outages::<AcOpfModel>(case, outages, options, start, false)
}

/// [`ac_opf_outages`] with the outages solved in parallel (feature `rayon`), on rayon's current
/// pool: bit for bit its results, in the same order, at any thread count.
#[cfg(feature = "rayon")]
pub fn par_ac_opf_outages(
    case: &Case,
    outages: &[OpfOutage],
    options: &AcOpfOptions,
    start: BatchStart,
) -> Result<OpfOutages<AcOpfSolution>, PowerflowError> {
    self::outages::<AcOpfModel>(case, outages, options, start, true)
}

/// The DC optimal power flow of `case` and of each outage in `outages`, serially, as
/// [`ac_opf_outages`] solves the AC one.
///
/// Returns the base's errors, as [`dc_opf_with_options`](super::dc::dc_opf_with_options)
/// returns them, and otherwise each outage's own result.
pub fn dc_opf_outages(
    case: &Case,
    outages: &[OpfOutage],
    options: &DcOpfOptions,
    start: BatchStart,
) -> Result<OpfOutages<DcOpfSolution>, PowerflowError> {
    self::outages::<DcOpfModel>(case, outages, options, start, false)
}

/// [`dc_opf_outages`] with the outages solved in parallel (feature `rayon`), on rayon's current
/// pool: bit for bit its results, in the same order, at any thread count.
#[cfg(feature = "rayon")]
pub fn par_dc_opf_outages(
    case: &Case,
    outages: &[OpfOutage],
    options: &DcOpfOptions,
    start: BatchStart,
) -> Result<OpfOutages<DcOpfSolution>, PowerflowError> {
    self::outages::<DcOpfModel>(case, outages, options, start, true)
}

// The batches' types cross threads
const _: () = {
    const fn send_sync<T: Send + Sync>() {}
    send_sync::<AcOpfModel>();
    send_sync::<DcOpfModel>();
    send_sync::<AcOpfSolution>();
    send_sync::<DcOpfSolution>();
    send_sync::<AcOpfOptions>();
    send_sync::<DcOpfOptions>();
    send_sync::<OpfOutage>();
    send_sync::<Case>();
};
