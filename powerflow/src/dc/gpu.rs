//! The DC power flow on the GPU (feature `cuda`): B′'s factorization, made on the CPU, uploaded
//! once through `ybus_cuda` and solved there for many right-hand sides at a time.
//!
//! The results agree with the CPU's to within rounding, not bit for bit: the GPU's triangular
//! solves order their arithmetic differently. Everything else (the right-hand sides, the checks
//! and the flows and slack of each solution) is computed on the host exactly as the CPU path
//! computes it.

use ybus_cuda::{CudaError, Gpu, GpuFactor};

use super::{DcSolution, DcStudy, Factor};
use crate::case::Case;
use crate::powerflow_error::PowerflowError;

/// `PowerflowError::Gpu` for a GPU failure. (Not a `From` impl: a second conversion into
/// `PowerflowError` would leave `?` in closures elsewhere without a type to infer.)
pub(super) fn gpu_failed(err: CudaError) -> PowerflowError {
    PowerflowError::Gpu {
        reason: err.to_string(),
    }
}

/// A [`DcStudy`]'s factorization of B′ on the GPU, from [`DcStudy::gpu`], for solving many
/// cases at once. It borrows the study and the [`Gpu`].
pub struct DcGpuStudy<'s, 'g> {
    study: &'s DcStudy<'s>,
    /// `None` when the study has nothing to solve (every connected bus is a reference bus)
    factor: Option<GpuFactor<'g>>,
}

impl DcStudy<'_> {
    /// Uploads B′'s factorization to `gpu`, for [`DcGpuStudy`]'s batched solves.
    ///
    /// Returns `PowerflowError::Gpu` if the upload fails.
    pub fn gpu<'s, 'g>(&'s self, gpu: &'g Gpu) -> Result<DcGpuStudy<'s, 'g>, PowerflowError> {
        let factor = match &self.factor {
            Some(Factor::Cholesky(cholesky)) => {
                Some(GpuFactor::from_cholesky(gpu, cholesky).map_err(gpu_failed)?)
            }
            Some(Factor::Lu(lu)) => Some(GpuFactor::from_lu(gpu, lu).map_err(gpu_failed)?),
            None => None,
        };
        Ok(DcGpuStudy {
            study: self,
            factor,
        })
    }
}

impl DcGpuStudy<'_, '_> {
    /// The study whose factorization this is.
    pub fn study(&self) -> &DcStudy<'_> {
        self.study
    }

    /// How many right-hand sides the GPU solves together (see `ybus_cuda`'s
    /// `GpuFactor::block_width`); 0 when there is nothing to solve.
    pub fn block_width(&self) -> usize {
        self.factor.as_ref().map_or(0, GpuFactor::block_width)
    }

    /// Sets how many right-hand sides the GPU solves together (at least 1). Each costs two
    /// vectors of B′'s order on the device.
    pub fn set_block_width(&mut self, width: usize) {
        if let Some(factor) = &mut self.factor {
            factor.set_block_width(width);
        }
    }

    /// [`DcStudy::solve`] for each of `cases`, with every right-hand side solved on the GPU,
    /// [`DcGpuStudy::block_width`] at a time: the same results to within rounding, in the
    /// same order. The same case may appear more than once.
    ///
    /// The outer result is `PowerflowError::Gpu` if the GPU fails; each inner one holds a
    /// case's own errors, those of [`DcStudy::solve`]. Every solution is held at once, so
    /// callers with many large cases pass them a block at a time.
    pub fn solve_many(
        &mut self,
        cases: &[&Case],
    ) -> Result<Vec<Result<DcSolution, PowerflowError>>, PowerflowError> {
        let study = self.study;
        let Some(factor) = &mut self.factor else {
            return Ok(cases.iter().map(|case| study.solve(case)).collect());
        };
        let b_ref = study.reference_coupling();
        let width = factor.block_width();
        let mut results = Vec::with_capacity(cases.len());
        for chunk in cases.chunks(width) {
            // The cases that match the model, with their injections, angles and right-hand
            //   sides; the others' errors go straight into their place
            let mut scheduled = Vec::with_capacity(chunk.len());
            let mut rhs = Vec::with_capacity(chunk.len());
            let mut chunk_results: Vec<Option<Result<DcSolution, PowerflowError>>> =
                Vec::with_capacity(chunk.len());
            for case in chunk {
                match study.model.topology.check_case(case) {
                    Ok(()) => {
                        let (pbus, theta) = study.scheduled(case);
                        rhs.push(study.rhs_with(&b_ref, &pbus, &theta));
                        scheduled.push((pbus, theta));
                        chunk_results.push(None);
                    }
                    Err(err) => chunk_results.push(Some(Err(err))),
                }
            }
            let solved = factor.solve_many(&rhs).map_err(gpu_failed)?;
            drop(rhs);
            let mut solved = scheduled.into_iter().zip(solved);
            for (case, result) in chunk.iter().zip(&mut chunk_results) {
                if result.is_some() {
                    continue;
                }
                let ((pbus, mut theta), x) = solved.next().expect("one per matching case");
                *result = Some(
                    study
                        .place(case, &mut theta, &x)
                        .map(|()| study.finish(case, &theta, &pbus, &[])),
                );
            }
            results.extend(chunk_results.into_iter().map(|r| r.expect("filled above")));
        }
        Ok(results)
    }
}
