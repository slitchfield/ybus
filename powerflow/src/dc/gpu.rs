//! The DC power flow on the GPU (feature `cuda`): B′'s factorization, made on the CPU, uploaded
//! once through `ybus_cuda` and solved there for many right-hand sides at a time.
//!
//! - [`DcGpuStudy`] ([`DcStudy::gpu`]): batches of cases.
//! - [`DcGpuSensitivity`] ([`DcSensitivity::gpu`]): PTDF columns and transfers, whose right-hand
//!   sides are built on the device and whose flows are computed there.
//!
//! The results agree with the CPU's to within rounding, not bit for bit: the GPU's triangular
//! solves order their arithmetic differently. Everything else (the right-hand sides, the checks
//! and the flows and slack of each solution) is computed on the host exactly as the CPU path
//! computes it.

use sparse_mat::compressed_matrix::CsrMatrix;
use ybus_cuda::{CudaError, Gpu, GpuFactor, GpuOperator, Pair};

use super::{DcModel, DcSensitivity, DcSolution, DcStudy, Factor};
use crate::case::Case;
use crate::powerflow_error::PowerflowError;

/// `PowerflowError::Gpu` for a GPU failure. (Not a `From` impl: a second conversion into
/// `PowerflowError` would leave `?` in closures elsewhere without a type to infer.)
pub(super) fn gpu_failed(err: CudaError) -> PowerflowError {
    PowerflowError::Gpu {
        reason: err.to_string(),
    }
}

/// Uploads a factorization of B′ (or of the reduced Bbus a sensitivity factors).
fn upload<'g>(gpu: &'g Gpu, factor: &Factor) -> Result<GpuFactor<'g>, PowerflowError> {
    match factor {
        Factor::Cholesky(cholesky) => GpuFactor::from_cholesky(gpu, cholesky),
        Factor::Lu(lu) => GpuFactor::from_lu(gpu, lu),
    }
    .map_err(gpu_failed)
}

/// Bf over the unknowns, `position[bus]` being an internal bus's unknown (`None` for a
/// reference bus, whose columns are left out), uploaded as an operator on `factor`'s solutions:
/// its product with the unknowns is every in-service branch's flow less the reference buses'
/// part and the phase shifters'.
fn flow_operator(
    factor: &GpuFactor<'_>,
    model: &DcModel,
    position: &[Option<usize>],
) -> Result<GpuOperator, PowerflowError> {
    let unknowns = position.iter().flatten().count();
    let triplets = model
        .bf
        .iter()
        .filter_map(|(k, bus, value)| Some((k, position[bus]?, value)));
    let bf = CsrMatrix::from_triplets(model.bf.shape().0, unknowns, triplets)
        .expect("branches and unknowns in range");
    factor.operator(&bf).map_err(gpu_failed)
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
        let factor = self.factor.as_ref().map(|f| upload(gpu, f)).transpose()?;
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

/// A [`DcSensitivity`]'s factorization on the GPU, from [`DcSensitivity::gpu`], with Bf as an
/// operator on its solutions, so PTDF columns are built, solved and turned into flows on the
/// device. It borrows the sensitivity and the [`Gpu`].
pub struct DcGpuSensitivity<'s, 'g> {
    sensitivity: &'s DcSensitivity<'s>,
    /// The factorization and Bf; `None` when there is nothing to solve (a single bus)
    device: Option<(GpuFactor<'g>, GpuOperator)>,
    /// The distributed slack's shift `H w` over the in-service branches, negated: the
    /// operator's offset for PTDF columns
    negated_shift: Option<Vec<f64>>,
    /// Whether the operator's offset is now `negated_shift` (otherwise zero)
    shifted: bool,
}

impl DcSensitivity<'_> {
    /// Uploads the sensitivity's factorization and Bf to `gpu`, for [`DcGpuSensitivity`]'s
    /// PTDF columns and transfers. With a distributed slack, its shift is computed here, on
    /// the CPU, as [`DcSensitivity::ptdf_columns`] computes it.
    ///
    /// Returns `PowerflowError::Gpu` if the upload fails.
    pub fn gpu<'s, 'g>(&'s self, gpu: &'g Gpu) -> Result<DcGpuSensitivity<'s, 'g>, PowerflowError> {
        let device = match &self.factor {
            Some(factor) => {
                let factor = upload(gpu, factor)?;
                let operator = flow_operator(&factor, self.model, &self.position)?;
                Some((factor, operator))
            }
            None => None,
        };
        let negated_shift = self.distributed_shift().map(|shift| {
            self.model
                .topology
                .branch_rows
                .iter()
                .map(|&row| -shift[row])
                .collect()
        });
        Ok(DcGpuSensitivity {
            sensitivity: self,
            device,
            negated_shift,
            shifted: false,
        })
    }
}

impl DcGpuSensitivity<'_, '_> {
    /// The sensitivity whose factorization this is.
    pub fn sensitivity(&self) -> &DcSensitivity<'_> {
        self.sensitivity
    }

    /// How many columns the GPU solves together; 0 when there is nothing to solve.
    pub fn block_width(&self) -> usize {
        self.device
            .as_ref()
            .map_or(0, |(factor, _)| factor.block_width())
    }

    /// Sets how many columns the GPU solves together (at least 1). Each costs two vectors of
    /// the factorization's order and one of the in-service branches on the device.
    pub fn set_block_width(&mut self, width: usize) {
        if let Some((factor, _)) = &mut self.device {
            factor.set_block_width(width);
        }
    }

    /// [`DcSensitivity::ptdf_columns`] on the GPU: the same columns to within rounding, with
    /// the same errors.
    pub fn ptdf_columns(&mut self, buses: &[usize]) -> Result<Vec<Vec<f64>>, PowerflowError> {
        let sensitivity = self.sensitivity;
        if self.device.is_none() {
            return sensitivity.ptdf_columns(buses);
        }
        let internal = sensitivity.internal_buses(buses)?;
        let pairs: Vec<Pair> = internal
            .iter()
            .map(|bus| Pair {
                plus: bus.and_then(|bus| sensitivity.position[bus]),
                minus: None,
            })
            .collect();
        self.set_shifted(true)?;
        let mut columns = self.flows(&pairs)?;
        // As on the CPU, a bus that isn't connected gets zeros, without the shift
        for (column, bus) in columns.iter_mut().zip(&internal) {
            if bus.is_none() {
                column.fill(0.0);
            }
        }
        Ok(columns)
    }

    /// [`DcSensitivity::ptdf_transfers`] on the GPU: the same columns to within rounding, with
    /// the same errors.
    pub fn ptdf_transfers(
        &mut self,
        transfers: &[(usize, usize)],
    ) -> Result<Vec<Vec<f64>>, PowerflowError> {
        let sensitivity = self.sensitivity;
        if self.device.is_none() {
            return sensitivity.ptdf_transfers(transfers);
        }
        let from: Vec<usize> = transfers.iter().map(|&(from, _)| from).collect();
        let to: Vec<usize> = transfers.iter().map(|&(_, to)| to).collect();
        let (from, to) = (
            sensitivity.internal_buses(&from)?,
            sensitivity.internal_buses(&to)?,
        );
        let position = |bus: Option<usize>| bus.and_then(|bus| sensitivity.position[bus]);
        let pairs: Vec<Pair> = from
            .into_iter()
            .zip(to)
            .map(|(from, to)| {
                let (plus, minus) = (position(from), position(to));
                // A transfer from a bus to itself injects nothing
                if plus.is_some() && plus == minus {
                    Pair {
                        plus: None,
                        minus: None,
                    }
                } else {
                    Pair { plus, minus }
                }
            })
            .collect();
        self.set_shifted(false)?;
        self.flows(&pairs)
    }

    /// Sets the operator's offset to the negated distributed shift, or to zero.
    fn set_shifted(&mut self, shifted: bool) -> Result<(), PowerflowError> {
        let Some(negated_shift) = &self.negated_shift else {
            return Ok(());
        };
        if self.shifted != shifted {
            let (_, operator) = self.device.as_mut().expect("checked by the caller");
            if shifted {
                operator.set_offset(negated_shift)
            } else {
                operator.set_offset(&vec![0.0; negated_shift.len()])
            }
            .map_err(gpu_failed)?;
            self.shifted = shifted;
        }
        Ok(())
    }

    /// Every case branch's flow (zero for one out of service) for each pair's injections, a
    /// block at a time.
    fn flows(&mut self, pairs: &[Pair]) -> Result<Vec<Vec<f64>>, PowerflowError> {
        let topology = &self.sensitivity.model.topology;
        let (factor, operator) = self.device.as_mut().expect("checked by the caller");
        let width = factor.block_width();
        let mut columns = Vec::with_capacity(pairs.len());
        for chunk in pairs.chunks(width) {
            let mut block = factor.pair_solves(chunk).map_err(gpu_failed)?;
            let flows = block
                .product(operator, &vec![None; chunk.len()])
                .map_err(gpu_failed)?;
            if topology.branch_rows.len() == topology.case_shape.2 {
                // Every branch is in service, so the in-service rows are the case's
                columns.extend(flows);
            } else {
                for in_service in flows {
                    let mut column = vec![0.0; topology.case_shape.2];
                    for (&row, flow) in topology.branch_rows.iter().zip(in_service) {
                        column[row] = flow;
                    }
                    columns.push(column);
                }
            }
        }
        Ok(columns)
    }
}
