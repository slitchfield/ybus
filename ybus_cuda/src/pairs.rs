//! Blocks whose right-hand sides are built on the device, `c = e_plus − e_minus`, solved and
//! optionally updated there, and kept there for products and per-column maxima, so that only
//! what a caller asks for is downloaded.

use std::sync::Arc;

use cudarc::driver::{CudaSlice, CudaStream};
use sparse_mat::compressed_matrix::CsrMatrix;

use crate::factor::GpuFactor;
use crate::sparse::{Csr, Product};
use crate::CudaError;

/// A right-hand side `e_plus − e_minus` over the factorization's unknowns (`None` for no
/// entry). `plus` and `minus` must differ when both are given.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Pair {
    /// The unknown with +1, if any.
    pub plus: Option<usize>,
    /// The unknown with −1, if any.
    pub minus: Option<usize>,
}

/// One column of [`GpuFactor::pair_updates`]: with `y = A⁻¹ c` for the pair's `c`, the column's
/// result is `x + z y`, where `x = base + alpha y`, `w = inverse_weight − cᵀy` and
/// `z = cᵀx / w`: the solution of `(A − c cᵀ / inverse_weight) θ = A base + alpha c`, by
/// Sherman–Morrison, in the order the CPU's outage update computes it. A column with `w = 0`
/// is [singular](GpuBlock::singular) and gets `x`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PairUpdate {
    /// The column's right-hand side.
    pub pair: Pair,
    /// The multiple of `y` added to the base before the update.
    pub alpha: f64,
    /// `1 / b` for the update's weight `b`.
    pub inverse_weight: f64,
}

/// A column's largest `|value| · weight`, from [`GpuBlock::maxima`] or
/// [`GpuBlock::product_maxima`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ColumnMax {
    /// The row of the largest product, the smallest on a tie; `None` if every product is zero.
    pub row: Option<usize>,
    /// The largest product (0.0 if every product is zero). Values that aren't finite take no
    /// part.
    pub value: f64,
    /// Whether any value in the column isn't finite.
    pub not_finite: bool,
}

/// A sparse matrix whose columns are a factorization's unknowns, uploaded once, to apply to the
/// columns of a [`GpuBlock`]: `M θ + offset` for each column `θ`, optionally with one entry per
/// column set to zero, downloaded whole or reduced to its largest `|value| · weight`.
pub struct GpuOperator {
    stream: Arc<CudaStream>,
    // Declared before the matrix, so it's dropped before it
    product: Option<Product>,
    matrix: Csr,
    rows: usize,
    columns: usize,
    offset: CudaSlice<f64>,
    weight: CudaSlice<f64>,
}

impl GpuOperator {
    /// The number of rows of the matrix.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Sets the offset added to every product (initially zero): one value per row.
    pub fn set_offset(&mut self, offset: &[f64]) -> Result<(), CudaError> {
        assert_eq!(offset.len(), self.rows, "one offset per row");
        self.stream.memcpy_htod(offset, &mut self.offset)?;
        Ok(())
    }

    /// Sets the weights [`GpuBlock::product_maxima`] multiplies each row's `|value|` by
    /// (initially one): one value per row.
    pub fn set_weights(&mut self, weights: &[f64]) -> Result<(), CudaError> {
        assert_eq!(weights.len(), self.rows, "one weight per row");
        self.stream.memcpy_htod(weights, &mut self.weight)?;
        Ok(())
    }
}

/// A block of solutions on the device, from [`GpuFactor::pair_solves`] or
/// [`GpuFactor::pair_updates`], one column per pair. It borrows the factor until dropped.
pub struct GpuBlock<'f, 'g> {
    factor: &'f mut GpuFactor<'g>,
    columns: usize,
    /// Whether the result is in the block's second buffer (after an update)
    updated: bool,
    singular: Vec<bool>,
}

/// A pair's entries as rows of a block (`-1` for none), padded to `width`
fn rows(pairs: impl Iterator<Item = Pair>, map: &[usize], width: usize) -> (Vec<i32>, Vec<i32>) {
    let row = |unknown: Option<usize>| unknown.map_or(-1, |u| map[u] as i32);
    let (mut plus, mut minus): (Vec<i32>, Vec<i32>) =
        pairs.map(|pair| (row(pair.plus), row(pair.minus))).unzip();
    plus.resize(width, -1);
    minus.resize(width, -1);
    (plus, minus)
}

impl<'g> GpuFactor<'g> {
    /// Uploads `matrix` as an operator on this factorization's solutions: it must have one
    /// column per unknown ([`GpuFactor::dim`]). Its columns are renumbered into the order the
    /// device's solutions come in, so products need no permuting.
    pub fn operator(&self, matrix: &CsrMatrix<f64>) -> Result<GpuOperator, CudaError> {
        let (rows, columns) = matrix.shape();
        assert_eq!(columns, self.n, "one column per unknown");
        let (pointers, indices, values) = (matrix.pointers(), matrix.indices(), matrix.values());
        let mut new_indices = Vec::with_capacity(indices.len());
        let mut new_values = Vec::with_capacity(values.len());
        let mut entries = vec![];
        for row in 0..rows {
            entries.clear();
            entries.extend(
                (pointers[row]..pointers[row + 1])
                    .map(|p| (self.solution_row[indices[p]], values[p])),
            );
            entries.sort_unstable_by_key(|&(column, _)| column);
            new_indices.extend(entries.iter().map(|&(column, _)| column));
            new_values.extend(entries.iter().map(|&(_, value)| value));
        }
        let gpu = self.gpu;
        let matrix = Csr::upload(
            &gpu.context,
            &gpu.stream,
            rows,
            columns,
            pointers,
            &new_indices,
            &new_values,
            None,
        )?;
        Ok(GpuOperator {
            stream: gpu.stream.clone(),
            product: None,
            matrix,
            rows,
            columns,
            offset: gpu.stream.alloc_zeros(rows.max(1))?,
            weight: gpu.stream.clone_htod(&vec![1.0; rows.max(1)])?,
        })
    }

    /// Checks a block's pairs and makes the block, with the right-hand sides scattered into its
    /// first buffer and solved.
    fn solve_pairs(
        &mut self,
        pairs: impl Iterator<Item = Pair> + Clone,
    ) -> Result<usize, CudaError> {
        let n = self.n;
        let count = pairs.clone().count();
        assert!(
            count <= self.block_width(),
            "at most block_width() ({}) pairs, not {count}",
            self.block_width()
        );
        for pair in pairs.clone() {
            assert!(
                pair.plus.is_none_or(|u| u < n) && pair.minus.is_none_or(|u| u < n),
                "unknowns below dim() ({n})"
            );
            assert!(
                pair.plus.is_none() || pair.plus != pair.minus,
                "a pair's two unknowns must differ"
            );
        }
        self.ensure_block()?;
        let width = self.block_width();
        let gpu = self.gpu;
        let stream = &gpu.stream;
        let (plus, minus) = rows(pairs, &self.rhs_row, width);
        let (plus, minus) = (stream.clone_htod(&plus)?, stream.clone_htod(&minus)?);
        let block = self.block.as_mut().expect("made above");
        let first = &mut block.first().data;
        stream.memset_zeros(first)?;
        gpu.kernels
            .scatter_pairs(stream, first, n, width, &plus, &minus)?;
        block.solve(&gpu.sparse, stream, &self.lower, &self.upper)?;
        Ok(count)
    }

    /// Solves `A y = e_plus − e_minus` for each pair, at most [`GpuFactor::block_width`] of
    /// them, keeping the solutions on the device.
    ///
    /// Panics for more pairs than the block width, an unknown past [`GpuFactor::dim`], or a
    /// pair whose two unknowns are the same.
    pub fn pair_solves(&mut self, pairs: &[Pair]) -> Result<GpuBlock<'_, 'g>, CudaError> {
        let columns = self.solve_pairs(pairs.iter().copied())?;
        Ok(GpuBlock {
            factor: self,
            columns,
            updated: false,
            singular: vec![false; columns],
        })
    }

    /// [`GpuFactor::pair_solves`] for the updates' pairs, then each column's rank-one update of
    /// `base` (one value per unknown), as [`PairUpdate`] describes, keeping the results on the
    /// device.
    ///
    /// Panics as `pair_solves` does, or if `base` doesn't have one value per unknown.
    pub fn pair_updates(
        &mut self,
        base: &[f64],
        updates: &[PairUpdate],
    ) -> Result<GpuBlock<'_, 'g>, CudaError> {
        let n = self.n;
        assert_eq!(base.len(), n, "one base value per unknown");
        let columns = self.solve_pairs(updates.iter().map(|update| update.pair))?;
        let width = self.block_width();
        let gpu = self.gpu;
        let stream = &gpu.stream;
        let base: Vec<f64> = self.scatter.iter().map(|&u| base[u]).collect();
        let (plus, minus) = rows(updates.iter().map(|u| u.pair), &self.solution_row, width);
        let mut alpha: Vec<f64> = updates.iter().map(|u| u.alpha).collect();
        let mut inverse_weight: Vec<f64> = updates.iter().map(|u| u.inverse_weight).collect();
        alpha.resize(width, 0.0);
        inverse_weight.resize(width, 1.0);
        let base = stream.clone_htod(&base)?;
        let (plus, minus) = (stream.clone_htod(&plus)?, stream.clone_htod(&minus)?);
        let alpha = stream.clone_htod(&alpha)?;
        let inverse_weight = stream.clone_htod(&inverse_weight)?;
        let mut singular = stream.alloc_zeros::<i32>(width)?;
        let block = self.block.as_mut().expect("made by solve_pairs");
        let (first, second) = block.both();
        gpu.kernels.pair_update(
            stream,
            &first.data,
            &mut second.data,
            n,
            width,
            &base,
            &plus,
            &minus,
            &alpha,
            &inverse_weight,
            &mut singular,
        )?;
        let singular = stream.clone_dtoh(&singular)?;
        Ok(GpuBlock {
            factor: self,
            columns,
            updated: true,
            singular: singular[..columns].iter().map(|&flag| flag != 0).collect(),
        })
    }
}

impl GpuBlock<'_, '_> {
    /// The number of columns.
    pub fn columns(&self) -> usize {
        self.columns
    }

    /// For each column of [`GpuFactor::pair_updates`], whether its update was singular
    /// (`w = 0`); all false for [`GpuFactor::pair_solves`].
    pub fn singular(&self) -> &[bool] {
        &self.singular
    }

    /// Downloads the solutions, one per column, indexed by unknown.
    pub fn solutions(&mut self) -> Result<Vec<Vec<f64>>, CudaError> {
        let factor = &mut *self.factor;
        let (n, columns) = (factor.n, self.columns);
        let gpu = factor.gpu;
        let block = factor.block.as_mut().expect("made by the block's solve");
        let (first, second) = block.both();
        let result = if self.updated { second } else { first };
        let host = factor.staging.get(gpu, n * columns)?;
        gpu.stream
            .memcpy_dtoh(&result.data.slice(..n * columns), host)?;
        gpu.stream.synchronize()?;
        Ok(host
            .chunks_exact(n)
            .map(|column| {
                let mut x = vec![0.0; n];
                for (&to, &value) in factor.scatter.iter().zip(column) {
                    x[to] = value;
                }
                x
            })
            .collect())
    }

    /// Each column's largest `|value|`, with its unknown and whether any value isn't finite.
    pub fn maxima(&mut self) -> Result<Vec<ColumnMax>, CudaError> {
        let factor = &mut *self.factor;
        let (n, width, columns) = (factor.n, factor.block_width(), self.columns);
        let gpu = factor.gpu;
        let stream = &gpu.stream;
        let block = factor.block.as_mut().expect("made by the block's solve");
        let (first, second) = block.both();
        let result = if self.updated { second } else { first };
        let ones = stream.clone_htod(&vec![1.0; n])?;
        let mut maxima = column_maxima(gpu, &result.data, n, width, &ones, columns)?;
        // The rows are the device's order; report unknowns
        for max in &mut maxima {
            max.row = max.row.map(|row| factor.scatter[row]);
        }
        Ok(maxima)
    }

    /// Applies `operator` (`M θ + offset`) to every column, then sets column `k`'s entry at
    /// `zero[k]` (if given) to zero, leaving the products on the device.
    fn apply(
        &mut self,
        operator: &mut GpuOperator,
        zero: &[Option<usize>],
    ) -> Result<(), CudaError> {
        let factor = &mut *self.factor;
        assert_eq!(
            operator.columns, factor.n,
            "an operator of this factorization"
        );
        assert_eq!(
            zero.len(),
            self.columns,
            "one entry to zero (or none) per column"
        );
        assert!(
            zero.iter().flatten().all(|&row| row < operator.rows),
            "rows to zero below the operator's rows"
        );
        let (width, rows) = (factor.block_width(), operator.rows);
        let gpu = factor.gpu;
        let stream = &gpu.stream;
        let block = factor.block.as_mut().expect("made by the block's solve");
        let (first, second) = block.both();
        let input = if self.updated { second } else { first };
        if operator
            .product
            .as_ref()
            .is_none_or(|product| product.width() != width)
        {
            operator.product = None;
            operator.product = Some(Product::new(
                &gpu.context,
                stream,
                &gpu.sparse,
                &operator.matrix,
                input,
                rows,
                width,
            )?);
        }
        let product = operator.product.as_mut().expect("made above");
        gpu.kernels
            .broadcast(stream, &mut product.out.data, rows, width, &operator.offset)?;
        product.run(&gpu.sparse, stream, &operator.matrix, input)?;
        let mut zero_rows: Vec<i32> = zero
            .iter()
            .map(|row| row.map_or(-1, |r| r as i32))
            .collect();
        zero_rows.resize(width, -1);
        let zero_rows = stream.clone_htod(&zero_rows)?;
        gpu.kernels
            .zero_entries(stream, &mut product.out.data, rows, width, &zero_rows)?;
        Ok(())
    }

    /// Downloads `operator`'s product with every column (see [`GpuOperator`]), with column `k`'s
    /// entry at `zero[k]` (if given) set to zero: one vector of the operator's rows per column.
    ///
    /// Panics unless `operator` was made by this block's factor and `zero` has one entry per
    /// column, each below the operator's rows.
    pub fn product(
        &mut self,
        operator: &mut GpuOperator,
        zero: &[Option<usize>],
    ) -> Result<Vec<Vec<f64>>, CudaError> {
        self.apply(operator, zero)?;
        let factor = &mut *self.factor;
        let (rows, columns) = (operator.rows, self.columns);
        let gpu = factor.gpu;
        let out = &operator.product.as_ref().expect("made by apply").out;
        let host = factor.staging.get(gpu, rows * columns)?;
        gpu.stream
            .memcpy_dtoh(&out.data.slice(..rows * columns), host)?;
        gpu.stream.synchronize()?;
        Ok(host
            .chunks_exact(rows.max(1))
            .take(columns)
            .map(<[f64]>::to_vec)
            .collect())
    }

    /// [`GpuBlock::product`] reduced on the device: each column's largest `|value| · weight`
    /// (see [`GpuOperator::set_weights`]) and its row.
    pub fn product_maxima(
        &mut self,
        operator: &mut GpuOperator,
        zero: &[Option<usize>],
    ) -> Result<Vec<ColumnMax>, CudaError> {
        self.apply(operator, zero)?;
        let factor = &*self.factor;
        let width = factor.block_width();
        let out = &operator.product.as_ref().expect("made by apply").out;
        column_maxima(
            factor.gpu,
            &out.data,
            operator.rows,
            width,
            &operator.weight,
            self.columns,
        )
    }
}

/// The first `columns` columns' maxima of a `rows`×`width` block on the device.
fn column_maxima(
    gpu: &crate::Gpu,
    block: &CudaSlice<f64>,
    rows: usize,
    width: usize,
    weight: &CudaSlice<f64>,
    columns: usize,
) -> Result<Vec<ColumnMax>, CudaError> {
    let stream = &gpu.stream;
    let mut best = stream.alloc_zeros::<f64>(width)?;
    let mut best_row = stream.alloc_zeros::<i32>(width)?;
    let mut not_finite = stream.alloc_zeros::<i32>(width)?;
    gpu.kernels.column_max(
        stream,
        block,
        rows,
        width,
        weight,
        &mut best,
        &mut best_row,
        &mut not_finite,
    )?;
    let best = stream.clone_dtoh(&best)?;
    let best_row = stream.clone_dtoh(&best_row)?;
    let not_finite = stream.clone_dtoh(&not_finite)?;
    Ok((0..columns)
        .map(|k| ColumnMax {
            row: usize::try_from(best_row[k]).ok(),
            value: best[k],
            not_finite: not_finite[k] != 0,
        })
        .collect())
}
