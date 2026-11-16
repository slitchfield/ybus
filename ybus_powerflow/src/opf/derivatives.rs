//! First and second derivatives of AC power injections and branch flows with respect to the
//! bus voltages in polar coordinates (angles `Va`, magnitudes `Vm`), which the AC optimal power
//! flow uses, in the complex matrix notation of MATPOWER Technical Note 2 (R. D. Zimmerman, "AC
//! power flows, generalized OPF costs and their derivatives using complex matrix notation",
//! 2010). The matrices come from [`AcModel`](crate::ac::AcModel): Ybus, Yf and Yt, over its
//! connected buses (in `bus_rows` order) and in-service branches.
//!
//! Every quantity here has the form `S = [C V] conj(Y V)`: the bus injections with `C = I` and
//! `Y = Ybus`, and the flows into a branch at its from (to) end with `C = Cf` (`Ct`) and
//! `Y = Yf` (`Yt`). Writing `[x]` for the diagonal matrix of a vector `x`, and `dₓV` for the
//! derivative of `V` with respect to one of its coordinates (`d_Va V = jV`, `d_Vm V = V/|V|`):
//!
//! - **First derivatives:** `∂S/∂x = [conj(Y V)] C [dₓV] + [C V] conj(Y) [conj(dₓV)]`, Technical
//!   Note 2's `dSbus_dV` and `dSbr_dV`.
//! - **Second derivatives:** a weighted sum `wᵀ S` (weights `w` real) is the bilinear form
//!   `Vᵀ M conj(V)` with `M = Cᵀ [w] conj(Y)`, whose Hessian block for coordinates `x` and `y`
//!   is
//!
//!   `[dₓV] M [conj(d_yV)] + [conj(dₓV)] Mᵀ [d_yV] + [dₓ_yV ∘ (M conj V) + conj(dₓ_yV) ∘ (Mᵀ V)]`
//!
//!   with `d_Va_Va V = −V`, `d_Va_Vm V = jV/|V|` and `d_Vm_Vm V = 0`. Its real part is the
//!   Hessian of `wᵀ Re(S)`, and its imaginary part that of `wᵀ Im(S)`: Technical Note 2's
//!   `d2Sbus_dV2` and `d2Sbr_dV2`, organized as one formula.
//! - **Squared apparent power** `A = |S|²`: `∂A/∂x = 2 Re([conj S] ∂S/∂x)`, and the Hessian of
//!   `μᵀ A` is `2 Re((∂S/∂V)ᴴ [μ] ∂S/∂V)` plus `2 Re` of the bilinear Hessian with weights
//!   `μ ∘ Re(S)` and `2 Im` of the one with weights `μ ∘ Im(S)`, Technical Note 2's `dAbr_dV`
//!   and `d2Abr_dV2`.
//!
//! Each Hessian is returned as the real `2nb × 2nb` matrix over `[Va, Vm]`.
//!
//! The functions panic if the matrices, voltages and weights don't agree in size: `C` and `Y`
//! must both be `rows × nb`, with one voltage per bus and one weight per row.
//!
//! **Plans.** [`WeightedHessianPlan`] and [`SquaredMagnitudeHessianPlan`] compute the same
//! Hessians as [`weighted_hessian`] and [`squared_magnitude_hessian`], bit for bit, from plans
//! made once for `C` and `Y`. Their patterns depend only on `C` and `Y`, so an OPF that
//! evaluates them at every iteration can keep the plans for the solve, and skip the products'
//! hashing and sorting, the merges and the stacking.

use ybus_sparse::compressed_matrix::{CsrMatrix, ProductPlan, SumPlan};
use ybus_sparse::{Complex, Scalar};

type C64 = Complex<f64>;

/// The voltages' derivatives with respect to their polar coordinates, one entry per bus.
struct PolarDerivatives {
    /// `∂V/∂Va = jV`
    va: Vec<C64>,
    /// `∂V/∂Vm = V/|V|`
    vm: Vec<C64>,
    /// `∂²V/∂Va² = −V`
    va_va: Vec<C64>,
    /// `∂²V/∂Va∂Vm = jV/|V|`
    va_vm: Vec<C64>,
}

impl PolarDerivatives {
    fn new(v: &[C64]) -> Self {
        let j = C64::new(0.0, 1.0);
        let unit: Vec<C64> = v.iter().map(|v| v / v.norm()).collect();
        PolarDerivatives {
            va: v.iter().map(|v| j * v).collect(),
            va_va: v.iter().map(|v| -v).collect(),
            va_vm: unit.iter().map(|u| j * u).collect(),
            vm: unit,
        }
    }
}

fn conj(values: &[C64]) -> Vec<C64> {
    values.iter().map(|v| v.conj()).collect()
}

/// `S = [C V] conj(Y V)` for `c` (rows × nb) and `y` (rows × nb), with its first derivatives
/// with respect to `Va` and `Vm` (complex, rows × nb).
#[derive(Clone, Debug, PartialEq)]
pub struct Flows {
    /// `S`, one per row
    pub s: Vec<C64>,
    /// `∂S/∂Va`
    pub ds_dva: CsrMatrix<C64>,
    /// `∂S/∂Vm`
    pub ds_dvm: CsrMatrix<C64>,
}

/// [`Flows`] for `S = [C V] conj(Y V)` at voltages `v`.
pub fn flows(c: &CsrMatrix<C64>, y: &CsrMatrix<C64>, v: &[C64]) -> Flows {
    let d = PolarDerivatives::new(v);
    let cv = c.mul_vec(v).expect("C has a column per bus");
    let i_conj = conj(&y.mul_vec(v).expect("Y has a column per bus"));
    let y_conj = y.conj();
    let derivative = |dv: &[C64]| {
        let first = c
            .scale_rows(&i_conj)
            .and_then(|m| m.scale_cols(dv))
            .expect("one entry per row and per bus");
        let second = y_conj
            .scale_rows(&cv)
            .and_then(|m| m.scale_cols(&conj(dv)))
            .expect("one entry per row and per bus");
        first.try_add(&second).expect("both rows × nb")
    };
    Flows {
        s: cv.iter().zip(&i_conj).map(|(cv, i)| cv * i).collect(),
        ds_dva: derivative(&d.va),
        ds_dvm: derivative(&d.vm),
    }
}

/// The complex Hessian blocks `(Va Va, Va Vm, Vm Va, Vm Vm)` of the bilinear form
/// `Vᵀ M conj(V)` at voltages `v` (see the module documentation).
fn bilinear_hessian(m: &CsrMatrix<C64>, v: &[C64]) -> [CsrMatrix<C64>; 4] {
    let d = PolarDerivatives::new(v);
    let mt = m.clone().transpose().to_csr();
    let m_conj_v = m.mul_vec(&conj(v)).expect("M is nb × nb");
    let mt_v = mt.mul_vec(v).expect("M is nb × nb");
    let block = |dx: &[C64], dy: &[C64], dxy: Option<&[C64]>| {
        let first = m
            .scale_rows(dx)
            .and_then(|b| b.scale_cols(&conj(dy)))
            .expect("one entry per bus");
        let second = mt
            .scale_rows(&conj(dx))
            .and_then(|b| b.scale_cols(dy))
            .expect("one entry per bus");
        let sum = first.try_add(&second).expect("both nb × nb");
        match dxy {
            Some(dxy) => {
                let diagonal = CsrMatrix::from_diagonal(
                    dxy.iter()
                        .zip(&m_conj_v)
                        .zip(&mt_v)
                        .map(|((d, a), b)| d * a + d.conj() * b),
                );
                sum.try_add(&diagonal).expect("both nb × nb")
            }
            None => sum,
        }
    };
    let va_va = block(&d.va, &d.va, Some(&d.va_va));
    let va_vm = block(&d.va, &d.vm, Some(&d.va_vm));
    let vm_va = va_vm.clone().transpose().to_csr();
    let vm_vm = block(&d.vm, &d.vm, None);
    [va_va, va_vm, vm_va, vm_vm]
}

/// Joins four nb × nb blocks into the 2nb × 2nb matrix `[[aa, av], [va, vv]]`.
fn join<T: Scalar>(blocks: [CsrMatrix<T>; 4]) -> CsrMatrix<T> {
    let [aa, av, va, vv] = blocks;
    let top = CsrMatrix::hstack([aa, av]).expect("both nb rows");
    let bottom = CsrMatrix::hstack([va, vv]).expect("both nb rows");
    CsrMatrix::vstack([top, bottom]).expect("both 2nb columns")
}

/// `M = Cᵀ [w] conj(Y)` for weights `w`, one per row of `c` and `y`.
fn weighted(c: &CsrMatrix<C64>, y: &CsrMatrix<C64>, w: &[f64]) -> CsrMatrix<C64> {
    let w: Vec<C64> = w.iter().map(|&w| C64::new(w, 0.0)).collect();
    let weighted_y = y.conj().scale_rows(&w).expect("one weight per row");
    c.clone()
        .transpose()
        .to_csr()
        .matmul(&weighted_y)
        .expect("Cᵀ is nb × rows, [w] conj(Y) rows × nb")
}

/// The Hessian (real, 2nb × 2nb over `[Va, Vm]`) of `pᵀ Re(S) + qᵀ Im(S)` for
/// `S = [C V] conj(Y V)` at voltages `v`, with weights `p` and `q`, one per row.
pub fn weighted_hessian(
    c: &CsrMatrix<C64>,
    y: &CsrMatrix<C64>,
    v: &[C64],
    p: &[f64],
    q: &[f64],
) -> CsrMatrix<f64> {
    let real = bilinear_hessian(&weighted(c, y, p), v);
    let imaginary = bilinear_hessian(&weighted(c, y, q), v);
    join([0, 1, 2, 3].map(|k| {
        real[k]
            .real()
            .try_add(&imaginary[k].imag())
            .expect("nb × nb")
    }))
}

/// The squared apparent power `|S|²` of each row of `flows`' `S`, and its derivatives with
/// respect to `[Va, Vm]` (real, rows × 2nb).
pub fn squared_magnitudes(flows: &Flows) -> (Vec<f64>, CsrMatrix<f64>) {
    let squared = flows.s.iter().map(|s| s.norm_sqr()).collect();
    let twice_conj: Vec<C64> = flows.s.iter().map(|s| 2.0 * s.conj()).collect();
    let derivative = |ds: &CsrMatrix<C64>| {
        ds.scale_rows(&twice_conj)
            .expect("one entry per row")
            .real()
    };
    let jacobian = CsrMatrix::hstack([derivative(&flows.ds_dva), derivative(&flows.ds_dvm)])
        .expect("both rows × nb");
    (squared, jacobian)
}

/// The Hessian (real, 2nb × 2nb over `[Va, Vm]`) of `μᵀ |S|²` for `S = [C V] conj(Y V)` at
/// voltages `v`, with `flows` the [`Flows`] there.
pub fn squared_magnitude_hessian(
    c: &CsrMatrix<C64>,
    y: &CsrMatrix<C64>,
    v: &[C64],
    flows: &Flows,
    mu: &[f64],
) -> CsrMatrix<f64> {
    // 2 Re(dSᴴ [μ] dS), from the first derivatives
    let ds =
        CsrMatrix::hstack([flows.ds_dva.clone(), flows.ds_dvm.clone()]).expect("both rows × nb");
    let mu_complex: Vec<C64> = mu.iter().map(|&m| C64::new(2.0 * m, 0.0)).collect();
    let outer = ds
        .conj()
        .transpose()
        .to_csr()
        .matmul(&ds.scale_rows(&mu_complex).expect("one weight per row"))
        .expect("2nb × rows times rows × 2nb")
        .real();
    // 2 Σ μ (Re S ∇² Re S + Im S ∇² Im S)
    let p: Vec<f64> = mu
        .iter()
        .zip(&flows.s)
        .map(|(m, s)| 2.0 * m * s.re)
        .collect();
    let q: Vec<f64> = mu
        .iter()
        .zip(&flows.s)
        .map(|(m, s)| 2.0 * m * s.im)
        .collect();
    outer
        .try_add(&weighted_hessian(c, y, v, &p, &q))
        .expect("both 2nb × 2nb")
}

// ---------------------------------------------------------------------------------------------
// Plans: the same Hessians, bit for bit, from plans made once for C and Y
// ---------------------------------------------------------------------------------------------

const ZERO: C64 = C64::new(0.0, 0.0);

/// `m`'s pattern, with each entry's value its slot plus `offset + 1`, so that the operations
/// that only move entries (transposes and stacks) can be traced: an entry with value `t` came
/// from slot `t - 1` of the values they were applied to.
fn tagged<T: Scalar>(m: &CsrMatrix<T>, offset: usize) -> CsrMatrix<usize> {
    let (rows, cols) = m.shape();
    CsrMatrix::try_from_raw_parts(
        rows,
        cols,
        m.pointers().to_vec(),
        m.indices().to_vec(),
        (offset + 1..=offset + m.num_stored()).collect(),
    )
    .expect("the same valid pattern")
}

/// The source slot of each entry of a traced matrix (see `tagged`)
fn sources(traced: &CsrMatrix<usize>) -> Vec<usize> {
    traced.values().iter().map(|&t| t - 1).collect()
}

/// `values` moved as a traced operation moved them: entry `k` is `values[sources[k]]`
fn gather<T: Copy>(values: &[T], sources: &[usize]) -> Vec<T> {
    sources.iter().map(|&s| values[s]).collect()
}

/// The row of each stored entry of `m`
fn entry_rows<T: Scalar>(m: &CsrMatrix<T>) -> Vec<usize> {
    m.pointers()
        .windows(2)
        .enumerate()
        .flat_map(|(row, bounds)| std::iter::repeat_n(row, bounds[1] - bounds[0]))
        .collect()
}

/// `m`'s pattern with every value zero, as a template for a plan's result
fn pattern<T: Scalar, U: Scalar>(m: &CsrMatrix<T>) -> CsrMatrix<U> {
    m.map(|_| U::zero())
}

/// `m`'s pattern holding `values`
fn with_values<T: Scalar>(m: &CsrMatrix<T>, values: Vec<T>) -> CsrMatrix<T> {
    let mut m = m.clone();
    m.values_mut().copy_from_slice(&values);
    m
}

/// The slot of entry `(row, col)` of `m`, if it's stored
fn slot<T: Scalar>(m: &CsrMatrix<T>, row: usize, col: usize) -> Option<usize> {
    let start = m.pointers()[row];
    m.indices()[start..m.pointers()[row + 1]]
        .binary_search(&col)
        .ok()
        .map(|offset| start + offset)
}

/// `mul_vec` on `pattern` holding `values`: each entry `0 + v₁x₁ + v₂x₂ + …` in the pattern's
/// order, as `mul_vec` computes it
fn mul_vec_values(pattern: &CsrMatrix<C64>, values: &[C64], x: &[C64]) -> Vec<C64> {
    let idx = pattern.indices();
    pattern
        .pointers()
        .windows(2)
        .map(|bounds| (bounds[0]..bounds[1]).fold(ZERO, |sum, p| sum + values[p] * x[idx[p]]))
        .collect()
}

/// How a bus's diagonal term in `bilinear_hessian` reaches its block, where `from_diagonal`
/// drops a zero
#[derive(Clone, Copy, Debug)]
enum Diagonal {
    /// `M + Mᵀ` stores `(i, i)`: the block adds the term, or keeps `M + Mᵀ`'s value where the
    /// term is zero and dropped
    Covered { block: usize, sum: usize },
    /// Only the term can store `(i, i)`: the block's pattern has it, so a zero term, which
    /// `from_diagonal` would drop, can't be planned
    Stored,
    /// Neither `M`'s nor `Mᵀ`'s row `i` stores anything, so the term is `d·0 + conj(d)·0`,
    /// which `from_diagonal` drops; only a NaN or infinite term would be kept
    Empty,
}

/// The plan of `bilinear_hessian` for one pattern of `M`
#[derive(Clone, Debug)]
struct BilinearPlan {
    m: CsrMatrix<C64>,
    m_rows: Vec<usize>,
    mt: CsrMatrix<C64>,
    mt_rows: Vec<usize>,
    /// `Mᵀ`'s entries from `M`'s
    mt_sources: Vec<usize>,
    /// `first + second`, over `M`'s and `Mᵀ`'s patterns
    sum: SumPlan,
    /// `first + second` plus the diagonal terms it can hold (`Covered` and `Stored`)
    with_diagonal: SumPlan,
    diagonal: Vec<Diagonal>,
    /// The `vm_va` block's entries from the `va_vm` block's
    vm_va_sources: Vec<usize>,
    /// The patterns of the blocks with a diagonal term (`va_va`, `va_vm`), of `vm_va`, and of
    /// `vm_vm`, which is `M + Mᵀ`'s
    block: CsrMatrix<C64>,
    vm_va: CsrMatrix<C64>,
    vm_vm: CsrMatrix<C64>,
}

impl BilinearPlan {
    fn new(m: CsrMatrix<C64>) -> Self {
        let nb = m.shape().0;
        let mt_traced = tagged(&m, 0).transpose().to_csr();
        let mt: CsrMatrix<C64> = pattern(&mt_traced);
        let sum = m.plan_sum(&mt).expect("both nb × nb");
        let sum_pattern = sum.add(&m, &mt).expect("the planned operands");
        let nonempty = |a: &CsrMatrix<C64>, i: usize| a.pointers()[i] < a.pointers()[i + 1];
        let mut kinds: Vec<Option<Diagonal>> = (0..nb)
            .map(|i| match slot(&sum_pattern, i, i) {
                Some(_) => None, // covered: its slots come from the joined pattern below
                None if nonempty(&m, i) || nonempty(&mt, i) => Some(Diagonal::Stored),
                None => Some(Diagonal::Empty),
            })
            .collect();
        let held: Vec<usize> = (0..nb)
            .filter(|&i| !matches!(kinds[i], Some(Diagonal::Empty)))
            .collect();
        let diagonal_pattern =
            CsrMatrix::<C64>::from_triplets(nb, nb, held.iter().map(|&i| (i, i, ZERO)))
                .expect("diagonal in range");
        let with_diagonal = sum_pattern
            .plan_sum(&diagonal_pattern)
            .expect("both nb × nb");
        let block = with_diagonal
            .add(&sum_pattern, &diagonal_pattern)
            .expect("the planned operands");
        for (i, kind) in kinds.iter_mut().enumerate() {
            if kind.is_none() {
                *kind = Some(Diagonal::Covered {
                    block: slot(&block, i, i).expect("held"),
                    sum: slot(&sum_pattern, i, i).expect("covered"),
                });
            }
        }
        let vm_va_traced = tagged(&block, 0).transpose().to_csr();
        BilinearPlan {
            m_rows: entry_rows(&m),
            mt_rows: entry_rows(&mt),
            mt_sources: sources(&mt_traced),
            m,
            mt,
            sum,
            with_diagonal,
            diagonal: kinds.into_iter().map(|k| k.expect("classified")).collect(),
            vm_va_sources: sources(&vm_va_traced),
            vm_va: pattern(&vm_va_traced),
            vm_vm: sum_pattern,
            block,
        }
    }

    /// `first + second` for one block, as `bilinear_hessian`'s `block` computes it
    fn sum(&self, m: &[C64], mt: &[C64], dx: &[C64], dy: &[C64]) -> Vec<C64> {
        let (dx_conj, dy_conj) = (conj(dx), conj(dy));
        // scale_rows multiplies factor · value, then scale_cols value · factor
        let first: Vec<C64> = (m.iter().zip(&self.m_rows))
            .zip(self.m.indices())
            .map(|((&value, &row), &col)| dx[row] * value * dy_conj[col])
            .collect();
        let second: Vec<C64> = (mt.iter().zip(&self.mt_rows))
            .zip(self.mt.indices())
            .map(|((&value, &row), &col)| dx_conj[row] * value * dy[col])
            .collect();
        let mut sum = vec![ZERO; self.sum.num_stored()];
        self.sum
            .add_values(&first, &second, &mut sum)
            .expect("the planned lengths");
        sum
    }

    /// One block with its diagonal terms `terms`, as `bilinear_hessian` adds them through
    /// `from_diagonal` and `try_add`, or `None` where that would store a different pattern
    fn add_diagonal(&self, sum: &[C64], terms: &[C64]) -> Option<Vec<C64>> {
        let mut held = Vec::with_capacity(terms.len());
        for (&kind, &term) in self.diagonal.iter().zip(terms) {
            match kind {
                Diagonal::Covered { .. } => held.push(term),
                Diagonal::Stored if term == ZERO => return None,
                Diagonal::Stored => held.push(term),
                Diagonal::Empty if term != ZERO => return None,
                Diagonal::Empty => {}
            }
        }
        let mut block = vec![ZERO; self.with_diagonal.num_stored()];
        self.with_diagonal
            .add_values(sum, &held, &mut block)
            .expect("the planned lengths");
        // A zero term is dropped, so try_add copies M + Mᵀ's value rather than adding 0
        for (&kind, &term) in self.diagonal.iter().zip(terms) {
            if let Diagonal::Covered { block: b, sum: s } = kind {
                if term == ZERO {
                    block[b] = sum[s];
                }
            }
        }
        Some(block)
    }

    /// The four blocks `(Va Va, Va Vm, Vm Va, Vm Vm)` of `bilinear_hessian(M, v)` for `M`'s
    /// values `m`: the first three over the block pattern with its diagonal, the last over
    /// `M + Mᵀ`'s. `None` where `bilinear_hessian` would store a different pattern.
    fn evaluate(&self, m: &[C64], v: &[C64], d: &PolarDerivatives) -> Option<[Vec<C64>; 4]> {
        let mt = gather(m, &self.mt_sources);
        let m_conj_v = mul_vec_values(&self.m, m, &conj(v));
        let mt_v = mul_vec_values(&self.mt, &mt, v);
        let terms = |dxy: &[C64]| -> Vec<C64> {
            dxy.iter()
                .zip(&m_conj_v)
                .zip(&mt_v)
                .map(|((d, a), b)| d * a + d.conj() * b)
                .collect()
        };
        let va_va = self.add_diagonal(&self.sum(m, &mt, &d.va, &d.va), &terms(&d.va_va))?;
        let va_vm = self.add_diagonal(&self.sum(m, &mt, &d.va, &d.vm), &terms(&d.va_vm))?;
        let vm_va = gather(&va_vm, &self.vm_va_sources);
        let vm_vm = self.sum(m, &mt, &d.vm, &d.vm);
        Some([va_va, va_vm, vm_va, vm_vm])
    }
}

/// [`weighted_hessian`] for one `C` and `Y`, planned once: the same Hessian, pattern and
/// values bit for bit, for any voltages and weights.
///
/// `evaluate` returns `None` in the one case the plan doesn't cover: a bus whose diagonal term
/// in the Hessian (see the module documentation) is exactly zero, where no other term stores
/// that position, so that `weighted_hessian`'s pattern lacks it. The caller then uses
/// `weighted_hessian` itself.
#[derive(Clone, Debug)]
pub struct WeightedHessianPlan {
    nb: usize,
    /// `Cᵀ`, whose values never change
    ct: CsrMatrix<C64>,
    /// `conj(Y)`'s values, and the row of each
    y_conj: CsrMatrix<C64>,
    y_rows: Vec<usize>,
    /// `M = Cᵀ ([w] conj Y)`
    product: ProductPlan,
    bilinear: BilinearPlan,
    /// The joined matrix's entries from the four blocks' values, concatenated
    join_sources: Vec<usize>,
    pattern: CsrMatrix<f64>,
}

impl WeightedHessianPlan {
    /// Plans `weighted_hessian(c, y, ·, ·, ·)`. Panics unless `c` and `y` are both `rows × nb`.
    pub fn new(c: &CsrMatrix<C64>, y: &CsrMatrix<C64>) -> Self {
        assert_eq!(c.shape(), y.shape(), "C and Y are both rows × nb");
        let nb = c.shape().1;
        let ct = c.clone().transpose().to_csr();
        let y_conj = y.conj();
        let product = ct.plan_product(&y_conj).expect("Cᵀ is nb × rows");
        let m = product
            .multiply(&ct, &y_conj)
            .expect("the planned operands");
        let bilinear = BilinearPlan::new(pattern(&m));
        // The four blocks' patterns, traced through join, each tagged with its offset in the
        //   blocks' values concatenated
        let (block, vm_va, vm_vm) = (&bilinear.block, &bilinear.vm_va, &bilinear.vm_vm);
        let n = block.num_stored();
        let traced = join([
            tagged(block, 0),
            tagged(block, n),
            tagged(vm_va, 2 * n),
            tagged(vm_vm, 3 * n),
        ]);
        WeightedHessianPlan {
            nb,
            y_rows: entry_rows(&y_conj),
            ct,
            y_conj,
            product,
            bilinear,
            join_sources: sources(&traced),
            pattern: pattern(&traced),
        }
    }

    /// The Hessian's pattern, with every value zero: what `evaluate` returns holds its values.
    pub fn pattern(&self) -> &CsrMatrix<f64> {
        &self.pattern
    }

    /// `M = Cᵀ [w] conj(Y)`'s values, as `weighted` computes them
    fn weighted(&self, w: &[f64]) -> Vec<C64> {
        let weighted_y: Vec<C64> = self
            .y_conj
            .values()
            .iter()
            .zip(&self.y_rows)
            .map(|(&value, &row)| C64::new(w[row], 0.0) * value)
            .collect();
        let mut m = vec![ZERO; self.product.num_stored()];
        self.product
            .multiply_values(self.ct.values(), &weighted_y, &mut m)
            .expect("the planned lengths");
        m
    }

    /// The Hessian's values at voltages `v` for weights `p` and `q`, in `pattern()`'s order,
    /// bit for bit `weighted_hessian(c, y, v, p, q)`'s, or `None` where its pattern differs.
    fn evaluate_values(&self, v: &[C64], p: &[f64], q: &[f64]) -> Option<Vec<f64>> {
        assert_eq!(v.len(), self.nb, "one voltage per bus");
        assert_eq!(p.len(), self.y_conj.shape().0, "one weight per row");
        assert_eq!(q.len(), p.len(), "one weight per row");
        // bilinear_hessian computes the voltages' derivatives afresh for each of its calls,
        //   always to the same values
        let d = PolarDerivatives::new(v);
        let real = self.bilinear.evaluate(&self.weighted(p), v, &d)?;
        let imaginary = self.bilinear.evaluate(&self.weighted(q), v, &d)?;
        // real[k].real() + imaginary[k].imag(), over the same pattern, so each entry's sum
        let mut concatenated = Vec::with_capacity(self.join_sources.len());
        for (re, im) in real.iter().zip(&imaginary) {
            concatenated.extend(re.iter().zip(im).map(|(re, im)| re.re + im.im));
        }
        Some(gather(&concatenated, &self.join_sources))
    }

    /// The Hessian (real, 2nb × 2nb over `[Va, Vm]`) of `pᵀ Re(S) + qᵀ Im(S)` at voltages `v`,
    /// bit for bit `weighted_hessian(c, y, v, p, q)`, or `None` where `weighted_hessian`'s
    /// pattern differs from the plan's (see the type's documentation).
    ///
    /// Panics unless there is one voltage per bus and one weight per row.
    pub fn evaluate(&self, v: &[C64], p: &[f64], q: &[f64]) -> Option<CsrMatrix<f64>> {
        self.evaluate_values(v, p, q)
            .map(|values| with_values(&self.pattern, values))
    }
}

/// [`squared_magnitude_hessian`] for one `C` and `Y`, planned once: the same Hessian, pattern
/// and values bit for bit. `evaluate` returns `None` where [`WeightedHessianPlan`] would.
#[derive(Clone, Debug)]
pub struct SquaredMagnitudeHessianPlan {
    rows: usize,
    /// The number of entries of `flows`' two derivatives, which `dS` stacks
    ds_dva_stored: usize,
    ds_dvm_stored: usize,
    /// `dS = [∂S/∂Va, ∂S/∂Vm]`'s entries from the two derivatives' values, concatenated
    ds_sources: Vec<usize>,
    ds_rows: Vec<usize>,
    /// `dSᴴ`'s entries from `dS`'s
    ds_adjoint_sources: Vec<usize>,
    /// `dSᴴ ([2μ] dS)`
    product: ProductPlan,
    weighted: WeightedHessianPlan,
    /// `Re(dSᴴ [2μ] dS)` plus the weighted Hessian
    sum: SumPlan,
    pattern: CsrMatrix<f64>,
}

impl SquaredMagnitudeHessianPlan {
    /// Plans `squared_magnitude_hessian(c, y, ·, ·, ·)`. Panics unless `c` and `y` are both
    /// `rows × nb`.
    pub fn new(c: &CsrMatrix<C64>, y: &CsrMatrix<C64>) -> Self {
        let (rows, nb) = c.shape();
        // flows' patterns depend only on C's and Y's, so any voltages give them
        let flows = flows(c, y, &vec![C64::new(1.0, 0.0); nb]);
        let (dva, dvm) = (flows.ds_dva.num_stored(), flows.ds_dvm.num_stored());
        let ds_traced = CsrMatrix::hstack([tagged(&flows.ds_dva, 0), tagged(&flows.ds_dvm, dva)])
            .expect("both rows × nb");
        let ds: CsrMatrix<C64> = pattern(&ds_traced);
        let adjoint_traced = tagged(&ds, 0).transpose().to_csr();
        let ds_adjoint: CsrMatrix<C64> = pattern(&adjoint_traced);
        let product = ds_adjoint
            .plan_product(&ds)
            .expect("2nb × rows times rows × 2nb");
        let outer: CsrMatrix<f64> = pattern(&product.multiply(&ds_adjoint, &ds).expect("planned"));
        let weighted = WeightedHessianPlan::new(c, y);
        let sum = outer.plan_sum(weighted.pattern()).expect("both 2nb × 2nb");
        let pattern = sum.add(&outer, weighted.pattern()).expect("planned");
        SquaredMagnitudeHessianPlan {
            rows,
            ds_dva_stored: dva,
            ds_dvm_stored: dvm,
            ds_sources: sources(&ds_traced),
            ds_rows: entry_rows(&ds),
            ds_adjoint_sources: sources(&adjoint_traced),
            product,
            weighted,
            sum,
            pattern,
        }
    }

    /// The Hessian's pattern, with every value zero.
    pub fn pattern(&self) -> &CsrMatrix<f64> {
        &self.pattern
    }

    /// The Hessian (real, 2nb × 2nb over `[Va, Vm]`) of `μᵀ |S|²` at voltages `v`, with `flows`
    /// the [`Flows`] there, bit for bit `squared_magnitude_hessian(c, y, v, flows, mu)`, or
    /// `None` where its pattern differs from the plan's.
    ///
    /// Panics unless `flows` is for this plan's `C` and `Y`, with one voltage per bus and one
    /// multiplier per row.
    pub fn evaluate(&self, v: &[C64], flows: &Flows, mu: &[f64]) -> Option<CsrMatrix<f64>> {
        assert_eq!(mu.len(), self.rows, "one multiplier per row");
        assert_eq!(flows.s.len(), self.rows, "flows for this plan's C and Y");
        assert_eq!(
            flows.ds_dva.num_stored(),
            self.ds_dva_stored,
            "flows for this plan"
        );
        assert_eq!(
            flows.ds_dvm.num_stored(),
            self.ds_dvm_stored,
            "flows for this plan"
        );
        // 2 Re(dSᴴ [μ] dS), from the first derivatives
        let concatenated: Vec<C64> = flows
            .ds_dva
            .values()
            .iter()
            .chain(flows.ds_dvm.values())
            .copied()
            .collect();
        let ds = gather(&concatenated, &self.ds_sources);
        let mu_complex: Vec<C64> = mu.iter().map(|&m| C64::new(2.0 * m, 0.0)).collect();
        let scaled: Vec<C64> = ds
            .iter()
            .zip(&self.ds_rows)
            .map(|(&value, &row)| mu_complex[row] * value)
            .collect();
        let adjoint: Vec<C64> = gather(&ds, &self.ds_adjoint_sources)
            .iter()
            .map(|value| value.conj())
            .collect();
        let mut outer = vec![ZERO; self.product.num_stored()];
        self.product
            .multiply_values(&adjoint, &scaled, &mut outer)
            .expect("the planned lengths");
        let outer: Vec<f64> = outer.iter().map(|value| value.re).collect();
        // 2 Σ μ (Re S ∇² Re S + Im S ∇² Im S)
        let p: Vec<f64> = mu
            .iter()
            .zip(&flows.s)
            .map(|(m, s)| 2.0 * m * s.re)
            .collect();
        let q: Vec<f64> = mu
            .iter()
            .zip(&flows.s)
            .map(|(m, s)| 2.0 * m * s.im)
            .collect();
        let weighted = self.weighted.evaluate_values(v, &p, &q)?;
        let mut values = vec![0.0; self.sum.num_stored()];
        self.sum
            .add_values(&outer, &weighted, &mut values)
            .expect("the planned lengths");
        Some(with_values(&self.pattern, values))
    }
}

/// [`flows`] for one `C` and `Y`, planned once: the same [`Flows`], patterns and values bit for
/// bit, at any voltages.
#[derive(Clone, Debug)]
pub struct FlowsPlan {
    c: CsrMatrix<C64>,
    c_rows: Vec<usize>,
    y: CsrMatrix<C64>,
    y_conj: CsrMatrix<C64>,
    y_rows: Vec<usize>,
    /// `first + second`, over `C`'s and `Y`'s patterns
    sum: SumPlan,
    pattern: CsrMatrix<C64>,
}

impl FlowsPlan {
    /// Plans `flows(c, y, ·)`. Panics unless `c` and `y` are both `rows × nb`.
    pub fn new(c: &CsrMatrix<C64>, y: &CsrMatrix<C64>) -> Self {
        assert_eq!(c.shape(), y.shape(), "C and Y are both rows × nb");
        let y_conj = y.conj();
        let sum = c.plan_sum(&y_conj).expect("both rows × nb");
        let pattern = sum.add(c, &y_conj).expect("planned");
        FlowsPlan {
            c_rows: entry_rows(c),
            c: c.clone(),
            y_rows: entry_rows(y),
            y: y.clone(),
            y_conj,
            sum,
            pattern: pattern.map(|_| ZERO),
        }
    }

    /// [`Flows`] at voltages `v`, bit for bit `flows(c, y, v)`. Panics unless there is one
    /// voltage per bus.
    pub fn evaluate(&self, v: &[C64]) -> Flows {
        let d = PolarDerivatives::new(v);
        let cv = self.c.mul_vec(v).expect("C has a column per bus");
        let i_conj = conj(&self.y.mul_vec(v).expect("Y has a column per bus"));
        let derivative = |dv: &[C64]| {
            let dv_conj = conj(dv);
            // scale_rows multiplies factor · value, then scale_cols value · factor
            let first: Vec<C64> = (self.c.values().iter().zip(&self.c_rows))
                .zip(self.c.indices())
                .map(|((&value, &row), &col)| i_conj[row] * value * dv[col])
                .collect();
            let second: Vec<C64> = (self.y_conj.values().iter().zip(&self.y_rows))
                .zip(self.y_conj.indices())
                .map(|((&value, &row), &col)| cv[row] * value * dv_conj[col])
                .collect();
            let mut values = vec![ZERO; self.sum.num_stored()];
            self.sum
                .add_values(&first, &second, &mut values)
                .expect("the planned lengths");
            with_values(&self.pattern, values)
        };
        Flows {
            s: cv.iter().zip(&i_conj).map(|(cv, i)| cv * i).collect(),
            ds_dva: derivative(&d.va),
            ds_dvm: derivative(&d.vm),
        }
    }
}

/// [`squared_magnitudes`] for the [`Flows`] of one `C` and `Y`, planned once: the same squared
/// magnitudes and Jacobian, bit for bit.
#[derive(Clone, Debug)]
pub struct SquaredMagnitudesPlan {
    ds_dva_stored: usize,
    ds_dvm_stored: usize,
    /// The row of each entry of the two derivatives, concatenated
    rows: Vec<usize>,
    /// The Jacobian's entries from the two derivatives' values, concatenated
    sources: Vec<usize>,
    pattern: CsrMatrix<f64>,
}

impl SquaredMagnitudesPlan {
    /// Plans `squared_magnitudes` for the flows of `c` and `y`. Panics unless `c` and `y` are
    /// both `rows × nb`.
    pub fn new(c: &CsrMatrix<C64>, y: &CsrMatrix<C64>) -> Self {
        let flows = FlowsPlan::new(c, y);
        let (dva, dvm) = (&flows.pattern, &flows.pattern);
        let traced = CsrMatrix::hstack([tagged(dva, 0), tagged(dvm, dva.num_stored())])
            .expect("both rows × nb");
        let mut rows = entry_rows(dva);
        rows.extend(entry_rows(dvm));
        SquaredMagnitudesPlan {
            ds_dva_stored: dva.num_stored(),
            ds_dvm_stored: dvm.num_stored(),
            rows,
            sources: sources(&traced),
            pattern: pattern(&traced),
        }
    }

    /// The Jacobian's pattern, with every value zero.
    pub fn pattern(&self) -> &CsrMatrix<f64> {
        &self.pattern
    }

    /// The squared apparent power of each row of `flows`' `S`, and its derivatives with respect
    /// to `[Va, Vm]` (real, rows × 2nb), bit for bit `squared_magnitudes(flows)`. Panics unless
    /// `flows` is for this plan's `C` and `Y`.
    pub fn evaluate(&self, flows: &Flows) -> (Vec<f64>, CsrMatrix<f64>) {
        assert_eq!(
            flows.ds_dva.num_stored(),
            self.ds_dva_stored,
            "flows for this plan"
        );
        assert_eq!(
            flows.ds_dvm.num_stored(),
            self.ds_dvm_stored,
            "flows for this plan"
        );
        let squared = flows.s.iter().map(|s| s.norm_sqr()).collect();
        let twice_conj: Vec<C64> = flows.s.iter().map(|s| 2.0 * s.conj()).collect();
        // ds.scale_rows(twice_conj).real(), for both derivatives
        let scaled: Vec<f64> = (flows.ds_dva.values().iter())
            .chain(flows.ds_dvm.values())
            .zip(&self.rows)
            .map(|(&value, &row)| (twice_conj[row] * value).re)
            .collect();
        let jacobian = with_values(&self.pattern, gather(&scaled, &self.sources));
        (squared, jacobian)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `add_diagonal` against what `bilinear_hessian` does, `sum.try_add(&from_diagonal(terms))`,
    /// with a `-0.0` in the sum where the term is zero: try_add copies it, while adding the
    /// term would make it `+0.0`. Searches of whole Hessians didn't reach this case, as the
    /// products start their sums at `+0.0`, but the plan must still copy as try_add does.
    #[test]
    fn a_zero_diagonal_term_leaves_the_sum_as_it_is() {
        // M's pattern covers (0, 0) and (1, 1); bus 2 is only reached through Mᵀ
        let m =
            CsrMatrix::<C64>::from_triplets(3, 3, vec![(0, 0, ZERO), (1, 1, ZERO), (0, 2, ZERO)])
                .unwrap();
        let plan = BilinearPlan::new(m);
        let sum_pattern = plan.vm_vm.clone();
        let sum: Vec<C64> = [
            C64::new(-0.0, -0.0),
            C64::new(1.5, -0.0),
            C64::new(-0.0, 2.0),
        ]
        .into_iter()
        .cycle()
        .take(sum_pattern.num_stored())
        .collect();
        let sum_matrix = with_values(&sum_pattern, sum.clone());
        for terms in [
            [ZERO, ZERO, C64::new(1.0, 1.0)],
            [C64::new(-0.0, 0.0), C64::new(0.5, 0.0), C64::new(-2.0, 0.0)],
            [C64::new(3.0, 0.0), ZERO, C64::new(1.0, -1.0)],
        ] {
            let planned = plan
                .add_diagonal(&sum, &terms)
                .expect("bus 2's term isn't zero");
            let direct = sum_matrix
                .try_add(&CsrMatrix::from_diagonal(terms))
                .unwrap();
            let bits = |values: &[C64]| -> Vec<(u64, u64)> {
                values
                    .iter()
                    .map(|v| (v.re.to_bits(), v.im.to_bits()))
                    .collect()
            };
            assert_eq!(plan.block.indices(), direct.indices());
            assert_eq!(bits(&planned), bits(direct.values()), "terms {terms:?}");
        }
        // Bus 2's term is the only entry at (2, 2), so a zero one changes the pattern
        assert!(plan.add_diagonal(&sum, &[ZERO, ZERO, ZERO]).is_none());
    }
}
