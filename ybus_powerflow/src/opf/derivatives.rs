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

use ybus_sparse::compressed_matrix::CsrMatrix;
use ybus_sparse::Complex;

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

/// Joins four nb × nb blocks into the real 2nb × 2nb matrix `[[aa, av], [va, vv]]`.
fn join(blocks: [CsrMatrix<f64>; 4]) -> CsrMatrix<f64> {
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
