//! First and second derivatives of AC power injections and branch flows with respect to the
//! bus voltages in polar coordinates (angles `Va`, magnitudes `Vm`), for the AC optimal power
//! flow, in the complex matrix notation of MATPOWER Technical Note 2 (R. D. Zimmerman, "AC power
//! flows, generalized OPF costs and their derivatives using complex matrix notation", 2010).
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

use sparse_mat::compressed_matrix::CsrMatrix;
use sparse_mat::Complex;

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
pub(crate) struct Flows {
    /// `S`, one per row
    pub(crate) s: Vec<C64>,
    /// `∂S/∂Va`
    pub(crate) ds_dva: CsrMatrix<C64>,
    /// `∂S/∂Vm`
    pub(crate) ds_dvm: CsrMatrix<C64>,
}

/// [`Flows`] for `S = [C V] conj(Y V)` at voltages `v`.
pub(crate) fn flows(c: &CsrMatrix<C64>, y: &CsrMatrix<C64>, v: &[C64]) -> Flows {
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
pub(crate) fn weighted_hessian(
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

/// The squared apparent power `|S|²` of each row of `S = [C V] conj(Y V)` at voltages `v`, and
/// its derivatives with respect to `[Va, Vm]` (real, rows × 2nb).
pub(crate) fn squared_magnitudes(flows: &Flows) -> (Vec<f64>, CsrMatrix<f64>) {
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
pub(crate) fn squared_magnitude_hessian(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ac::Admittances;
    use crate::case::Case;
    use crate::topology::Network;

    /// Four buses with a transformer (tap and phase shift), line charging and bus shunts.
    const CASE: &str = "function mpc = derivatives
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.02 0 135 1 1.1 0.9;
2 1 50 20 2 5 1 0.98 -3 135 1 1.1 0.9;
3 2 30 10 0 -3 1 1.01 4 135 1 1.1 0.9;
4 1 40 15 0 0 1 0.97 -6 135 1 1.1 0.9;
];
mpc.gen = [
1 0 0 100 -100 1 100 1 200 0 0 0 0 0 0 0 0 0 0 0 0;
3 40 0 100 -100 1 100 1 200 0 0 0 0 0 0 0 0 0 0 0 0;
];
mpc.branch = [
1 2 0.01 0.1 0.02 0 0 0 0 0 1 -360 360;
2 3 0.02 0.15 0.01 0 0 0 0.98 5 1 -360 360;
3 4 0.015 0.12 0.03 0 0 0 0 0 1 -360 360;
1 4 0.03 0.2 0 0 0 0 1.03 -2 1 -360 360;
2 4 0.02 0.18 0.01 0 0 0 0 0 1 -360 360;
];
";

    struct Fixture {
        v: Vec<C64>,
        ybus: CsrMatrix<C64>,
        yf: CsrMatrix<C64>,
        cf: CsrMatrix<C64>,
        identity: CsrMatrix<C64>,
    }

    fn fixture() -> Fixture {
        let case = Case::parse(CASE).unwrap();
        let network = Network::new(&case, |_, _| Ok(())).unwrap();
        let y = Admittances::new(
            &case,
            &network.bus_rows,
            &network.branch_rows,
            &network.ends,
        );
        let nb = network.bus_rows.len();
        let v = network
            .bus_rows
            .iter()
            .map(|&row| C64::from_polar(case.bus[row].vm_pu, case.bus[row].va_deg.to_radians()))
            .collect();
        let cf = CsrMatrix::<C64>::from_triplets(
            network.ends.len(),
            nb,
            network
                .ends
                .iter()
                .enumerate()
                .map(|(k, &(f, _))| (k, f, C64::new(1.0, 0.0))),
        )
        .unwrap();
        Fixture {
            v,
            ybus: y.ybus,
            yf: y.yf,
            cf,
            identity: CsrMatrix::identity(nb),
        }
    }

    /// The voltages with the polar coordinate `k` (angles first, then magnitudes) moved by `h`.
    fn moved(v: &[C64], k: usize, h: f64) -> Vec<C64> {
        let nb = v.len();
        let mut v = v.to_vec();
        let (m, a) = (v[k % nb].norm(), v[k % nb].arg());
        v[k % nb] = if k < nb {
            C64::from_polar(m, a + h)
        } else {
            C64::from_polar(m + h, a)
        };
        v
    }

    fn dense(m: &CsrMatrix<f64>) -> Vec<Vec<f64>> {
        let (rows, cols) = m.shape();
        let mut out = vec![vec![0.0; cols]; rows];
        for (i, j, value) in m.iter() {
            out[i][j] += value;
        }
        out
    }

    fn assert_close(actual: f64, expected: f64, what: &str) {
        let tolerance = 1e-8 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "{what}: {actual} vs {expected}"
        );
    }

    const H: f64 = 1e-6;

    #[test]
    fn first_derivatives_match_finite_differences() {
        let f = fixture();
        let nb = f.v.len();
        for (name, c, y) in [("bus", &f.identity, &f.ybus), ("from", &f.cf, &f.yf)] {
            let flows = flows(c, y, &f.v);
            let jacobian = CsrMatrix::hstack([flows.ds_dva.clone(), flows.ds_dvm.clone()]).unwrap();
            let (re, im) = (dense(&jacobian.real()), dense(&jacobian.imag()));
            let (_, squared_jacobian) = squared_magnitudes(&flows);
            let squared_jacobian = dense(&squared_jacobian);
            for k in 0..2 * nb {
                let plus = super::flows(c, y, &moved(&f.v, k, H)).s;
                let minus = super::flows(c, y, &moved(&f.v, k, -H)).s;
                for row in 0..plus.len() {
                    let ds = (plus[row] - minus[row]) / (2.0 * H);
                    assert_close(re[row][k], ds.re, &format!("{name} Re dS[{row}][{k}]"));
                    assert_close(im[row][k], ds.im, &format!("{name} Im dS[{row}][{k}]"));
                    let da = (plus[row].norm_sqr() - minus[row].norm_sqr()) / (2.0 * H);
                    assert_close(
                        squared_jacobian[row][k],
                        da,
                        &format!("{name} dA[{row}][{k}]"),
                    );
                }
            }
        }
    }

    #[test]
    fn bus_derivatives_agree_with_the_newton_jacobian() {
        // The power flow's Jacobian holds the P rows at PV and PQ buses and the Q rows at PQ
        //   buses, against the angles at PV and PQ buses and the magnitudes at PQ buses
        let case = Case::parse(CASE).unwrap();
        let model = crate::ac::AcModel::new(&case).unwrap();
        let f = fixture();
        let jacobian = model.jacobian(&f.v).to_csr();
        let jacobian = dense(&jacobian);
        let flows = flows(&f.identity, &f.ybus, &f.v);
        let (dp_dva, dq_dva) = (dense(&flows.ds_dva.real()), dense(&flows.ds_dva.imag()));
        let (dp_dvm, dq_dvm) = (dense(&flows.ds_dvm.real()), dense(&flows.ds_dvm.imag()));
        let pv_pq: Vec<usize> = model
            .pv_buses()
            .iter()
            .chain(model.pq_buses())
            .copied()
            .collect();
        let pq = model.pq_buses();
        let rows: Vec<(usize, bool)> = pv_pq
            .iter()
            .map(|&b| (b, true))
            .chain(pq.iter().map(|&b| (b, false)))
            .collect();
        let cols: Vec<(usize, bool)> = pv_pq
            .iter()
            .map(|&b| (b, true))
            .chain(pq.iter().map(|&b| (b, false)))
            .collect();
        for (i, &(bus, real)) in rows.iter().enumerate() {
            for (k, &(other, angle)) in cols.iter().enumerate() {
                let expected = match (real, angle) {
                    (true, true) => dp_dva[bus][other],
                    (true, false) => dp_dvm[bus][other],
                    (false, true) => dq_dva[bus][other],
                    (false, false) => dq_dvm[bus][other],
                };
                assert!(
                    (jacobian[i][k] - expected).abs() <= 1e-12 * expected.abs().max(1.0),
                    "J[{i}][{k}]: {} vs {expected}",
                    jacobian[i][k]
                );
            }
        }
    }

    #[test]
    fn second_derivatives_match_finite_differences_of_the_first() {
        let f = fixture();
        let nb = f.v.len();
        for (name, c, y) in [("bus", &f.identity, &f.ybus), ("from", &f.cf, &f.yf)] {
            let rows = c.shape().0;
            // Distinct weights, so that no term can hide behind another
            let p: Vec<f64> = (0..rows).map(|r| 0.7 + 0.3 * r as f64).collect();
            let q: Vec<f64> = (0..rows).map(|r| -0.4 + 0.25 * r as f64).collect();
            let mu: Vec<f64> = (0..rows).map(|r| 1.1 + 0.2 * r as f64).collect();
            let hessian = dense(&weighted_hessian(c, y, &f.v, &p, &q));
            let flows_at = flows(c, y, &f.v);
            let squared_hessian = dense(&squared_magnitude_hessian(c, y, &f.v, &flows_at, &mu));
            // The gradient of pᵀ Re S + qᵀ Im S, and of μᵀ |S|², at voltages v
            let gradients = |v: &[C64]| {
                let flows = flows(c, y, v);
                let jacobian =
                    CsrMatrix::hstack([flows.ds_dva.clone(), flows.ds_dvm.clone()]).unwrap();
                let weighted: Vec<f64> = (0..2 * nb)
                    .map(|k| {
                        jacobian.real().vec_mul(&p).unwrap()[k]
                            + jacobian.imag().vec_mul(&q).unwrap()[k]
                    })
                    .collect();
                let squared = squared_magnitudes(&flows).1.vec_mul(&mu).unwrap();
                (weighted, squared)
            };
            for k in 0..2 * nb {
                let (wp, sp) = gradients(&moved(&f.v, k, H));
                let (wm, sm) = gradients(&moved(&f.v, k, -H));
                for i in 0..2 * nb {
                    let expected = (wp[i] - wm[i]) / (2.0 * H);
                    assert_close(hessian[i][k], expected, &format!("{name} H[{i}][{k}]"));
                    let expected = (sp[i] - sm[i]) / (2.0 * H);
                    assert_close(
                        squared_hessian[i][k],
                        expected,
                        &format!("{name} HA[{i}][{k}]"),
                    );
                }
            }
            // Symmetric, as a Hessian must be
            for (i, row) in hessian.iter().enumerate() {
                for (k, &value) in row.iter().enumerate() {
                    assert_close(value, hessian[k][i], &format!("{name} symmetry"));
                }
            }
        }
    }
}
