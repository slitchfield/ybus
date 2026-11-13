//! The AC derivatives of `powerflow::opf::derivatives` against central finite differences, for
//! symmetry, and against the power flow's Newton Jacobian, on a small case with a transformer
//! (tap and phase shift), line charging and bus shunts.

use powerflow::ac::AcModel;
use powerflow::case::Case;
use powerflow::opf::derivatives::{
    flows, squared_magnitude_hessian, squared_magnitudes, weighted_hessian,
};
use ybus_sparse::compressed_matrix::CsrMatrix;
use ybus_sparse::Complex;

type C64 = Complex<f64>;

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
    let model = AcModel::new(&case).unwrap();
    let nb = model.bus_rows().len();
    let v = model
        .bus_rows()
        .iter()
        .map(|&row| C64::from_polar(case.bus[row].vm_pu, case.bus[row].va_deg.to_radians()))
        .collect();
    // Each in-service branch's from bus, as an internal index
    let internal = |label: usize| {
        model
            .bus_rows()
            .iter()
            .position(|&row| case.bus[row].bus_label == label)
            .unwrap()
    };
    let cf = CsrMatrix::<C64>::from_triplets(
        model.branch_rows().len(),
        nb,
        model
            .branch_rows()
            .iter()
            .enumerate()
            .map(|(k, &row)| (k, internal(case.branch[row].from_label), C64::new(1.0, 0.0))),
    )
    .unwrap();
    Fixture {
        v,
        ybus: model.ybus().clone(),
        yf: model.yf().clone(),
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
            let plus = powerflow::opf::derivatives::flows(c, y, &moved(&f.v, k, H)).s;
            let minus = powerflow::opf::derivatives::flows(c, y, &moved(&f.v, k, -H)).s;
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
    let model = AcModel::new(&case).unwrap();
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
            let jacobian = CsrMatrix::hstack([flows.ds_dva.clone(), flows.ds_dvm.clone()]).unwrap();
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
