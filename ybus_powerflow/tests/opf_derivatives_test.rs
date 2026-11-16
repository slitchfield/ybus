//! The AC derivatives of `ybus_powerflow::opf::derivatives` against central finite differences, for
//! symmetry, and against the power flow's Newton Jacobian, on a small case with a transformer
//! (tap and phase shift), line charging and bus shunts. Then the Hessians' plans against the
//! functions they replay, bit for bit, on that case and on MATPOWER's case14, case30 and
//! case118.

mod common;

use ybus_powerflow::ac::AcModel;
use ybus_powerflow::case::Case;
use ybus_powerflow::opf::derivatives::{
    flows, squared_magnitude_hessian, squared_magnitudes, weighted_hessian, FlowsPlan,
    SquaredMagnitudeHessianPlan, SquaredMagnitudesPlan, WeightedHessianPlan,
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
            let plus = ybus_powerflow::opf::derivatives::flows(c, y, &moved(&f.v, k, H)).s;
            let minus = ybus_powerflow::opf::derivatives::flows(c, y, &moved(&f.v, k, -H)).s;
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

// ---------------------------------------------------------------------------------------------
// The plans against the functions, bit for bit
// ---------------------------------------------------------------------------------------------

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    /// Uniform in [-1, 1)
    fn unit(&mut self) -> f64 {
        (self.next() % 2_000_001) as f64 / 1_000_000.0 - 1.0
    }
}

/// A case's three pairs `(C, Y)`: the buses `(I, Ybus)`, and every in-service branch's from and
/// to ends `(Cf, Yf)` and `(Ct, Yt)`
struct Network {
    name: String,
    nb: usize,
    pairs: Vec<(&'static str, CsrMatrix<C64>, CsrMatrix<C64>)>,
}

fn network(name: &str, case: &Case) -> Network {
    let model = AcModel::new(case).unwrap();
    let nb = model.bus_rows().len();
    let internal = |label: usize| {
        model
            .bus_rows()
            .iter()
            .position(|&row| case.bus[row].bus_label == label)
            .unwrap()
    };
    let incidence = |end: &dyn Fn(usize) -> usize| {
        CsrMatrix::<C64>::from_triplets(
            model.branch_rows().len(),
            nb,
            model
                .branch_rows()
                .iter()
                .enumerate()
                .map(|(k, &row)| (k, internal(end(row)), C64::new(1.0, 0.0))),
        )
        .unwrap()
    };
    let cf = incidence(&|row| case.branch[row].from_label);
    let ct = incidence(&|row| case.branch[row].to_label);
    Network {
        name: name.to_string(),
        nb,
        pairs: vec![
            ("bus", CsrMatrix::identity(nb), model.ybus().clone()),
            ("from", cf, model.yf().clone()),
            ("to", ct, model.yt().clone()),
        ],
    }
}

/// The inline case, and MATPOWER's case14, case30 and case118 where they've been fetched
fn networks() -> Vec<Network> {
    let mut networks = vec![network("inline", &Case::parse(CASE).unwrap())];
    for name in ["case14", "case30", "case118"] {
        if let Some(path) = common::matpower_case(name) {
            networks.push(network(name, &Case::from_path(&path).unwrap()));
        }
    }
    networks
}

/// Voltages near 1 p.u. at random angles
fn random_voltages(nb: usize, rng: &mut Lcg) -> Vec<C64> {
    (0..nb)
        .map(|_| C64::from_polar(1.0 + 0.1 * rng.unit(), 0.5 * rng.unit()))
        .collect()
}

/// The same pattern, and the same values bit for bit
fn assert_identical(planned: &CsrMatrix<f64>, direct: &CsrMatrix<f64>, what: &str) {
    assert_eq!(planned.pointers(), direct.pointers(), "{what}: pointers");
    assert_eq!(planned.indices(), direct.indices(), "{what}: indices");
    let bits = |m: &CsrMatrix<f64>| m.values().iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(planned), bits(direct), "{what}: values");
}

/// The weight vectors tried: random, zero (as MIPS starts its equality multipliers), and
/// random with some zeros
fn weight_sets(rows: usize, rng: &mut Lcg) -> Vec<(Vec<f64>, Vec<f64>)> {
    let random = |rng: &mut Lcg| (0..rows).map(|_| rng.unit()).collect::<Vec<f64>>();
    let some_zero = |rng: &mut Lcg| {
        (0..rows)
            .map(|_| {
                if rng.next().is_multiple_of(3) {
                    0.0
                } else {
                    rng.unit()
                }
            })
            .collect::<Vec<f64>>()
    };
    // Signed zeros, so that a value the functions copy and a plan adds 0 to would differ in
    //   its sign: -0.0 + 0.0 is +0.0
    let signed_zero = |rng: &mut Lcg| {
        (0..rows)
            .map(|_| match rng.next() % 3 {
                0 => 0.0,
                1 => -0.0,
                _ => rng.unit(),
            })
            .collect::<Vec<f64>>()
    };
    vec![
        (random(rng), random(rng)),
        (vec![0.0; rows], vec![0.0; rows]),
        (some_zero(rng), some_zero(rng)),
        (vec![-0.0; rows], vec![-0.0; rows]),
        (vec![-0.0; rows], vec![0.0; rows]),
        (signed_zero(rng), signed_zero(rng)),
    ]
}

#[test]
fn weighted_hessian_plans_match_the_function_bit_for_bit() {
    let mut rng = Lcg(7);
    let (mut compared, mut declined) = (0, 0);
    for network in networks() {
        for (end, c, y) in &network.pairs {
            let plan = WeightedHessianPlan::new(c, y);
            for _ in 0..4 {
                let v = random_voltages(network.nb, &mut rng);
                for (set, (p, q)) in weight_sets(c.shape().0, &mut rng).into_iter().enumerate() {
                    let direct = weighted_hessian(c, y, &v, &p, &q);
                    let what = format!("{} {end}", network.name);
                    match plan.evaluate(&v, &p, &q) {
                        Some(planned) => {
                            assert_identical(&planned, &direct, &what);
                            compared += 1;
                        }
                        // A plan declines where a diagonal term of the real or the imaginary
                        //   part is zero with nothing else at its position. With every weight
                        //   zero, both are, and the function's pattern lacks the position;
                        //   with some zero, the other part may still store it, and the decline
                        //   only costs a call to the function
                        None => {
                            if set == 1 || set == 3 {
                                assert!(
                                    direct.num_stored() < plan.pattern().num_stored(),
                                    "{what}"
                                );
                            }
                            declined += 1;
                        }
                    }
                }
            }
        }
    }
    assert!(compared > 0);
    // Zero weights on a branch end leave a to-end bus's diagonal term zero and dropped
    assert!(declined > 0, "zero weights on the branch ends are declined");
}

#[test]
fn squared_magnitude_hessian_plans_match_the_function_bit_for_bit() {
    let mut rng = Lcg(11);
    let mut compared = 0;
    for network in networks() {
        for (end, c, y) in network.pairs.iter().skip(1) {
            let plan = SquaredMagnitudeHessianPlan::new(c, y);
            for _ in 0..4 {
                let v = random_voltages(network.nb, &mut rng);
                let flows = flows(c, y, &v);
                // Positive multipliers, as the interior point keeps them
                let mu: Vec<f64> = (0..c.shape().0).map(|_| 1.0 + rng.unit()).collect();
                let direct = squared_magnitude_hessian(c, y, &v, &flows, &mu);
                let planned = plan
                    .evaluate(&v, &flows, &mu)
                    .expect("positive multipliers are planned");
                assert_identical(&planned, &direct, &format!("{} {end}", network.name));
                compared += 1;
            }
        }
    }
    assert!(compared > 0);
}

#[test]
fn a_plan_declines_only_a_dropped_diagonal() {
    // One branch, limited at its from end: bus 1 is only a to end, so the to-end diagonal term
    //   is the only entry at (1, 1), and a zero weight drops it
    let c = CsrMatrix::<C64>::from_triplets(1, 2, vec![(0, 0, C64::new(1.0, 0.0))]).unwrap();
    let y = CsrMatrix::<C64>::from_triplets(
        1,
        2,
        vec![(0, 0, C64::new(1.0, -10.0)), (0, 1, C64::new(-1.0, 10.0))],
    )
    .unwrap();
    let v = [C64::new(1.0, 0.0), C64::from_polar(0.98, -0.1)];
    let plan = WeightedHessianPlan::new(&c, &y);
    assert!(plan.evaluate(&v, &[0.0], &[0.0]).is_none());
    let direct = weighted_hessian(&c, &y, &v, &[0.0], &[0.0]);
    assert!(direct.num_stored() < plan.pattern().num_stored());
    // A nonzero weight stores it, and the plan matches
    let planned = plan.evaluate(&v, &[0.5], &[-0.25]).unwrap();
    assert_identical(
        &planned,
        &weighted_hessian(&c, &y, &v, &[0.5], &[-0.25]),
        "one branch",
    );
}

#[test]
fn flows_and_squared_magnitude_plans_match_the_functions_bit_for_bit() {
    let mut rng = Lcg(13);
    let complex_bits = |m: &CsrMatrix<C64>| -> Vec<(u64, u64)> {
        m.values()
            .iter()
            .map(|v| (v.re.to_bits(), v.im.to_bits()))
            .collect()
    };
    let mut compared = 0;
    for network in networks() {
        for (end, c, y) in &network.pairs {
            let what = format!("{} {end}", network.name);
            let (flows_plan, magnitudes_plan) =
                (FlowsPlan::new(c, y), SquaredMagnitudesPlan::new(c, y));
            for _ in 0..4 {
                let v = random_voltages(network.nb, &mut rng);
                let (planned, direct) = (flows_plan.evaluate(&v), flows(c, y, &v));
                let s_bits = |f: &ybus_powerflow::opf::derivatives::Flows| -> Vec<(u64, u64)> {
                    f.s.iter()
                        .map(|s| (s.re.to_bits(), s.im.to_bits()))
                        .collect()
                };
                assert_eq!(s_bits(&planned), s_bits(&direct), "{what}: S");
                for (p, d) in [
                    (&planned.ds_dva, &direct.ds_dva),
                    (&planned.ds_dvm, &direct.ds_dvm),
                ] {
                    assert_eq!(p.pointers(), d.pointers(), "{what}");
                    assert_eq!(p.indices(), d.indices(), "{what}");
                    assert_eq!(complex_bits(p), complex_bits(d), "{what}");
                }
                let (squared, jacobian) = magnitudes_plan.evaluate(&direct);
                let (direct_squared, direct_jacobian) = squared_magnitudes(&direct);
                let bits = |x: &[f64]| x.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
                assert_eq!(bits(&squared), bits(&direct_squared), "{what}: |S|²");
                assert_identical(&jacobian, &direct_jacobian, &what);
                compared += 1;
            }
        }
    }
    assert!(compared > 0);
}
