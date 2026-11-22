//! cuDSS's `L D Lᵀ` (`GpuLdlt`) against `ybus_sparse`'s pivoted LDLᵀ, on symmetric indefinite
//! and saddle-point (KKT-shaped) matrices: the same inertia, solves accurate once refined,
//! refactoring, repeatability, cuDSS's perturbation of small pivots, and the errors. Each test skips, saying so, on a machine without
//! CUDA or cuDSS, unless `YBUS_REQUIRE_CUDA=1` makes that a failure.

mod common;

use common::{gpu, Lcg};
use ybus_cuda::{CudaError, Gpu, GpuLdlt};
use ybus_sparse::compressed_matrix::CsrMatrix;
use ybus_sparse::Factorization;

/// A solve's relative residual, `‖b − A x‖∞ / (‖A‖∞ ‖x‖∞ + ‖b‖∞)`
fn residual(a: &CsrMatrix<f64>, x: &[f64], b: &[f64]) -> f64 {
    let (ptr, idx, values) = (a.pointers(), a.indices(), a.values());
    let (mut worst, mut norm) = (0.0f64, 0.0f64);
    for i in 0..b.len() {
        let (mut ax, mut row) = (0.0, 0.0);
        for p in ptr[i]..ptr[i + 1] {
            ax += values[p] * x[idx[p]];
            row += values[p].abs();
        }
        worst = worst.max((b[i] - ax).abs());
        norm = norm.max(row);
    }
    let inf = |v: &[f64]| v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    worst / (norm * inf(x) + inf(b))
}

/// `x` refined against `a` by `steps` corrections, each a solve of the residual, as the KKT
///   solver refines on the host
fn refine(ldlt: &mut GpuLdlt<'_>, a: &CsrMatrix<f64>, b: &[f64], steps: usize) -> Vec<f64> {
    let mut x = ldlt.solve(b).unwrap();
    for _ in 0..steps {
        let ax = a.mul_vec(&x).unwrap();
        let r: Vec<f64> = b.iter().zip(&ax).map(|(b, ax)| b - ax).collect();
        let d = ldlt.solve(&r).unwrap();
        for (x, d) in x.iter_mut().zip(&d) {
            *x += d;
        }
    }
    x
}

/// A symmetric indefinite matrix: random off-diagonals, and diagonals of both signs
fn indefinite(n: usize, fill_percent: u64, seed: u64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed);
    let mut triplets = vec![];
    for i in 0..n {
        for j in 0..i {
            if rng.next() % 100 < fill_percent {
                let value = rng.value();
                triplets.extend([(i, j, value), (j, i, value)]);
            }
        }
        let sign = if rng.next().is_multiple_of(3) {
            -1.0
        } else {
            1.0
        };
        triplets.push((i, i, sign * (2.0 + rng.value())));
    }
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// `[[H, Jᵀ], [J, 0]]`: H a grid Laplacian plus a diagonal (`g × g` nodes), and J one row per
///   `stride` nodes coupling neighbours, with its zero block stored, as a KKT matrix's is.
///   `shift` scales H's diagonal, for a second matrix on the same pattern
fn saddle(g: usize, stride: usize, seed: u64, shift: f64) -> CsrMatrix<f64> {
    let mut rng = Lcg(seed);
    let nh = g * g;
    let rows: Vec<usize> = (0..nh).step_by(stride).collect();
    let n = nh + rows.len();
    let mut triplets = vec![];
    for r in 0..g {
        for c in 0..g {
            let i = r * g + c;
            triplets.push((i, i, shift * (4.5 + rng.value())));
            for (dr, dc) in [(0, 1), (1, 0)] {
                let (r2, c2) = (r + dr, c + dc);
                if r2 < g && c2 < g {
                    let j = r2 * g + c2;
                    triplets.extend([(i, j, -1.0), (j, i, -1.0)]);
                }
            }
        }
    }
    for (k, &node) in rows.iter().enumerate() {
        let row = nh + k;
        triplets.push((row, row, 0.0));
        for j in [node, (node + 1) % nh] {
            let v = 1.5 + rng.value();
            triplets.extend([(row, j, v), (j, row, v)]);
        }
    }
    CsrMatrix::from_triplets(n, n, triplets).unwrap()
}

/// A right-hand side
fn rhs(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed);
    (0..n).map(|_| rng.value()).collect()
}

/// cuDSS on `gpu` for `a`, or `None` (after printing why) where the test should skip
fn analyze<'g>(gpu: &'g Gpu, a: &CsrMatrix<f64>, test: &str) -> Option<GpuLdlt<'g>> {
    let ldlt = GpuLdlt::new(gpu, a).expect("cuDSS is present but failed to analyze");
    if ldlt.is_none() {
        if std::env::var("YBUS_REQUIRE_CUDA").is_ok_and(|v| v == "1") {
            panic!("{test}: no cuDSS, and YBUS_REQUIRE_CUDA=1");
        }
        eprintln!("{test}: skipped, no cuDSS");
    }
    ldlt
}

#[test]
fn solves_and_inertia_match_the_cpu_ldlt() {
    let test = "solves_and_inertia_match_the_cpu_ldlt";
    let Some(gpu) = gpu(test) else { return };
    for (label, a) in [
        ("indefinite", indefinite(400, 2, 3)),
        ("saddle point", saddle(20, 3, 5, 1.0)),
        ("larger saddle point", saddle(60, 2, 7, 1.0)),
    ] {
        let Some(mut ldlt) = analyze(&gpu, &a, test) else {
            return;
        };
        let cpu = a.ldlt().unwrap();
        ldlt.factor(&a).unwrap();
        assert_eq!(ldlt.inertia().unwrap(), cpu.inertia(), "{label}: inertia");
        let b = rhs(a.shape().0, 11);
        // cuDSS perturbs small pivots (its default), so an unrefined solve is loose on a saddle
        //   point (about 1e-4 here), and three corrections bring it to the CPU's accuracy
        let unrefined = residual(&a, &ldlt.solve(&b).unwrap(), &b);
        assert!(unrefined < 1e-3, "{label}: unrefined {unrefined:e}");
        let x = refine(&mut ldlt, &a, &b, 3);
        let r = residual(&a, &x, &b);
        assert!(r < 1e-11, "{label}: refined relative residual {r:e}");
        let cpu_x = cpu.solve(&b).unwrap();
        let worst = x
            .iter()
            .zip(&cpu_x)
            .map(|(g, c)| (g - c).abs() / c.abs().max(1.0))
            .fold(0.0f64, f64::max);
        assert!(worst < 1e-8, "{label}: against the CPU {worst:e}");
    }
}

#[test]
fn refactoring_repeats_exactly_and_solves_as_a_fresh_factorization() {
    let test = "refactoring_repeats_exactly_and_solves_as_a_fresh_factorization";
    let Some(gpu) = gpu(test) else { return };
    let (first, second) = (saddle(30, 3, 5, 1.0), saddle(30, 3, 5, 0.25));
    assert_eq!(first.pointers(), second.pointers());
    assert_eq!(first.indices(), second.indices());
    let b = rhs(first.shape().0, 13);
    // Analyzed on `first`, refactored for `second`: the matching and scaling come from the
    //   analysis's values, so this isn't bit for bit a fresh analysis of `second`, but the same
    //   sequence repeats exactly (deterministic mode)
    let run = || {
        let mut ldlt = analyze(&gpu, &first, test).unwrap();
        ldlt.factor(&first).unwrap();
        ldlt.factor(&second).unwrap();
        let x = ldlt.solve(&b).unwrap();
        (ldlt, x)
    };
    let Some(_) = analyze(&gpu, &first, test) else {
        return;
    };
    let (mut refactored, x) = run();
    let (_, again) = run();
    assert_eq!(again, x, "the same sequence, the same solution");
    assert_eq!(refactored.solve(&b).unwrap(), x, "solved again");
    let inertia = refactored.inertia().unwrap();
    assert_eq!(inertia, second.ldlt().unwrap().inertia());
    let mut fresh = analyze(&gpu, &second, test).unwrap();
    fresh.factor(&second).unwrap();
    assert_eq!(fresh.inertia().unwrap(), inertia);
    // H indefinite here (its diagonal scaled down), so cuDSS's perturbed pivots cost more: each
    //   correction gains only about two digits; both factorizations converge
    for (label, ldlt) in [("refactored", &mut refactored), ("fresh", &mut fresh)] {
        let steps: Vec<f64> = (0..8)
            .map(|k| residual(&second, &refine(ldlt, &second, &b, k), &b))
            .collect();
        assert!(
            steps
                .windows(2)
                .all(|pair| pair[1] < pair[0] / 10.0 || pair[1] < 1e-15)
                && steps[7] < 1e-14,
            "{label}: {steps:?}"
        );
    }
    let timings = refactored.timings();
    assert_eq!((timings.factorizations, timings.solves), (2, 2 + 36));
}

#[test]
fn the_errors() {
    let test = "the_errors";
    let Some(gpu) = gpu(test) else { return };
    let a = saddle(10, 3, 5, 1.0);
    let n = a.shape().0;
    let Some(mut ldlt) = analyze(&gpu, &a, test) else {
        return;
    };
    assert!(matches!(
        ldlt.solve(&rhs(n, 1)),
        Err(CudaError::NotFactored)
    ));
    assert!(matches!(ldlt.inertia(), Err(CudaError::NotFactored)));
    assert!(matches!(
        ldlt.perturbed_pivots(),
        Err(CudaError::NotFactored)
    ));
    let other = indefinite(n, 2, 9);
    assert!(matches!(
        ldlt.factor(&other),
        Err(CudaError::PatternMismatch)
    ));
    ldlt.factor(&a).unwrap();
    assert!(matches!(
        ldlt.solve(&rhs(n + 1, 1)),
        Err(CudaError::ShapeMismatch { .. })
    ));
    let wide = CsrMatrix::from_triplets(2, 3, vec![(0, 0, 1.0), (1, 2, 1.0)]).unwrap();
    assert!(matches!(
        GpuLdlt::new(&gpu, &wide),
        Err(CudaError::ShapeMismatch { .. })
    ));
}

/// cuDSS perturbs a zero pivot rather than failing, so a singular matrix factors, with the
///   inertia of the perturbed matrix (no zero eigenvalue). Pinned here so a change in cuDSS's
///   behavior is noticed; `GpuLdlt`'s docs say callers must not rely on it to detect singularity
#[test]
fn a_singular_matrix_is_perturbed_not_reported() {
    let test = "a_singular_matrix_is_perturbed_not_reported";
    let Some(gpu) = gpu(test) else { return };
    // Row and column 5 stored but zero, so the matrix is exactly singular
    let base = indefinite(40, 5, 3);
    let mut triplets = vec![];
    for i in 0..40 {
        for p in base.pointers()[i]..base.pointers()[i + 1] {
            let j = base.indices()[p];
            let zero = i == 5 || j == 5;
            triplets.push((i, j, if zero { 0.0 } else { base.values()[p] }));
        }
    }
    let a = CsrMatrix::from_triplets(40, 40, triplets).unwrap();
    assert!(a.ldlt().is_err(), "singular on the CPU");
    let Some(mut ldlt) = analyze(&gpu, &a, test) else {
        return;
    };
    ldlt.factor(&a).unwrap();
    let inertia = ldlt.inertia().unwrap();
    assert_eq!(inertia.positive + inertia.negative, 40, "{inertia:?}");
    assert!(
        ldlt.perturbed_pivots().unwrap() >= 1,
        "the zero pivot was perturbed"
    );
}
