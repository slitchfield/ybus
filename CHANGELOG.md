# Changelog

The notable changes to the four published crates of the `ybus` workspace: `ybus_sparse`,
`ybus_ipm`, `ybus_cuda` and `ybus_powerflow`. They share one version, so each release has one
entry here, with a subsection for each crate it changes. `ybus_ffi`, the C interface, isn't
published and isn't covered.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the crates
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html). In the 0.x releases a minor
version bump (0.1 to 0.2) may break the API, and a patch bump doesn't.

## [Unreleased]

## [0.1.0] - YYYY-MM-DD

The first release on crates.io. The API is young, so expect breaking changes in the 0.x
releases.

### Added

#### `ybus_sparse`

- `SparseMatrix` for building and editing a matrix, and `CompressedMatrix` (`CsrMatrix`,
  `CscMatrix`) for the numerical work, generic over the element type, `f64` by default.
- Arithmetic, sparse products, norms and selections, and refilling values into a fixed pattern
  through an assembly map.
- Cholesky and LDLᵀ analyses with a permutation of the caller's
  (`analyze_cholesky_with_permutation`, `analyze_ldlt_with_permutation`), for orderings that
  need more than the pattern.
- Product and sum plans (`ProductPlan`, `SumPlan`): a product or sum planned once for fixed
  patterns, then recomputed from new values bit for bit as `matmul` and `try_add` compute it,
  with no hashing, sorting or allocation.
- Factorizations: sparse LU with threshold partial pivoting, Cholesky, and LDLᵀ with
  Bunch–Kaufman pivoting, with orderings (AMD, COLAMD, reverse Cuthill–McKee), condition estimates,
  and solves with many right-hand sides.
- Conjugate gradient, and Matrix Market reading and writing.
- Optional features: `serde`, and `rayon` for parallel products and solves.

#### `ybus_ipm`

- MIPS, the MATPOWER Interior Point Solver, reimplemented from MATPOWER's `mips.m` and checked
  against MATPOWER's own iterations and results (`mips::mips`, `mips::qps_mips`).
- A primal-dual method with Mehrotra's predictor-corrector (`mehrotra::mehrotra`,
  `mehrotra::qps_mehrotra`).
- The `Problem` trait for general nonlinear programs, and `QuadraticProblem` for quadratic and
  linear ones. Either solver factors its KKT systems by LU or by LDLᵀ. The predictor-corrector
  method corrects an indefinite Hessian using LDLᵀ's inertia.
- The matched ordering for the pivoted LDLᵀ (`kkt::KktOrdering::Matched`, Duff and Pralet,
  2005): each constraint paired with a variable it touches, so that its pivot isn't delayed.
  AMD stays the default.
- The KKT factorizations' statistics on `Solution` (`kkt: KktStatistics`): factorizations,
  delayed pivots, 2×2 pivots and factor entries, totalled over a solve.
- Buffer-taking evaluations on `Problem` (`objective_hessian_into`, `constraints_into`,
  `lagrangian_hessian_into`), which the solvers keep across iterations, so a problem whose
  patterns don't change can refill values in place. Each defaults to its allocating method.
- The solvers keep the Newton system's product and sum, and the KKT matrix's assembly, planned
  for the solve while their patterns hold, with the same results bit for bit.

#### `ybus_cuda` (experimental)

- `Gpu` and `GpuFactor`: sparse triangular solves with many right-hand sides on CUDA through
  cuSPARSE, against a Cholesky or LU factorization made on the CPU.
- Blocks of right-hand sides built on the device, with rank-one updates, for the DC studies'
  PTDF columns and outages.
- CUDA is loaded at run time, so the crate builds without the CUDA toolkit.

#### `ybus_powerflow`

- `Case`, read from MATPOWER case files.
- The DC power flow, and DC studies that reuse one factorization: batched injections,
  single-branch outages, and PTDF and LODF sensitivities.
- The AC power flow, by Newton's method or the fast-decoupled XB and BX methods, optionally
  enforcing generators' reactive limits.
- The DC and AC optimal power flows, solved with `ybus_ipm`'s MIPS or predictor-corrector
  method. The AC OPF keeps plans of its derivatives' sparse products and sums for the solve
  (`opf::derivatives::WeightedHessianPlan`, `SquaredMagnitudeHessianPlan`, `FlowsPlan`,
  `SquaredMagnitudesPlan`), with the same results bit for bit.
- Islands and bridges of a network.
- Each solver is tested against results computed by MATPOWER.
- Optional features: `rayon` for parallel outages and sensitivities, and `cuda`
  (experimental) for the DC studies on a GPU, through `ybus_cuda`.

[Unreleased]: https://github.com/slitchfield/ybus/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/slitchfield/ybus/releases/tag/v0.1.0
