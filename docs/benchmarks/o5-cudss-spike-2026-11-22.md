# A GPU KKT factorization with cuDSS: the spike, 2026-11-22

This is phase O5's target 4, PR 4a (`OptimizationChecklist.md`, "Phase O5"): NVIDIA's cuDSS
sparse direct solver, bound in `ybus_cuda`, measured on real KKT matrices before the solvers
change. It has two gates, each written before it was measured.

**Outcome: the first gate passes, the second fails.**
- **cuDSS would be fast:** a KKT iteration's work (a refactorization and two refined solves)
  runs ×2.3–6.2 faster than the CPU LDLᵀ, nominally.
- **Its results are wrong:**
  - **The inertia:** wrong on all 9 dumped KKT matrices.
  - **The solves:** relative residuals of 5 × 10⁴ to 10⁷, which refinement never brought
    toward the 1e-12 tolerance.
  - **No setting tried fixed it:** matching, Bunch–Kaufman pivoting, the scaled epsilon, AMD
    ordering, symmetric equilibration, and cuDSS's general mode.

So the solvers don't get cuDSS. Target 4 closes, and phase O5 with it.
- **The bindings stay** as `ybus_cuda::GpuLdlt`, tested on synthetic matrices, with this
  finding in their docs.

## The setup

- **The machine:** an Intel Core i7-14700K (28 threads, 15 GB under WSL2) and an NVIDIA GeForce
  RTX 4070 Ti SUPER (16 GB), Linux 6.6.
- **The software:** CUDA toolkit 13.5 on a driver reporting CUDA 13.3 (KMD 610.88), and cuDSS
  1.0.0 for CUDA 13.

## The first gate: is the factorization still the largest stage?

Same-day samply profiles of the predictor-corrector method's AC solves (`ldlt-uncorrected`,
AMD). The KKT factorization is the dense kernel, the front assembly and the numeric
factorization combined, which `profile_summary.py` now prints. It's set against the largest
other stage, the evaluation:

| Case | 1 thread | 14 threads |
|---|---|---|
| case9241pegase | 36.5% against 23.7% | 32.6% against 28.3% |
| case_ACTIVSg10k | 42.2% against 27.4% | 36.2% against 30.9% |
| case_ACTIVSg25k | 38.1% against 25.6% | 32.0% against 35.7% |

- **It passes at one thread on every case,** and at 14 on all but case_ACTIVSg25k. Either
  thread count was enough (decided 2026-11-22).
- **The transfer room:** each iteration's upload of the KKT values is at most its factor's
  entries times 8 bytes, 6.6–21.0 MB. At O2's pinned rate (1.07 GB in 55 ms,
  `cuda-spike-2026-11-07.md`), that's at most 0.34–1.08 ms, against 17.7–55.1 ms of
  factorization at one thread: at most 2%.

## The bindings

`ybus_cuda::GpuLdlt` (`src/cudss.rs`):
- **What it does:** analyzes a pattern once, refactors for new values, solves, and reports
  the inertia.
- **How it loads cuDSS:** at run time (`libcudss.so.1`), so the crate builds without it, and
  `GpuLdlt::new` returns `Ok(None)` if it's missing.
- **The FFI:** hand-written against the installed 1.0 header. Its constants were checked
  against the header with a C compiler, and 1.0's `cudssMatrixCreateCsr` takes separate offset
  and index types.
- **Deterministic mode is on:** refactored, fresh and repeated solves were exactly identical.

**What its tests found on synthetic matrices** (`tests/gpu_ldlt_test.rs`):
- **cuDSS pivots within supernodes and perturbs small pivots,** rather than delaying them as
  `ybus_sparse`'s LDLᵀ does.
- **On saddle points,** the inertia matched the CPU's, but an unrefined solve's relative
  residual was about 1e-4. Three host corrections brought it to about 1e-13, and to 2e-9 when
  H itself was indefinite, which needed about five.
- **Turning the perturbation off** (epsilon 0) made the inertia wrong, under every pivot type.
- **A singular matrix factors,** perturbed, with no error.

**Not run:** `compute-sanitizer`. Under WSL it needs NVIDIA's debugger interface enabled on
Windows (`EnableDebuggerInterface.bat`, as administrator).

## The second gate: the dumped KKT matrices

- **The matrices:** a throwaway patch, never committed, wrote `KktSolver`'s matrix as Matrix
  Market at each case's first, middle and last factorization:
  - case9241pegase (order 39,855)
  - case_ACTIVSg10k (order 45,179)
  - case_ACTIVSg25k (order 108,012)
- **The benchmark:** `scripts/benchmark-cudss-kkt.sh` and the `cudss_kkt` example time one
  iteration's work, a refactorization and two refined solves:
  - **The CPU:** as the KKT solver runs it: AMD, threshold 1e-4, refactoring, and normwise
    refinement to 1e-12 in at most 3 steps.
  - **cuDSS:** the same refinement, with its residuals computed on the host.
  - **Reported:** each stage's median over 5 runs, and the corrections cuDSS needs to reach
    1e-12 at all, up to 10.
- **The rows:** [`o5-cudss-spike-2026-11-22.csv`](o5-cudss-spike-2026-11-22.csv).

| Case | Iteration, CPU 1 thread | 14 threads | cuDSS | Nominal speedup | Inertia | cuDSS's residual | CPU's residual |
|---|---|---|---|---|---|---|---|
| case9241pegase | 23–28 ms | 16–21 ms | 6.7–7.7 ms | ×3.3–4.0, ×2.3–3.1 | wrong, 3 of 3 | 6e6–1e7 | 3e-8–3e-3 |
| case_ACTIVSg10k | 29–30 ms | 18–21 ms | 7.0–7.6 ms | ×4.0–4.2, ×2.4–3.0 | wrong, 3 of 3 | 5e4–6e5 | 1e-4–3e-3 |
| case_ACTIVSg25k | 97–102 ms | 50–58 ms | 16–18 ms | ×5.8–6.2, ×2.9–3.2 | wrong, 3 of 3 | 3e6–5e6 | 1e-10–3e-6 |

- **The inertia was wrong by 5 to 960 negative eigenvalues.** For example, case9241pegase's
  middle matrix: 18,485 negative on the CPU, 19,445 on cuDSS.
- **Refinement didn't help:** cuDSS's solutions never reached the tolerance in 10 corrections.
- **Not an artifact of refactoring:** a single fresh factorization was as wrong (residuals
  3–5 × 10⁶, the same inertia).
- **The device memory used:** about 1.4 GB.
- **The CPU's residuals are loose too** on case_ACTIVSg10k and on the solves' last iterations.
  That's the known ill-conditioning of late interior-point iterations, and case_ACTIVSg10k's
  `ldlt-uncorrected` AC solve is a known numerical failure. The CPU still solves to within
  3e-3; cuDSS's solutions are off by orders of magnitude more than their right-hand sides.

**The settings tried,** on case9241pegase's three matrices (a temporary, uncommitted override):

| Setting | Effect |
|---|---|
| matching (`CUDSS_MATCHING_ALG_AUTO`, `MAX_DIAG_PRODUCT`) | none: identical results |
| Bunch–Kaufman pivoting | residuals 3e6–8e7 |
| the scaled pivot epsilon | the best: residuals 1.5e4–4e6, inertia still wrong |
| AMD ordering | residuals 5e6–2e7, the factorization ten times slower |
| symmetric equilibration `D A D` on the host (`D = 1/√(row max)`) | residuals 4e3–2e5 |
| equilibration with Bunch–Kaufman or the scaled epsilon | residuals 2e4–2e5 |
| general mode (LU) on the full matrix, scaled or not | NaN solutions |

**The likely cause:**
- **The suspected mechanism:** cuDSS's static pivoting perturbs small pivots instead of
  delaying them. An interior-point KKT matrix has a zero block and barrier terms spanning many
  orders of magnitude, so many pivots are small, and perturbing them changes the matrix too
  much for refinement to recover. The synthetic saddle points have neither the scale nor the
  near-singularity, and factored well.
- **Supporting it:** the synthetic tests showed the perturbation directly, and `ybus_sparse`'s
  LDLᵀ, which delays such pivots to the parent instead, solves the same matrices.
- **Not measured:** the mechanism was inferred, not traced, and cuDSS's documentation wasn't
  checked for it.

**Not tried:** regularizing the KKT matrix before factoring (adding δ to the Hessian block and
−δ to the zero block), as some GPU interior-point codes do so cuDSS sees a quasi-definite
matrix. It changes the solver's arithmetic and its inertia correction. It would be a solver
change judged by O3's bar, not a drop-in factorization, so it's left for later.

## The gate

Written into the checklist before measuring:

| Condition | Required | Measured | |
|---|---|---|---|
| Speed, case_ACTIVSg10k and case_ACTIVSg25k, at 1 or 14 threads | at least ×1.5 | ×2.4–6.2, nominally | met |
| The inertia, on every dumped matrix | matching | wrong on all 9 | not met |

**No go.** The solvers don't get cuDSS, so PRs 4b and 4c aren't built. Target 4 closes, and
with it phase O5. Next is the publishing plan, P5.
