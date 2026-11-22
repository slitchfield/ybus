# A GPU KKT factorization with cuDSS: the spike, 2026-11-22

This is phase O5's target 4, PR 4a (`OptimizationChecklist.md`, "Phase O5"): NVIDIA's cuDSS
sparse direct solver, bound in `ybus_cuda`, measured on real KKT matrices before the solvers
change. It has two gates, each written before it was measured.

**Outcome: the first gate passes, the second fails.**
- **cuDSS would be fast:** a KKT iteration's work (a refactorization and two refined solves)
  runs ×2.1–5.7 faster than the CPU LDLᵀ, nominally.
- **Its results are wrong:**
  - **The inertia:** wrong on all 9 dumped KKT matrices.
  - **The solves:** relative residuals of 7 × 10² to 10¹¹, which refinement never brought
    toward the 1e-12 tolerance.
  - **Why:** the KKT matrix's zero constraint block gives cuDSS exact zero pivots it can't
    pivot around, and it replaces a fifth to a third of all pivots with ±1e-13 (see "Why cuDSS
    disagrees" below).
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
| case9241pegase | 24–28 ms | 16–24 ms | 6.9–8.4 ms | ×3.1–3.8, ×2.1–3.1 | wrong, 3 of 3 | 2e4–3e6 | 3e-8–3e-3 |
| case_ACTIVSg10k | 30 ms | 19 ms | 7.0–8.1 ms | ×3.8–4.3, ×2.4–2.7 | wrong, 3 of 3 | 2e5–5e7 | 1e-4–3e-3 |
| case_ACTIVSg25k | 101–111 ms | 53–54 ms | 18–20 ms | ×5.5–5.7, ×2.9–3.0 | wrong, 3 of 3 | 7e2–1e11 | 1e-10–3e-6 |

- **The inertia was wrong by 6 to 491 negative eigenvalues.**
- **Refinement didn't help:** cuDSS's solutions never reached the tolerance in 10 corrections.
- **Not an artifact of refactoring:** a single fresh factorization was as wrong.
- **The device memory used:** about 1.4 GB.
- **The CPU's residuals are loose too** on case_ACTIVSg10k and on the solves' last iterations.
  That's the known ill-conditioning of late interior-point iterations, and case_ACTIVSg10k's
  `ldlt-uncorrected` AC solve is a known numerical failure. The CPU still solves to within
  3e-3.

**These numbers are after a fix to the bindings,** found while explaining the failure (below).
- **The bug:** the first version ran cuDSS's analysis before uploading any values. cuDSS 1.0's
  default matching (`AUTO`) computes its permutation and scaling from the values at analysis,
  so it worked on zeros.
- **The fix:** `GpuLdlt::new` now uploads them first.
- **The effect:** it changed the residuals, better on some matrices and worse on others, but
  not the verdict. Before the fix, the residuals were 5e4–1e7, the inertia was off by 5–960,
  and the matching settings tried below did nothing.

**The settings tried,** on case9241pegase's three matrices (a temporary, uncommitted override),
before the fix:

| Setting | Effect |
|---|---|
| matching (`CUDSS_MATCHING_ALG_AUTO`, `MAX_DIAG_PRODUCT`) | none, since the analysis saw zeros |
| Bunch–Kaufman pivoting | residuals 3e6–8e7 |
| the scaled pivot epsilon | residuals 1.5e4–4e6, inertia still wrong |
| AMD ordering | residuals 5e6–2e7, the factorization ten times slower |
| symmetric equilibration `D A D` on the host (`D = 1/√(row max)`) | residuals 4e3–2e5 |
| equilibration with Bunch–Kaufman or the scaled epsilon | residuals 2e4–2e5 |
| general mode (LU) on the full matrix, scaled or not | NaN solutions |

## Why cuDSS disagrees with the CPU

**The short answer:** the KKT matrix's zero constraint block needs pivots that cuDSS can't
choose.
- **Exactly zero pivots:** cuDSS's ordering eliminates most constraint rows before the
  variables they couple to, so their pivots are exactly zero.
- **What could take them:** a 2×2 pivot pairing such a row with a coupled variable. When that
  variable is in another front, the zero pivot must be delayed to the parent, as our LDLᵀ does.
- **What cuDSS does instead:** it pivots only within a supernode and doesn't delay. It replaces
  each such pivot with ±1e-13.
- **The damage:** dividing by those makes the factors grow until rounding swamps the matrix.

The evidence, step by step, all measured on the dumped matrices unless said otherwise:

1. **The CPU's inertia is the matrix's own.** Sylvester's law makes the inertia independent of
   the ordering, for a stable factorization.
   - **It's invariant:** on case9241pegase's and case_ACTIVSg25k's matrices, the CPU gave the
     same inertia under AMD with pivot thresholds 1e-4, 0.1 and 0.5, and under RCM.
   - **It's the theory's:** it equals (n_x, n_eq), for example (+21,372, −18,483) on
     case9241pegase's first matrix. That's the inertia of a KKT matrix whose Hessian is
     positive definite on the constraints' null space.

   So cuDSS's is the wrong one.
2. **The structure.** Each matrix is `[W Jᵀ; J 0]`, with n_eq constraint rows whose diagonal
   isn't even stored (0 of 18,483 on case9241pegase).
   - **W's diagonal has zeros too:** 4,692 zero entries on case9241pegase's first matrix.
   - **The entries span 24–26 orders of magnitude,** for example 2.5e-18 to 2.0e8.
3. **cuDSS perturbs a fifth to a third of all pivots.** `CUDSS_DATA_NPIVOTS`, which cuDSS's
   documentation defines as the number of perturbed 1×1 pivots, is:
   - 7,755–13,090 of case9241pegase's 39,855
   - 11,510–14,402 of case_ACTIVSg10k's 45,179
   - 25,565–36,786 of case_ACTIVSg25k's 108,012
4. **They're the zero-block pivots its ordering creates.** In cuDSS's own elimination order
   (`CUDSS_DATA_PERM_REORDER_ROW`, with matching off), a constraint row eliminated before all
   of the W variables it couples to still has its zero diagonal when its turn comes. Such rows:
   - 9,752 and 11,952 on two of case9241pegase's matrices, against 11,872 and 12,513 perturbed
     pivots
   - 31,521 on case_ACTIVSg25k's middle matrix, against 32,230
   
   Our AMD order creates such zero pivots too, 3,956–21,606 of them. Our LDLᵀ resolves them
   with 2×2 pivots (7,423–26,930) and delays (2,971–10,650).
5. **cuDSS can't take them otherwise.** Its documentation describes pivoting only "within the
   diagonal sub-block" of a supernode, Bunch–Kaufman included, and describes no delayed pivots.
   An unacceptable pivot is replaced by an "appropriately signed" epsilon, 1e-13 in double
   precision, an absolute value.
6. **The same rule fails the same way on the CPU.** `ybus_sparse`'s static LDLᵀ
   (`factor_static_ldlt_replacing`), with pivots below 1e-13 replaced by ±1e-13 on the AMD
   ordering:

   | Matrix | Pivots replaced | Max \|D\| | Residual, unrefined and refined | Inertia |
   |---|---|---|---|---|
   | case9241pegase, first | 7,969 | 1e21 | 6e6 | wrong |
   | case9241pegase, middle | 11,097 | 2e21 | 7e6 | wrong |
   | case_ACTIVSg25k, middle | 25,128 | 4e22 | 2e5 | wrong |

   Refinement over 10 steps didn't reduce the residual. Phase O4 saw the same with this
   factorization in the solvers (about 9,700 pivots replaced on case9241pegase's AC, D up to
   7e12), and set static pivots aside for the delayed-pivot LDLᵀ.
7. **The perturbation is small, but dividing by it isn't.** Changing a zero to 1e-13 is a tiny
   backward error. But the multipliers become 1e13 times the coupling entries, and D grows to
   1e21 and beyond.
   - **The rounding that follows:** in later updates it's about the unit roundoff times that
     growth, 1e-16 × 1e21 = 1e5 relative to the matrix. That's the size of the residuals seen
     (1e2–1e7).
   - **Why refinement can't recover:** the factorization isn't backward stable, so each
     correction is as wrong as the first solve.
   - **The inertia:** it counts the signs of these corrupted and arbitrarily signed pivots, so
     it's wrong too.
8. **With the zero block regularized, cuDSS is exact.** On case9241pegase's middle matrix, with
   δ added to W's diagonal and −δ to the constraint block (`[W + δI, Jᵀ; J, −δI]`, every
   diagonal stored):

   | δ | Perturbed pivots | Inertia | Residual on the regularized matrix |
   |---|---|---|---|
   | 0 (explicit zeros) | 13,090 | wrong | 1e4 |
   | 1e-8 | 302 | wrong | 3e3 |
   | 1e-4 | 0 | (21,372, 18,483), the CPU's on it | 1e-1 |
   | 1e-2 | 0 | the CPU's | 5e-5 |
   | 1 | 0 | the CPU's | 9e-9 |

   - **The zero block is the problem:** with no zero pivots left, cuDSS's factorization is
     clean and its inertia is right.
   - **But regularization has a cost:** its accuracy improves only as δ grows, the known
     stability bound for quasi-definite factorizations, about 1/δ. At δ = 1e-4, the regularized
     matrix's solution is far from the original's (residual 1e0 against the original).
   - **Where it can't help:** a regularization large enough for cuDSS changes the problem the
     interior-point step solves. Where W itself is indefinite (the 7 zero and 299 negative
     diagonal entries of that matrix), W + δI isn't positive definite at small δ either.

**What this doesn't settle:**
- **The growth was measured on our CPU static LDLᵀ under cuDSS's rule,** not read from cuDSS's
  own factors. cuDSS's residuals are of the same size, which is consistent with it.
- **Whether a better ordering would fix cuDSS:** one that eliminates each constraint row after
  a coupled variable, as `ybus_ipm`'s matched ordering does, supplied through
  `CUDSS_DATA_USER_PERM`. It wasn't tried.
  - **For it:** it would remove the zero pivots at elimination.
  - **Against it:** cuDSS's matching, its best remaining tool, is disabled with a user
    permutation.
  - **In supernodes:** whether a variable and its constraint would land in the same supernode
    for a 2×2 pivot is up to cuDSS's symbolic analysis.

## The gate

Written into the checklist before measuring:

| Condition | Required | Measured | |
|---|---|---|---|
| Speed, case_ACTIVSg10k and case_ACTIVSg25k, at 1 or 14 threads | at least ×1.5 | ×2.4–5.7, nominally | met |
| The inertia, on every dumped matrix | matching | wrong on all 9 | not met |

**No go.** The solvers don't get cuDSS, so PRs 4b and 4c aren't built. Target 4 closes, and
with it phase O5. Next is the publishing plan, P5.
