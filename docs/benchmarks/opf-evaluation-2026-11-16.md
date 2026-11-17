# The AC OPF's evaluation through kept plans, 2026-11-16

This is O4's second target, from
[`ipm-profile-2026-11-08.md`](ipm-profile-2026-11-08.md): the AC OPF's evaluation. After the
LDLᵀ work (#61, #62), the evaluation was the largest stage on AC, at 32–45% of the
predictor-corrector method's time. Most of that was the Lagrangian Hessian, rebuilt every
iteration by sparse algebra whose patterns never change.

The change keeps a plan of each product, sum, transpose and stack for the solve, and replays the
same floating-point operations in the same order. **The results are the same, bit for bit, on
both solvers.**

## The change

- **`ybus_sparse`:**
  - **`ProductPlan`** (`plan_product`) records each term of a product with its result slot, in
    `matmul`'s order.
  - **`SumPlan`** (`plan_sum`) records each result slot's sources, as `try_add` merges them.
  - **Replaying either** on new values gives `matmul`'s or `try_add`'s pattern and values, with
    no hashing, sorting or merging.
- **`ybus_powerflow`'s `opf::derivatives`:** `WeightedHessianPlan`,
  `SquaredMagnitudeHessianPlan`, `FlowsPlan` and `SquaredMagnitudesPlan` replay the functions of
  the same names.
  - **Transposes and stacks** become gathers, traced by running the same operations once on a
    matrix whose values are its slot numbers plus 1.
  - **Each per-entry operation** is the function's own, in the same order: for example,
    `(dx[i] · m) · conj(dy[j])` for a row scaling followed by a column scaling.
- **The AC OPF's problem** keeps the plans, the Lagrangian's two sums, and each `widen` as a
  kept `AssemblyMap`.

### The one pattern that depends on values

`bilinear_hessian` adds its diagonal terms through `from_diagonal`, which drops zeros. Where
nothing else stores `(i, i)`, the function's pattern then lacks the position.

- **The first version** declined to plan such an evaluation, and fell back to the functions. I
  expected that to be rare.
- **The surveys showed otherwise.**
  - **What they found:**
    - Every result column matched `main`'s. But the cases above about 2,000 buses ran
      0.95–1.00×, against 1.13–1.28× on the small cases.
    - On case2383wp, every one of the 32 Hessian evaluations declined, and each computed the
      Hessian twice.
  - **The cause:** near-zero-impedance ties. Branch row 82 of case2383wp has `r = 0` and
    `x = 1e-4`.
    - The solver keeps the tie's two buses at exactly equal voltages, so the tie's flow is
      exactly `0` at every iteration.
    - The flow Hessian's weights, `2μ Re S` and `2μ Im S`, are then `0`, and so is its diagonal
      term at the far bus.
- **The fix, decided with the maintainer:** hold a dropped term as `-0.0`.
  - **Why it's exact:** `x + (-0.0)` is `x` bit for bit for every `x`, `-0.0` included. So every
    later sum gives `try_add`'s copy of the other operand.
  - **Where a position is really absent:** the plans report it where both the real and the
    imaginary part dropped the term and nothing else stores it.
  - **Why the Lagrangian's pattern never changes:** the bus Hessian stores every diagonal
    position, since Ybus stores its diagonal. The plans check this when they're made.
  - **The fallback remains** only if that check fails, or for a NaN term.

## The same results, checked

- **The 72-case OPF snapshot.** Both solvers on DC and AC hash their solutions, iterations,
  points and multipliers. Every case's hashes match `main`'s, after each of the four commits.
- **The AC surveys.** All 73 cases with costs ran on each solver: MIPS on LU, and the
  predictor-corrector method on LDLᵀ, `AcOpfOptions::mehrotra()`'s default. Of the 1,825 result
  columns on each solver, none differs from `main`'s survey. The columns compared are:
  - the outcome and iterations
  - the agreement with MATPOWER's MP-Core, legacy and core results
  - the KKT residuals
  - the analyses, refinement steps and inertia corrections
- **The tests:**
  - the plans against the functions, bit for bit (`to_bits`), on case14, case30, case118 and a
    four-bus case
  - random voltages, with random, zero and signed-zero weights
  - the sparse plans on random real and complex matrices, in both orientations
- **Planted bugs,** each caught by a test:
  - a product's terms in reverse order
  - a copied entry turned into an add of zero
  - an output not zeroed first
  - a reassociated scaling
  - a `+0.0` sentinel, and a covered zero term added as is
  - the absent positions taken as the union of the two parts' drops, not the intersection
  - two of `dg`'s blocks swapped

  Three changes that were expected to be harmless were tried as well. Each changed nothing, as
  expected:
  - **Swapping `a + b`:** IEEE addition commutes bit for bit.
  - **Reassociating `C`'s scaling:** `C`'s entries are all `1`.
  - **The diagonal copy step, before the sentinel:** a search of 200,000 small random Hessians
    found no input where it shows, since the products' sums start at `+0.0`. Only a unit test of
    that step reached it.

## Speed

**Mean time per run** of `profile_pf`'s whole pipeline (parse, model, solve):
- **Both binaries:** `main` at `d8e9ed1`, and `03e3219`, built at the same time.
- **Runs:** interleaved before, after, before, after, two rounds, on the maintainer's laptop
  with nothing else running.
- **Repeats:** 20 s per case, and 60 s on case_ACTIVSg25k.

| Workload | Case | Before | After | Speedup |
|---|---|--:|--:|--:|
| MIPS (LU) | case300 | 52.6 ms | 45.8 ms | 1.15× |
| MIPS (LU) | case2383wp | 1,583.3 ms | 1,300.3 ms | 1.22× |
| MIPS (LU) | case9241pegase | 9,896.8 ms | 8,550.3 ms | 1.16× |
| MIPS (LU) | case_ACTIVSg25k | 71,551.4 ms | 66,355.6 ms | 1.08× |
| Predictor-corrector (LDLᵀ) | case300 | 33.8 ms | 26.4 ms | 1.28× |
| Predictor-corrector (LDLᵀ) | case2383wp | 748.6 ms | 504.9 ms | 1.48× |
| Predictor-corrector (LDLᵀ) | case9241pegase | 3,783.2 ms | 2,642.2 ms | 1.43× |
| Predictor-corrector (LDLᵀ) | case_ACTIVSg25k | 23,918.6 ms | 17,377.1 ms | 1.38× |

The profile's bound was 1.4–1.6× on AC. The predictor-corrector method reaches 1.38–1.48× on the
large cases. MIPS gains less, since its LU factorization is 49–89% of its time.

**The surveys' times,** from each case's `opf_ms`. The before and after surveys ran some hours
apart on the same machine, so differences of a few percent are noise:

| Solver | Median | Under 1,000 buses | 1,000–10,000 | Over 10,000 | Total |
|---|--:|--:|--:|--:|--:|
| MIPS (LU) | 1.22× | 1.29× | 1.11× | 1.03× | 0.99× |
| Predictor-corrector (LDLᵀ) | 1.26× | 1.23× | 1.32× | 1.28× | 1.28× |

- **The largest cases on the predictor-corrector method:**
  - case_SyntheticUSA: 160.1 s to 124.9 s
  - case_ACTIVSg70k: 81.4 s to 64.7 s
  - case_ACTIVSg25k: 23.3 s to 17.6 s
- **MIPS's total** is set by case_SyntheticUSA and case_ACTIVSg70k. Those run 0.96–0.98× in the
  surveys, within their noise, and LU dominates them. `profile_pf`'s interleaved timing on
  case_ACTIVSg25k gives MIPS 1.07–1.09×.

**Where it came from:** samply at 1 kHz, each pipeline repeated for 10 s, folded by
`profile_summary.py`. The OPF evaluation stage's share, and its time per run (share × time):

| Workload | Case | Share | Evaluation per run | Faster |
|---|---|--:|--:|--:|
| MIPS (LU) | case2383wp | 23.5% → 10.3% | 371 → 135 ms | 2.8× |
| MIPS (LU) | case9241pegase | 17.0% → 6.8% | 1,668 → 596 ms | 2.8× |
| MIPS (LU) | case_ACTIVSg25k | 10.3% → 4.2% | 7,278 → 2,767 ms | 2.6× |
| Predictor-corrector (LDLᵀ) | case2383wp | 45.2% → 24.1% | 337 → 121 ms | 2.8× |
| Predictor-corrector (LDLᵀ) | case9241pegase | 39.3% → 19.3% | 1,443 → 525 ms | 2.7× |
| Predictor-corrector (LDLᵀ) | case_ACTIVSg25k | 36.6% → 18.8% | 9,254 → 3,398 ms | 2.7× |

**After the change,** on the predictor-corrector method on case9241pegase:
- **The largest stages:** the LDLᵀ's dense kernel at 20.7%, the OPF's evaluation at 19.3%, and
  the KKT assembly at 18.1%.
- **The top frames of their own:**
  - `bunch_kaufman::factor_front`, 20.4%
  - `CompressedMatrix::from_triplets`, 10.6%, mostly the KKT assembly's
  - `Ldlt::solve`, 9.5%
  - the plans' `SumPlan::add_values`, 5.1% across its two instances

## What's left

- **The KKT assembly** is now as large as the evaluation: 18.1% on case9241pegase.
  `from_triplets` rebuilds the KKT matrix from triplets at every iteration, so a kept
  `AssemblyMap` would be the analogue there. That's O4's next target, allocation and KKT
  assembly.
- **The evaluation's own allocation:**
  - The plans return fresh matrices, cloning the patterns, and build fresh value vectors at each
    step. Writing into kept buffers would need the `Problem` trait to take buffers, or interior
    mutability.
  - The remaining 19% is now mostly the per-entry arithmetic and these copies.
- **Fusing steps:** a block's scalings, its sum and its diagonal are separate passes. Fusing them
  is possible, as long as the operations keep their order.
