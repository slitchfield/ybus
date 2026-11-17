# The Newton system and the KKT assembly kept planned, 2026-11-17

This is O4's third target: allocation and the KKT assembly. After #79, the KKT assembly was 18.1%
of case9241pegase's AC predictor-corrector solve. Both solvers rebuilt the Newton system from
scratch at every iteration:
- `Jhᵀ`, transposed, then `Jhᵀ (D Jh)` by `matmul`
- `W = H + …` by `try_add`
- `W` collected as triplets
- the KKT matrix by `from_triplets`: once for LU, and twice for LDLᵀ

The patterns of all of these hold for a whole solve. **The results are the same, bit for bit, on
both solvers and both models.**

## The change

Three commits.

1. **The Newton system planned** (`ybus_ipm`'s `newton` and `kkt` modules).
   - **`WeightedGram`** keeps `W = H + Jhᵀ (D Jh)` planned: `Jhᵀ` as a traced gather, then
     #79's `ProductPlan` and `SumPlan`.
   - **`KktSolver`** keeps the KKT matrix's assembly as `AssemblyMap`s: one in CSC for LU, and
     for LDLᵀ the lower triangle's, the mirror's gather and the full matrix's.
   - **Both are keyed by their inputs' patterns,** and planned afresh when one differs, as the
     kept analysis is. That covers a problem whose patterns change, such as #79's fallback.
2. **Buffer-taking `Problem` methods:** `objective_hessian_into`, `constraints_into` and
   `lagrangian_hessian_into`.
   - **Each defaults to the allocating method,** so implementors are unchanged.
   - **The solvers keep one buffer of each** across iterations.
   - **`AcOpfProblem` refills its kept patterns in place.**
3. **The LU factor's arrays reserved:** `SymbolicLu` keeps the last factorization's sizes of `L`
   and `U`, and the next factorization reserves them.
   - **Why:** an allocation profile of case9241pegase after commits 1 and 2 put almost all of
     MIPS's allocation in `gilbert_peierls`, growing those arrays one push at a time: 5.9% of
     the solve.
   - **What the same profile ruled out:**
     - #79's intermediate vectors came to about 0.1%, so they stay as they are.
     - The LDLᵀ's allocation is spread over sites of 0.4% or less.

## The same results, checked

- **The 72-case OPF snapshot** matches `main`'s (`e94f898`) after each commit. It hashes both
  solvers on DC and AC: their solutions, iterations, points and multipliers.
- **The surveys:** all 73 cases with costs, on four solver and model pairs:
  - MIPS on LU, AC and DC
  - the predictor-corrector method on LDLᵀ, AC and DC, the defaults of
    `AcOpfOptions::mehrotra()` and `DcOpfOptions::mehrotra()`

  Of the 1,825 result columns in each of the four surveys, none differs from `main`'s. The
  columns are the outcome, iterations, agreement with MATPOWER, KKT residuals, analyses,
  refinement steps and inertia corrections.
- **The tests:**
  - the kept assembly against the triplet assembly, for LU and LDLᵀ, with signed zeros,
    magnitudes far apart, reused patterns and changed ones
  - `WeightedGram` against `try_add` and `matmul`, through phases of reused patterns, then only
    `Jh`'s changing, then only `H`'s
  - a kept LU analysis whose fill changes between factorizations, against a fresh one
  - `QuadraticProblem::objective_hessian_into` from an empty buffer, a foreign one, and a stale
    one
- **Planted bugs,** each caught:
  - an assembly that's never planned afresh
  - plans that ignore `Jh`'s pattern, or `H`'s

  The first version of the `W` test drew new patterns every round. It never reused a plan, and
  missed the first of those; it now runs in phases.

## Speed

**Mean time per run** of `profile_pf`'s whole pipeline (parse, model, solve):
- **Both binaries:** `main` at `e94f898`, and `19d3b86`, built together.
- **Runs:** interleaved before, after, two rounds, on the maintainer's laptop with nothing else
  running.
- **Repeats:** 15 s per case, and 60 s for AC on case_ACTIVSg25k.

| Workload | case300 | case2383wp | case9241pegase | case_ACTIVSg25k |
|---|--:|--:|--:|--:|
| AC, MIPS (LU) | 1.03× | 1.10× | 1.07× | 1.06× |
| AC, predictor-corrector (LDLᵀ) | 1.34× | 1.21× | 1.26× | 1.20× |
| DC, MIPS (LU) | 1.16× | 1.15× | 1.13× | 1.12× |
| DC, predictor-corrector (LDLᵀ) | 1.09× | 1.26× | 1.09× | 1.12× |

The times themselves:

| Workload | Case | Before | After |
|---|---|--:|--:|
| AC, predictor-corrector | case9241pegase | 2,485.0 ms | 1,969.1 ms |
| AC, predictor-corrector | case_ACTIVSg25k | 16,426.8 ms | 13,741.6 ms |
| AC, MIPS | case_ACTIVSg25k | 62,301.0 ms | 58,774.4 ms |
| DC, MIPS | case_ACTIVSg25k | 20,455.6 ms | 18,187.5 ms |
| DC, predictor-corrector | case_ACTIVSg25k | 1,417.0 ms | 1,262.3 ms |

**The surveys' times,** from each case's `opf_ms`. The before and after surveys ran back to back,
on the same morning:

| Pair | Median | Under 1,000 buses | 1,000–10,000 | Over 10,000 | Total |
|---|--:|--:|--:|--:|--:|
| AC, MIPS | 1.09× | 1.12× | 1.06× | 1.03× | 1.02× |
| AC, predictor-corrector | 1.20× | 1.15× | 1.31× | 1.26× | 1.26× |
| DC, MIPS | 1.36× | 1.47× | 1.06× | 1.07× | 1.06× |
| DC, predictor-corrector | 1.09× | 1.05× | 1.19× | 1.15× | 1.15× |

The largest case, case_SyntheticUSA, takes 122.2 s → 97.0 s on AC with the predictor-corrector
method.

**Where it came from:** samply at 1 kHz, each pipeline repeated for 10 s, folded by
`profile_summary.py`. Each stage's time per run (share × time):

| Workload | Case | KKT assembly | Allocation |
|---|---|--:|--:|
| AC, predictor-corrector | case2383wp | 84 → 7 ms | 31 → 23 ms |
| AC, predictor-corrector | case9241pegase | 449 → 43 ms | 147 → 115 ms |
| AC, predictor-corrector | case_ACTIVSg25k | 2,762 → 204 ms | 891 → 655 ms |
| AC, MIPS | case9241pegase | 471 → 30 ms | 543 → 365 ms |
| AC, MIPS | case_ACTIVSg25k | 1,814 → 117 ms | 2,439 → 1,170 ms |
| DC, MIPS | case_ACTIVSg25k | 783 → 70 ms | 1,847 → 650 ms |
| DC, predictor-corrector | case_ACTIVSg25k | 195 → 19 ms | 101 → 91 ms |

- **The KKT assembly** is 10–20× smaller on every workload: 0.2–3.2% of each solve, from
  2.9–18.7%.
- **MIPS's allocation** halves, from the LU's reserved arrays.

## What's left

**On case9241pegase's AC predictor-corrector solve,** the largest stages are now:

| Stage | Share |
|---|--:|
| LDLᵀ's dense kernel | 27.6% |
| OPF evaluation | 22.6% |
| Triangular solves | 10.5% |
| Front assembly and extend-add | 8.9% |

The dense kernel is target 4's first half, the next PR.

**MIPS on LU** is 83.4% numeric factorization. That's for O5's exploration (rayon, SIMD), since
the rule on results keeps its LU's arithmetic as it is.
