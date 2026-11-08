# The LDLᵀ's bookkeeping kept with its analysis, 2026-11-08

This is the first half of the target that [`ipm-profile-2026-11-08.md`](ipm-profile-2026-11-08.md)
proposed first: the pivoted LDLᵀ's bookkeeping. That profile found the factorization spending
more on bookkeeping than on its arithmetic:

- **Front assembly:** 11–21%.
- **Per-pattern work redone at every factorization:** the symmetry check, the permutation and the
  supernode partition, 7.5–11.7%.
- **Allocation:** 11–27%.

This half keeps the per-pattern work with the analysis and reuses the buffers, with the results
bit for bit the same. The second half, relaxed amalgamation of the tiny fronts, changes the
pivots, so it comes separately.

## The change

- **`SymbolicIndefinite`** (`CsrMatrix::analyze_indefinite`, or `SymbolicIndefinite::new` on a
  Cholesky analysis) holds:
  - the Cholesky analysis;
  - the fundamental supernodes, with each one's parent;
  - the gather map from `A`'s stored values to `P A Pᵀ`'s, made by permuting, once, the matrix
    whose values are their own positions;
  - each stored entry's mirror.
- **`SymbolicIndefinite::factor`:**
  - **Gathering:** it gathers the values instead of permuting the matrix.
  - **The symmetry check:** one pass through the mirror map, refusing exactly what
    `is_hermitian` refuses, instead of building the adjoint. `SymbolicCholesky::factor_indefinite`
    builds one per call and runs it.
- **Buffers:**
  - **one front buffer,** zeroed for each front;
  - **L's columns** copied out at their exact size (each front used to keep its whole square's
    capacity);
  - **spent contribution blocks and variable lists,** handed out again. Blocks are still consumed
    in the order they were made, so every sum is as before;
  - **the kernel's per-pivot vectors** come from one scratch buffer.
- **The solvers' `KktSolver`** keeps a `SymbolicIndefinite`.

## The same results, checked

- **The snapshot:** a temporary snapshot of every solver's DC and AC OPF solutions on the 72
  normalized cases up to 3 MB hashes identically before and after. The AC predictor-corrector
  method factors by LDLᵀ by default, so this covers the factorization in the solvers.
- **Repeated factorizations:** with one kept analysis they give the same results, bit for bit,
  as fresh analyses and as `factor_indefinite`.
- **Planted bugs,** each failing a test:
  - the front buffer not cleared;
  - a missing mirror ignored;
  - the gather map shifted;
  - a contribution block sent to the wrong parent.

## Speed

#60's LDLᵀ profiles, rerun on the same cases: samply at 1 kHz, each pipeline repeated for 5 s
and at least once. Times per run:

| Workload | Case | Before | After | Speedup |
|---|---|--:|--:|--:|
| DC, predictor-corrector, LDLᵀ | case300 | 2.6 ms | 2.0 ms | 1.25× |
| | case2383wp | 144 ms | 109 ms | 1.32× |
| | case9241pegase | 110 ms | 95 ms | 1.16× |
| | case_ACTIVSg10k | 420 ms | 324 ms | 1.30× |
| | case_ACTIVSg70k | 8.52 s | 6.48 s | 1.31× |
| AC, predictor-corrector (LDLᵀ) | case2383wp | 834 ms | 698 ms | 1.20× |
| | case9241pegase | 4.28 s | 3.44 s | 1.25× |
| | case_ACTIVSg25k | 26.9 s | 25.3 s | 1.06× |
| | case_ACTIVSg70k | 95.4 s | 87.2 s | 1.09× |
| | case_SyntheticUSA | 196 s | 178 s | 1.10× |

**Where it came from,** each stage's time per run (share × time), before and after:

| Stage | DC case_ACTIVSg70k | AC case9241pegase | AC case_ACTIVSg70k |
|---|--:|--:|--:|
| Symmetry check and permutation | 707 → 39 ms | 437 → 34 ms | 6.68 → 0.52 s |
| Supernode partition | 204 → 0 ms | 64 → 0 ms | 1.53 → 0 s |
| Front assembly and extend-add | 1,763 → 1,121 ms | 480 → 334 ms | 11.6 → 8.6 s |
| Allocation (all of it) | 1,550 → 1,127 ms | 557 → 409 ms | 10.8 → 9.4 s |
| Dense kernel (unchanged) | 1,243 → 1,218 ms | 360 → 354 ms | 16.3 → 16.8 s |

- **The per-pattern work is gone.** What's left of the symmetry check and the permutation
  (0.4–1.0% of the run) is the one pass through the mirror map, and the gather.
- **Front assembly** takes 22–36% less, and allocation 10–38% less.
- **The largest AC cases gain least:** 1.06–1.10×.
  - **Their factorization's savings** are about the size of mid-size AC's.
  - **But the untouched stages read slower:** this profile measured the OPF's evaluation and KKT
    assembly 2–24% slower than #60's on the same cases. Neither changed, and each case_ACTIVSg25k,
    case_ACTIVSg70k and case_SyntheticUSA profile is a single run, so that difference is
    run-to-run noise. It hides part of the gain.

## What's left

- **On DC, allocation is still 17–28% of the run.** Its largest remaining sources are outside the
  factorization: the KKT assembly's triplets (`from_triplets` is 5.5–9.1% of samples by itself on
  every case) and the solves.
- **The fronts are still tiny** (a median of 4 rows on DC), so front assembly and the kernel's cost
  per front remain. That is the second half: relaxed amalgamation.
- **On AC, the OPF's evaluation** (32–45% now) is the largest stage. It is #60's second target.
