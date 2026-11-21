# Allocation, and the wider-SIMD build, 2026-11-21

This is phase O5's target 3, PR 3a (`OptimizationChecklist.md`, "Phase O5"). It has two parts,
both bit for bit:
- **Fewer allocations** in the LDLᵀ's factorizations.
- **A measurement of the `x86-64-v3` build,** which compiles the existing code for AVX2.

**Outcome:** both are bit for bit, and neither is faster.
- **The allocation work** halves the allocator's share of the profiles, but gains nothing
  measurable against `main`: medians ×0.98–1.03, within run-to-run noise. It's kept (decided
  2026-11-21). It costs little, and `SymbolicLdlt::refactor` is useful to any caller factoring
  a sequence.
- **The `x86-64-v3` build is slower:** the predictor-corrector method's AC solves are ×0.80–0.91
  of the default build's, in two repeated rounds, because the dense kernel takes 1.83× as long.
  So the flag is documented as bit for bit, but not recommended for these OPFs.

## Finding the allocations

The predictor-corrector method's DC and AC solves on case9241pegase and case_ACTIVSg25k were
profiled with samply at 1 and 8 threads. The archived `allocation-callers.py` attributed the
"allocation and copying" stage (5.8–15.2% of the samples) to the nearest ybus caller. The
largest sites were all in the LDLᵀ:
- **`factor_node`'s copies** (1.1–2.2%): each front's `L` into a fresh `Vec`, and its variables
  collected anew.
- **`RawVec::grow_one` in `factor_front`** (0.8–2.2%): `D`'s blocks grown pivot by pivot.
- **`drop_glue` in the solver loop** (1.0–2.0%): the previous iteration's factorization dropped,
  every front's buffers freed.
- **The rest:** smaller sites, and on DC some once-per-solve setup (`Split::new`, the AMD
  ordering, parsing) that shows only because DC solves are short.

## The fix

- **`SymbolicLdlt::refactor`,** a new public method: a factorization into an earlier `Ldlt`'s
  storage. The old fronts' `L`, variables and blocks are handed out again through the
  workspace's spare lists, and `factor_front` takes its order and blocks buffers, with the
  blocks reserved.
- **The KKT solver keeps the last factor its solver recycles** (`KktSolver::recycle`, from both
  solvers, including the inertia correction's rejected factors), and refactors into it.
- **Only the allocations change,** so the results are bit for bit.

**After,** at one thread, the allocation stage falls:
- case9241pegase's AC, 7.4% → 4.0%
- case_ACTIVSg25k's DC, 9.8% → 4.9%
- case_ACTIVSg25k's AC, 5.8% → 2.9%

case9241pegase's DC falls less, 15.1% → 13.0%, since its allocations are mostly the
once-per-solve setup.

## The survey

- **The machine:** an Intel Core i7-14700K (28 threads, 15 GB under WSL2), Linux 6.6.
- **The runs:** `scripts/benchmark-opf.sh`, on case2383wp, case3120sp, case9241pegase,
  case_ACTIVSg10k and case_ACTIVSg25k, DC and AC, the predictor-corrector method
  (`ldlt-uncorrected`, matched on DC, AMD on AC) and MIPS (LU).
- **The builds,** all with `rayon`, at 1 and 8 threads:
  - `main`
  - this branch
  - this branch with `RUSTFLAGS="-C target-cpu=x86-64-v3"`
  
  Also `main`'s serial build, as the reference. All the same day.
- **The rows:** [`o5-allocation-2026-11-21.csv`](o5-allocation-2026-11-21.csv).

### Bit for bit

- **The survey:** every result column of the three builds, at both thread counts, matches
  `main`'s serial build: 2,280 compared, none differing.
- **The 72-case snapshot** is unchanged in the default, `rayon` and `x86-64-v3` builds.
- **No fused instructions:** the `x86-64-v3` binary has none (`vfmadd` and the rest). Rust
  doesn't let LLVM contract a multiply and an add unasked, so wider vectors don't change the
  arithmetic.

### The allocation work: no measurable gain

This branch against `main`, medians over the five cases:

| | DC, predictor-corrector | DC, MIPS | AC, predictor-corrector | AC, MIPS |
|---|---|---|---|---|
| 1 thread | ×0.99 | ×1.01 | ×1.00 | ×0.98 |
| 8 threads | ×1.02 | ×1.03 | ×0.99 | ×1.01 |

Single cases range from ×0.88 to ×1.11 either way, about the run-to-run noise found earlier. A
quick same-day check before the survey also put it at about ×1.02–1.03. The allocator's fast
path is cheap, so halving its samples saves little time.

### The `x86-64-v3` build: slower

| | DC, predictor-corrector | DC, MIPS | AC, predictor-corrector | AC, MIPS |
|---|---|---|---|---|
| 1 thread | ×0.95 | ×1.01 | ×0.88 | ×0.96 |
| 8 threads | ×0.96 | ×0.97 | ×0.95 | ×0.98 |

These are medians against the same branch built normally. Repeated twice more, alternating
the builds, the predictor-corrector method's AC solves at one thread were consistently slower:

| Case | Round 1 | Round 2 |
|---|---|---|
| case3120sp | ×0.90 | ×0.91 |
| case_ACTIVSg10k | ×0.80 | ×0.83 |
| case_ACTIVSg25k | ×0.82 | ×0.82 |

**The cause, by profile** (case_ACTIVSg10k's AC, one thread): the dense kernel takes 991 ms in
the default build and 1,812 ms in `x86-64-v3`, while every other stage is unchanged
(evaluation 1,041 against 1,095 ms, solves 403 against 454 ms).
- **The likely reason,** inferred, not measured instruction by instruction: the kernel updates
  short columns, and the median AC front has 7 rows. AVX2's wider loops spend their time in
  their setup and scalar remainders, where SSE2's two lanes, which the default build already
  uses, fit better.
- **What it means for later:** a wider kernel would need restructuring for short columns before
  the flag could pay. PR 3c's fused arithmetic, which needs the flag for FMA, stays planned but
  unstarted (decided 2026-11-21).

## Documented

The README's Development section notes the flag: it's bit for bit, and slower on these OPFs,
so not recommended.
