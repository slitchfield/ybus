# The LU refactorization kernel, 2026-11-05

Phase O2's second CPU target, chosen from [`profile-2026-11-04.md`](profile-2026-11-04.md):
`Lu::refactor` takes 41–64% of Newton's method with reactive limits enforced one at a time, and
25% of a large Newton solve. The profile bounded the gain at 1.3–1.5× for reactive limits if
the kernel's cost halved. This record measures what bit-identical changes reach, and what they
don't.

## Setup

| | |
|---|---|
| Code | The `dev/lu-refactor` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K (28 hardware threads), Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Kernel | A temporary benchmark (not committed): each case's Newton Jacobian at its case voltages, factored with AMD, then refactored repeatedly with two alternating sets of values on the same pattern; with a hash of every factor value to check bit-identity |
| Workloads | The `profile_pf` example on `main` and on this branch, built alike, run back to back |

Run-to-run noise on the kernel timings was about 10%; each figure is from two or three runs.

## What the kernel does, and what limits it

For each column, `refactor` scatters the matrix's column into a dense vector `x`, subtracts the
earlier columns of `L` in the order of `U`'s pattern, checks the reused pivot, and writes the new
values. Its cost:

| Jacobian | n | nnz(L) | Multiply-subtracts | Refactorization | Per multiply-subtract |
|---|--:|--:|--:|--:|--:|
| case2383wp | 4,438 | 24,787 | 157,603 | 0.19 ms | 1.3 ns |
| case9241pegase | 17,036 | 116,487 | 1,155,245 | 1.39 ms | 1.2 ns |
| case_ACTIVSg70k | 134,104 | 1,213,983 | 47,358,403 | 35.7 ms | 0.72 ns |
| case_SyntheticUSA | 156,255 | 1,389,814 | 50,104,825 | 39.1 ms | 0.79 ns |

At under 1.5 ns per multiply-subtract, the time goes to indirect gathers and scatters through `x`,
not to arithmetic. Four ways to speed it up were measured:

- **Keeping its buffers** between calls (the dense column and the values of `L` and `U`, which it
  used to allocate and zero each time) made no measurable difference: allocation was a small
  part of its cost. It's kept, since it costs nothing.
- **Level-scheduled parallelism** can't pay. Column `k` depends on the columns in `U`'s column
  `k`, and in that dependency graph most of the work sits in a long serial chain of separator
  columns near the top of the elimination tree. On case_ACTIVSg70k there are 1,422 levels, and
  84% of the work is in levels narrower than 28 columns. With ideal, zero-overhead threads, the
  speedup would be at most 1.4× (2.3× on case9241pegase), before paying a barrier per level.
- **Fusing nested column pairs** works. On these Jacobians, 85–98% of the work comes in pairs of
  columns of `L` where column `j + 1`'s rows below its diagonal are column `j`'s after row `j + 1`
  (each bus's angle and magnitude). Such a pair now updates together, in one pass over the
  shared rows, each row taking column `j`'s term and then column `j + 1`'s: the same operations in
  the same order, so the factors are bit for bit unchanged, with half the scatter traffic.
- **Bounds checks** on `x` cost about 10% more (25.6 ms against 28.2 ms on case_ACTIVSg70k, with
  unchecked indexing tried and reverted). Removing them would need `unsafe`, which this
  repository doesn't use, so they stay.

## Results

| Jacobian | Before | After | Speedup |
|---|--:|--:|--:|
| case2383wp | 0.19 ms | 0.19 ms | none measurable |
| case9241pegase | 1.39 ms | 1.25–1.30 ms | 1.1× |
| case_ACTIVSg70k | 35.7 ms | 28.0–28.5 ms | 1.27× |
| case_SyntheticUSA | 39.1 ms | 32.9–33.1 ms | 1.19× |

Every factor value hashes the same before and after, and the power flow tests (MATPOWER's
iteration counts and agreement, with and without reactive limits) pass unchanged.

| Workload | `main` | This branch | Speedup |
|---|--:|--:|--:|
| Newton, case_ACTIVSg70k | 711 ms | 685 ms | 1.04× |
| Newton, case_SyntheticUSA | 823 ms | 779 ms | 1.06× |
| Reactive limits one at a time, case2383wp | 463 ms | 455 ms | 1.02× |
| Reactive limits one at a time, case9241pegase | 1.79 s | 1.76 s | 1.02× |
| Reactive limits one at a time, case_ACTIVSg70k | 565 s | 514 s | 1.10× |

- **Below the profile's bound.** The bound assumed halving the kernel's cost; bit-identical
  changes reach 1.19–1.27× on the largest Jacobians, which is 10% on the workload this target was
  for.
- **What would go further** changes the arithmetic or the algorithm: unchecked indexing (about
  10% more), parallelism within each large column (splitting its rows among threads keeps each
  row's order, but needs a column large enough to pay for it), or, for one-at-a-time reactive
  limits, fewer refactorizations, by low-rank updates of the factorization as buses change type,
  which would no longer follow MATPOWER's sequence bit for bit. These were set aside with the
  maintainer (2026-11-05) in favour of moving on to parsing.
