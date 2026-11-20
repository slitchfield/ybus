# A parallel LDLᵀ factorization, 2026-11-20

This is phase O5's target 2, PR 2a (`OptimizationChecklist.md`, "Phase O5"). With the
`rayon` feature, the multifrontal LDLᵀ now factors its front tree's independent subtrees in
parallel.
- **Bit for bit:** each front assembles its children's contribution blocks in their index
  order, as the serial loop did, so the results are the serial factorization's at any thread
  count.
- **Also in the code:** a front's trailing update can split across columns, bit for bit too,
  but it's off by default (see the calibration).
- **The serial path** runs on a one-thread pool, and inside `opf::batch`'s `par_` forms.

O5's profile (`o5-profile-2026-11-19.md`) put the factorization at 32–42% of the
predictor-corrector method's time, and its front trees bound subtree parallelism at ×4–9.

**Outcome:**
- **Bit for bit everywhere:** 1,740 survey result columns at 1–28 threads, and the 72-case
  snapshot in both builds.
- **The gains are modest,** on the predictor-corrector method's single solves:
  - **AC:** a median ×1.22 at 4 and 8 threads, up to ×1.34 (case_ACTIVSg10k).
  - **DC:** about ×1.0, but ×1.36 on case_ACTIVSg25k.
  - **More threads don't help:** speedups peak at 4–8 threads, and fall by 28.
- **Utilization is poor:** at 14 threads, a solve keeps only 1.3–1.6 threads busy.
- **Batches don't use it:** inside `opf::batch`'s parallel forms it gave no gain, so it's off
  there.

## The survey

- **The machine:** an Intel Core i7-14700K (28 threads, 15 GB under WSL2), Linux 6.6, release
  builds.
- **The runs:** `scripts/benchmark-opf.sh`, with the new `YBUS_BENCH_FEATURES=rayon` and
  `RAYON_NUM_THREADS` at 1, 2, 4, 8, 14 and 28, against the default build, serial, timed the
  same day.
- **The workloads:** the predictor-corrector method at its defaults (`ldlt-uncorrected`,
  matched on DC, AMD on AC), on case2383wp, case3120sp, case9241pegase, case_ACTIVSg10k and
  case_ACTIVSg25k.
- **Also measured:** samply profiles at 14 threads, and `opf_batch` on case9241pegase's AC
  outages.
- **The rows:** [`o5-parallel-ldlt-2026-11-20.csv`](o5-parallel-ldlt-2026-11-20.csv) has the
  survey, and [`o5-parallel-ldlt-2026-11-20-calibration.csv`](o5-parallel-ldlt-2026-11-20-calibration.csv)
  the calibration.
- **The one-thread rows** are from the final code, which takes the serial path on a one-thread
  pool. The rest are from the code before that change, which doesn't touch them.

## Bit for bit

- **The survey:** every result column (outcomes, iterations, agreement with MATPOWER, KKT
  residuals, the factorizations' statistics), 870 per model, matches the serial build's at
  every thread count.
- **The 72-case snapshot** is unchanged in the default build, which runs the restructured serial
  path, and in the `rayon` build.
- **The tests** run on pools of 1, 2, 7 and 28 threads, with every subtree and update forced
  onto the parallel paths:
  - random indefinite and saddle-point matrices, with 2×2 pivots and delays
  - a larger saddle point at the default cutoffs
  - a matrix singular in two independent blocks, which gives the serial order's first failure
- **Planted bugs, each caught:**
  - children assembled in reverse
  - fronts placed in arrival order
  - a column skipped by the parallel split
  - the last failure returned instead of the first

  Two needed the tests strengthened first. The singular case needed two independent failures,
  and the round trip had compared only one.

## The calibration

The two cutoffs are a subtree's work, past which it runs in parallel, and a front's trailing
order, past which its update splits by columns. They were calibrated at 14 threads, before the
survey.
- **On case2383wp alone, as planned,** the parallel factorization gains nothing on DC (the
  serial 61 ms is fastest), and ×1.14 on AC.
- **The large cases gain:** case_ACTIVSg25k ×1.44 (AC) and ×1.42 (DC), and case9241pegase
  ×1.24 (AC) and ×1.11 (DC).
- **The choice, decided 2026-11-20, over all six workloads:** subtrees of 2²⁰ multiply-adds or
  more in parallel, and no column split. That's a geometric-mean ×1.15, with the best worst
  case, ×0.97.
- **The column split** never beat the noise. From trailing order 256 it made AC solves 10–20×
  slower: dispatching rayon work on every pivot of many mid-sized fronts costs far more than it
  saves.
- **The survey's cases include the calibration's,** so their speedups are at tuned cutoffs.

## The speedups

Each solve's time against the serial build's (its medians are over repeated runs):

| Case | Serial | 1 thread | 2 | 4 | 8 | 14 | 28 |
|---|---|---|---|---|---|---|---|
| **DC** case2383wp | 66 ms | ×1.02 | ×1.03 | ×1.01 | ×0.97 | ×1.02 | ×0.96 |
| case3120sp | 45 ms | ×0.99 | ×1.00 | ×0.98 | ×0.96 | ×0.96 | ×0.95 |
| case9241pegase | 63 ms | ×0.97 | ×1.07 | ×1.05 | ×1.07 | ×1.01 | ×0.98 |
| case_ACTIVSg10k | 196 ms | ×0.91 | ×0.92 | ×0.98 | ×1.05 | ×1.01 | ×0.88 |
| case_ACTIVSg25k | 889 ms | ×1.01 | ×1.16 | ×1.17 | ×1.36 | ×1.25 | ×1.03 |
| **AC** case2383wp | 364 ms | ×0.99 | ×1.11 | ×1.06 | ×1.07 | ×1.03 | ×0.94 |
| case3120sp | 595 ms | ×1.00 | ×1.06 | ×1.08 | ×1.11 | ×1.02 | ×0.91 |
| case9241pegase | 1.91 s | ×0.98 | ×1.14 | ×1.22 | ×1.22 | ×1.17 | ×1.06 |
| case_ACTIVSg10k | 3.92 s | ×0.99 | ×1.13 | ×1.24 | ×1.34 | ×1.25 | ×1.19 |
| case_ACTIVSg25k | 11.86 s | ×0.98 | ×1.16 | ×1.26 | ×1.30 | ×1.30 | ×1.25 |

**Peak memory** rises with the threads, since more contribution blocks are alive at once. On
case_ACTIVSg25k's AC it goes from 816 MiB serially to 1,014 MiB at 28 threads, and on its DC
from 123 to 206 MiB.

## Against the bounds

The trees bounded the factorization at ×4–9, and O5's profile bounded the whole solve at
×1.5–2.1 (DC) and ×2.6–3.0 (AC). The whole-solve bound counts the AC evaluation (20–27%) and
the triangular solves (11–16%), which PR 2b and target 3 take. **What's left between the
factorization's share and these results:**
- **The factorization itself falls short of the trees' bound.** In the 14-thread profiles of
  case9241pegase's and case_ACTIVSg25k's AC, the dense kernel and front assembly still take
  32–38%, against 28–32% serially, now with only 1.3–1.6 threads busy. The root fronts are the
  largest, and their parents wait for them, so most of the factorization's time is on the
  critical path, which is serial.
- **Overhead:** on DC, rayon's own runtime takes 5–10% of the profile at 14 threads, since the
  DC factorizations are short (milliseconds).
- **Past 8 threads,** the efficiency cores and memory bandwidth add little, as O5's profile
  found for independent work.

The column split was meant to parallelize those root fronts. At the cutoffs tried, its
per-pivot dispatch cost more than it saved. A blocked kernel would split each front's update
into a few large tasks per block of pivots, rather than one per pivot, and could be what the
root fronts need. That's target 3's (SIMD and layouts).

## Inside a batch

The first measurement compared a batch with the nested parallelism on against #89's record,
from another day. It looked 10% slower. Against `main` on the same day it was at most 0–4%
slower, within the drift.

Either way the nested parallelism gave batches nothing, so `opf::batch`'s `par_` forms now run
each solve's factorizations serially (decided 2026-11-20). They do so through `ybus_sparse`'s
`serially`, a scope kept in a thread-local count. Same-day, case9241pegase's AC outages:

| | Serial batch loop | `par_` form, 14 threads |
|---|---|---|
| `main` | 55.2 s cold, 53.6 s warm | 15.2–15.5 s cold, 15.4–15.8 s warm |
| this branch | 48.7 s cold, 48.1 s warm (×1.13) | 15.1–15.2 s cold, 15.1–15.6 s warm |

- **The serial loop gains ×1.13,** since each of its solves now factors in parallel.
- **The parallel form is unchanged.**
- **On a one-thread pool,** the factorization takes its serial path, so the parallel build there
  is the serial build's speed (×0.98–1.02, the one-thread column above).
