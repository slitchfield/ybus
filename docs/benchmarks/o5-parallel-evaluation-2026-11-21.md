# The AC evaluation in parallel, 2026-11-21

This is phase O5's target 2, PR 2b (`OptimizationChecklist.md`, "Phase O5"), after PR 2a's
parallel LDLᵀ (`o5-parallel-ldlt-2026-11-20.md`).
- **The change:** with the `rayon` feature, on networks of 1,000 buses or more, the AC OPF's
  evaluation runs its independent parts concurrently with `rayon::join`:
  - the bus flows beside the branch ends' flows and magnitudes
  - the Lagrangian Hessian's bus part beside its branch ends'
  - each branch end's outer product beside its weighted Hessian
  - the weighted Hessian's real and imaginary parts
  - the bilinear Hessian's three blocks
- **Bit for bit:** each part is the existing code on the same inputs, combined in the same
  order.
- **Where it runs serially:** as PR 2a's LDLᵀ does, on a one-thread pool and inside
  `opf::batch`'s `par_` forms.

**Outcome:** with PR 2a, the predictor-corrector method's AC solves reach a median ×1.50 at 8
threads, up to ×1.67 (case_ACTIVSg25k at 14). That's from PR 2a's ×1.15 alone, measured the
same day.
- **This PR's own gain is about ×1.3.** That's above its planned bound of ×1.1–1.25, which
  assumed it would share the threads with the factorization rather than add to them.
- **MIPS gains ×1.03,** since its time is LU.
- **Bit for bit throughout:** 2,320 result columns, and the 72-case snapshot in both builds.

## The survey

- **The machine:** an Intel Core i7-14700K (28 threads, 15 GB under WSL2), Linux 6.6, release
  builds, one run at a time.
- **The runs:** `scripts/benchmark-opf.sh` (`YBUS_BENCH_FEATURES=rayon`, `RAYON_NUM_THREADS` at
  1, 2, 4, 8, 14 and 28) on the AC OPF, both solvers at their defaults. The predictor-corrector
  method used `ldlt-uncorrected` with AMD, and MIPS used LU.
- **The cases:** case2383wp, case3120sp, case9241pegase, case_ACTIVSg10k and case_ACTIVSg25k.
- **The comparisons, all the same day:**
  - the default build, serial
  - PR 2a's branch (`069a040`) at 8 and 14 threads
  - `main`'s 14-thread batch
- **The rows:** [`o5-parallel-evaluation-2026-11-21.csv`](o5-parallel-evaluation-2026-11-21.csv)
  has the survey, and [`…-calibration.csv`](o5-parallel-evaluation-2026-11-21-calibration.csv)
  the cutoff's calibration.

## Bit for bit

- **The survey:** every result column, 1,160 per solver across all thread counts and PR 2a's
  runs, matches the serial build's.
- **The 72-case snapshot** is unchanged in the default and the `rayon` build.
- **The tests:**
  - each plan built `.parallel(true)` is its serial self, declines included, on pools of 1, 2,
    7 and 28 threads
  - case2383wp's AC solve is bit for bit at 2, 7 and 28 threads against one
- **Planted bugs, each caught:**
  - the real and imaginary parts swapped
  - the blocks swapped
  - the from and to ends' Hessians swapped, also in 20 of the snapshot's hash lines
- **Not plantable:** ignoring one part's `None`. With finite voltages, the real and imaginary
  parts are `None` together, so one part alone failing isn't reachable. The joined code returns
  `None` if any part does, the serial result by construction.

## The cutoff

The evaluation goes concurrent from `EVALUATION_BUSES`, 1,000 buses, calibrated before the
survey at 14 threads, each case on and off:

| Case | MIPS, on/off | Predictor-corrector, on/off |
|---|---|---|
| case300 | ×0.87 | ×0.75 |
| case1354pegase | ×1.02 | ×1.05 |
| case2383wp | ×1.02 | ×1.10 |
| case9241pegase | ×1.02 | ×1.13 |
| case_ACTIVSg25k | ×1.01 | ×1.18 |

By PR 2a's rule, the cutoff lies above 300 buses and at most 1,354.

## The speedups

These are against the serial build. The last two columns are PR 2a's branch alone.

**The predictor-corrector method:**

| Case | Serial | 1 thread | 2 | 4 | 8 | 14 | 28 | PR 2a, 8 | PR 2a, 14 |
|---|---|---|---|---|---|---|---|---|---|
| case2383wp | 362 ms | ×0.99 | ×1.25 | ×1.21 | ×1.32 | ×1.07 | ×1.17 | ×0.99 | ×0.97 |
| case3120sp | 580 ms | ×0.98 | ×1.29 | ×1.29 | ×1.41 | ×1.03 | ×1.21 | ×1.02 | ×0.97 |
| case9241pegase | 1.87 s | ×1.02 | ×1.32 | ×1.34 | ×1.50 | ×1.35 | ×1.27 | ×1.15 | ×1.13 |
| case_ACTIVSg10k | 3.93 s | ×1.04 | ×1.36 | ×1.56 | ×1.63 | ×1.65 | ×1.55 | ×1.27 | ×1.25 |
| case_ACTIVSg25k | 11.51 s | ×1.01 | ×1.26 | ×1.56 | ×1.55 | ×1.67 | ×1.56 | ×1.35 | ×1.37 |
| **Median** | | ×1.01 | ×1.29 | ×1.34 | ×1.50 | ×1.35 | ×1.27 | ×1.15 | ×1.13 |

**MIPS:** a median ×1.03 at every thread count from 2 to 28, and ×1.00 at one. PR 2a alone is
×0.98–0.99 on MIPS, since MIPS's LU doesn't use the parallel LDLᵀ.

- **The best thread count is 8:** 14 and 28 give less on the smaller cases. That matches PR 2a's
  and O5's profile's finding that the efficiency cores and memory bandwidth add little past 8.
- **PR 2a's gains on the smaller cases** were ×1.06–1.11 at 8 threads two days ago, and ×0.97–1.02
  today. The day-to-day spread is as large as those gains, so they're noise. This record compares
  only same-day runs.

## Where the time goes now

These are the 14-thread profiles of the predictor-corrector method's AC solves:

| Case | Threads busy (PR 2a → 2b) | Largest stages |
|---|---|---|
| case9241pegase | 1.5 → 2.0 | evaluation 29.9%, dense kernel 16.2%, front assembly 11.1%, triangular solves 8.9%, allocation 8.1%, parallel runtime 6.7% |
| case_ACTIVSg25k | 1.6 → 2.4 | evaluation 37.6%, dense kernel 20.3%, front assembly 9.9%, triangular solves 9.8%, allocation 6.3% |

A share in a parallel profile counts the samples on every thread, so the evaluation's share
rises as its work spreads across threads, while its wall time falls.
- **What limits it now:** the Hessian splits into at most about a dozen parts at once, with the
  bus part the largest, so the evaluation's critical path is that part's own serial chain. A
  finer split, by element ranges, is what PR 2a's lesson advised against.
- **What's left for target 3:** the triangular solves (9–10%), serial throughout, and the
  root fronts' dense kernel.

## Inside a batch

case9241pegase's AC outage batch, 28 outages, the same day:

| | Serial batch loop | `par_` form, 14 threads |
|---|---|---|
| `main` | 55.0 s cold, 53.4 s warm | 14.9 s cold, 14.9 s warm |
| this branch | 40.3 s cold, 41.9 s warm (×1.37, ×1.27) | 15.1 s cold, 14.7 s warm |

- **The parallel form is unchanged:** its solves run their factorizations and evaluations
  serially.
- **The serial loop gains ×1.27–1.37,** since each of its solves now runs both parallel parts.

## Against the bound

O5's profile bounded parallelism within one solve at ×2.6–3.0 on the predictor-corrector
method's AC. Targets 2a and 2b together reach ×1.50 at 8 threads, about half of it.
- **The evaluation's critical path:** its largest part, the bus Hessian.
- **The factorization's critical path:** the root fronts.
- **The serial triangular solves,** which are target 3's.
- **Rayon's runtime,** 5–7% at 14 threads.

**Target 2 closes with this PR.**
