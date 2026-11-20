# Batched OPF solves, 2026-11-20

This is phase O5's target 1 (`OptimizationChecklist.md`, "Phase O5"). Many OPF solves
sharing a base now run in parallel, behind the optional `rayon` feature. There are two
batches:
- **Outages** (`ac_opf_outages`, `dc_opf_outages`, and their `par_` forms): a base case and its
  outages (`OpfOutage`, any set of branches and generators). Each outage gets its own model,
  cold or warm from the base solution.
- **Cases on one network** (`solve_many` and `par_solve_many` on both models): one kept model
  solving many cases with new values from one start.

O5's profile (`o5-profile-2026-11-19.md`) bounded this target at about ×10, the machine's
ceiling for independent work.

**Outcome:**
- **Every parallel batch is bit for bit its serial one,** at every thread count surveyed (321
  rows, 1 to 28 threads).
- **MIPS reaches the ceiling:** a median ×9.35 (DC) and ×7.69 (AC) at 28 threads, up to ×11.6.
- **The predictor-corrector method stops near ×4–5.5,** at 14 threads, and falls back at 28.
  That's not load imbalance: its cold batches could reach ×22–27 by that measure.
- **Warm outage batches** reach ×2–3.5 on several cases, held back by their slowest outage.

## The survey

- **The machine:** an Intel Core i7-14700K (8 performance cores with two threads each, and 12
  efficiency cores: 28 threads, 15 GB under WSL2), Linux 6.6, the release build.
- **The tool:** `scripts/benchmark-opf-batch.sh`, the `opf_batch` example.
  - **The serial batch** is timed once per process.
  - **The `par_` form** is then timed on rayon pools of 1, 2, 4, 8, 14 and 28 threads, once each.
  - **Every parallel result** is compared with the serial one through its `Debug` form, and a
    difference fails the run.
- **The batches:**
  - `outages`: a seeded sample of 14 branch outages (not bridges) and 14 generator outages
  - `hours`: the walk's 24 load levels through `solve_many`
- **The starts:** cold, and warm at each solver's own `warm_shift` (1e-4).
- **The cases:**
  - case2383wp, case3120sp and case9241pegase: DC and AC, both solvers, both batches
  - case_ACTIVSg25k: the predictor-corrector method, with 7 outages of each kind
- **The memory limit:** 12 GiB per run, below this machine's 15 GB. Every thread's solve holds
  its own model and factorization.
- **The rows** are in [`o5-batched-solves-2026-11-20.csv`](o5-batched-solves-2026-11-20.csv).

## Bit for bit

All 321 rows are `bit_for_bit: true`, so every batch, at every thread count, gave the serial
batch's results. The tests check the same at 1, 2, 7 and 28 threads, and caught three planted
bugs:
- results collected out of order
- a mapped start shared across outages
- an outage started from the previous outage's solution

## The speedups

These are medians over each model and solver's batches (outages and hours, cold and warm, all
cases):

| Model, solver | 2 threads | 4 | 8 | 14 | 28 |
|---|---|---|---|---|---|
| DC, MIPS | ×1.91 | ×3.57 | ×5.78 | ×7.75 | ×9.35 |
| AC, MIPS | ×1.95 | ×3.54 | ×5.55 | ×6.98 | ×7.69 |
| DC, predictor-corrector | ×1.87 | ×3.21 | ×4.82 | ×5.56 | ×5.00 |
| AC, predictor-corrector | ×1.85 | ×3.09 | ×4.01 | ×4.13 | ×3.68 |

Some batches, with their serial times and their best speedups:

| Batch | Serial | Best |
|---|---|---|
| case9241pegase AC, MIPS, 24 hours cold | 181 s | ×9.06 at 28 (20 s) |
| case9241pegase AC, MIPS, 28 outages warm | 261 s | ×7.09 at 28 (37 s) |
| case3120sp AC, MIPS, 24 hours cold | 41.7 s | ×9.05 at 28 |
| case2383wp DC, MIPS, 24 hours warm | 6.9 s | ×11.56 at 28 |
| case9241pegase AC, predictor-corrector, 28 outages cold | 52.8 s | ×3.83 at 14 |
| case_ACTIVSg25k AC, predictor-corrector, 24 hours cold | 265 s | ×4.04 at 14 (66 s) |
| case_ACTIVSg25k DC, predictor-corrector, 14 outages cold | 13.1 s | ×3.93 at 14 |

The full table, per case, batch and start, is in the CSV.

## Why some batches fall short of ×10

### Warm outage batches: their slowest outage

A batch can finish no sooner than its slowest solve. Cold solves take similar times, but warm
ones don't: an outage that moves the solution far takes many more iterations than one that
barely moves it.
- **How it was measured:** each outage of the survey's samples was timed alone (a temporary
  probe, not committed).
- **The bound:** the batch's serial time over its slowest outage's.

| Batch | Slowest outage / median | Bound | Measured |
|---|---|---|---|
| case2383wp AC, MIPS, warm | 4.15 s / 0.33 s | ×3.92 | ×2.92 |
| case2383wp AC, predictor-corrector, warm | 1.60 s / 0.16 s | ×4.33 | ×3.05 |
| case3120sp DC, MIPS, warm | 0.29 s / 0.06 s | ×6.79 | ×2.15 |
| case2383wp DC, predictor-corrector, warm | 0.04 s / 0.01 s | ×8.33 | ×3.30 |
| case2383wp AC, MIPS, cold | 1.37 s / 1.09 s | ×22.7 | ×7.24 |
| case9241pegase AC, predictor-corrector, cold | 1.84 s / 1.74 s | ×26.5 | ×3.83 |

- **For the warm batches,** the slowest outage explains most of the shortfall.
- **For the small ones,** such as case3120sp's DC on MIPS, warm (2.7 s serially), rayon's overhead
  and the base solve, which runs before the outages and serially, take the rest.
- **A batch of more outages would balance better,** since the slow ones would overlap with
  more fast ones.

### The predictor-corrector method: not the balance

Its cold batches could reach ×22–27 by the balance bound, but stop at ×4–5.5, with no gain or
a loss past 14 threads. MIPS's cold batches on the same cases reach ×7–11. **This isn't
measured here, but the likeliest causes are:**
- **Memory bandwidth.** The predictor-corrector method's time is dense-front kernels and
  triangular solves over large factors, which stream memory where MIPS's LU spends longer in
  each column. The parallel DC N-1 stopped at ×5 on case_ACTIVSg70k for the same reason
  (O5's profile).
- **Allocation.** Allocation and copying are 5–12% of its time (O5's profile), and many
  threads allocating at once contend in the allocator.

Target 3's layouts and allocation work is where either would show. A profile of a parallel
batch could tell them apart before then.

### Memory

case_ACTIVSg25k's AC batches went over the 12 GiB limit: the outages beyond 8 threads, and the
hours beyond 14. Each thread holds a model and an LDLᵀ factorization of a 25,000-bus KKT
system. On this 15 GB machine that caps how many large AC solves run at once. A caller with a
large case should size the pool to its memory, for example with `rayon::ThreadPoolBuilder`.

## Against the bound

O5's profile bounded batched solves at about ×10 on 14 threads.
- **MIPS** reaches ×7–11 at 28 threads, the bound, where its batches are balanced.
- **The predictor-corrector method** reaches about half of it, for reasons the balance doesn't
  explain.
- **Warm outage batches** are held to ×2–4 by their slowest outages.

Nothing was chosen: the rule recorded the speedups, since target 1 adds APIs and changes no
default.
