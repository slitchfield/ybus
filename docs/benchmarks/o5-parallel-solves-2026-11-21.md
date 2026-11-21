# The LDLᵀ's solves in parallel, 2026-11-21

This is phase O5's target 3, PR 3b (`OptimizationChecklist.md`, "Phase O5"). With the `rayon`
feature, the multifrontal LDLᵀ's triangular solves now run the front tree's independent
subtrees in parallel, as PR 2a's factorization does.
- **Bit for bit:** every entry of the vector sees the serial solve's operations in the serial
  order.
- **Where it runs serially:** as the factorization does, on a one-thread pool, inside
  `opf::batch`'s `par_` forms, and for a tree not in postorder. Also when the whole tree's
  sweep is below the cutoff, which is the case below case_ACTIVSg25k's size on the OPF.

**Outcome:** bit for bit throughout (2,320 result columns, and the 72-case snapshot in both
builds).
- **The gain is on the largest case only:** case_ACTIVSg25k's AC solve is ×1.12 at 8 threads
  and ×1.11 at 14, against `main` the same day. That's at the plan's bound of ×1.1–1.15.
- **Everywhere else, ×0.96–1.02,** within the noise. At the calibrated cutoff, the smaller
  cases' solves don't go parallel.

## How it stays bit for bit

The fronts follow a postorder, so each subtree is a run of fronts and of columns, and sibling
subtrees never touch each other's columns. A pivot delayed out of a subtree stays in its
columns, and an ancestor eliminates it after the subtree has finished.
- **The forward sweep:** a task updates the entries in its own columns in place. Each update
  of an ancestor's entry is recorded, `(row, term)`, with the term computed as the serial
  sweep computes it. A parent applies its children's records in index order before its own
  front, which is the serial order for every entry, and passes up those above it.
- **The diagonal solve:** each front's `D` solve follows its own forward step, which is its
  entries' last update, rather than a separate pass. The serial path shares the code.
- **The backward sweep:** top down. A subtree's front goes first, then its children's subtrees
  in parallel, reading their ancestors' entries, which are final by then.
  - **A change from the plan:** a subtree's ancestors' columns aren't contiguous, since its
    parent's siblings lie between them, and they're being written concurrently. So each task
    reads them through a chain of read-only slices, one per parallel level.
- **Safe Rust throughout:** disjoint mutable slices, rayon's ordered collects, and no unsafe
  code.

**Tested** against the serial solve, on pools of 1, 2, 7 and 28 threads, with every subtree
forced onto the parallel path and at the default cutoff. The tests cover single, many and
refined right-hand sides:
- random indefinite and saddle-point matrices, with 2×2 pivots and delays
- a larger saddle point
- two independent blocks, for several roots
- a Hermitian matrix, for the backward sweep's conjugates

**Planted bugs, each caught:**
- children's records applied in reverse
- the outgoing records dropped
- the `D` solve before a front's forward step
- the backward sweep's children run before their parent

## The cutoff

Subtrees of 2²⁰ sweep multiply-adds or more go parallel, the factorization's own cutoff.
- **The calibration:** before the survey, at 8 and 14 threads, two rounds each, on the
  predictor-corrector method's DC and AC solves of case2383wp, case9241pegase and
  case_ACTIVSg25k. A second calibration extended it past 2²⁰. The rows are in
  [`o5-parallel-solves-2026-11-21-calibration.csv`](o5-parallel-solves-2026-11-21-calibration.csv).
- **Small cutoffs cost far more than they save:** at 2¹², case2383wp's DC solve took ×2.6 as
  long (×0.38), and below 2¹⁶ every small case lost.
- **The one steady gain** is case_ACTIVSg25k's AC at 2²⁰, ×1.09–1.14 in all four rounds. It's
  gone from 2²².
- **The choice:**

  | Cutoff | Geometric mean, calibration 1 | Worst | Calibration 2 | Worst |
  |---|---|---|---|---|
  | 2¹⁶ | ×0.97 | ×0.90 | | |
  | 2¹⁸ | ×1.03 | ×0.97 | | |
  | 2²⁰ | ×1.03 | ×0.99 | ×1.02 | ×0.95 |
  | 2²² | | | ×1.00 | ×0.97 |
  | 2²⁴ | | | ×1.01 | ×0.98 |

  By PR 2a's rule, the setting with the best worst case within 2% of the best geometric mean,
  the first calibration chooses 2²⁰ and the second 2²⁴. The second's worst cases are the small
  cases' noise (±5%): 2²⁰'s ×0.95 was case2383wp's 56 ms DC, at ×1.00 in the first. 2²⁰ was
  taken, for the steady gain.

## The survey

- **The machine:** an Intel Core i7-14700K (28 threads, 15 GB under WSL2), Linux 6.6.
- **The runs:** `scripts/benchmark-opf.sh`, the predictor-corrector method (`ldlt-uncorrected`,
  matched on DC, AMD on AC), on case2383wp, case3120sp, case9241pegase, case_ACTIVSg10k and
  case_ACTIVSg25k. MIPS factors by LU, so it's unaffected and wasn't surveyed.
- **The builds, all the same day:**
  - this branch, serial
  - this branch with `rayon` at 1, 2, 4, 8, 14 and 28 threads
  - `main` (target 2 and PR 3a) with `rayon` at 8 and 14
- **The rows:** [`o5-parallel-solves-2026-11-21.csv`](o5-parallel-solves-2026-11-21.csv).

**Bit for bit:** every result column, 1,160 per model across all eight parallel runs,
`main`'s included, matches the serial build's.

### The speedups

These are against the serial build. The last two columns are this branch over `main` at the
same thread count.

**AC:**

| Case | Serial | 1 thread | 2 | 4 | 8 | 14 | 28 | `main`, 8 | `main`, 14 | Over `main`, 8 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| case2383wp | 347 ms | ×1.03 | ×1.17 | ×1.28 | ×1.23 | ×1.14 | ×1.17 | ×1.22 | ×1.18 | ×1.01 | ×0.97 |
| case3120sp | 573 ms | ×1.04 | ×1.17 | ×1.37 | ×1.37 | ×1.23 | ×1.19 | ×1.35 | ×1.28 | ×1.01 | ×0.97 |
| case9241pegase | 1.73 s | ×0.95 | ×1.11 | ×1.26 | ×1.33 | ×1.30 | ×1.16 | ×1.32 | ×1.32 | ×1.01 | ×0.99 |
| case_ACTIVSg10k | 3.60 s | ×0.98 | ×1.15 | ×1.39 | ×1.56 | ×1.54 | ×1.42 | ×1.53 | ×1.53 | ×1.02 | ×1.01 |
| case_ACTIVSg25k | 11.41 s | ×1.00 | ×1.29 | ×1.55 | ×1.74 | ×1.74 | ×1.66 | ×1.56 | ×1.57 | ×1.12 | ×1.11 |
| **Median** | | ×1.00 | ×1.17 | ×1.37 | ×1.37 | ×1.30 | ×1.19 | ×1.35 | ×1.32 | ×1.01 | ×0.99 |

**DC:**

| Case | Serial | 1 thread | 2 | 4 | 8 | 14 | 28 | `main`, 8 | `main`, 14 | Over `main`, 8 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| case2383wp | 55 ms | ×0.98 | ×0.97 | ×0.97 | ×0.96 | ×0.95 | ×0.96 | ×0.97 | ×0.99 | ×0.99 | ×0.96 |
| case3120sp | 40 ms | ×0.99 | ×0.98 | ×0.95 | ×0.94 | ×0.97 | ×0.92 | ×0.98 | ×0.98 | ×0.96 | ×0.99 |
| case9241pegase | 59 ms | ×0.99 | ×1.00 | ×1.04 | ×1.04 | ×1.02 | ×1.00 | ×1.02 | ×1.02 | ×1.02 | ×1.00 |
| case_ACTIVSg10k | 189 ms | ×1.03 | ×0.92 | ×1.03 | ×1.14 | ×1.02 | ×0.91 | ×1.12 | ×1.05 | ×1.02 | ×0.97 |
| case_ACTIVSg25k | 832 ms | ×1.01 | ×0.97 | ×1.16 | ×1.26 | ×1.21 | ×1.11 | ×1.25 | ×1.25 | ×1.01 | ×0.97 |
| **Median** | | ×0.99 | ×0.97 | ×1.03 | ×1.04 | ×1.02 | ×0.96 | ×1.02 | ×1.02 | ×1.01 | ×0.97 |

- **case_ACTIVSg25k's AC now reaches ×1.74** at 8 and 14 threads, from `main`'s ×1.56–1.57.
- **Its DC solve gains nothing over `main`** (×0.97–1.01). Its sweeps are a smaller share of a
  DC solve, and its tree is smaller.
- **The other cases are `main`'s speed,** within the noise.

## Where the time goes now

These are the 14-thread profiles of the predictor-corrector method's AC solves, `main` against
this branch, the same day.

| Case | Wall time | Threads busy | Triangular solves' share of the samples |
|---|---|---|---|
| case9241pegase | 1.33 → 1.31 s | 2.0 → 2.0 | 8.9% → 9.3%, all serial |
| case_ACTIVSg25k | 7.82 → 6.88 s | 2.4 → 2.6 | 10.8% → 11.9%, now on several threads |

- **case9241pegase's solves never take the parallel path:** none of the parallel sweeps' frames
  appear in its profile. Its tree's whole sweep is below the cutoff.
- **On case_ACTIVSg25k,** the forward and backward sweeps take 4.7% and 5.7% of the samples.
  A share in a parallel profile counts the samples on every thread, so it can rise while the
  solves' wall time falls.
- **A fix to the profile's stage rules:** `scripts/profile_summary.py` first missed the
  parallel sweeps. Their frames (`multifrontal::Ldlt>::forward_run`) carry no generics, and on
  rayon's threads the `Factorization` trait's frames that catch the serial solve aren't on the
  stack. It counted them as rayon's runtime (14.3%), and now counts them as solves (rayon's
  runtime 2.8%).

## Inside a batch

case9241pegase's AC outage batch, 28 outages, the same day:

| | Serial batch loop | `par_` form, 14 threads |
|---|---|---|
| `main` | 46.5 s cold, 45.8 s warm | 14.8 s cold, 15.5 s warm |
| this branch | 42.9 s cold, 44.9 s warm | 14.7 s cold, 14.9 s warm |

- **The parallel form is unchanged:** its solves run serially.
- **The serial loop's difference is noise:** case9241pegase's solves don't go parallel.
- **The warm rows are failing solves:** the batch reports none of the 28 converging. That's
  #87's recorded finding, not this PR's: at the default shift, every one of case9241pegase's
  warm AC predictor-corrector outages fails (`opf-outage-warm-start-2026-11-19.md`). The
  earlier batch records time the same failing solves.

## Against the bound

The plan bounded this PR at ×1.1–1.15, from the solves' 9–15% share and their trees' ×9–20
bound on a sweep.
- **On case_ACTIVSg25k's AC it reached ×1.11–1.12.**
- **Elsewhere it gives nothing,** since the smaller trees' sweeps are below the cutoff. Below
  it, the tasks' overhead cost more than the sweeps saved, as the calibration found.

Target 3's PR 3c, opt-in fused arithmetic, stays planned but unstarted (decided 2026-11-21).
