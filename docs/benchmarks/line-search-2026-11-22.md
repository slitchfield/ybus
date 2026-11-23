# The predictor-corrector method's line search, 2026-11-22

This is goal 5's survey ([`LineSearchChecklist.md`](../LineSearchChecklist.md), PR 4). The
predictor-corrector method has three line search settings (`MehrotraOptions::line_search`):
- **`Off`:** the default. Every step is taken whole.
- **`MonotoneFilter`** (#101, #104): Wächter and Biegler's monotone barrier, with their filter
  line search and feasibility restoration.
- **`AdaptiveFilter`** (#105): Nocedal, Wächter and Waltz's adaptive barrier. Free mode is
  Mehrotra's steps; fixed mode is the monotone filter.

This survey measures both filters on every MATPOWER case with costs, on MIPS's test problems,
and on warm-started outage batches. It applies the rule written before it ran (the checklist,
commit d18f45b) to decide the OPF presets' default.

**The verdict: neither filter qualifies, so both presets keep `Off`.**
- **DC:** `adaptive` met every condition but the iteration cost. It took 1.106× the
  iterations, against the 1.10× limit.
- **AC:** `adaptive` lost O3's bar on one case, case2737sop. Its LAM_Q moved from 9.7e-4 to
  1.05e-3, just over the 1e-3 price bar.
- **`monotone`** failed several conditions on both models.

**The other finding:** the filters make warm-started outage batches far more reliable. With
`Off`, case9241pegase's AC outages converge on none of 28 warm starts. With either filter, all
28 do. The rule measures only cold starts, so this doesn't change a default here; it's
recorded for a rule of its own (below).

## Setup

| | |
|---|---|
| Code | The `dev/line-search-survey` branch at `d18f45b`: `main` after #105, with the survey's tooling (313a7e3) and the rule (d18f45b) |
| Machine and system | As in the earlier records: Intel Core i7-14700K, 15 GB of memory available to WSL2, Ubuntu 24.04.3 under WSL2, Rust 1.99.1, release profile. All runs were on the same evening, 18:17 to 19:26 (UTC−5), one after another, with nothing else running |
| OPF runs | `scripts/benchmark-opf.sh` for DC and AC, all 73 cases with costs, at the default limits (3600 s and 10 GiB per case). Five settings per model, in this order: MIPS-LU; then the predictor-corrector method on LDLᵀ, `ldlt-uncorrected` with `Off` (today's preset), `ldlt` with `Off`, `ldlt` with `monotone`, and `ldlt` with `adaptive`. The results are in `line-search-2026-11-22-<model>-<solver>-<kkt>-<line search>.csv` beside this file |
| Batches | `scripts/benchmark-opf-batch.sh` at its default cases (case2383wp, case3120sp, case9241pegase): the predictor-corrector method, DC and AC, 28 outages each (14 branches, 14 generators), cold and warm, serially (`YBUS_BATCH_THREADS=1`), once per line search. The results are in `line-search-2026-11-22-batch-<line search>.csv` |
| MIPS's test problems | `every_line_search_meets_the_bar_on_mipss_problems` in `ybus_ipm/tests/mehrotra_test.rs` |

**On the filters' `ldlt` runs:** either filter turns the inertia correction on, so they run as
`ldlt`. The `ldlt`-with-`Off` run separates the correction's effect from the line search's.
For `adaptive`, it is also what free mode is, bit for bit.

## The default path is unchanged

- **Today's preset** (`ldlt-uncorrected`, `Off`) matches the last survey of it
  (`ldlt-amalgamation-2026-11-08-*-mehrotra-ldlt-uncorrected.csv`, #62) in every outcome and
  iteration count, on all 146 runs.
- **MIPS-LU** matches its own runs in #59's survey (`ldlt-opf-2026-11-08-*-mips-lu.csv`) in
  every outcome and iteration count.
- **MIPS-LU converges** on 59 DC and 46 AC cases, and **the preset** on 64 DC and 49 AC.

**The tooling's `--line-search off`** matched the unchanged example on every result column, in
36 runs (313a7e3).

## Outcomes

| | DC converged | DC: O3's bar | AC converged | AC: O3's bar |
|---|--:|--:|--:|--:|
| MIPS-LU | 59 | | 46 | |
| Today's preset (`ldlt-uncorrected`, `Off`) | 64 | 58 | 49 | 30 |
| `ldlt`, `Off` (free mode's equivalent) | 64 | 58 | 48 | 30 |
| `ldlt`, `monotone` | 63 | 58 | 50 | 29 |
| `ldlt`, `adaptive` | 64 | 58 | 49 | 29 |

**O3's bar,** where MATPOWER converged:
- the objective within 1e-6 of MATPOWER's
- the prices within 1e-3
- each KKT residual within 100× MIPS-LU's, or below 1e-6

**Every outcome that differs from today's preset:**

| Model | Setting | Case | Today's preset | This setting |
|---|---|---|---|---|
| DC | `monotone` | case_ACTIVSg70k | converged | max iterations |
| DC | `monotone` | case_SyntheticUSA | numerical failure | max iterations |
| AC | `ldlt`, `Off` | case8387pegase | converged | max iterations |
| AC | `ldlt`, `Off` | case_SyntheticUSA | max iterations | numerical failure |
| AC | `monotone` | case13659pegase | max iterations | **converged** |
| AC | `monotone` | case145 | numerical failure | max iterations |
| AC | `monotone` | case_ACTIVSg10k | numerical failure | max iterations |
| AC | `adaptive` | case_ACTIVSg10k | numerical failure | max iterations |
| AC | `adaptive` | case_SyntheticUSA | max iterations | numerical failure |

- **`adaptive` keeps every case today's preset converges on,** on both models. That includes
  AC's case8387pegase, which the inertia correction alone loses.
- **`monotone` gains case13659pegase,** the AC failure goal 5 was recorded for. It's the only
  case either filter gains, and `monotone` loses DC's case_ACTIVSg70k.

## The rule

| Condition | DC, `monotone` | DC, `adaptive` | AC, `monotone` | AC, `adaptive` |
|---|---|---|---|---|
| 1. Every converged case still converges | **no:** case_ACTIVSg70k | yes | yes | yes |
| 1. No converged objective further from MATPOWER's (beyond 1e-6) | yes | yes | **no:** case30pwl (5.1e-7 to 3.8e-5), case9241pegase (2.3e-6 to 2.8e-6) | **no:** case9241pegase (2.3e-6 to 2.8e-6) |
| 2. O3's bar kept | yes | yes | **no:** case30, case30pwl, case2383wp | **no:** case2737sop |
| 3. Total iterations (cases both converge on) | **no:** 2.59× (1203 / 464) | **no:** 1.106× (552 / 499) | **no:** 1.27× (1166 / 915) | yes: 1.060× (970 / 915) |
| 3. Median time ratio | **no:** 1.24 | yes: 1.02 | **no:** 1.12 | yes: 1.00 |
| **Qualifies** | no | no | no | no |

**`adaptive`'s two AC misses are slight, and partly not the line search's:**
- **case9241pegase:** the inertia correction alone (`ldlt` with `Off`) moves the objective the
  same way, to 2.78e-6, so the line search isn't the cause. Both values are already over O3's
  1e-6.
- **case2737sop:** one switch to fixed mode moves LAM_Q from 9.7e-4 to 1.05e-3 of MATPOWER's.
  The objective (6.4e-10) and the residuals stay within the bar.

**DC `adaptive`'s 1.106×:** 46 of its 64 converged cases never switch to fixed mode, and take
exactly today's iterations. The extra iterations all come from the cases that switch: 14 of the
18 take more iterations than today's preset. On AC, 24 of 49 converged cases switch, and only 6
of those take more.

## The line search's work

Over all 73 cases, each candidate's counts from the method's record (`opf_survey`'s columns):

| | Backtracks | Second-order corrections | Restorations | Their iterations | Switches to fixed | Returns | Fixed iterations |
|---|--:|--:|--:|--:|--:|--:|--:|
| DC, `monotone` | 14 | 7 | 20 | 85 | – | – | 1542 |
| DC, `adaptive` | 0 | 0 | 17 | 162 | 27 | 18 | 201 |
| AC, `monotone` | 252 | 15 | 68 | 1570 | – | – | 1743 |
| AC, `adaptive` | 333 | 9 | 79 | 1479 | 66 | 44 | 500 |

- **Most restorations are in cases that fail anyway.** They're why the filters' total solve time
  over all 73 cases is about twice today's: DC 44 s (`adaptive`) and 40 s (`monotone`) against
  16 s; AC 371 s and 463 s against 166 s. The failing large cases spend their iterations in
  fixed mode and restoration before failing.
- **On the cases both converge on,** `adaptive`'s median time ratio is 1.00 (AC) and 1.02 (DC).
- **Fewer or more iterations than today's preset,** on the cases both converge on:
  - **`adaptive`:** fewer on 7 and more on 6 of 49 (AC); fewer on 1 and more on 14 of 64 (DC).
  - **`monotone`:** more on 47 of 49 (AC) and on 60 of 63 (DC).

## The objective-gap cases

The checklist's "Why" listed four AC cases where the method converges to a different objective
than MATPOWER's while MIPS-LU matches it. **This survey corrects that for three of them:**
MATPOWER itself fails on case60nordic, case_ACTIVSg25k and case_ACTIVSg70k (AC), and so does
MIPS-LU here. Their "gap" is to MATPOWER's failed final iterate, not to a solution.

| Case (AC) | MATPOWER | MIPS-LU | Today's preset | `monotone` | `adaptive` |
|---|---|---|---|---|---|
| case33mg | converged | converged, 5.3e-14 | 2.0e-5 | 2.0e-5 | 2.0e-5 |
| case60nordic | failed | failed | 3.2e-1 | 3.2e-1 | 3.2e-1 |
| case_ACTIVSg25k | failed | failed | 1.1e-1 | 1.1e-1 | 1.1e-1 |
| case_ACTIVSg70k | failed | failed | 9.7e-2 | 9.7e-2 | 9.7e-2 |

- **Neither filter changes any of them.**
- **case33mg's gap is the one real difference:** 2.0e-5 of MATPOWER's objective, unchanged by
  either line search.
- **On DC,** case33mg and case60nordic agree with MATPOWER under every setting, and MATPOWER
  fails on case_ACTIVSg25k and case_ACTIVSg70k.

## MIPS's test problems

Every problem MIPS converges on meets `check`'s bar under every setting, on LDLᵀ. The infeasible
LP fails under each.

| Problem | MIPS | `Off` | `monotone` | `adaptive` |
|---|--:|--:|--:|--:|
| banana | 6 | 6 | 24 | 6 |
| quadratic3 | 2 | 2 | 2 | 2 |
| qp4 | 10 | 13 | 10 | 13 |
| lp3 | 11 | 10 | 10 | 8 |
| qps_quadratic3 | 2 | 2 | 2 | 2 |
| qp2 | 11 | 8 | 9 | 8 |
| qps_qp4 | 12 | 13 | 10 | 13 |
| nonlinear2 | 8 | 8 | 9 | 8 |
| nonlinear3 | 8 | 6 | 6 | 6 |
| hs71 | 8 | 9 | 9 | 9 |
| infeasible LP (fails) | | 5 | 12 | 14 |

Where `adaptive`'s count is `Off`'s, it stays in free mode, which is `Off`'s step. It differs
only on lp3, 8 iterations against 10. The tests of PRs 1 to 3 cover the starts where the filters matter: hs71 from 1 and 2,
and the Wächter–Biegler counterexample.

## Warm-started outage batches

Each line is 28 outages solved serially. "Converged" counts the outages that converged.

| Case | Model | Start | `Off` | `monotone` | `adaptive` |
|---|---|---|---|---|---|
| case2383wp | DC | cold | 27, 1.72 s | 26, 6.21 s | 27, 2.43 s |
| case2383wp | DC | warm | **22**, 0.45 s | 27, 1.15 s | 26, 1.04 s |
| case2383wp | AC | cold | 26, 9.54 s | 26, 16.4 s | 26, 15.4 s |
| case2383wp | AC | warm | 26, 6.87 s | 26, 13.0 s | 26, 9.28 s |
| case3120sp | DC | cold | 28, 1.35 s | 28, 6.58 s | 28, 4.13 s |
| case3120sp | DC | warm | **25**, 0.62 s | 28, 1.06 s | 28, 0.78 s |
| case3120sp | AC | cold | 28, 14.5 s | 28, 14.9 s | 28, 14.7 s |
| case3120sp | AC | warm | 28, 6.01 s | 28, 12.7 s | 28, 6.00 s |
| case9241pegase | DC | cold | 28, 1.73 s | 28, 4.91 s | 28, 1.77 s |
| case9241pegase | DC | warm | 28, 2.52 s | 28, 4.19 s | 28, 2.53 s |
| case9241pegase | AC | cold | 28, 42.2 s | 28, 70.7 s | 28, 49.9 s |
| case9241pegase | AC | warm | **0**, 42.0 s | 28, 75.3 s | 28, 51.0 s |

- **With `Off`, warm starts lose outages a cold start solves:**
  - case9241pegase's AC outages: none of 28
  - case2383wp's DC outages: 22 of 28
  - case3120sp's DC outages: 25 of 28
  
  O5's record has the same counts (`o5-batched-solves-2026-11-20.csv`), unremarked there.
- **`adaptive` converges on every warm outage its cold start converges on,** except one of
  case2383wp DC's 27.
- **Its warm times** are today's on case3120sp and case9241pegase DC. They're longer on
  case2383wp: 1.04 s against 0.45 s (DC), and 9.28 s against 6.87 s (AC).
- **This is the filters' clearest gain.** The rule written for this survey judged only cold
  starts on the 73 cases, so it changes no default. **A rule of its own, for the warm-start
  default, is left as a follow-up** (in `LongTermGoals.md`, goal 5's closing note).
