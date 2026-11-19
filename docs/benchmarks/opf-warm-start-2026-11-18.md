# Warm-started OPF sequences, 2026-11-18

This is the third stage of the warm-start plan in `OptimizationChecklist.md` ("Warm starts for
OPF sequences"). A kept model (`AcOpfModel`, `DcOpfModel`) now solves a case from an earlier
solution of the same model, with `solve_from` and an `OpfStart`:
- `Cold`: MATPOWER's starting point, bit for bit `solve`'s path
- `Primal`: the previous point, with cold multipliers
- `Warm`: the previous point and multipliers, shifted off the boundary by the solver's
  `warm_shift`

The survey below measured what warm starts save on load sequences, and applied the rule written
before it to choose `warm_shift`'s default.

**Outcome:**
- **No shift qualifies** under the rule, for either solver, so `warm_shift` stays at the
  provisional 1e-4.
- **Warm starts still pay where they work.** MIPS at 1e-6 never ends worse than a cold start,
  rescues 38 DC steps where the cold solve fails, and takes a median ×2.11 less time on DC and
  ×3.40 on AC. It misses the bar on 3 AC steps of 616.
- **The predictor-corrector method** has no shift without regressions or many misses. It fails
  most on case9241pegase's AC, where an absolute shift fits badly.

## Found and fixed on the way

The OPF sequences found three defects in the solvers' warm start (#85), each fixed and tested
before the survey's final run.
- **The multipliers' scale.** A solution reports its multipliers divided by `cost_mult`, and the
  warm start read them back unscaled. At the AC OPF's `cost_mult`, 1e-4, they started 10⁴ times
  too large. Now they're scaled back (`9b8521a`).
- **The predictor-corrector method's bound slacks.** The warm start raised them to the shift.
  But they are `x − xmin` and `xmax − x` by definition, and between two close bounds the push
  leaves less than the shift. The raised slack overstated the room, a step crossed the bound,
  and on case2383wp's DC every shift failed in 2–3 iterations. Now they stay the point's
  (`0fe3f19`).
- **The thresholded report.** Each solver zeroes an inactive inequality's multiplier below
  `mu_threshold`, as MATPOWER reports it, and the warm start read that report. On an OPF's many
  inactive flow limits the zeros move the start off the solution. Warm from its own solution at
  a 1e-12 shift, case9241pegase's AC started at gradcond 0.27 against the solution's 2.8e-7.
  `Solution::iterate_multipliers` now keeps the final iterate's multipliers for warm starts,
  and the reported multipliers are unchanged (`180345d`).

## The survey

- **The machine:** an Intel Core i7-14700K (28 threads) under WSL2, Linux 6.6. The release
  build ran one solve at a time.
- **The cases:** tier 1, case9, case30, case118, case300, case2383wp, case3120sp and
  case9241pegase. Each ran DC and AC, both solvers at MATPOWER's OPF settings with their own KKT
  defaults, and both sequences.
- **The sequences** (`opf_sequence`, deterministic):
  - `sweep`: every load ×0.90 to ×1.10, 21 steps
  - `walk`: 25 steps, each bus's load moved ±2% a step, uniformly, by a seeded SplitMix64
    stream
- **Each step is solved by independent chains:** cold, primal, and warm at 1, 1e-1, 1e-2, 1e-4
  and 1e-6. Each chain starts from its own previous solution, and cold after a failure. A
  predictor-corrector run also solves each step cold with MIPS-LU, the residuals' reference.
- **The run:** `YBUS_SEQUENCE_CASES="case9 case30 case118 case300 case2383wp case3120sp
  case9241pegase" YBUS_SEQUENCE_SHIFTS="1,1e-1,1e-2,1e-4,1e-6"
  scripts/benchmark-opf-sequence.sh`. Its 8,680 rows are in
  [`opf-warm-start-2026-11-18.csv.gz`](opf-warm-start-2026-11-18.csv.gz).
- **Times are single runs** (a chain can't repeat a step). A speedup is the median, over the
  steps where both solves converged, of the cold time over the warm time.
- **Tier 2 wasn't run.** It was to run the chosen shift on the larger cases, and no shift was
  chosen.

## The rule

The rule was written before the survey and corrected once, after the first run and before the
rerun, both in the checklist (PR C):
- **Unacceptable:** a regression (the cold step converged and the warm one didn't), or a miss
  where both converged. A miss means the objective beyond 1e-6 or `LAM_P` beyond 1e-3 of the
  cold step's, or a KKT residual over 100× MIPS-LU's cold solve and over 1e-6.
- **Counted, not barred:**
  - rescues: warm converged where cold failed
  - AC local optima: the objective's miss alone, with the residuals within the bar
  - `LAM_Q`'s difference, since two cold solvers already differ in it by up to 9.4e-3
- **The choice:** each solver's default is the shift with the fewest median warm iterations
  among those with no unacceptable step, ties to the larger.

The first run found no shift for reasons outside the warm start. It counted rescues as
unacceptable, and compared residuals with a cold solve that had over-converged. It also exposed
the thresholded report. The rule was corrected and the defect fixed, then the survey was rerun.
The results here are the rerun's.

## The rule's tables

MIPS:

| Shift | Median iterations | DC regressions | DC misses | AC regressions | AC misses | Rescues | AC local optima | Qualifies |
|---|---|---|---|---|---|---|---|---|
| 1 | 13 | 0 | 0 | 1 | 4 | 39 | 14 | no |
| 1e-1 | 11 | 0 | 0 | 0 | 6 | 39 | 14 | no |
| 1e-2 | 9 | 0 | 0 | 5 | 5 | 39 | 9 | no |
| 1e-4 | 6 | 0 | 0 | 3 | 3 | 39 | 11 | no |
| 1e-6 | 4 | 0 | 0 | 0 | 3 | 38 | 12 | no |

The predictor-corrector method:

| Shift | Median iterations | DC regressions | DC misses | AC regressions | AC misses | Rescues | AC local optima | Qualifies |
|---|---|---|---|---|---|---|---|---|
| 1 | 7 | 0 | 12 | 19 | 22 | 1 | 0 | no |
| 1e-1 | 6 | 0 | 17 | 19 | 22 | 1 | 1 | no |
| 1e-2 | 4 | 0 | 11 | 19 | 22 | 1 | 1 | no |
| 1e-4 | 3 | 10 | 10 | 19 | 20 | 1 | 0 | no |
| 1e-6 | 3 | 21 | 13 | 1 | 36 | 1 | 2 | no |

Each shift has 616 warm steps per solver. The median iterations are over all of them.

## What warm starts save

These are medians over the steps where both solves converged: iterations warm/cold, then the
speedup, then regressions and rescues.

| Model, solver | Primal | Warm 1 | Warm 1e-2 | Warm 1e-4 | Warm 1e-6 |
|---|---|---|---|---|---|
| DC, MIPS | 9/9 ×1.03 (−6, +3) | 9/9 ×1.10 (+39) | 6/9 ×1.36 (+39) | 4/9 ×1.67 (+39) | 2/9 ×2.11 (+38) |
| DC, predictor-corrector | 3/4 ×1.09 (−24) | 3/4 ×1.30 (+1) | 3/4 ×1.35 (+1) | 2/4 ×1.36 (−10, +1) | 2/4 ×1.35 (−21, +1) |
| AC, MIPS | 17/19 ×1.10 | 15/19 ×1.14 (−1) | 10/19 ×1.59 (−5) | 6/19 ×2.45 (−3) | 4/19 ×3.40 |
| AC, predictor-corrector | 8/19 ×1.53 | 11/17 ×1.11 (−19) | 7/17 ×1.42 (−19) | 6/17 ×2.08 (−19) | 5/19 ×2.22 (−1) |

- **MIPS gains most at the smallest shift,** and more on larger cases. At 1e-6 it is ×4.91 on
  case300's AC, ×4.41 on case2383wp's and ×7.16 on case3120sp's. On case2383wp's and
  case3120sp's DC, where MIPS's cold solves take about 145 iterations and often fail, it takes
  9 and 8.5, ×16.9 and ×18.4.
- **The predictor-corrector method's cold solves are already short,** 4 iterations on DC, so
  warm starts save less. On DC its best shifts (1 to 1e-2) take 3 iterations. Below 1e-2,
  case2383wp and case3120sp regress, as the probe that widened the grid found: their
  multipliers reach 1e5.
- **A primal-only start saves little,** ×1.03–1.10, apart from the predictor-corrector method
  on AC (×1.53). It regresses 24 times on the predictor-corrector method's DC.

## Why no shift qualifies

- **MIPS:** 3 AC misses at 1e-6, against a cold step that converged:
  - case30: one stationarity miss, and one `LAM_P` and stationarity miss
  - case2383wp: one `LAM_P` miss
  
  At the larger shifts, case9241pegase's AC also regresses or misses on 3–7 steps.
- **The predictor-corrector method:**
  - **case9241pegase's AC:** 19 regressions and 17 misses at every shift from 1 to 1e-4, and 35
    misses at 1e-6. With exact multipliers, a shift of 1e-4 already starts this case at
    gradcond 17, and 1e-2 at 1.7e3. Its scaled multipliers are about 0.4, and it has tens of
    thousands of constraints, so a uniform absolute shift moves the start far from stationary.
  - **case300's DC:** misses on complementarity at every shift (9–16 steps).
  - **case2383wp's and case3120sp's DC:** regress below 1e-2.
- **No absolute shift suits both ends.** One shift has to serve multipliers of about 0.4 on
  case9241pegase's AC and of 1e5 on case2383wp's DC. A shift relative to the problem, such as
  to the previous solution's complementarity (Gondzio and Grothey, 2008), is left for later.

## AC local optima

These are steps that converged with residuals within the bar, but an objective 1.0e-6 to 5.0e-6
from the cold step's:
- **MIPS on case9241pegase:** 9–14 steps per shift, and 14 primal steps, so the warm start
  isn't their cause. The objectives differ by 1.8e-6 to 5.0e-6.
- **The predictor-corrector method:**
  - case9: one step at 1e-1 and one at 1e-2, at 1.4e-6 and 1.0e-6
  - case9241pegase: 2 steps at 1e-6, and 2 primal steps

The full list comes from the CSV's `verdict` column (`local optimum`).
