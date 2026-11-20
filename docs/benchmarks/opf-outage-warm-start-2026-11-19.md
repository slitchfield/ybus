# Warm starts across outages, 2026-11-19

This is the last stage of the warm-start plan in `OptimizationChecklist.md` ("Warm starts for
OPF sequences"). A kept model refuses a case whose topology differs, so a warm start can't
cross an outage through `solve_from` alone. Each model now keys its variables and constraints
by what they belong to, by case row:
- a bus's voltage and power balances
- a generator's outputs
- a branch's flow and angle limits
- a piecewise cost's helper and segments
- a dispatchable load's power factor

`map_start` carries a solution of one model onto another model of the same case's tables, and
`OpfStart::Mapped` and `MappedPrimal` start a solve from it.
- **What carries over:** a matched variable or multiplier takes the previous value, from the
  final iterate's multipliers.
- **What starts fresh:** a new variable takes its cold start's value, and a new multiplier 0.

The survey measured what this saves on single outages, solved as an N-1 study solves them:
each outage starts from the base case's solution. As decided for this PR, `warm_shift`'s
default isn't revisited. It stays at the provisional 1e-4 until a relative shift is tried.

**Outcome:** warm starts across outages save what they save within one network (PR C's
[record](opf-warm-start-2026-11-18.md)), with the same pattern by solver.
- **MIPS at 1e-6** takes a median ×2.24 less time on DC branch outages, ×2.05 on DC generator
  outages, and ×3.38 and ×2.99 on AC. It rescues 49 outages where the cold solve fails. It
  regresses on 4 outages of 1,106 and misses the bar on 6.
- **The predictor-corrector method** gains less, ×1.3–1.5 on DC and ×2.2–2.3 on AC at the small
  shifts, and fails where it failed in PR C:
  - on case9241pegase's AC at every shift from 1 to 1e-4
  - on case2383wp's and case3120sp's DC below 1e-2
- **A new model per outage** costs a median 2% of its cold solve on DC, and 1.6% on AC.

## The survey

- **The machine:** an Intel Core i7-14700K (28 threads) under WSL2, Linux 6.6. The release
  build ran one solve at a time.
- **The cases:** PR C's tier 1: case9, case30, case118, case300, case2383wp, case3120sp and
  case9241pegase, DC and AC, both solvers at MATPOWER's OPF settings with their own KKT defaults.
- **The outages** (`opf_sequence`):
  - `branch-outages`: the in-service branches that aren't bridges (`find_bridges`)
  - `generator-outages`: the in-service generators
  
  Each is a seeded sample of at most 50 per case, so 2,212 outage steps in all.
- **Each outage step:** a new model of the outage case (`model_ms`), solved cold, primal from
  the mapped base point, and warm at 1, 1e-1, 1e-2, 1e-4 and 1e-6 from the mapped point and
  multipliers. A predictor-corrector run also solves each step cold with MIPS-LU, the
  residuals' reference.
- **The verdicts** are PR C's corrected rule's:
  - `regression`: cold converged and warm didn't
  - `rescue`: warm converged where cold failed
  - `miss`: O3's bar without `LAM_Q`, against the cold step, with residuals against MIPS-LU's
  - `local optimum`: AC only
- **The run:** `YBUS_SEQUENCE_CASES="case9 case30 case118 case300 case2383wp case3120sp
  case9241pegase" YBUS_SEQUENCE_SEQUENCES="branch-outages generator-outages"
  YBUS_SEQUENCE_SHIFTS="1,1e-1,1e-2,1e-4,1e-6" scripts/benchmark-opf-sequence.sh`. Its rows
  are in [`opf-outage-warm-start-2026-11-19.csv.gz`](opf-outage-warm-start-2026-11-19.csv.gz).
- **Times are single runs.** A speedup is the median, over the steps where both solves
  converged, of the cold solve's time over the warm one's.
- **case2383wp's DC on MIPS** has no warm steps. Its base solve fails, as in PR C, so every
  outage there started cold.

## What warm starts save across outages

These are medians over the steps where both solves converged: iterations warm/cold, then the
speedup, then regressions, misses and rescues.

| Model, solver, outages | Primal | Warm 1 | Warm 1e-2 | Warm 1e-4 | Warm 1e-6 |
|---|---|---|---|---|---|
| DC, MIPS, branches | 10/11 ×1.07 (−33) | 9/11 ×1.15 (+17) | 7/11 ×1.42 (+17) | 5/11 ×1.74 (+17) | 3/11 ×2.24 (−1, +15) |
| DC, MIPS, generators | 10/11 ×1.06 (−14) | 9/11 ×1.11 (+35) | 7/11 ×1.30 (+35) | 5/11 ×1.58 (+35) | 3/11 ×2.05 (+34) |
| DC, predictor-corrector, branches | 3/4 ×1.12 (−48, ✗45) | 3/5 ×1.23 (✗1) | 2/4 ×1.42 (−5, ✗90) | 2/4 ×1.54 (−13, ✗115) | 2/4 ×1.55 (−22, ✗119) |
| DC, predictor-corrector, generators | 4/4 ×1.00 (−51) | 4/5 ×1.23 (✗4, +1) | 3/5 ×1.35 (✗1, +1) | 2/5 ×1.27 (−16, ✗5, +1) | 2/4 ×1.32 (−29, ✗2, +1) |
| AC, MIPS, branches | 30/34 ×1.11 (✗4) | 28/34 ×1.20 (✗1) | 15.5/26.5 ×1.69 (−10, ✗2) | 11/34 ×2.43 (−2, ✗5) | 6/34 ×3.38 (−1, ✗3) |
| AC, MIPS, generators | 31/35 ×1.14 (✗1) | 28/35 ×1.22 (✗2) | 22/35 ×1.73 (−14) | 11/35 ×2.56 (−4, ✗1) | 6/35 ×2.99 (−2, ✗3) |
| AC, predictor-corrector, branches | 16/32 ×1.50 (✗52, +7) | 12/23 ×1.10 (−50, ✗7, +7) | 10/23 ×1.33 (−50, ✗4, +7) | 7/23 ×2.13 (−50, ✗3, +7) | 5/32 ×2.30 (−2, ✗34, +7) |
| AC, predictor-corrector, generators | 20/32 ×1.46 (✗50, +6) | 34/32 ×1.08 (−50, ✗2, +6) | 26/32 ×1.29 (−50, ✗1, +6) | 9/32 ×2.04 (−50, ✗2, +6) | 7/32 ×2.22 (−2, ✗49, +6) |

- **MIPS gains most at the smallest shift,** and more on larger cases, as in PR C. At 1e-6 it
  is ×6.45 on case2383wp's AC branch outages, ×7.47 on case3120sp's, and ×21.9 on case3120sp's
  DC. On case3120sp's DC, MIPS's cold solves take about 148 iterations and fail on 17 branch
  and 36 generator outages; warm, it takes 7.
- **The predictor-corrector method's 50 AC regressions** at shifts 1 to 1e-4 are all
  case9241pegase's, every outage there. That's PR C's finding: a uniform absolute shift moves
  this case's start far from stationary. At 1e-6 the regressions fall to 4, on case2383wp and
  case3120sp, but case9241pegase then misses on 34 branch and 48 generator outages.
- **Its DC misses** are mostly complementarity on case118's and case300's branch outages, where
  warm solves stop at the tolerance in 1–2 iterations. On case118, the objectives agree within
  3e-11 and the complementarity is 3–7e-6, just over the bar's 1e-6 floor. Below 1e-2,
  case2383wp's and case3120sp's DC regress, as in PR C.
- **A primal-only start saves little** (×1.00–1.14), except on AC on the predictor-corrector
  method (×1.46–1.50). It regresses 33–51 times on DC.

## The rule's table, for the record

This applies PR C's corrected rule, which chooses no default here.

| Shift | MIPS: median iterations | Regressions | Misses | Rescues | Predictor-corrector: median iterations | Regressions | Misses | Rescues |
|---|---|---|---|---|---|---|---|---|
| 1 | 15 | 0 | 3 | 52 | 8.5 | 100 | 14 | 14 |
| 1e-1 | 12 | 0 | 9 | 52 | 8 | 100 | 112 | 14 |
| 1e-2 | 11 | 24 | 2 | 52 | 6 | 105 | 96 | 14 |
| 1e-4 | 7 | 6 | 6 | 52 | 4 | 129 | 125 | 14 |
| 1e-6 | 5 | 4 | 6 | 49 | 3 | 55 | 204 | 14 |

There are 1,106 warm steps per shift and solver. No shift is free of unacceptable steps on
either solver, which would have stopped the rule choosing a default, as in PR C.
- **MIPS's 24 regressions at 1e-2** are case9241pegase's AC: 10 on branch outages and 14 on
  generator outages, where its warm solves take a median 112 iterations against 44 cold.
- **MIPS's AC local optima** are counted, not barred: 61–75 at the larger shifts, falling to
  2 at 1e-6. All are case9241pegase's, as in PR C.

## The cost of a model per outage

An outage changes the topology, so each needs its own model. Building one costs a median 2.1%
of the outage's cold solve on DC (at most 11%, on the smallest cases), and 1.6% on AC (at most
31%, on case9's fraction-of-a-millisecond solves). On the large cases it's negligible: 27 ms
against 1.7–7.3 s on case9241pegase's AC, and 7 ms against 0.33–1.1 s on case2383wp's.

## The mapping, tested

In `ybus_powerflow/tests/opf_outage_test.rs`:
- **The identity:** a model's own solution maps to `WarmStart::from_solution`, bit for bit.
- **What an outage loses:** a branch outage loses its flow rows (AC: `Sf` and `St`; DC: one
  row of `A`). A generator outage loses its outputs and their bounds' multipliers, and on
  case30pwl its piecewise helper and segments. A round trip changes nothing else, bit for bit.
- **Values on their own rows:** a mapped start, reported by case row (solved with `max_it` 0),
  has the base solution's angles, outputs, prices and binding flow multipliers, row by row.
  case30pwl's objective is the base's less the outaged generator's cost. The test asserts it
  compared at least one flow multiplier, since a vacuous comparison would pass.
- **End to end:** on case118, a branch and a generator outage meet the bar in fewer iterations.
- **Planted bugs, each caught by the values test:**
  - AC's flow rows keyed by internal index
  - DC's flow rows keyed by internal index
  - the piecewise helpers keyed by internal index
  - the thresholded multipliers

  The round trip alone missed the first three: a key misnamed the same way in both directions
  maps back to itself.
