# AC power flow with reactive limits, reusing one Jacobian across rounds, 2026-11-03

The reactive-limit survey of [`ac-reactive-limits-2026-11-01.md`](ac-reactive-limits-2026-11-01.md)
rerun after a limit loop's later rounds began solving against one superset Jacobian: an angle and
a magnitude unknown for every bus, laid out, ordered and analyzed once, with the unknowns a round
doesn't solve for pinned, and its LU refactored across rounds instead of an AMD and a fresh
factorization each round. The first round still solves as plain Newton does. The survey's output,
one row per case and mode, is in
[`ac-reactive-limits-reuse-2026-11-03.csv`](ac-reactive-limits-reuse-2026-11-03.csv).

## Setup

| | |
|---|---|
| Code | The `dev/qlim-setup-reuse` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Reference | MATPOWER 8.1's legacy `runpf` with `pf.enforce_q_lims` 1 and 2, as before |
| Run | `scripts/benchmark-ac.sh` with `YBUS_BENCH_QLIM` 1 (06:28–06:29 EST) and 2 (06:29–06:56 EST), then the default Newton survey (06:56–06:57), with the default limits |

## Results

| | All at once | One at a time |
|---|--:|--:|
| Converged, with MATPOWER's iteration count and held generators | all 74 | all 74 |
| Same outcome as the earlier record, every case | yes | yes |
| Worst relative difference from MATPOWER | 1.7e-10 (unchanged) | 1.7e-10 (unchanged) |
| Worst absolute difference, any quantity | 3.3e-9 (2.7e-9 before) | 8.4e-9 (9.1e-9 before) |
| Solve time, all converged cases | 4.12 s (6.54 s before) | 1,519 s (4,397 s before) |
| Speedup on the 36 cases that hold a generator, median (range) | 1.45× (0.79–1.80×) | 3.38× (0.79–5.77×) |
| Peak memory, any case | 416 MiB (370 MiB before) | 374 MiB (328 MiB before) |

- **The results still match MATPOWER.** Every case converges where MATPOWER does, in its number of
  iterations, holding the same generators, in both modes. The values change only by rounding
  from a loop's second round, since the factorization's ordering and pivots differ; the worst
  relative difference (case141, which holds nothing) is unchanged.
- **A run that holds nothing is unchanged,** bit for bit: its one round solves as plain Newton
  does. That's 38 of the 74 converged cases.
- **The default Newton survey is unchanged,** matching
  [`ac-power-flow-2026-11-01.csv`](ac-power-flow-2026-11-01.csv) in every column but time and
  memory.
- **One at a time gains most:** the loop now refactors one LU, and on every case the pivots
  stayed acceptable from round to round, so a run takes two or three full factorizations
  however many rounds it has, against one per round before.
- **The smallest cases are slightly slower** (case5, case14, case39, case118, case145: 0.79–0.93×),
  where a solve takes well under a millisecond and building the superset outweighs what it saves.
- **Memory rises:** the superset carries every bus's magnitude unknown, and the factor
  stores zeros for the pinned ones, adding 4–29% to the peak on the cases of 9,000 buses or more
  (the most, case_ACTIVSg25k one at a time, 100 to 129 MiB). The fill itself wasn't counted
  separately.

| Case (one at a time) | Before | After | Speedup | Full factorizations |
|---|--:|--:|--:|--:|
| case2746wop | 1.87 s | 0.56 s | 3.4× | 300 → 2 |
| case3375wp | 1.47 s | 0.26 s | 5.8× | 212 → 2 |
| case9241pegase | 5.94 s | 1.75 s | 3.4× | 197 → 2 |
| case_ACTIVSg10k | 36.6 s | 10.7 s | 3.4× | 1,245 → 2 |
| case_ACTIVSg25k | 185 s | 53.4 s | 3.5× | 2,171 → 3 |
| case_ACTIVSg70k | 28.3 min | 9.6 min | 2.9× | 4,969 → 2 |
| case_SyntheticUSA | 40.8 min | 14.5 min | 2.8× | 6,337 → 2 |

**What remains per round** is the Newton work itself: about 2.2 refills, refactors and solves,
and the mismatches and `pfsoln`'s results, each O(n) or O(nnz).

## What a record may hold

Case names and sizes, and the time, memory and accuracy of this crate's code on them: no case data
or power flow results.
