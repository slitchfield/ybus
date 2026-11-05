# N-1 outages without per-outage allocation, 2026-11-04

Phase O2's first CPU target, chosen from [`profile-2026-11-04.md`](profile-2026-11-04.md): every
single-branch outage of every MATPOWER case by rank-k updates of B′'s kept factorization, as in
[`dc-contingency-2026-11-02.md`](dc-contingency-2026-11-02.md) (serial) and
[`dc-contingency-parallel-2026-11-02.md`](dc-contingency-parallel-2026-11-02.md) (parallel), after
three changes that leave every result bit for bit the same. The survey's output is in
[`dc-contingency-allocation-2026-11-04.csv`](dc-contingency-allocation-2026-11-04.csv) and
[`dc-contingency-allocation-parallel-2026-11-04.csv`](dc-contingency-allocation-parallel-2026-11-04.csv).

## The changes

1. **Reachability without a search per outage.** Each outage used to rebuild the network's
   adjacency (about one allocation per bus) and search every bus breadth-first to check that no
   bus was stranded. A single branch can only strand buses if it's a bridge, so a bridge forest
   built once per model answers every single-branch outage, error included.
2. **Sparse solves into caller buffers** (`sparse_mat`'s `solve_sparse_into` with a reusable
   `SparseSolveWorkspace`), instead of fresh marks, stacks and dense vectors per solve.
3. **Per-thread outage workspaces** for the copies of the injections, angles and base solution
   and the columns of `A⁻¹C`, and `finish` computing each flow straight into the result and
   `Bbus θ` only at the reference buses.

Every result is the same operations in the same order:
- **A snapshot** of every outage's results, single and paired, serial and parallel, on case14
  through case_ACTIVSg10k, is bit for bit unchanged after each change.
- **This survey's outcome columns** (outages solved and islanding, mismatches against refactoring,
  worst deviation) are identical to the earlier records' on all 78 cases, serial and parallel.

## Setup

| | |
|---|---|
| Code | The `dev/n1-allocation` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K (28 hardware threads), Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Run | `scripts/benchmark-n1.sh`, then with `YBUS_BENCH_PARALLEL=1`, finishing at 22:20 and 22:24 EST, with the default limits |

## Results

| | Serial | Parallel |
|---|--:|--:|
| All outages of all 78 cases | 4.8 min, from 21.5 (4.5×) | 46 s, from 3.5 min (4.6×) |
| Median speedup, cases of 1,000 buses or more (27) | 5.2× | 2.8× |
| case2383wp | 0.40 → 0.080 s (5.2×) | 0.040 → 0.016 s (2.5×) |
| case9241pegase | 12.8 → 2.85 s (4.5×) | 0.95 → 0.25 s (3.8×) |
| case_ACTIVSg10k | 8.30 → 1.66 s (5.0×) | 0.62 → 0.16 s (3.8×) |
| case_ACTIVSg70k | 477 → 105 s (4.6×) | 81 → 17.2 s (4.7×) |
| case_SyntheticUSA | 680 → 157 s (4.3×) | 118 → 25.9 s (4.6×) |
| Peak memory, any case | 462 MiB (was 456) | 688 MiB (was 787) |

- **Most of the gain is the reachability check.** Timed at each change on case2383wp,
  case9241pegase and case_ACTIVSg10k (serial): the bridge forest alone gives 4.7×, 3.5× and 4.0×,
  and the workspaces 7%, 26% and 27% more. In parallel, 2.3–3.4× and 3–11% more.
- **The gain beats the profile's bound (1.5–2.2×).** The bound covered the allocation it could
  see; the breadth-first search's own work, O(buses + branches) per outage, was inlined into
  `DcStudy::outage` and counted as the outage update.
- **Radial networks gain most.** case1197, where every branch is a bridge and every outage
  strands buses, takes 0.11 ms serially, from 32.5 ms (290×). In parallel it takes 7.8 ms (it was
  7.1): with so little work, handing it to 28 threads costs far more than doing it, so the serial
  pass is 70 times faster there.
- **Parallel speedups narrowed** on mid-sized cases (to 5–11× over serial, from 11–13×), since
  the serial path improved more; on the largest cases, parallel is 6× serial, as before.

## Where the time goes now

Profiled again with `scripts/profile-pf.sh` (workloads `n1` and `n1-parallel`):

| Case | Mode | Outage updates | Triangular solves | Allocation and copying | Threads busy |
|---|---|--:|--:|--:|--:|
| case2383wp | serial | 58% | (with the factorization) | 12% | 1.0 |
| case9241pegase | serial | 65% | 28% | 6% | 1.0 |
| case9241pegase | parallel | 66% | 26% | 4% | 15.7 |
| case_ACTIVSg70k | serial | 62% | 35% | 2.5% | 1.0 |
| case_ACTIVSg70k | parallel | 61% | 38% | 0.7% | 25.3 |
| case_SyntheticUSA | serial | 63% | 34% | 2.4% | 1.0 |

Allocation fell from 35–55% to 1–12%. What's left is arithmetic every outage's full result
needs: the dense update of the angles, flows and outputs (the outage updates) and the dense back
substitution of each sparse solve, since `A⁻¹ c_k` is dense. Reducing either means computing less
per outage, such as only the flows a caller asks for, through the streaming API this PR deferred.
