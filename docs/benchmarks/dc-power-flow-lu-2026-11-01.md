# DC power flow with the Gilbert–Peierls LU, 2026-11-01

Phase 3's third measurement: the 19 MATPOWER 8.1 cases whose B′ isn't positive definite, so
the DC power flow solves them by LU. They're solved with the new Gilbert–Peierls LU (threshold
partial pivoting, 0.1), ordered by AMD as `DcModel::solve` does, and unordered to separate the
algorithm from the ordering. The baseline is the old LU in
[`dc-power-flow-2026-10-31.md`](dc-power-flow-2026-10-31.md); the 59 Cholesky cases are in
[`dc-power-flow-amd-2026-11-01.md`](dc-power-flow-amd-2026-11-01.md). The AMD survey's output
is in [`dc-power-flow-lu-2026-11-01.csv`](dc-power-flow-lu-2026-11-01.csv).

## Setup

| | |
|---|---|
| Code | The `dev/gilbert-peierls` branch, at the commit that adds this record |
| Machine and system | As before: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Run | `scripts/benchmark-dc.sh` over the 19 LU cases with `YBUS_BENCH_ORDERING=amd` (09:05 to 09:06 EST), then `natural` (09:06 to 09:07 EST) |

The baseline was recorded before Phase 2, so its times and memory include costs Phase 2
removed; its fill and backward errors don't depend on that.

## Results

| The 19 LU cases | Old LU (unordered) | New LU, unordered | New LU, AMD |
|---|--:|--:|--:|
| Solves, together | 478.5 s | 13.2 s | **0.33 s** |
| Slowest case | 83.3 s (case13659pegase) | 2.2 s (case13659pegase) | 0.12 s (case_SyntheticUSA) |
| Worst fill (factor ÷ B′) | 108.9× | about the same | **2.5×** (case_ACTIVSg70k) |
| Median ratio of our backward error to MATPOWER's | 2.32 | 2.57 | **1.23** |
| Largest ratio | 30.0 | 28.0 | **2.7** (case2848rte) |
| Cases with a ratio above 3 | 8 | 9 | **0** |
| Largest backward error | 6.5 ε | 6.0 ε | **0.49 ε** (case60nordic) |
| Peak memory | 928 MiB | | **128 MiB** |

- **The algorithm and the ordering each matter.**
  - Unordered, the new LU fills the same as the old one (a median of 1.00 times its factor
    entries, from 0.94 to 1.11, since its threshold pivoting picks different pivots), but runs
    26 times faster at the median (3.5 to 49 times). Its work follows the arithmetic, where
    the old LU scanned every remaining row for each pivot.
  - AMD then cuts the fill from up to 109 times B′'s entries to at most 2.5 times, and the
    time by another factor of 40. case13659pegase solves in 15 ms instead of 83 s, and
    case6515rte 6,500 times faster than before.
- **All 78 cases now solve in 0.35 s together** (0.33 s of LU and 19 ms of Cholesky), against
  496 s in the Phase 1.5 baseline.
- **Every case agrees with MATPOWER more closely.** On the LU cases, no difference exceeds
  2.5e-12 of each quantity's largest value (case13659pegase's angles), against 1.1e-10
  (case_SyntheticUSA) with the old LU. That was the worst case overall; with AMD-ordered
  Cholesky and LU, every one of the 78 cases is now within 2.5e-12.
- **The criterion benchmark** (case2383wp, a Cholesky case) is unchanged, as it should be.

## The fill-in hypothesis: the LU verdict

The unordered run is the clean test: the same new algorithm, at the old fill, gives the old
backward errors back (median ratio 2.57, largest 28), and at AMD's fill it removes them
(median 1.23, largest 2.7, none above 3). case_SyntheticUSA falls from 5.4 ε to 0.18 ε, equal
to MATPOWER's; case_ACTIVSg70k from 6.5 ε to 0.37 ε. So the excess followed the fill, not the
algorithm.

Together with the Cholesky verdict, the hypothesis is supported for both factorizations:
ordered to MATPOWER's level of fill, our backward error matches MATPOWER's, and every case's is
below 0.8 ε. Threshold pivoting at 0.1 hasn't cost stability here: the largest backward error
of an LU case is 0.49 ε.

## LU cases

Fill divides the factor's entries (L and U, their shared diagonal counted once) by B′'s. ε is
`f64::EPSILON`.

| Case | Buses | Solve, old LU (ms) | Solve, new LU unordered (ms) | Solve, new LU AMD (ms) | Fill, old LU | Fill, new LU AMD | Backward error, old LU (ε) | unordered (ε) | AMD (ε) | MATPOWER's (ε) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| case60nordic | 60 | 0.056 | 0.016 | 0.037 | 1.6 | 1.16 | 0.23 | 1.11 | 0.49 | 0.28 |
| case300 | 300 | 2.959 | 0.381 | 0.314 | 7.1 | 1.46 | 0.07 | 0.05 | 0.06 | 0.04 |
| case1888rte | 1888 | 315.7 | 19.5 | 1.658 | 23.1 | 1.36 | 1.05 | 0.54 | 0.47 | 0.73 |
| case1951rte | 1951 | 382.0 | 25.0 | 1.725 | 25.5 | 1.35 | 0.22 | 0.30 | 0.19 | 0.13 |
| case2848rte | 2848 | 1,900 | 81.0 | 2.646 | 42.2 | 1.41 | 0.23 | 0.31 | 0.20 | 0.07 |
| case2868rte | 2868 | 1,946 | 80.1 | 2.575 | 41.9 | 1.40 | 0.20 | 0.24 | 0.33 | 0.16 |
| case3012wp | 3012 | 1,912 | 73.5 | 2.980 | 53.2 | 1.87 | 0.56 | 0.43 | 0.35 | 0.28 |
| case3120sp | 3120 | 1,873 | 75.6 | 3.075 | 52.1 | 1.87 | 1.40 | 0.42 | 0.41 | 0.23 |
| case3375wp | 3374 | 1,981 | 78.3 | 3.330 | 49.5 | 1.87 | 0.87 | 0.42 | 0.19 | 0.19 |
| case6468rte | 6468 | 41,688 | 1,263 | 6.438 | 108.9 | 1.58 | 0.23 | 0.41 | 0.14 | 0.10 |
| case6470rte | 6470 | 39,225 | 1,226 | 6.445 | 104.1 | 1.58 | 0.23 | 0.41 | 0.12 | 0.12 |
| case6495rte | 6495 | 40,377 | 1,231 | 6.425 | 105.1 | 1.56 | 0.68 | 0.52 | 0.15 | 0.14 |
| case6515rte | 6515 | 41,386 | 1,355 | 6.364 | 106.3 | 1.56 | 0.75 | 0.46 | 0.18 | 0.20 |
| case9241pegase | 9241 | 68,696 | 1,893 | 11.9 | 98.6 | 1.73 | 0.52 | 0.69 | 0.19 | 0.18 |
| case_ACTIVSg10k | 10000 | 1,165 | 49.8 | 11.0 | 17.4 | 2.02 | 0.10 | 0.09 | 0.02 | 0.04 |
| case13659pegase | 13659 | 83,337 | 2,222 | 15.2 | 79.9 | 1.53 | 1.55 | 1.72 | 0.33 | 0.67 |
| case_ACTIVSg25k | 25000 | 24,387 | 699.3 | 29.7 | 41.2 | 2.22 | 3.97 | 2.60 | 0.24 | 0.16 |
| case_ACTIVSg70k | 70000 | 59,361 | 1,395 | 98.4 | 37.3 | 2.52 | 6.48 | 5.98 | 0.37 | 0.25 |
| case_SyntheticUSA | 82000 | 68,518 | 1,406 | 116.3 | 34.3 | 2.46 | 5.42 | 5.06 | 0.18 | 0.18 |
