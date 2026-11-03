# DC outages on every MATPOWER case, 2026-11-02

Phase 5's N-1 survey: every in-service branch of all 78 MATPOWER 8.1 cases taken out alone, by
rank-k (here rank-1) updates of B′'s kept factorization (`DcStudy::solve_outages`), against
refactoring: a new model of the case with that branch out of service, solved afresh. The
survey's output, one row per case, is in
[`dc-contingency-2026-11-02.csv`](dc-contingency-2026-11-02.csv).

## Setup

| | |
|---|---|
| Code | The `dev/woodbury` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Comparison | `powerflow`'s own refactored DC power flow (`dc_power_flow` on the outaged case), which agrees with MATPOWER's `rundcpf` on every case (earlier records); no MATPOWER results needed |
| Run | `scripts/benchmark-n1.sh`, from 21:02 to 21:28 EST, with the default limits (1 hour and 10 GiB per case) |

The whole pass over a case's outages is timed by its median over repeated runs, for up to 1 s or
101 runs; it runs in chunks of 256 outages, sharing a base solve per chunk. Every outage is
compared with refactoring on the 64 cases of up to 3,000 buses, and about 300 evenly spaced ones
on the 14 larger cases. Refactoring's time is the median over the compared outages, building the
model and solving it.

## Results

| | |
|---|--:|
| Cases | 78, all completed |
| Single outages | 380,629: 275,897 solved, 104,732 islanding (27.5%) |
| Compared with refactoring | 48,945 |
| Mismatches (one solves and the other doesn't, or different errors) | 0 |
| Worst relative difference | 4.9e-12 (case13659pegase) |
| Speedup over refactoring, median (range) | 10.1× (2.8–26.9×) |
| Speedup, the 27 cases of 1,000 buses or more, median | 18.9× |
| All outages of all cases, by rank-1 updates | 21.5 minutes |
| Peak memory, any case | 456 MiB |

- **The outages match refactoring.** Every compared outage either solves in both, within
  4.9e-12 relative in angles, flows and generator outputs, or islands in both with the same
  `IslandWithoutReference`, naming the same bus. That covers the 59 Cholesky cases and the 19
  whose B′ needs the LU fallback.
- **An islanding outage costs only the graph check.** 27.5% of single outages leave a bus without
  a reference: radial feeders' branches, and pendant buses in transmission cases.
- **The speedup grows with the network.** Per outage, a rank-1 update costs one sparse solve,
  whose back substitution is dense, plus building a full solution (flows and injections are
  O(n)); refactoring costs a new model, an analysis and a factorization. The distribution
  feeders case533mt_lo, case533mt_hi, case141 and case1197 gain least (2.8–3.0×): their B′
  fills in so little that refactoring is already cheap.

| Case | Buses | Outages | B′ | Per outage, rank-1 update | Refactoring | Speedup |
|---|--:|--:|---|--:|--:|--:|
| case2383wp | 2,383 | 2,896 | Cholesky | 137 µs | 1.93 ms | 14.2× |
| case9241pegase | 9,241 | 16,049 | LU | 795 µs | 16.6 ms | 20.8× |
| case13659pegase | 13,659 | 20,467 | LU | 998 µs | 20.2 ms | 20.2× |
| case_ACTIVSg25k | 25,000 | 32,229 | LU | 1.72 ms | 35.8 ms | 20.8× |
| case_ACTIVSg70k | 70,000 | 88,207 | LU | 5.41 ms | 124 ms | 22.9× |
| case_SyntheticUSA | 82,000 | 104,121 | LU | 6.53 ms | 147 ms | 22.5× |

On case_SyntheticUSA, all 104,121 single outages take 11.3 minutes, against about 4.3 hours by
refactoring at its median.

## What a record may hold

Case names and sizes, and the time, memory and accuracy of this crate's code on them, as in the
other records: no case data or power flow results.
