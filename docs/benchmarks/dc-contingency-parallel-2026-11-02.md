# DC outages in parallel on every MATPOWER case, 2026-11-02

The N-1 survey of [`dc-contingency-2026-11-02.md`](dc-contingency-2026-11-02.md) rerun with the
outages solved in parallel (`DcStudy::par_solve_outages`, feature `rayon`), on rayon's default
pool of 28 threads. The survey's output is in
[`dc-contingency-parallel-2026-11-02.csv`](dc-contingency-parallel-2026-11-02.csv).

## Setup

| | |
|---|---|
| Code | The `dev/refinement-parallel` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K (28 hardware threads), Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Run | `YBUS_BENCH_PARALLEL=1 scripts/benchmark-n1.sh`, from 21:53 to 22:00 EST, with the default limits |

As before, the outages run in chunks of 256 (now each chunk in parallel), and the pass is timed
by its median over repeated runs.

## Results

- **The same results as the serial survey.** Every case's outcomes (outages solved and
  islanding, the comparison with refactoring, no mismatch) and its worst relative difference
  (4.9e-12 at most) match the serial record exactly, as `par_solve_outages`' bit-for-bit
  equality with `solve_outages` implies.
- **All outages of all 78 cases take 3.5 minutes,** against 21.5 serially.
- **On cases of 1,000 buses or more, the median speedup is 11.0×** (best 13.5×, case6468rte).
  Smaller cases are slower in parallel than serially: their outages take microseconds, less
  than handing them to threads.
- **The largest cases scale less,** 5.9× on case_ACTIVSg70k and 5.8× on case_SyntheticUSA,
  probably because of the survey's chunking rather than the method: each chunk of 256 repeats
  the serial base solve, and case_SyntheticUSA's 407 chunks would spend about 50 s of the
  parallel 118 s in it (estimated from its DC solve time, 123 ms, in the condition record). A
  caller solving all outages in one call shares one base solve.
- **Memory rises with the threads' working sets:** 787 MiB at most (case_SyntheticUSA), against
  456 MiB serially.

| Case | Buses | Outages | Serial | Parallel | Speedup | Against refactoring |
|---|--:|--:|--:|--:|--:|--:|
| case2383wp | 2,383 | 2,896 | 0.40 s | 0.040 s | 9.9× | 136× |
| case9241pegase | 9,241 | 16,049 | 12.8 s | 0.95 s | 13.5× | 282× |
| case13659pegase | 13,659 | 20,467 | 20.4 s | 1.71 s | 11.9× | 255× |
| case_ACTIVSg25k | 25,000 | 32,229 | 55.4 s | 5.95 s | 9.3× | 205× |
| case_ACTIVSg70k | 70,000 | 88,207 | 477 s | 81.2 s | 5.9× | 135× |
| case_SyntheticUSA | 82,000 | 104,121 | 679 s | 118 s | 5.8× | 144× |

"Against refactoring" compares the parallel pass's time per outage with refactoring one outage
serially (the survey's `speedup` column).

## What a record may hold

Case names and sizes, and the time, memory and accuracy of this crate's code on them: no case
data or power flow results.
