# DC power flow on every MATPOWER case, 2026-10-31

The Phase 1.5 baseline: `powerflow`'s DC power flow on all 78 cases in MATPOWER 8.1, before
Phase 2's fixed-pattern matrices and Phase 3's ordered LU. The survey's output, with more
columns than the table below, is in [`dc-power-flow-2026-10-31.csv`](dc-power-flow-2026-10-31.csv).

## Setup

| | |
|---|---|
| Code | `powerflow` at `1d6ac9b`, whose library is unchanged since `7036a35`, with the survey tools as of this record's latest update |
| Machine | Intel Core i7-14700K (8 performance and 12 efficiency cores, 28 threads), 15 GiB of memory visible to WSL2 |
| System | Ubuntu 24.04.3 under WSL2 (Linux 6.6.87.2), Rust 1.97.1, release profile |
| Cases | MATPOWER 8.1's 78 cases, written as literals by `scripts/normalize-cases.sh` |
| Reference | MATPOWER's `rundcpf` under Octave 8.4.0, run locally by `scripts/generate-reference.sh`; its results aren't committed |
| Run | `scripts/benchmark-dc.sh` with its defaults, from 10:20 to 10:38 EDT: each case in its own process, limited to 1 hour and 10 GiB of address space |

Each step is timed by its median over repeated runs, for up to 1 s or 101 runs, so the larger
cases ran once, without warm-up. Expect a few percent of noise: case2383wp's solve took
334–360 ms across today's runs of the survey tool, and 317–330 ms in two runs of the criterion
benchmark, which warms up first. Peak memory is the whole process's, including a second
factorization after the timed runs that counts the factor's entries; a run without that step
peaked up to 50 MiB lower (case13659pegase: 384 MiB).

## Results

- **All 78 cases solve, and agree with MATPOWER.** Relative to the largest value of each
  quantity in its case, no difference exceeds 1.1e-10 (case_SyntheticUSA's angles). The
  largest absolute differences are 7.6e-8° (case13659pegase) and 1.7e-6 MW (case_SyntheticUSA,
  whose flows reach 16,148 MW). Every case of up to 300 buses is within the 1e-9 that the
  committed reference tests allow.
- **Every solve is backward stable.** Our angles solve B′θ = P′ to within 10.7 ε
  (case2383wp), where ε is `f64::EPSILON`, with a median of 0.31 ε. MATPOWER's angles, measured
  against the same system, are within 1.4 ε (case145), with a median of 0.23 ε.
- **59 cases solve by Cholesky and 19 by LU.** All 19 have in-service branches with a negative
  series susceptance, which can make B′ indefinite; case145 has 24 and stays positive definite.
- **The solve dominates.** Parsing and building the model take at most 163 ms and 35 ms
  (case_SyntheticUSA). The solve takes up to 83 s (case13659pegase), and 496 s over all 78
  cases.
- **Run time follows fill-in, not size.** Without a fill-reducing ordering, a factor holds up to
  109 times the entries it would without fill-in (case6468rte). How the buses are numbered
  decides it: case_ACTIVSg10k's LU factor fills 17× and takes 1.2 s, while case9241pegase, of
  similar size, fills 99× and takes 69 s.
- **Memory stays modest.** The largest peak is 928 MiB (case_SyntheticUSA).

## For Phase 3

- **Order before factoring.** The PEGASE and RTE cases above 6,000 buses fill 70–109×, and even
  the large synthetic cases (case_ACTIVSg25k, case_ACTIVSg70k, case_SyntheticUSA) fill 34–41×.
  A fill-reducing ordering attacks the main cost.
- **Make the LU fallback cheaper.** At similar size and fill, LU costs about four times as much
  as Cholesky: case2848rte (LU, 42×) solves in 1.9 s and case2869pegase (Cholesky, 35×) in
  0.51 s. LU's pivot search also visits every remaining row at every step, however sparse the
  matrix: 3.4 billion row lookups for case_SyntheticUSA. Phase 3's Gilbert–Peierls LU looks only
  at a column's nonzeros.
- **The failed Cholesky attempt costs nothing.** In every LU case it stops at an early pivot
  within 13 ms (measured once, with a scratch program), so an LU case's time is the LU's.

## Fill-in and backward error: a hypothesis to track

Our backward error exceeds MATPOWER's on the cases whose factors fill in most, and MATPOWER's
solve reorders B′ to limit fill-in. The guess is that fill-in drives it, since every fill entry
adds rounding to the entries it updates. This record is the baseline for testing that once
Phase 3's orderings cut the fill, which a Phase 3 item in `docs/PowerFlowChecklist.md` asks for.

- **Ours against MATPOWER's:** the median ratio is 1.5, and 15 cases exceed 3×, up to 30×
  (case_SyntheticUSA). The largest ratios are on high-fill cases: case_SyntheticUSA,
  case_ACTIVSg70k and case_ACTIVSg25k (LU, 116–140 factor entries per row) and case2383wp
  (Cholesky, 88).
- **By fill:** over the 42 cases with fewer than 10 factor entries per row, the medians are
  0.25 ε (ours) and 0.23 ε (MATPOWER's). Over the 25 with 50 or more, they're 0.75 ε and 0.19 ε.
- **Rank correlations** (Spearman, over all 78 cases): factor entries per row against our
  backward error, 0.46; against MATPOWER's, −0.03; against the ratio of the two, 0.56. Size
  against ours, 0.41; against MATPOWER's, −0.06.

So our backward error rises with fill while MATPOWER's, for the same systems, doesn't, which
points at our factorization rather than at the cases' size or conditioning. Without an
ordering, though, size and fill rise together, so this record can't separate them. Phase 3 can:
an ordering changes the fill of the same systems, most cleanly for Cholesky, whose algorithm
stays the same. If the guess is right, our backward error will fall toward MATPOWER's on the
high-fill cases, and its correlations with fill and size will weaken. If it doesn't, look at
the summation order in the triangular solves, or at pivot growth in LU.

## Columns

- **Solve** is `DcModel::solve`: building B′, factoring it, the triangular solves, the flows and
  the slack. For an LU case it includes the failed Cholesky attempt.
- **Factor entries** are those of L for Cholesky, or of L and U, with their diagonal counted
  once, for LU. **Fill** divides them by the entries the factor would hold without fill-in: B′'s
  lower triangle for Cholesky, or all of B′ for LU.
- **Backward error** is the normwise backward error of our angles,
  ‖P′ − B′θ‖∞ / (‖B′‖∞ ‖θ‖∞ + ‖P′‖∞): the smallest relative change to B′ and P′ that makes
  them an exact solution (`powerflow/tests/common/accuracy.rs`). **MATPOWER's** is the same
  measure of MATPOWER's angles, for the same system.
- **Angle**, **Flow** and **Output** are the largest absolute differences from MATPOWER's bus
  angles (degrees), branch flows (MW) and generator outputs (MW).

## All cases

Sorted by bus count.

| Case | Buses | Branches | Solver | Parse (ms) | Model (ms) | Solve (ms) | Peak (MiB) | Factor entries | Fill | Backward error | MATPOWER's | Angle (°) | Flow (MW) | Output (MW) |
|---|--:|--:|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| case4_dist | 4 | 3 | Cholesky | 0.006 | 0.001 | 0.002 | 3 | 4 | 1.0 | 8.3e-17 | 2.8e-17 | 5.6e-17 | 1.1e-16 | 0 |
| case4gs | 4 | 4 | Cholesky | 0.006 | 0.001 | 0.002 | 3 | 6 | 1.2 | 4.1e-17 | 8.2e-17 | 4.4e-16 | 2.8e-14 | 2.8e-14 |
| case5 | 5 | 6 | Cholesky | 0.010 | 0.001 | 0.003 | 3 | 8 | 1.1 | 7.8e-17 | 2.6e-17 | 8.9e-16 | 2.8e-14 | 6.7e-14 |
| case6ww | 6 | 11 | Cholesky | 0.011 | 0.002 | 0.004 | 3 | 15 | 1.2 | 3.9e-17 | 5.9e-17 | 4.4e-16 | 7.1e-15 | 0 |
| case9 | 9 | 9 | Cholesky | 0.011 | 0.002 | 0.005 | 3 | 19 | 1.2 | 3.2e-17 | 4.0e-17 | 7.1e-15 | 2.8e-14 | 2.8e-14 |
| case9Q | 9 | 9 | Cholesky | 0.012 | 0.002 | 0.004 | 3 | 19 | 1.2 | 3.2e-17 | 4.0e-17 | 7.1e-15 | 2.8e-14 | 2.8e-14 |
| case9target | 9 | 9 | Cholesky | 0.012 | 0.002 | 0.005 | 3 | 19 | 1.2 | 9.1e-17 | 5.2e-17 | 1.1e-14 | 9.2e-14 | 1.1e-13 |
| case10ba | 10 | 9 | Cholesky | 0.011 | 0.002 | 0.005 | 3 | 17 | 1.0 | 4.0e-17 | 2.1e-17 | 7.1e-15 | 1.6e-14 | 1.2e-14 |
| case12da | 12 | 11 | Cholesky | 0.012 | 0.002 | 0.006 | 3 | 21 | 1.0 | 7.3e-17 | 6.5e-17 | 3.9e-15 | 1.3e-15 | 1.1e-15 |
| case14 | 14 | 20 | Cholesky | 0.021 | 0.003 | 0.010 | 3 | 50 | 1.6 | 4.9e-17 | 9.8e-17 | 1.4e-14 | 1.7e-13 | 2.8e-13 |
| case15da | 15 | 14 | Cholesky | 0.014 | 0.003 | 0.011 | 3 | 63 | 2.3 | 1.1e-16 | 4.3e-17 | 2.2e-15 | 1.8e-15 | 1.8e-15 |
| case15nbr | 15 | 14 | Cholesky | 0.014 | 0.002 | 0.011 | 3 | 63 | 2.3 | 1.5e-16 | 6.3e-17 | 3.6e-15 | 3.6e-15 | 3.6e-15 |
| case16am | 15 | 14 | Cholesky | 0.014 | 0.003 | 0.009 | 3 | 46 | 1.7 | 2.4e-23 | 2.5e-23 | 8.9e-16 | 1.3e-14 | 3.6e-15 |
| case16ci | 16 | 16 | Cholesky | 0.016 | 0.003 | 0.008 | 3 | 29 | 1.3 | 3.7e-17 | 7.8e-17 | 1.6e-15 | 2.1e-14 | 2.1e-14 |
| case17me | 17 | 16 | Cholesky | 0.016 | 0.003 | 0.009 | 3 | 42 | 1.4 | 7.8e-17 | 1.3e-16 | 1.3e-14 | 5.3e-14 | 5.3e-15 |
| case18 | 18 | 17 | Cholesky | 0.017 | 0.003 | 0.012 | 3 | 64 | 1.9 | 9.2e-17 | 1.2e-16 | 7.3e-14 | 1.2e-13 | 1.2e-13 |
| case18nbr | 18 | 17 | Cholesky | 0.017 | 0.003 | 0.013 | 3 | 79 | 2.4 | 6.3e-17 | 1.0e-16 | 1.8e-15 | 2.8e-15 | 4.4e-16 |
| case22 | 22 | 21 | Cholesky | 0.021 | 0.004 | 0.012 | 3 | 52 | 1.3 | 5.4e-17 | 2.6e-17 | 1.8e-14 | 2.0e-14 | 1.8e-14 |
| case24_ieee_rts | 24 | 38 | Cholesky | 0.050 | 0.006 | 0.017 | 3 | 93 | 1.7 | 5.9e-17 | 5.9e-17 | 2.8e-14 | 7.4e-13 | 7.0e-13 |
| case28da | 28 | 27 | Cholesky | 0.023 | 0.005 | 0.019 | 3 | 111 | 2.1 | 7.1e-17 | 7.1e-17 | 9.1e-14 | 2.9e-14 | 2.8e-14 |
| case30 | 30 | 41 | Cholesky | 0.033 | 0.006 | 0.024 | 3 | 140 | 2.1 | 5.5e-17 | 2.1e-17 | 1.2e-14 | 5.5e-14 | 5.0e-14 |
| case30Q | 30 | 41 | Cholesky | 0.033 | 0.006 | 0.025 | 3 | 140 | 2.1 | 5.5e-17 | 2.1e-17 | 1.2e-14 | 5.5e-14 | 5.0e-14 |
| case30pwl | 30 | 41 | Cholesky | 0.035 | 0.007 | 0.025 | 3 | 140 | 2.1 | 5.5e-17 | 2.1e-17 | 1.2e-14 | 5.5e-14 | 5.0e-14 |
| case_ieee30 | 30 | 41 | Cholesky | 0.039 | 0.006 | 0.024 | 3 | 133 | 2.0 | 7.9e-17 | 5.2e-17 | 3.2e-14 | 1.8e-13 | 2.3e-13 |
| case33bw | 33 | 37 | Cholesky | 0.031 | 0.006 | 0.021 | 3 | 117 | 1.9 | 6.5e-17 | 5.2e-17 | 1.3e-14 | 1.4e-14 | 7.5e-15 |
| case33mg | 33 | 37 | Cholesky | 0.032 | 0.006 | 0.021 | 3 | 117 | 1.9 | 4.7e-17 | 1.6e-16 | 1.9e-14 | 2.8e-14 | 1.2e-14 |
| case34sa | 34 | 33 | Cholesky | 0.030 | 0.006 | 0.022 | 3 | 124 | 1.9 | 5.5e-17 | 6.1e-17 | 1.1e-14 | 7.3e-14 | 7.3e-14 |
| case38si | 38 | 37 | Cholesky | 0.033 | 0.006 | 0.034 | 3 | 230 | 3.2 | 3.9e-17 | 5.2e-17 | 6.3e-14 | 5.7e-14 | 5.6e-14 |
| case39 | 39 | 46 | Cholesky | 0.050 | 0.007 | 0.029 | 3 | 157 | 1.9 | 5.8e-17 | 8.1e-17 | 8.3e-14 | 3.3e-12 | 3.3e-12 |
| case51ga | 51 | 50 | Cholesky | 0.044 | 0.008 | 0.035 | 3 | 221 | 2.2 | 6.5e-17 | 6.3e-17 | 2.4e-14 | 1.3e-14 | 1.2e-14 |
| case51he | 51 | 50 | Cholesky | 0.044 | 0.008 | 0.064 | 3 | 449 | 4.5 | 7.7e-17 | 4.3e-17 | 8.9e-15 | 2.5e-14 | 3.1e-15 |
| case57 | 57 | 80 | Cholesky | 0.070 | 0.012 | 0.077 | 3 | 522 | 4.0 | 8.5e-17 | 6.9e-17 | 1.2e-13 | 4.9e-13 | 7.4e-13 |
| case59 | 59 | 138 | Cholesky | 0.090 | 0.019 | 0.040 | 3 | 207 | 1.6 | 2.0e-17 | 4.8e-17 | 5.2e-13 | 1.1e-11 | 5.7e-13 |
| case60nordic | 60 | 88 | LU | 0.081 | 0.013 | 0.056 | 3 | 327 | 1.6 | 5.2e-17 | 6.3e-17 | 4.2e-13 | 9.3e-12 | 9.3e-12 |
| case69 | 69 | 68 | Cholesky | 0.058 | 0.010 | 0.058 | 3 | 427 | 3.2 | 2.3e-17 | 2.3e-17 | 2.1e-13 | 4.9e-13 | 4.9e-13 |
| case70da | 70 | 76 | Cholesky | 0.062 | 0.011 | 0.037 | 3 | 184 | 1.4 | 4.0e-17 | 4.8e-17 | 2.0e-14 | 8.2e-15 | 1.1e-14 |
| case_RTS_GMLC | 73 | 120 | Cholesky | 0.259 | 0.021 | 0.064 | 3 | 375 | 2.1 | 6.9e-17 | 5.6e-17 | 1.4e-13 | 1.1e-12 | 3.3e-12 |
| case74ds | 74 | 73 | Cholesky | 0.062 | 0.011 | 0.036 | 3 | 143 | 1.0 | 1.1e-16 | 1.4e-16 | 1.4e-14 | 8.5e-14 | 5.3e-15 |
| case85 | 85 | 84 | Cholesky | 0.071 | 0.013 | 0.084 | 3 | 603 | 3.6 | 1.3e-16 | 4.0e-17 | 3.0e-13 | 2.8e-13 | 2.8e-13 |
| case89pegase | 89 | 210 | Cholesky | 0.148 | 0.029 | 0.250 | 3 | 1,162 | 4.0 | 4.9e-17 | 4.0e-17 | 7.1e-13 | 4.9e-11 | 9.5e-11 |
| case94pi | 94 | 93 | Cholesky | 0.077 | 0.014 | 0.206 | 3 | 1,347 | 7.3 | 1.2e-16 | 9.3e-17 | 7.9e-13 | 4.5e-13 | 3.9e-13 |
| case118 | 118 | 186 | Cholesky | 0.182 | 0.027 | 0.218 | 3 | 1,276 | 4.4 | 8.4e-17 | 5.8e-17 | 2.3e-13 | 3.6e-12 | 2.4e-12 |
| case118zh | 118 | 132 | Cholesky | 0.104 | 0.018 | 0.067 | 3 | 377 | 1.6 | 1.6e-17 | 1.4e-17 | 9.3e-14 | 8.9e-14 | 1.7e-13 |
| case136ma | 136 | 156 | Cholesky | 0.127 | 0.020 | 0.072 | 3 | 336 | 1.3 | 4.0e-17 | 1.7e-17 | 6.0e-14 | 1.6e-13 | 2.2e-13 |
| case141 | 141 | 140 | Cholesky | 0.120 | 0.021 | 0.386 | 3 | 2,248 | 8.1 | 2.3e-16 | 1.0e-16 | 7.9e-11 | 2.9e-10 | 2.2e-10 |
| case145 | 145 | 453 | Cholesky | 0.311 | 0.059 | 0.504 | 3 | 2,166 | 3.9 | 2.1e-16 | 3.0e-16 | 1.8e-12 | 1.4e-10 | 1.8e-11 |
| case_ACTIVSg200 | 200 | 245 | Cholesky | 0.380 | 0.037 | 0.263 | 3 | 1,374 | 3.1 | 3.2e-17 | 9.6e-17 | 3.9e-13 | 4.0e-12 | 4.0e-12 |
| case300 | 300 | 411 | LU | 0.386 | 0.058 | 2.959 | 3 | 7,887 | 7.1 | 1.5e-17 | 9.6e-18 | 3.0e-13 | 7.1e-12 | 4.3e-12 |
| case_ACTIVSg500 | 500 | 597 | Cholesky | 0.911 | 0.097 | 1.344 | 3 | 5,805 | 5.4 | 1.1e-16 | 8.6e-17 | 3.3e-13 | 7.5e-12 | 1.1e-12 |
| case533mt_hi | 533 | 577 | Cholesky | 0.510 | 0.078 | 1.130 | 3 | 6,765 | 6.4 | 2.4e-17 | 9.0e-18 | 8.7e-14 | 6.4e-14 | 1.6e-14 |
| case533mt_lo | 533 | 577 | Cholesky | 0.518 | 0.080 | 1.115 | 3 | 6,765 | 6.4 | 2.0e-16 | 6.9e-17 | 4.2e-14 | 1.2e-13 | 1.1e-13 |
| case1197 | 1,197 | 1,196 | Cholesky | 1.035 | 0.179 | 3.298 | 5 | 16,876 | 7.1 | 3.7e-18 | 2.9e-18 | 2.0e-12 | 5.8e-13 | 5.8e-13 |
| case1354pegase | 1,354 | 1,991 | Cholesky | 1.906 | 0.330 | 101.3 | 13 | 74,834 | 24.5 | 3.9e-17 | 4.2e-17 | 7.8e-12 | 1.5e-10 | 6.8e-10 |
| case1888rte | 1,888 | 2,531 | LU | 2.627 | 0.425 | 315.7 | 17 | 149,940 | 23.1 | 2.3e-16 | 1.6e-16 | 3.7e-11 | 1.6e-09 | 3.0e-09 |
| case1951rte | 1,951 | 2,596 | LU | 2.454 | 0.417 | 382.0 | 18 | 170,090 | 25.5 | 4.9e-17 | 3.0e-17 | 5.1e-12 | 1.4e-10 | 6.6e-10 |
| case_ACTIVSg2000 | 2,000 | 3,206 | Cholesky | 4.381 | 0.505 | 33.8 | 15 | 59,868 | 12.8 | 6.6e-17 | 4.7e-17 | 1.7e-12 | 5.4e-11 | 5.4e-11 |
| case2383wp | 2,383 | 2,896 | Cholesky | 2.862 | 0.487 | 335.8 | 27 | 209,361 | 39.8 | 2.4e-15 | 1.4e-16 | 7.0e-10 | 9.4e-09 | 3.1e-08 |
| case2736sp | 2,736 | 3,504 | Cholesky | 5.170 | 0.538 | 279.2 | 25 | 211,223 | 35.2 | 2.2e-16 | 1.1e-16 | 3.6e-11 | 1.8e-09 | 2.8e-09 |
| case2737sop | 2,737 | 3,506 | Cholesky | 5.126 | 0.539 | 242.3 | 26 | 194,979 | 32.5 | 4.2e-16 | 1.1e-16 | 1.3e-11 | 4.6e-10 | 9.1e-10 |
| case2746wop | 2,746 | 3,514 | Cholesky | 5.027 | 0.554 | 377.2 | 43 | 237,376 | 39.3 | 3.3e-16 | 5.4e-17 | 5.7e-11 | 2.3e-09 | 4.5e-09 |
| case2746wp | 2,746 | 3,514 | Cholesky | 5.145 | 0.555 | 342.3 | 26 | 226,454 | 37.7 | 1.2e-16 | 2.5e-17 | 5.5e-11 | 1.0e-09 | 1.9e-09 |
| case2848rte | 2,848 | 3,776 | LU | 3.449 | 0.606 | 1,900 | 35 | 409,914 | 42.2 | 5.2e-17 | 1.6e-17 | 3.2e-12 | 1.1e-10 | 4.8e-11 |
| case2868rte | 2,868 | 3,808 | LU | 3.614 | 0.640 | 1,946 | 36 | 409,669 | 41.9 | 4.3e-17 | 3.5e-17 | 1.6e-11 | 1.4e-10 | 7.3e-10 |
| case2869pegase | 2,869 | 4,582 | Cholesky | 3.920 | 0.761 | 514.3 | 43 | 240,089 | 35.1 | 1.0e-16 | 6.5e-17 | 4.1e-11 | 6.4e-10 | 9.1e-10 |
| case3012wp | 3,012 | 3,572 | LU | 3.379 | 0.593 | 1,912 | 64 | 539,432 | 53.2 | 1.2e-16 | 6.3e-17 | 8.5e-12 | 4.5e-10 | 7.3e-11 |
| case3120sp | 3,120 | 3,693 | LU | 3.555 | 0.629 | 1,873 | 65 | 546,023 | 52.1 | 3.1e-16 | 5.2e-17 | 4.0e-11 | 1.3e-09 | 2.6e-09 |
| case3375wp | 3,374 | 4,161 | LU | 4.091 | 0.721 | 1,981 | 54 | 568,958 | 49.5 | 1.9e-16 | 4.1e-17 | 1.9e-11 | 6.4e-10 | 5.6e-10 |
| case6468rte | 6,468 | 9,000 | LU | 8.911 | 1.592 | 41,688 | 230 | 2,458,503 | 108.9 | 5.1e-17 | 2.3e-17 | 1.2e-11 | 1.8e-10 | 7.8e-10 |
| case6470rte | 6,470 | 9,005 | LU | 8.716 | 1.540 | 39,225 | 226 | 2,350,095 | 104.1 | 5.2e-17 | 2.7e-17 | 3.8e-11 | 3.6e-10 | 2.0e-09 |
| case6495rte | 6,495 | 9,019 | LU | 8.525 | 1.595 | 40,377 | 228 | 2,380,132 | 105.1 | 1.5e-16 | 3.1e-17 | 8.1e-11 | 7.8e-10 | 4.0e-09 |
| case6515rte | 6,515 | 9,037 | LU | 8.439 | 1.627 | 41,386 | 228 | 2,413,707 | 106.3 | 1.7e-16 | 4.4e-17 | 1.5e-10 | 1.5e-09 | 5.9e-09 |
| case8387pegase | 8,387 | 14,561 | Cholesky | 14.1 | 2.639 | 15,305 | 162 | 1,497,459 | 70.1 | 5.4e-16 | 2.1e-16 | 4.6e-10 | 8.7e-09 | 3.7e-08 |
| case9241pegase | 9,241 | 16,049 | LU | 14.3 | 2.954 | 68,696 | 421 | 3,710,049 | 98.6 | 1.2e-16 | 3.9e-17 | 5.2e-10 | 7.3e-09 | 3.1e-08 |
| case_ACTIVSg10k | 10,000 | 12,706 | LU | 19.9 | 2.303 | 1,165 | 73 | 600,482 | 17.4 | 2.2e-17 | 9.4e-18 | 1.7e-11 | 1.8e-09 | 1.8e-09 |
| case13659pegase | 13,659 | 20,467 | LU | 21.5 | 4.305 | 83,337 | 434 | 4,069,190 | 79.9 | 3.4e-16 | 1.5e-16 | 7.6e-08 | 8.5e-07 | 8.5e-07 |
| case_ACTIVSg25k | 25,000 | 32,230 | LU | 50.5 | 9.657 | 24,387 | 289 | 3,509,551 | 41.2 | 8.8e-16 | 3.5e-17 | 6.0e-10 | 2.2e-07 | 2.4e-07 |
| case_ACTIVSg70k | 70,000 | 88,207 | LU | 146.6 | 28.9 | 59,361 | 883 | 8,826,607 | 37.3 | 1.4e-15 | 5.5e-17 | 3.9e-08 | 1.5e-06 | 1.5e-06 |
| case_SyntheticUSA | 82,000 | 104,121 | LU | 163.2 | 35.4 | 68,518 | 928 | 9,548,931 | 34.3 | 1.2e-15 | 4.0e-17 | 4.4e-08 | 1.7e-06 | 1.7e-06 |
