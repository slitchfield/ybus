# DC power flow with RCM-ordered Cholesky, 2026-10-31

The first Phase 3 measurement: `powerflow`'s DC power flow on all 78 MATPOWER 8.1 cases, with
B′'s Cholesky factorization ordered by reverse Cuthill–McKee. LU, which 19 cases need, isn't
ordered yet. The baseline is the Phase 1.5 record,
[`dc-power-flow-2026-10-31.md`](dc-power-flow-2026-10-31.md). The survey's output is in
[`dc-power-flow-rcm-2026-10-31.csv`](dc-power-flow-rcm-2026-10-31.csv).

## Setup

| | |
|---|---|
| Code | The `dev/orderings-rcm` branch, at the commit that adds this record |
| Machine and system | As in the baseline: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Run | `YBUS_BENCH_ORDERING=rcm scripts/benchmark-dc.sh`, from 15:13 to 15:30 EDT, with its default limits (1 hour, 10 GiB per case) |
| Check | `YBUS_BENCH_ORDERING=natural` over the 70 cases that solved in under 30 s in the baseline, from 15:11 to 15:13 EDT |

The baseline was recorded before Phase 2, so its times and memory include costs Phase 2
removed. Its fill and backward errors don't depend on that: the check below reproduces them
exactly.

## The check: the new factorization reproduces the old one

With `--ordering natural`, the new Cholesky (symbolic analysis, then numeric factorization
with a dense work vector) does the same arithmetic in the same order as the baseline's. On all
70 cases, every column but time and memory matches the baseline exactly: solver, B′'s
dimension and entries, factor entries, our backward error, MATPOWER's, and the differences
from MATPOWER's results. So no case has a factor entry that comes out exactly zero, which
would count in the new symbolic pattern but not in the baseline's count.

## Results

- **All 78 cases solve, and agree with MATPOWER.** Relative to the largest value of each
  quantity in its case, no difference exceeds 1.1e-10 (case_SyntheticUSA's angles, an LU case,
  unchanged from the baseline).
- **The 59 Cholesky cases solve in 19 ms together, against 17.5 s.**
  - case2383wp's solve takes 1.45 ms instead of 336 ms.
  - case8387pegase, the slowest Cholesky case, takes 5.3 ms instead of 15.3 s.
  - The criterion benchmark (case2383wp, warmed up) measures 1.56 ms against 292 ms on `main`
    this afternoon: 187 times faster.
- **Fill falls everywhere.** No Cholesky case fills more than in the baseline. The median case
  holds 1.86 times fewer factor entries, case8387pegase 15.5 times fewer. The worst fill,
  relative to B′'s lower triangle, falls from 70× (case8387pegase) to 11× (case_ACTIVSg2000).
- **Every solve stays backward stable.** The largest backward error is 6.5 ε
  (case_ACTIVSg70k, an LU case). Among Cholesky cases it's 2.35 ε (case8387pegase), against
  10.7 ε (case2383wp) in the baseline.
- **The LU cases are unchanged** apart from noise (median time ratio 0.98). They now account
  for all but 19 ms of the survey's 476 s of solving. Phase 3's Gilbert–Peierls LU, with AMD
  ordering, is next for them.
- **Peak memory** of a Cholesky case is at most 17 MiB (case8387pegase), against 162 MiB.

## The fill-in hypothesis, on the Cholesky cases

The baseline guessed that fill-in drives our backward error, since ours exceeded MATPOWER's
on the cases that fill most, and MATPOWER's solve is ordered. RCM changes the fill of the same
systems with the same algorithm, so it tests that cleanly for Cholesky.

| Cholesky cases (59) | Baseline (unordered) | RCM |
|---|--:|--:|
| Median backward error, ours | 0.29 ε | 0.26 ε |
| Median backward error, MATPOWER's | 0.25 ε | 0.25 ε |
| Median ratio, ours to MATPOWER's | 1.23 | 1.00 |
| Largest ratio | 16.5 (case2383wp) | 3.6 (case22) |
| Cases with a ratio above 3 | 7 | 2 |
| Rank correlation: factor entries per row against our backward error | 0.45 | 0.20 |
| Rank correlation: factor entries per row against the ratio | 0.39 | 0.07 |

- **The worst cases support it.** case2383wp's backward error falls from 10.7 ε to 0.72 ε as its
  fill falls 7.8 times, to within 11% of MATPOWER's 0.65 ε; case2746wop's falls from 1.5 ε to
  0.22 ε. The correlation with fill all but disappears.
- **It isn't the whole story.** case8387pegase's fill falls 15.5 times, but its backward error
  only from 2.44 ε to 2.35 ε, still 2.5 times MATPOWER's 0.94 ε. Over all 59 cases, the
  backward error fell in 35, rose in 21 and stayed the same in 3; the rank correlation between
  the fall in fill and the fall in backward error is 0.35.
- **Verdict so far:** fill explains the large excesses, but not all of the gap. AMD, which cuts
  fill further than RCM, is the next test; if case8387pegase's excess survives it, look at the
  summation order of the triangular solves.

## Cholesky cases

Fill divides the factor's entries by those of B′'s lower triangle. ε is `f64::EPSILON`.

| Case | Buses | Solve, baseline (ms) | Solve, RCM (ms) | Fill, baseline | Fill, RCM | Backward error, baseline (ε) | Backward error, RCM (ε) | MATPOWER's (ε) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| case4_dist | 4 | 0.002 | 0.001 | 1.0 | 1.0 | 0.37 | 0.12 | 0.12 |
| case4gs | 4 | 0.002 | 0.001 | 1.2 | 1.0 | 0.19 | 0.00 | 0.37 |
| case5 | 5 | 0.003 | 0.002 | 1.1 | 1.0 | 0.35 | 0.35 | 0.12 |
| case6ww | 6 | 0.004 | 0.003 | 1.2 | 1.0 | 0.18 | 0.35 | 0.26 |
| case9 | 9 | 0.005 | 0.003 | 1.2 | 1.2 | 0.14 | 0.14 | 0.18 |
| case9Q | 9 | 0.004 | 0.003 | 1.2 | 1.2 | 0.14 | 0.14 | 0.18 |
| case9target | 9 | 0.005 | 0.003 | 1.2 | 1.2 | 0.41 | 0.18 | 0.23 |
| case10ba | 10 | 0.005 | 0.003 | 1.0 | 1.0 | 0.18 | 0.22 | 0.10 |
| case12da | 12 | 0.006 | 0.004 | 1.0 | 1.0 | 0.33 | 0.35 | 0.29 |
| case14 | 14 | 0.010 | 0.004 | 1.6 | 1.4 | 0.22 | 0.47 | 0.44 |
| case15da | 15 | 0.011 | 0.004 | 2.3 | 1.0 | 0.48 | 0.22 | 0.19 |
| case15nbr | 15 | 0.011 | 0.004 | 2.3 | 1.0 | 0.66 | 0.22 | 0.28 |
| case16am | 15 | 0.009 | 0.004 | 1.7 | 1.0 | 0.00 | 0.00 | 0.00 |
| case16ci | 16 | 0.008 | 0.004 | 1.3 | 1.0 | 0.16 | 0.18 | 0.35 |
| case17me | 17 | 0.009 | 0.005 | 1.4 | 1.0 | 0.35 | 0.51 | 0.61 |
| case18 | 18 | 0.012 | 0.005 | 1.9 | 1.0 | 0.41 | 0.43 | 0.53 |
| case18nbr | 18 | 0.013 | 0.005 | 2.4 | 1.0 | 0.28 | 0.35 | 0.47 |
| case22 | 22 | 0.012 | 0.006 | 1.3 | 1.0 | 0.24 | 0.42 | 0.12 |
| case24_ieee_rts | 24 | 0.017 | 0.007 | 1.7 | 1.4 | 0.27 | 0.40 | 0.27 |
| case28da | 28 | 0.019 | 0.007 | 2.1 | 1.0 | 0.32 | 0.63 | 0.32 |
| case30 | 30 | 0.024 | 0.009 | 2.1 | 1.4 | 0.25 | 0.15 | 0.09 |
| case30Q | 30 | 0.025 | 0.009 | 2.1 | 1.4 | 0.25 | 0.15 | 0.09 |
| case30pwl | 30 | 0.025 | 0.008 | 2.1 | 1.4 | 0.25 | 0.15 | 0.09 |
| case_ieee30 | 30 | 0.024 | 0.008 | 2.0 | 1.4 | 0.35 | 0.31 | 0.23 |
| case33bw | 33 | 0.021 | 0.008 | 1.9 | 1.0 | 0.29 | 0.23 | 0.23 |
| case33mg | 33 | 0.021 | 0.008 | 1.9 | 1.0 | 0.21 | 0.34 | 0.70 |
| case34sa | 34 | 0.022 | 0.008 | 1.9 | 1.0 | 0.25 | 0.30 | 0.27 |
| case38si | 38 | 0.034 | 0.010 | 3.2 | 1.0 | 0.17 | 0.34 | 0.23 |
| case39 | 39 | 0.029 | 0.010 | 1.9 | 1.5 | 0.26 | 0.21 | 0.37 |
| case51ga | 51 | 0.035 | 0.012 | 2.2 | 1.0 | 0.29 | 0.34 | 0.29 |
| case51he | 51 | 0.064 | 0.011 | 4.5 | 1.0 | 0.35 | 0.23 | 0.19 |
| case57 | 57 | 0.077 | 0.014 | 4.0 | 2.0 | 0.38 | 0.35 | 0.31 |
| case59 | 59 | 0.040 | 0.015 | 1.6 | 1.1 | 0.09 | 0.17 | 0.22 |
| case69 | 69 | 0.058 | 0.015 | 3.2 | 1.0 | 0.10 | 0.10 | 0.10 |
| case70da | 70 | 0.037 | 0.016 | 1.4 | 1.0 | 0.18 | 0.17 | 0.22 |
| case_RTS_GMLC | 73 | 0.064 | 0.020 | 2.1 | 1.9 | 0.31 | 0.21 | 0.25 |
| case74ds | 74 | 0.036 | 0.015 | 1.0 | 1.0 | 0.50 | 0.63 | 0.65 |
| case85 | 85 | 0.084 | 0.017 | 3.6 | 1.0 | 0.60 | 0.30 | 0.18 |
| case89pegase | 89 | 0.250 | 0.034 | 4.0 | 2.4 | 0.22 | 0.12 | 0.18 |
| case94pi | 94 | 0.206 | 0.019 | 7.3 | 1.0 | 0.52 | 0.33 | 0.42 |
| case118 | 118 | 0.218 | 0.029 | 4.4 | 1.6 | 0.38 | 0.17 | 0.26 |
| case118zh | 118 | 0.067 | 0.024 | 1.6 | 1.0 | 0.07 | 0.05 | 0.06 |
| case136ma | 136 | 0.072 | 0.026 | 1.3 | 1.0 | 0.18 | 0.06 | 0.08 |
| case141 | 141 | 0.386 | 0.026 | 8.1 | 1.0 | 1.01 | 0.45 | 0.45 |
| case145 | 145 | 0.504 | 0.070 | 3.9 | 1.9 | 0.93 | 0.79 | 1.36 |
| case_ACTIVSg200 | 200 | 0.263 | 0.055 | 3.1 | 2.2 | 0.14 | 0.21 | 0.43 |
| case_ACTIVSg500 | 500 | 1.344 | 0.171 | 5.4 | 2.0 | 0.48 | 0.26 | 0.39 |
| case533mt_hi | 533 | 1.130 | 0.119 | 6.4 | 1.0 | 0.11 | 0.04 | 0.04 |
| case533mt_lo | 533 | 1.115 | 0.119 | 6.4 | 1.0 | 0.90 | 0.31 | 0.31 |
| case1197 | 1197 | 3.298 | 0.206 | 7.1 | 1.0 | 0.02 | 0.01 | 0.01 |
| case1354pegase | 1354 | 101.3 | 0.578 | 24.5 | 2.4 | 0.18 | 0.26 | 0.19 |
| case_ACTIVSg2000 | 2000 | 33.8 | 2.573 | 12.8 | 11.2 | 0.30 | 0.43 | 0.21 |
| case2383wp | 2383 | 335.8 | 1.454 | 39.8 | 5.1 | 10.69 | 0.72 | 0.65 |
| case2736sp | 2736 | 279.2 | 1.743 | 35.2 | 6.4 | 0.97 | 0.74 | 0.49 |
| case2737sop | 2737 | 242.3 | 1.619 | 32.5 | 6.0 | 1.89 | 0.61 | 0.52 |
| case2746wop | 2746 | 377.2 | 1.502 | 39.3 | 5.3 | 1.48 | 0.22 | 0.24 |
| case2746wp | 2746 | 342.3 | 1.694 | 37.7 | 6.0 | 0.53 | 0.29 | 0.11 |
| case2869pegase | 2869 | 514.3 | 1.452 | 35.1 | 2.8 | 0.46 | 0.16 | 0.29 |
| case8387pegase | 8387 | 15,305 | 5.277 | 70.1 | 4.5 | 2.44 | 2.35 | 0.94 |
