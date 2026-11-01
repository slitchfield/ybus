# DC power flow with AMD-ordered Cholesky, 2026-11-01

Phase 3's second measurement: the 59 MATPOWER 8.1 cases that `powerflow`'s DC power flow solves
by Cholesky, ordered by approximate minimum degree (AMD), against reverse Cuthill–McKee (RCM) in
the same session. The 19 LU cases don't use this ordering, so they're unchanged from
[`dc-power-flow-rcm-2026-10-31.md`](dc-power-flow-rcm-2026-10-31.md). The AMD survey's output is
in [`dc-power-flow-amd-2026-11-01.csv`](dc-power-flow-amd-2026-11-01.csv).

## Setup

| | |
|---|---|
| Code | The `dev/amd` branch, at the commit that adds this record |
| Machine and system | As before: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Run | `scripts/benchmark-dc.sh` over the 59 Cholesky cases, with `YBUS_BENCH_ORDERING=rcm` then `amd`, from 08:37 to 08:38 EST |

Today's RCM run matches the RCM record in every column but time and memory.

## Results

- **AMD cuts fill further, and never adds any.** No case fills more than with RCM. The large
  cases hold far fewer factor entries: case_ACTIVSg2000 4.7 times fewer, case8387pegase 2.9,
  case2383wp 3.2. The worst fill, relative to B′'s lower triangle, falls from 11.2×
  (case_ACTIVSg2000) to 2.4×. Small cases fill the same under both.
- **For one solve, the time is a wash.** The 59 solves take 19.0 ms with AMD and 18.6 ms with
  RCM; the median case is 1.33 times slower with AMD, the large ones between 0.6 and 1.5 times.
  The criterion benchmark measures case2383wp at 1.36 ms against 1.52 ms (−10%).
- **Why:** AMD's ordering costs about three times RCM's, and saves it back in the factor and the
  solve. Measured separately (medians of 51 runs):

  | Case | Ordering, RCM / AMD (ms) | Factor, RCM / AMD (ms) | Solve, RCM / AMD (ms) |
  |---|--:|--:|--:|
  | case8387pegase | 0.97 / 3.05 | 2.08 / 1.17 | 0.244 / 0.105 |
  | case2383wp | 0.22 / 0.69 | 0.56 / 0.25 | 0.064 / 0.024 |
  | case2869pegase | 0.32 / 0.93 | 0.45 / 0.35 | 0.057 / 0.031 |
  | case_ACTIVSg2000 | 0.19 / 0.76 | 1.10 / 0.27 | 0.095 / 0.025 |

  So AMD pays off wherever an analysis is reused: refactoring, several right-hand sides, and
  Newton's iterations, which is what Phase 3's LU and Phase 4.5 need.
- **Every case still agrees with MATPOWER**, to within 1.9e-12 of each quantity's largest value
  on these cases (case2736sp's generator outputs; 4.4e-12 with RCM). The worst case overall is
  still an LU case, at 1.1e-10.
- **Peak memory** is at most 14 MiB (case8387pegase; 17 MiB with RCM).

## The fill-in hypothesis: the Cholesky verdict

| Cholesky cases (59) | Baseline (unordered) | RCM | AMD |
|---|--:|--:|--:|
| Median backward error, ours | 0.29 ε | 0.26 ε | 0.25 ε |
| Median backward error, MATPOWER's | 0.25 ε | 0.25 ε | 0.25 ε |
| Median ratio, ours to MATPOWER's | 1.23 | 1.00 | 1.00 |
| Largest backward error, ours | 10.7 ε (case2383wp) | 2.35 ε (case8387pegase) | 0.78 ε (case2737sop) |
| Rank correlation: factor entries per row against our backward error | 0.45 | 0.20 | −0.03 |
| Rank correlation: factor entries per row against the ratio | 0.39 | 0.07 | −0.22 |

- **case8387pegase**, the case RCM left unexplained, falls from 2.35 ε to 0.71 ε as its fill
  falls another 2.9 times, below MATPOWER's 0.94 ε.
- **No Cholesky case is above 0.78 ε**, and our backward error no longer rises with fill.
- **The remaining ratios above 3** (case22 3.6, case4_dist and case51he 3.0) are ratios of
  rounding-level numbers on small radial networks with little or no fill: ours at most 0.58 ε,
  MATPOWER's at least 0.12 ε.
- **Verdict for Cholesky:** supported. Ordered to MATPOWER's level of fill, our backward error
  matches MATPOWER's. The LU cases remain to be tested with Phase 3's ordered LU.

## Cholesky cases

Fill divides the factor's entries by those of B′'s lower triangle. Times are from this
session's runs. ε is `f64::EPSILON`.

| Case | Buses | Fill, RCM | Fill, AMD | Solve, RCM (ms) | Solve, AMD (ms) | Backward error, RCM (ε) | Backward error, AMD (ε) | MATPOWER's (ε) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| case4_dist | 4 | 1.00 | 1.00 | 0.001 | 0.001 | 0.12 | 0.37 | 0.12 |
| case4gs | 4 | 1.00 | 1.00 | 0.001 | 0.002 | 0.00 | 0.37 | 0.37 |
| case5 | 5 | 1.00 | 1.00 | 0.002 | 0.002 | 0.35 | 0.12 | 0.12 |
| case6ww | 6 | 1.00 | 1.00 | 0.003 | 0.003 | 0.35 | 0.18 | 0.26 |
| case9 | 9 | 1.19 | 1.19 | 0.003 | 0.004 | 0.14 | 0.11 | 0.18 |
| case9Q | 9 | 1.19 | 1.19 | 0.003 | 0.004 | 0.14 | 0.11 | 0.18 |
| case9target | 9 | 1.19 | 1.19 | 0.003 | 0.004 | 0.18 | 0.18 | 0.23 |
| case10ba | 10 | 1.00 | 1.00 | 0.003 | 0.004 | 0.22 | 0.18 | 0.10 |
| case12da | 12 | 1.00 | 1.00 | 0.004 | 0.004 | 0.35 | 0.33 | 0.29 |
| case14 | 14 | 1.42 | 1.13 | 0.004 | 0.006 | 0.47 | 0.33 | 0.44 |
| case15da | 15 | 1.00 | 1.00 | 0.004 | 0.005 | 0.22 | 0.38 | 0.19 |
| case15nbr | 15 | 1.00 | 1.00 | 0.004 | 0.005 | 0.22 | 0.47 | 0.28 |
| case16am | 15 | 1.00 | 1.00 | 0.004 | 0.005 | 0.00 | 0.00 | 0.00 |
| case16ci | 16 | 1.00 | 1.00 | 0.004 | 0.004 | 0.18 | 0.18 | 0.35 |
| case17me | 17 | 1.00 | 1.00 | 0.004 | 0.005 | 0.51 | 0.44 | 0.61 |
| case18 | 18 | 1.00 | 1.00 | 0.005 | 0.006 | 0.43 | 0.46 | 0.53 |
| case18nbr | 18 | 1.00 | 1.00 | 0.005 | 0.006 | 0.35 | 0.35 | 0.47 |
| case22 | 22 | 1.00 | 1.00 | 0.006 | 0.007 | 0.42 | 0.42 | 0.12 |
| case24_ieee_rts | 24 | 1.43 | 1.31 | 0.007 | 0.010 | 0.40 | 0.32 | 0.27 |
| case28da | 28 | 1.00 | 1.00 | 0.007 | 0.009 | 0.63 | 0.63 | 0.32 |
| case30 | 30 | 1.37 | 1.19 | 0.009 | 0.012 | 0.15 | 0.09 | 0.09 |
| case30Q | 30 | 1.37 | 1.19 | 0.009 | 0.013 | 0.15 | 0.09 | 0.09 |
| case30pwl | 30 | 1.37 | 1.19 | 0.009 | 0.012 | 0.15 | 0.09 | 0.09 |
| case_ieee30 | 30 | 1.37 | 1.19 | 0.008 | 0.013 | 0.31 | 0.14 | 0.23 |
| case33bw | 33 | 1.00 | 1.00 | 0.008 | 0.011 | 0.23 | 0.28 | 0.23 |
| case33mg | 33 | 1.00 | 1.00 | 0.008 | 0.010 | 0.34 | 0.58 | 0.70 |
| case34sa | 34 | 1.00 | 1.00 | 0.009 | 0.011 | 0.30 | 0.40 | 0.27 |
| case38si | 38 | 1.00 | 1.00 | 0.009 | 0.012 | 0.34 | 0.45 | 0.23 |
| case39 | 39 | 1.46 | 1.25 | 0.010 | 0.014 | 0.21 | 0.37 | 0.37 |
| case51ga | 51 | 1.00 | 1.00 | 0.012 | 0.016 | 0.34 | 0.28 | 0.29 |
| case51he | 51 | 1.00 | 1.00 | 0.012 | 0.015 | 0.23 | 0.58 | 0.19 |
| case57 | 57 | 1.98 | 1.41 | 0.015 | 0.023 | 0.35 | 0.20 | 0.31 |
| case59 | 59 | 1.15 | 1.10 | 0.015 | 0.022 | 0.17 | 0.17 | 0.22 |
| case69 | 69 | 1.00 | 1.00 | 0.016 | 0.020 | 0.10 | 0.11 | 0.10 |
| case70da | 70 | 1.00 | 1.00 | 0.016 | 0.020 | 0.17 | 0.15 | 0.22 |
| case_RTS_GMLC | 73 | 1.87 | 1.43 | 0.019 | 0.032 | 0.21 | 0.19 | 0.25 |
| case74ds | 74 | 1.00 | 1.00 | 0.016 | 0.020 | 0.63 | 0.50 | 0.65 |
| case85 | 85 | 1.00 | 1.00 | 0.018 | 0.025 | 0.30 | 0.30 | 0.18 |
| case89pegase | 89 | 2.39 | 1.39 | 0.034 | 0.061 | 0.12 | 0.24 | 0.18 |
| case94pi | 94 | 1.00 | 1.00 | 0.019 | 0.026 | 0.33 | 0.60 | 0.42 |
| case118 | 118 | 1.63 | 1.28 | 0.029 | 0.057 | 0.17 | 0.24 | 0.26 |
| case118zh | 118 | 1.00 | 1.00 | 0.024 | 0.032 | 0.05 | 0.06 | 0.06 |
| case136ma | 136 | 1.00 | 1.00 | 0.026 | 0.036 | 0.06 | 0.07 | 0.08 |
| case141 | 141 | 1.00 | 1.00 | 0.026 | 0.037 | 0.45 | 0.24 | 0.45 |
| case145 | 145 | 1.88 | 1.20 | 0.071 | 0.126 | 0.79 | 0.66 | 1.36 |
| case_ACTIVSg200 | 200 | 2.17 | 1.42 | 0.057 | 0.098 | 0.21 | 0.21 | 0.43 |
| case_ACTIVSg500 | 500 | 2.02 | 1.30 | 0.168 | 0.237 | 0.26 | 0.27 | 0.39 |
| case533mt_hi | 533 | 1.00 | 1.00 | 0.118 | 0.166 | 0.04 | 0.04 | 0.04 |
| case533mt_lo | 533 | 1.00 | 1.00 | 0.118 | 0.167 | 0.31 | 0.31 | 0.31 |
| case1197 | 1197 | 1.00 | 1.00 | 0.203 | 0.301 | 0.01 | 0.01 | 0.01 |
| case1354pegase | 1354 | 2.43 | 1.33 | 0.587 | 0.757 | 0.26 | 0.12 | 0.19 |
| case_ACTIVSg2000 | 2000 | 11.22 | 2.41 | 2.334 | 1.395 | 0.43 | 0.24 | 0.21 |
| case2383wp | 2383 | 5.12 | 1.61 | 1.499 | 1.324 | 0.72 | 0.68 | 0.65 |
| case2736sp | 2736 | 6.37 | 1.64 | 1.716 | 1.502 | 0.74 | 0.49 | 0.49 |
| case2737sop | 2737 | 6.04 | 1.64 | 1.552 | 1.504 | 0.61 | 0.78 | 0.52 |
| case2746wop | 2746 | 5.33 | 1.68 | 1.472 | 1.528 | 0.22 | 0.25 | 0.24 |
| case2746wp | 2746 | 6.04 | 1.64 | 1.666 | 1.504 | 0.29 | 0.13 | 0.11 |
| case2869pegase | 2869 | 2.82 | 1.46 | 1.370 | 1.783 | 0.16 | 0.20 | 0.29 |
| case8387pegase | 8387 | 4.51 | 1.58 | 5.249 | 5.971 | 2.35 | 0.71 | 0.94 |
