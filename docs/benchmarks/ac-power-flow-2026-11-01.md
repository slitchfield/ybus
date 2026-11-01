# AC power flow on every MATPOWER case, 2026-11-01

Phase 4.5's survey: `powerflow`'s AC power flow on all 78 MATPOWER 8.1 cases, by Newton's
method (NR) and fast-decoupled XB and BX, against MATPOWER's own `runpf` results for each method.
The survey's output, one row per case and method, is in
[`ac-power-flow-2026-11-01.csv`](ac-power-flow-2026-11-01.csv).

## Setup

| | |
|---|---|
| Code | The `dev/ac-survey` branch, at the commit that adds this record |
| Machine and system | As in the DC records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Reference | MATPOWER 8.1 under Octave 8.4.0, run locally by `scripts/generate-reference.sh` (results not committed): `runpf` with its defaults for NR (MATPOWER's MP-Core implementation), and the legacy `fdpf` for fast-decoupled, which `powerflow` follows |
| Run | `scripts/benchmark-ac.sh` with `YBUS_BENCH_METHOD` set to `nr`, `fdxb` and `fdbx` in turn, from 11:37 to 11:40 EST, with the default limits (1 hour and 10 GiB per case) |

Each step is timed by its median over repeated runs, for up to 1 s or 101 runs. A rerun gave the
same results in every column but time and memory.

## Results

| | NR | FDXB | FDBX |
|---|--:|--:|--:|
| Cases that converge | 77 | 67 | 66 |
| Cases where MATPOWER converges | 77 | 67 | 66 |
| Converged with MATPOWER's iteration count | all 77 | all 67 | all 66 |
| Iterations, median (most) | 3 (6) | 11 (30) | 12 (26) |
| Worst relative difference from MATPOWER | 2.3e-9 | 4.0e-9 | 1.9e-10 |
| Solve time, all converged cases | 2.17 s | 0.39 s | 0.39 s |
| Slowest converged case | 860 ms (case_SyntheticUSA) | 203 ms (case_ACTIVSg70k) | 211 ms (case_ACTIVSg70k) |
| Peak memory, any case | 294 MiB | 160 MiB | 160 MiB |

- **Convergence matches MATPOWER exactly.** We converge on exactly the cases MATPOWER
  converges on, for every method, and in exactly MATPOWER's number of iterations.
  - **NR:** case16am fails in both.
  - **Fast-decoupled:** fails in both on case9target, case16am, the RTE cases case1888rte,
    case1951rte, case2868rte, case6468rte, case6470rte, case6495rte and case6515rte,
    case_ACTIVSg10k and case_SyntheticUSA, and for BX also case94pi.
- **The results agree with MATPOWER's.** The relative difference is the largest absolute
  difference over the largest magnitude of that quantity in MATPOWER's results for the case, as
  in the DC records.
  - **Worst case:** case141, a distribution feeder whose angles are all under 0.3°. Its NR
    angles differ from MATPOWER's by 6.8e-10° absolute.
  - **Worst on cases of 1,000 buses or more, any method:** voltage magnitudes within 1.2e-12
    p.u. and flows within 2.7e-9 MW/MVAr, absolute.
- **Fast-decoupled is about three times faster on the cases it solves.** On those 66–67
  cases, NR takes 1.15 s together against 0.39 s, though fast-decoupled needs about four times
  the iterations (medians 11–12 against 3).

## Generator outputs: MATPOWER's two implementations differ

In 10 cases, NR's voltages and flows agree with MATPOWER's within the figures above, but
generator outputs differ by up to 185 MVAr, and one reference value is NaN. The 10 cases are
case2848rte, case2868rte, case3012wp, case3120sp, case3375wp, case6470rte, case_ACTIVSg10k,
case_ACTIVSg25k, case_ACTIVSg70k and case_SyntheticUSA. They're excluded from the NR row above.

- **The cause is MATPOWER, not `powerflow`.** The NR reference comes from MATPOWER's default
  MP-Core implementation, whose post-processing shares reactive output among generators
  differently from the legacy `pfsoln`, which `powerflow` follows. MP-Core and legacy differ in:
  - buses whose generators all have zero reactive range (Qmin = Qmax = 0): on case6470rte,
    MP-Core reports NaN for three generators where the legacy code gives finite values
  - some single generators, for example case3012wp's generator 1: MP-Core −17.4 MVAr against
    the legacy 84.1 MVAr
- **Against MATPOWER's legacy NR** (`exp.use_legacy_core`), all 10 cases agree within 2.9e-9
  MW/MVAr in generator output, with the same iteration counts.
- **The same 10 cases under fast-decoupled,** where the reference is the legacy code, show no
  generator differences.

## Benchmarks: case2383wp

Criterion (`cargo bench -p powerflow --bench ac_power_flow`):

| Step | Time |
|---|--:|
| Building the AC model | 575 µs |
| NR solve (6 iterations) | 12.9 ms |
| FDXB solve (18 iterations) | 4.15 ms |
| FDBX solve (14 iterations) | 3.87 ms |
| One Newton Jacobian, built at the solution | 1.12 ms |

- **The Jacobian build is about half of a Newton solve:** six builds take 6.7 ms of the
  12.9 ms. It's built from Ybus with matrix operations each iteration, as decided for Phase 4.5,
  so the build allocates every intermediate matrix. An assembly-map Jacobian, which refills
  values in place, is now a follow-up item.
- **The DC power flow is unaffected.** Its solve on the same case measured 1.474 and 1.460 ms on
  the code from before the AC work, and 1.477 ms on this branch, back to back.

## All cases

Sorted by bus count.

| Case | Buses | NR iterations | NR solve (ms) | FDXB iterations | FDXB solve (ms) | FDBX iterations | FDBX solve (ms) |
|---|--:|--:|--:|--:|--:|--:|--:|
| case4_dist | 4 | 3 | 0.012 | 11 | 0.006 | 11 | 0.006 |
| case4gs | 4 | 3 | 0.011 | 6 | 0.005 | 6 | 0.005 |
| case5 | 5 | 3 | 0.013 | 5 | 0.005 | 5 | 0.005 |
| case6ww | 6 | 3 | 0.016 | 9 | 0.009 | 8 | 0.008 |
| case9 | 9 | 4 | 0.032 | 6 | 0.010 | 6 | 0.011 |
| case9Q | 9 | 4 | 0.031 | 8 | 0.011 | 7 | 0.011 |
| case9target | 9 | 6 | 0.043 | no conv. | — | no conv. | — |
| case10ba | 10 | 4 | 0.033 | 12 | 0.013 | 17 | 0.015 |
| case12da | 12 | 3 | 0.029 | 14 | 0.016 | 13 | 0.015 |
| case14 | 14 | 2 | 0.028 | 6 | 0.016 | 8 | 0.017 |
| case15da | 15 | 3 | 0.037 | 11 | 0.016 | 12 | 0.017 |
| case15nbr | 15 | 3 | 0.038 | 8 | 0.015 | 8 | 0.015 |
| case16am | 15 | no conv. | — | no conv. | — | no conv. | — |
| case16ci | 16 | 3 | 0.036 | 7 | 0.014 | 8 | 0.015 |
| case17me | 17 | 4 | 0.048 | 22 | 0.025 | 23 | 0.025 |
| case18 | 18 | 4 | 0.050 | 8 | 0.017 | 10 | 0.018 |
| case18nbr | 18 | 3 | 0.041 | 8 | 0.018 | 9 | 0.018 |
| case22 | 22 | 3 | 0.048 | 11 | 0.022 | 11 | 0.022 |
| case24_ieee_rts | 24 | 4 | 0.073 | 8 | 0.027 | 9 | 0.027 |
| case28da | 28 | 3 | 0.061 | 15 | 0.032 | 14 | 0.030 |
| case30 | 30 | 3 | 0.084 | 11 | 0.035 | 8 | 0.032 |
| case30Q | 30 | 3 | 0.087 | 11 | 0.036 | 8 | 0.034 |
| case30pwl | 30 | 3 | 0.085 | 11 | 0.036 | 8 | 0.033 |
| case_ieee30 | 30 | 2 | 0.064 | 7 | 0.034 | 8 | 0.034 |
| case33bw | 33 | 3 | 0.069 | 14 | 0.035 | 13 | 0.034 |
| case33mg | 33 | 4 | 0.087 | 15 | 0.036 | 15 | 0.036 |
| case34sa | 34 | 3 | 0.072 | 19 | 0.043 | 16 | 0.039 |
| case38si | 38 | 4 | 0.096 | 15 | 0.043 | 15 | 0.043 |
| case39 | 39 | 1 | 0.051 | 4 | 0.035 | 4 | 0.035 |
| case51ga | 51 | 4 | 0.127 | 15 | 0.057 | 15 | 0.056 |
| case51he | 51 | 3 | 0.106 | 9 | 0.046 | 10 | 0.047 |
| case57 | 57 | 3 | 0.181 | 7 | 0.064 | 9 | 0.071 |
| case59 | 59 | 5 | 0.219 | 16 | 0.077 | 16 | 0.076 |
| case60nordic | 60 | 1 | 0.083 | 5 | 0.057 | 5 | 0.054 |
| case69 | 69 | 4 | 0.164 | 17 | 0.082 | 14 | 0.073 |
| case70da | 70 | 4 | 0.175 | 14 | 0.072 | 15 | 0.072 |
| case_RTS_GMLC | 73 | 4 | 0.240 | 10 | 0.086 | 10 | 0.088 |
| case74ds | 74 | 3 | 0.137 | 10 | 0.068 | 11 | 0.069 |
| case85 | 85 | 4 | 0.211 | 19 | 0.098 | 16 | 0.092 |
| case89pegase | 89 | 5 | 0.589 | 8 | 0.196 | 8 | 0.198 |
| case94pi | 94 | 4 | 0.230 | 21 | 0.111 | no conv. | — |
| case118 | 118 | 3 | 0.326 | 8 | 0.136 | 7 | 0.131 |
| case118zh | 118 | 4 | 0.280 | 18 | 0.129 | 16 | 0.118 |
| case136ma | 136 | 4 | 0.335 | 8 | 0.104 | 12 | 0.116 |
| case141 | 141 | 3 | 0.292 | 12 | 0.129 | 14 | 0.135 |
| case145 | 145 | 3 | 0.632 | 15 | 0.333 | 16 | 0.336 |
| case_ACTIVSg200 | 200 | 2 | 0.550 | 5 | 0.264 | 5 | 0.267 |
| case300 | 300 | 5 | 1.426 | 9 | 0.451 | 9 | 0.443 |
| case_ACTIVSg500 | 500 | 3 | 1.432 | 8 | 0.651 | 8 | 0.648 |
| case533mt_hi | 533 | 3 | 1.235 | 8 | 0.469 | 12 | 0.532 |
| case533mt_lo | 533 | 3 | 1.226 | 7 | 0.465 | 12 | 0.526 |
| case1197 | 1197 | 3 | 2.504 | 12 | 0.968 | 14 | 1.023 |
| case1354pegase | 1354 | 4 | 5.328 | 8 | 1.861 | 9 | 1.923 |
| case1888rte | 1888 | 2 | 4.681 | no conv. | — | no conv. | — |
| case1951rte | 1951 | 3 | 6.030 | no conv. | — | no conv. | — |
| case_ACTIVSg2000 | 2000 | 3 | 9.800 | 11 | 4.550 | 11 | 4.539 |
| case2383wp | 2383 | 6 | 12.4 | 18 | 4.066 | 14 | 3.758 |
| case2736sp | 2736 | 4 | 10.3 | 12 | 4.154 | 13 | 4.295 |
| case2737sop | 2737 | 5 | 12.5 | 12 | 4.286 | 14 | 4.569 |
| case2746wop | 2746 | 4 | 10.5 | 13 | 4.147 | 15 | 4.372 |
| case2746wp | 2746 | 4 | 10.5 | 15 | 4.261 | 15 | 4.302 |
| case2848rte | 2848 | 2 | 7.267 | 30 | 6.491 | 26 | 6.221 |
| case2868rte | 2868 | 5 | 13.3 | no conv. | — | no conv. | — |
| case2869pegase | 2869 | 6 | 16.5 | 9 | 4.792 | 11 | 5.088 |
| case3012wp | 3012 | 3 | 10.0 | 9 | 4.187 | 12 | 4.628 |
| case3120sp | 3120 | 6 | 16.0 | 13 | 4.766 | 18 | 5.493 |
| case3375wp | 3374 | 2 | 9.038 | 8 | 4.657 | 12 | 5.180 |
| case6468rte | 6468 | 3 | 24.1 | no conv. | — | no conv. | — |
| case6470rte | 6470 | 3 | 24.1 | no conv. | — | no conv. | — |
| case6495rte | 6495 | 2 | 19.5 | no conv. | — | no conv. | — |
| case6515rte | 6515 | 3 | 25.1 | no conv. | — | no conv. | — |
| case8387pegase | 8387 | 3 | 34.5 | 16 | 18.4 | 15 | 18.3 |
| case9241pegase | 9241 | 6 | 63.9 | 14 | 20.2 | 15 | 21.0 |
| case_ACTIVSg10k | 10000 | 4 | 45.5 | no conv. | — | no conv. | — |
| case13659pegase | 13659 | 5 | 77.0 | 14 | 28.9 | 15 | 29.9 |
| case_ACTIVSg25k | 25000 | 4 | 141.7 | 17 | 60.5 | 11 | 52.2 |
| case_ACTIVSg70k | 70000 | 6 | 684.1 | 16 | 207.7 | 17 | 210.3 |
| case_SyntheticUSA | 82000 | 6 | 860.3 | no conv. | — | no conv. | — |
