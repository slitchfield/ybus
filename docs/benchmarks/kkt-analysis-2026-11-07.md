# Keeping the KKT matrix's analysis, 2026-11-07

Phase O4's first target, chosen from [`ipm-profile-2026-11-07.md`](ipm-profile-2026-11-07.md).
Both interior-point solvers analyzed the KKT matrix (COLAMD's ordering and the symbolic LU)
afresh at every iteration, though its pattern doesn't change. They now analyze it once and keep
the analysis while the pattern holds, with the same results bit for bit, as decided for O4.

## Setup

| | |
|---|---|
| Code | The `dev/kkt-analysis` branch at `43b007f`, against O3's survey ([`mehrotra-opf-2026-11-07.md`](mehrotra-opf-2026-11-07.md)): MIPS at `c3826ea` and the predictor-corrector method at `b7a28ea`, whose code differs from this branch's only by this change |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.99.1, release profile; both surveys ran on the same machine and day |
| Run | `scripts/benchmark-opf.sh`, for DC and AC, MIPS and the predictor-corrector method, all 73 cases with costs, with the default limits; results in [`…-dc-mips.csv`](kkt-analysis-2026-11-07-dc-mips.csv), [`…-dc-mehrotra.csv`](kkt-analysis-2026-11-07-dc-mehrotra.csv), [`…-ac-mips.csv`](kkt-analysis-2026-11-07-ac-mips.csv) and [`…-ac-mehrotra.csv`](kkt-analysis-2026-11-07-ac-mehrotra.csv) |

## The change

`SymbolicLu` holds the column ordering, chosen from the pattern alone, and a copy of that
pattern. `SymbolicLu::factor` refuses a matrix with another pattern (`PatternMismatch`).

**How the solvers use it:** both now share a crate-private `KktAnalysis`. It factors each
iteration's KKT matrix with the kept analysis, and analyzes afresh only when there is none yet,
or when `factor` reports a different pattern.

**Why the results can't change:** on the same pattern, the same ordering gives the same
factorization, so the results are bit for bit the same.

**The count:** each solver's result reports the analyses it took (`MipsSolution::analyses`,
`MehrotraSolution::analyses`, and `SolverResult::analyses` in the OPF), and the survey has a
column for it.

## The same results, checked

- **A snapshot:** a temporary example (not committed) hashed both solvers' DC and AC OPF
  solutions on every normalized case up to 3 MB, 72 cases, at `main` and after the change:
  - the case quantities;
  - the solver's point and multipliers;
  - the status, the iterations and the KKT residuals.

  The hashes are identical.
- **The survey:** every result column matches O3's survey, on all 292 solves (73 cases, DC and
  AC, two solvers). That covers the outcome, the iterations, the differences from MATPOWER
  in the objective, prices, primal solution and multipliers, and the four KKT residuals.
- **The tests:** MIPS's tests against MATPOWER, the pinned iteration counts of both solvers,
  and the OPF's MATPOWER reference tests pass unchanged.

## How often the analysis is kept

**Every one of the 292 solves analyzed exactly once.** The OPF's matrices are built from
triplets and sparse products, which keep every position the operands reach, even where values
cancel, so the KKT pattern never changes.

**The fallback is exercised by the tests:**

- **Two of MIPS's test problems,** nonlinear2 and nonlinear3, analyze twice. They build their
  Jacobians with `from_dense`, which stores no zeros, and start with a coordinate at zero.
- **On the new solver,** only nonlinear3 does; it moves nonlinear2's coordinate off its bound
  first.
- **A test problem whose Hessian stores more entries** past a point solves correctly on both
  solvers, with two analyses.

## Speed

Time before over time after, per case:

| | Median | Range | Under 1,000 buses (46) | 1,000–10,000 (22) | Over 10,000 (5) | Total time |
|---|--:|--:|--:|--:|--:|--:|
| DC, MIPS | 2.15× | 1.37–3.31× | 2.21× | 2.09× | 1.47× | 327 → 229 s (1.43×) |
| DC, predictor-corrector | 1.66× | 1.28–2.46× | 1.63× | 1.88× | 1.48× | 98.7 → 70.6 s (1.40×) |
| AC, MIPS | 1.38× | 1.05–1.80× | 1.47× | 1.24× | 1.07× | 1,017 → 957 s (1.06×) |
| AC, predictor-corrector | 1.33× | 1.01–1.67× | 1.40× | 1.29× | 1.15× | 1,353 → 1,231 s (1.10×) |

The largest cases:

| Case | DC, MIPS | DC, new | AC, MIPS | AC, new |
|---|--:|--:|--:|--:|
| case9241pegase | 0.66 → 0.35 s (1.87×) | 0.23 → 0.15 s (1.57×) | 11.4 → 8.8 s (1.30×) | 10.3 → 8.0 s (1.28×) |
| case_ACTIVSg10k | 4.3 → 2.4 s (1.82×) | 1.11 → 0.69 s (1.61×) | 28.7 → 25.2 s (1.14×) | 24.7 → 21.5 s (1.15×) |
| case_ACTIVSg25k | 29.6 → 20.1 s (1.47×) | 5.6 → 3.8 s (1.48×) | 71.4 → 66.6 s (1.07×) | 75.6 → 65.2 s (1.16×) |
| case_ACTIVSg70k | 130 → 95 s (1.37×) | 30.0 → 21.8 s (1.37×) | 209 → 200 s (1.05×) | 322 → 296 s (1.09×) |
| case_SyntheticUSA | 148 → 104 s (1.42×) | 57.7 → 41.9 s (1.38×) | 626 → 599 s (1.05×) | 807 → 751 s (1.07×) |

- **On DC, more than the profile's bound.** The profile bounded the gain at 1.4–1.5× on mid-size
  DC, from an analysis share of 27–35%. Measured, case9241pegase gains 1.87× (MIPS) and
  case_ACTIVSg10k 1.82×.
  - **Why the bound understated it:** the profile counts a sample's innermost stage, so the
    analysis's own allocation (COLAMD's workspace, and `analyze_lu`'s copy of the pattern) fell
    under "allocation and copying", which was 16–18% there. Keeping the analysis removes those
    allocations too.
  - **Below 1,000 buses,** which the profile didn't cover, the analysis weighs more, and the gain
    reaches 3.31×.
- **On AC, as bounded.** The analysis's share is smaller (5–15%), so the gain is 1.3× on
  mid-size cases and 1.05–1.16× on the largest, where the factorization takes 66–80%.
- **Where the time goes now:** the factorization, even more than before. That is O4's next
  target, a symmetric indefinite LDLᵀ.
