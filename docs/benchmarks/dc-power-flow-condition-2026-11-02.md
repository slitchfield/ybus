# DC power flow: B′'s condition and iterative refinement, 2026-11-02

Phase 5's DC survey: the DC power flow on all 78 MATPOWER 8.1 cases, as in the AMD and LU
records, with three new columns: B′'s estimated 1-norm condition number (`condition_estimate_1`,
Hager's estimator with Higham's refinements), and the normwise backward error of the angles after
iterative refinement (`solve_refined`, at most 10 corrections), with the corrections taken. The
survey's output is in [`dc-power-flow-condition-2026-11-02.csv`](dc-power-flow-condition-2026-11-02.csv).

## Setup

| | |
|---|---|
| Code | The `dev/refinement-parallel` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Run | `scripts/benchmark-dc.sh`, AMD ordering, from 21:52 to 21:53 EST |

B′ and P′ for the new columns are built independently of the solve (`tests/common/accuracy.rs`),
factored as the solve factors them (Cholesky, or LU where B′ isn't positive definite), and the
refined angles are measured by the same normwise backward error as the `backward_error` column.

## Results

- **The rest is unchanged.** Every column the AMD and LU records hold matches them on all 78
  cases, backward errors and differences from MATPOWER included, apart from time and memory.
- **B′'s condition ranges over eight orders of magnitude.** κ₁ runs from 7.8 (case4gs) to 5.8e8
  (case_SyntheticUSA), with a median of 1.5e4. The 19 cases needing LU (some branch has a
  negative reactance) are the badly conditioned ones: their median is 4.8e6, against 4.1e3 for
  the 59 Cholesky cases. The largest: case_SyntheticUSA 5.8e8, case_ACTIVSg70k 4.8e8,
  case16am 1.6e8, case13659pegase 1.3e8, case_ACTIVSg10k 5.5e7.
- **Condition doesn't explain the backward error.** The rank correlation between log κ₁ and the
  backward error is −0.08 (for MATPOWER's angles, −0.15). That is what a backward-stable solve
  gives: the backward error stays at rounding level however badly B′ is conditioned, and it's
  the forward error that grows with κ₁. Since Phase 3's orderings, no case has a gap left for
  conditioning to explain: the worst is 0.78 ε (case2737sop).
- **Refinement has little left to do.** It took no correction on 24 cases, one on 48 and two
  on 6. The normwise backward error fell on 32 cases and rose on 22, always within rounding: at
  most 0.75 ε afterwards, against 0.78 ε before. Refinement stops on the componentwise backward
  error, which it minimizes, so the normwise figure can tick up at this level, as on
  case_ACTIVSg70k (0.37 ε to 0.75 ε). The largest gain is case2737sop's, 0.78 ε to 0.31 ε in
  one correction.

| Case | B′ | κ₁ (estimate) | Backward error | After refinement | Corrections |
|---|---|--:|--:|--:|--:|
| case4gs | Cholesky | 7.8 | 0.37 ε | 0.37 ε | 0 |
| case2383wp | Cholesky | 7.5e5 | 0.68 ε | 0.46 ε | 2 |
| case2737sop | Cholesky | 1.1e6 | 0.78 ε | 0.31 ε | 1 |
| case13659pegase | LU | 1.3e8 | 0.33 ε | 0.53 ε | 1 |
| case_ACTIVSg70k | LU | 4.8e8 | 0.37 ε | 0.75 ε | 1 |
| case_SyntheticUSA | LU | 5.8e8 | 0.18 ε | 0.18 ε | 2 |

So refinement isn't needed for the DC power flow on these cases, and `DcModel::solve` doesn't
use it; it's available in `sparse_mat` for systems that need it. The condition estimate says
how far a case's angles can move for a given backward error: on case_SyntheticUSA, a backward
error of 0.18 ε allows a relative forward error up to about κ₁ · 0.18 ε ≈ 2e-8.

## What a record may hold

Case names and sizes, and the time, memory and accuracy of this crate's code on them: no case
data or power flow results.
