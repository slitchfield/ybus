# Where the interior-point solvers spend their time, after the LDLᵀ, 2026-11-08

This is phase O4's second profile, taken after its first two targets:

- **the KKT matrix's analysis kept across iterations:** #54,
  [`kkt-analysis-2026-11-07.md`](kkt-analysis-2026-11-07.md);
- **the pivoted LDLᵀ:** #57–#59, [`ldlt-opf-2026-11-08.md`](ldlt-opf-2026-11-08.md), which the
  AC OPF's predictor-corrector settings now use.

It re-ranks the remaining targets. It uses the method of the first,
[`ipm-profile-2026-11-07.md`](ipm-profile-2026-11-07.md).

## Setup

| | |
|---|---|
| Code | The `dev/ipm-profile-2` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.99.1 |
| Build | `--profile profiling`: the release profile with debug line tables |
| Profiler | samply at 1 kHz, with `perf_event_paranoid` at 1 |
| Workloads | `profile_pf`'s `dcopf` and `acopf` (MIPS, LU), `dcopf-mehrotra` (the predictor-corrector method, LU), `dcopf-mehrotra-ldlt` (the same with LDLᵀ, new here) and `acopf-mehrotra` (LDLᵀ, its default). All at MATPOWER's OPF settings, each repeating the whole pipeline for 5 s and at least once |
| Cases | DC: case300 (new: a small case), case2383wp, case9241pegase, case_ACTIVSg10k and case_ACTIVSg70k. AC: case2383wp, case9241pegase, case_ACTIVSg25k and case_ACTIVSg70k, and case_SyntheticUSA, which fails on both solvers |
| Run | `scripts/profile-pf.sh` with `YBUS_PROFILE_CASES` and `YBUS_PROFILE_WORKLOADS`; `scripts/profile_summary.py` folds the samples into stages |

**The stages** are #53's. For the pivoted LDLᵀ, `profile_summary.py` now splits a sample under
`factor_indefinite` into four parts:

- **its symmetry check and permutation:** `is_hermitian` and `permute_symmetric`;
- **its supernode partition;**
- **its dense Bunch–Kaufman kernel;**
- **the rest:** front assembly and extend-add.

It also counts a residual's product inside refinement as "refinement", and the LDLᵀ's solves
as triangular solves. LU and Cholesky samples fold as before.

**Check:** the profiled times agree with the survey's for the same settings. For example, AC
case_ACTIVSg70k on the predictor-corrector method takes 95.4 s here against the survey's 96.5 s,
and AC case_SyntheticUSA on MIPS 593 s against 602 s.

## Stages: the solvers on LU

| Workload | Case | Per run | Factor | Analysis | OPF evaluation | Allocation | KKT assembly | Solves | Other |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| DC, MIPS | case300 | 3.2 ms | 48.5% | 11.1% | 0.3% | 16.8% | 6.7% | 2.6% | 14.0% |
| DC, new | case300 | 2.0 ms | 36.2% | 18.9% | 0.2% | 18.3% | 3.8% | 2.8% | 19.6% |
| DC, MIPS | case9241pegase | 0.36 s | 71.0% | 5.2% | 0.2% | 8.7% | 6.8% | 1.5% | 6.6% |
| DC, new | case9241pegase | 0.16 s | 58.1% | 11.8% | 0.1% | 11.4% | 5.4% | 2.8% | 10.4% |
| DC, MIPS | case_ACTIVSg70k | 92.6 s | 88.9% | 0.2% | 0.0% | 5.5% | 2.6% | 1.2% | 1.6% |
| DC, new | case_ACTIVSg70k | 21.4 s | 86.8% | 0.7% | 0.0% | 5.1% | 2.8% | 2.4% | 2.2% |
| AC, MIPS | case2383wp | 1.46 s | 58.9% | 0.6% | 23.7% | 7.6% | 6.3% | 1.1% | 1.8% |
| AC, MIPS | case9241pegase | 8.9 s | 67.5% | 0.5% | 17.0% | 7.0% | 5.0% | 1.3% | 1.7% |
| AC, MIPS | case_ACTIVSg25k | 65.8 s | 81.3% | 0.1% | 9.9% | 4.3% | 2.6% | 0.8% | 1.0% |
| AC, MIPS | case_ACTIVSg70k | 197 s | 87.3% | 0.1% | 7.2% | 2.5% | 1.7% | 0.6% | 0.6% |
| AC, MIPS | case_SyntheticUSA | 593 s | 86.7% | 0.1% | 7.8% | 2.4% | 1.8% | 0.7% | 0.5% |

"Other" includes parsing and the model build, which are up to 12.4% of the small case's 2.0 ms.

- **With the analysis kept (#54),** LU's numeric factorization takes 49–89%, more as cases
  grow.
- **MIPS stays on LU by default,** for MATPOWER parity, so its time is the LU's.

## Stages: the solvers on LDLᵀ

| Workload | Case | Per run | LDLᵀ kernel | Front assembly | Symmetry check | Supernodes | Solves | Refinement | OPF evaluation | Allocation | KKT assembly | Other |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| DC, new | case300 | 2.6 ms | 8.4% | 14.9% | 7.4% | 2.7% | 2.5% | 1.3% | 0.2% | 27.5% | 5.8% | 29.3% |
| DC, new | case2383wp | 0.14 s | 14.3% | 18.1% | 9.7% | 2.0% | 5.5% | 2.6% | 0.2% | 22.6% | 12.5% | 12.5% |
| DC, new | case9241pegase | 0.11 s | 8.1% | 13.2% | 8.4% | 1.9% | 9.3% | 3.0% | 0.1% | 19.1% | 8.8% | 28.1% |
| DC, new | case_ACTIVSg10k | 0.42 s | 13.8% | 17.9% | 8.8% | 2.6% | 8.4% | 3.3% | 0.1% | 19.7% | 11.3% | 14.1% |
| DC, new | case_ACTIVSg70k | 8.5 s | 14.6% | 20.7% | 8.3% | 2.4% | 16.1% | 2.6% | 0.1% | 18.2% | 8.5% | 8.5% |
| AC, new | case2383wp | 0.83 s | 9.6% | 10.7% | 6.1% | 1.4% | 4.0% | 1.5% | 37.5% | 12.2% | 10.3% | 6.7% |
| AC, new | case9241pegase | 4.3 s | 8.4% | 11.2% | 10.2% | 1.5% | 7.0% | 2.1% | 30.9% | 13.0% | 9.5% | 6.2% |
| AC, new | case_ACTIVSg25k | 26.9 s | 15.1% | 10.9% | 7.2% | 1.5% | 7.9% | 1.5% | 30.0% | 11.8% | 9.0% | 5.1% |
| AC, new | case_ACTIVSg70k | 95.4 s | 17.1% | 12.2% | 7.0% | 1.6% | 9.2% | 1.5% | 27.7% | 11.3% | 7.9% | 4.5% |
| AC, new | case_SyntheticUSA | 196 s | 18.5% | 12.2% | 7.6% | 1.7% | 5.0% | 0.8% | 30.1% | 11.7% | 8.1% | 4.3% |

"Other" includes the analysis (up to 10.8%, on the small DC cases), parsing (up to 10.0%), and
1.5–3.5% under `factor_indefinite`'s callers outside its four parts.

- **On AC, the OPF's evaluation is now the largest stage:** 28–38%, against 7–24% under LU. The
  factorization shrank about 5× around it.
- **The LDLᵀ factorization's arithmetic is a minority of its own time.**
  - **The dense kernel** takes 8–19%.
  - **Its bookkeeping** together, front assembly and extend-add, the symmetry check and
    permutation, and the supernode partition, takes 18–31%.
  - **The symmetry check and the supernode partition** are redone at every factorization,
    though neither changes while the pattern holds. They come to 7.5–11.7%.
- **On DC, allocation is the largest stage, 18–27%,** and with the bookkeeping it's most of the
  time. The kernel is 8–15%. That is the cost per front behind DC's small-case loss in #59.

## The fronts

At the first factorization, from a temporary measurement in `factor_indefinite`, not committed:

| | Order | Fronts | Front rows: median / p90 / p99 / max | Fronts of 4 rows or fewer | Delays | 2×2 pivots | Work in fronts over 64 rows |
|---|--:|--:|---|--:|--:|--:|--:|
| DC case300 | 670 | 612 | 4 / 9 / 14 / 18 | 54% | 331 | 256 | 0% |
| DC case9241pegase | 19,928 | 17,999 | 4 / 9 / 28 / 96 | 51% | 8,982 | 3,818 | 14% |
| DC case_ACTIVSg70k | 148,108 | 137,384 | 4 / 8 / 38 / 227 | 62% | 98,295 | 21,762 | 79% |
| AC case2383wp | 10,318 | 6,978 | 7 / 15 / 59 / 121 | 11% | 7,092 | 356 | 49% |
| AC case9241pegase | 39,855 | 25,662 | 7 / 16 / 54 / 196 | 11% | 19,326 | 3,311 | 53% |
| AC case_ACTIVSg25k | 108,012 | 67,628 | 7 / 14 / 60 / 274 | 12% | 65,216 | 11,430 | 86% |

"Work" is Σ m² e over the fronts, with m the front's rows and e its pivots.

- **Fundamental supernodes make almost one front per variable on DC** (0.9), with a median of
  4 rows. Every front allocates its dense square and walks its children, so the cost per front
  dominates.
- **AC's fronts are larger** (a median of 7 rows), and its large cases' work sits in fronts of
  over 64 rows.
- **Delays run to 0.5–0.9 per variable.** They come from the zero-diagonal constraint rows (#57),
  and each enlarges its parent front.

## What changed since the first profile

| Stage | First profile (#53) | Now |
|---|---|---|
| The KKT analysis | 27–35% of mid-size DC | 0.2–12% (kept; what's left is the one analysis on short solves) |
| The factorization, AC, predictor-corrector | 66–78% (LU) | 30–42% (LDLᵀ, all its parts) |
| The OPF's evaluation, AC, predictor-corrector | 7–19% | 28–38% |
| Allocation | 5–18% | 2–8% under LU on large cases; 11–27% under LDLᵀ |

## The candidates, with their bounds

Each bound is what removing the stage entirely would give, `1 / (1 − share)`, as #53 bounded its
targets. Each is for the setting where the stage matters.

| Candidate | Share | Bound | Where |
|---|--:|--:|---|
| **The LDLᵀ's bookkeeping:** the symmetry check and supernode partition kept with the analysis, no allocation per front and per pivot, relaxed amalgamation of the tiny fronts | DC: 24–31% bookkeeping, plus allocation at 18–27%, much of it the fronts'. AC: 18–23%, plus 11–13% | From the bookkeeping alone: DC 1.3–1.5×, AC 1.2–1.3×. With all the allocation too, at most DC 1.7–2.1×, AC 1.4–1.55× | the predictor-corrector method on LDLᵀ. It's the cost per front that keeps DC's default on LU |
| **The OPF's evaluation** | AC: 28–38% | 1.4–1.6× | AC, either solver; the largest single AC stage now |
| **Allocation outside the factorization, and KKT assembly** | 6–13% assembly, plus allocation | about 1.1–1.2× | both solvers |
| **The dense kernel** (blocked updates) | 8–19% | 1.1–1.2× | LDLᵀ, mostly the largest AC cases |
| **A KKT-aware ordering** for the delays | — | not bounded by a stage: it shrinks the fronts and their count | LDLᵀ |
| **The LU factorization** | 49–89% | — | MIPS, which keeps LU for MATPOWER parity; outside O4's targets |

**The proposed order** comes from these bounds, and is to be decided with the maintainer at
review:

1. **The LDLᵀ's bookkeeping.**
   - **What it gains:** the largest bound on DC and a sizeable one on AC. Its first steps are
     cheap: the symmetry check and the supernode partition only need to move into the analysis.
   - **What it decides:** whether DC's predictor-corrector default qualifies for LDLᵀ under
     #59's rule.
2. **The OPF's evaluation,** now the largest AC stage.
3. **Allocation and KKT assembly.**
4. **The kernel and the KKT ordering,** once the bookkeeping is gone and they're what's left.
