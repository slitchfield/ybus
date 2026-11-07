# Where the interior-point solvers spend their time, 2026-11-07

Phase O4's profile: where the optimal power flow spends its time on both interior-point solvers
(MIPS, and phase O3's predictor-corrector method), measured with a sampling profiler, to decide
what to accelerate first. It uses the same method as phase O2's
[`profile-2026-11-04.md`](profile-2026-11-04.md).

## Setup

| | |
|---|---|
| Code | The `dev/ipm-profile` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.99.1 |
| Build | `--profile profiling`: the release profile with debug line tables |
| Profiler | samply 0.13.1 at 1 kHz, with `perf_event_paranoid` at 1 |
| Workloads | `profile_pf`'s `dcopf` and `acopf` (MIPS) and `dcopf-mehrotra` and `acopf-mehrotra` (the predictor-corrector method), at MATPOWER's OPF settings, each repeating the whole pipeline for 5 s and at least once |
| Cases | DC: case9241pegase, case_ACTIVSg10k and case_ACTIVSg70k. AC: case2383wp, case9241pegase and case_ACTIVSg25k, which converge, and the large AC cases case_ACTIVSg70k and case_SyntheticUSA, which fail on both solvers within O3's limits |
| Run | `scripts/profile-pf.sh` with `YBUS_PROFILE_CASES` and `YBUS_PROFILE_WORKLOADS`; `scripts/profile_summary.py` folds the samples into stages |

A sample's stage is its innermost frame that belongs to one, as in O2, with two refinements for
the OPF. O2's 30 power-flow profiles fold into exactly the same stages with them.

- **OPF evaluation and derivatives:** any frame of the problem's own objective, constraints or
  Lagrangian Hessian anywhere on the stack, unless the sample is allocation.
- **KKT assembly:** sparse products, scalings, sums and assembly from triplets, when an
  interior-point solver called them.

"Factor" is the KKT matrix's numeric LU factorization (`gilbert_peierls`). "Analysis" is its
COLAMD ordering and symbolic LU (`analyze_lu`), which both solvers redo every iteration.
Allocation that libc leaves unsymbolized is counted with its caller, as in O2.

## Stages

| Workload | Case | Per run | Factor | Analysis | OPF evaluation | Allocation | KKT assembly | Solves | Other |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| DC, MIPS | case9241pegase | 0.68 s | 38.9% | 34.5% | 0.1% | 17.3% | 4.1% | 1.0% | 2.2% |
| DC, new | case9241pegase | 0.24 s | 38.3% | 32.7% | 0.1% | 16.7% | 3.1% | 1.8% | 2.4% |
| DC, MIPS | case_ACTIVSg10k | 4.4 s | 42.4% | 31.0% | 0.1% | 17.9% | 4.4% | 1.1% | 2.4% |
| DC, new | case_ACTIVSg10k | 1.1 s | 46.4% | 27.4% | 0.1% | 16.0% | 3.8% | 2.2% | 2.5% |
| DC, MIPS | case_ACTIVSg70k | 129 s | 65.8% | 17.8% | 0.0% | 11.9% | 2.2% | 0.9% | 1.2% |
| DC, new | case_ACTIVSg70k | 29 s | 65.1% | 17.6% | 0.0% | 11.4% | 2.1% | 1.8% | 1.2% |
| AC, MIPS | case2383wp | 1.8 s | 46.3% | 14.2% | 19.3% | 12.9% | 4.5% | 0.9% | 1.6% |
| AC, new | case2383wp | 1.7 s | 46.8% | 14.6% | 19.0% | 11.8% | 4.0% | 1.6% | 1.9% |
| AC, MIPS | case9241pegase | 11.6 s | 53.3% | 14.7% | 13.5% | 12.0% | 3.9% | 0.9% | 1.5% |
| AC, new | case9241pegase | 10.4 s | 53.5% | 14.2% | 13.2% | 11.7% | 3.6% | 2.0% | 1.5% |
| AC, MIPS | case_ACTIVSg25k | 76 s | 72.1% | 7.2% | 8.8% | 7.8% | 2.3% | 0.8% | 0.9% |
| AC, new | case_ACTIVSg25k | 77 s | 65.9% | 8.8% | 10.8% | 9.2% | 2.5% | 1.5% | 1.2% |
| AC, MIPS | case_ACTIVSg70k | 221 s | 79.9% | 5.3% | 6.7% | 5.3% | 1.5% | 0.6% | 0.6% |
| AC, new | case_ACTIVSg70k | 338 s | 74.8% | 6.5% | 8.2% | 6.4% | 1.9% | 1.1% | 0.9% |
| AC, MIPS | case_SyntheticUSA | 670 s | 78.7% | 5.8% | 7.1% | 5.5% | 1.6% | 0.6% | 0.6% |
| AC, new | case_SyntheticUSA | 842 s | 78.0% | 5.7% | 7.2% | 5.5% | 1.6% | 1.2% | 0.8% |

Parsing and the model build are at most 3.3% (DC case9241pegase on the new solver, whose whole
run is 0.24 s) and are left out. Every run spent 1.0 thread on average.

## The KKT matrices

At the first iteration, by MIPS, from a temporary measurement not committed:

| KKT at the first iteration | Order | Nonzeros | L | U | (L + U) / KKT |
|---|--:|--:|--:|--:|--:|
| DC case9241pegase | 19,928 | 95,477 | 204,954 | 315,443 | 5.5 |
| DC case_ACTIVSg25k | 53,780 | 242,727 | 659,084 | 1,082,289 | 7.2 |
| DC case_ACTIVSg70k | 148,108 | 680,277 | 2,303,813 | 3,884,171 | 9.1 |
| AC case9241pegase | 39,855 | 460,525 | 935,919 | 1,587,291 | 5.5 |
| AC case_ACTIVSg25k | 108,012 | 1,045,769 | 2,881,787 | 5,015,947 | 7.6 |

The KKT matrix `[W Jgᵀ; Jg 0]` is symmetric. The LU factors it as unsymmetric, with threshold
pivoting, which leaves U about 1.5–1.7 times L.

## What the profile shows

- **The numeric factorization dominates,** from 38% (DC, mid-size) to 80% (AC, the largest),
  and more so as cases grow, since fill grows faster than the matrix (5.5× to 9.1×).
- **The analysis comes second:** COLAMD and the symbolic LU, redone every iteration though the
  KKT matrix's pattern doesn't change. It's a third of the time on mid-size DC (27–35%), 18% on
  DC case_ACTIVSg70k, and 5–15% on AC.
- **Allocation** is 5–18%, the most on DC.
- **The OPF's own evaluation** (the costs, the constraints and their derivatives, the Hessian
  above all) is 7–19% on AC, and nothing on DC, whose problem is a fixed QP.
- **The second solve costs little:** the triangular solves are 0.6–2.2% in all. The
  predictor-corrector method's per-iteration cost matches MIPS's, as O3 measured.
- **The same hot spots in both solvers:** their stage shares agree within a few points on every
  case. Their per-run times differ by their iteration counts, so the targets serve both.

## Decision: what to accelerate first

Decided with the maintainer, 2026-11-07. Each bound is the speedup if the stage cost nothing, or
here, where noted, half as much (Amdahl's law). It's a bound, not a promise.

**The rule on results:** MIPS's default path stays bit for bit the same, keeping its MATPOWER
parity (MATPOWER's iteration counts, 1e-8 against MATPOWER). A change that alters the
arithmetic comes in as an option, or as the predictor-corrector solver's default, judged by
O3's bar.

1. **Analyze the KKT matrix once per solve.** Its pattern is fixed across the iterations, so
   COLAMD's ordering and the symbolic analysis can be computed at the first and reused. The
   ordering is the same, so the results stay bit for bit the same, on MIPS's default too.
   - **Bound:** 1.4–1.5× on mid-size DC, 1.2× on DC case_ACTIVSg70k, 1.17× on mid-size AC, and
     1.06× on the largest.
   - **First:** the KKT assembly must keep every position stored, even where a value is zero, so
     the pattern really is the same at every iteration.
2. **A symmetric indefinite LDLᵀ factorization,** with Bunch–Kaufman-style pivoting, in
   `sparse_mat`, written from the literature.
   - **What it saves:** it stores about one triangle where the LU stores two, with about half
     the work.
   - **What it adds:** the inertia, with which the predictor-corrector method can correct an
     indefinite Hessian (hs71's constrained local maximum, from O3's record).
   - **Bound:** if it halves the factorization, 1.3× on mid-size AC, and 1.5–1.7× on the
     largest cases, DC and AC.
   - **Where it goes:** it changes the arithmetic, so it's an option for MIPS, and a candidate
     default for the predictor-corrector method.
3. **Allocation in the KKT assembly and the solvers' loops:** buffers reused across iterations,
   as O2 did for the outages. Bound: 1.06–1.2×.
4. **The OPF's evaluation and derivatives,** last. Bound: about 1.1× on mid-size AC if halved,
   1.04× on the largest.

The checklist's "fuller study of the KKT ordering and factorization" follows LDLᵀ, since it
compares against it, and "parallel and GPU" after that.
