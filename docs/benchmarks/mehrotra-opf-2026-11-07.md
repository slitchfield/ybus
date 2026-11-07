# The predictor-corrector solver on every MATPOWER case with costs, 2026-11-07

Phase O3's last item: the primal-dual method with Mehrotra's predictor-corrector
(`interior_point::mehrotra`, PRs #49–#51) is verified and benchmarked against phase O1's MIPS
and against MATPOWER, on the DC and AC optimal power flow over all 73 cases with costs.

## Setup

| | |
|---|---|
| Code | MIPS at `c3826ea` (`main`, with #51); the predictor-corrector solver at `b7a28ea`, this branch's safeguard (below). MIPS's code is the same at both. |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.99.1, release profile |
| Reference | MATPOWER 8.1's `rundcopf` (MIPS forced) and `runopf` (MIPS), run locally as for [`opf-2026-11-04.md`](opf-2026-11-04.md) |
| Runs | `scripts/benchmark-opf.sh`, with `YBUS_BENCH_MODEL` set to `dc` or `ac` and `YBUS_BENCH_SOLVER` to `mips` or `mehrotra`, one case per process with the default limits (1 hour and 10 GiB). Both solvers run at MATPOWER's OPF settings (`DcOpfOptions::mehrotra()`, `AcOpfOptions::mehrotra()`). |
| Results | [`…-dc-mips.csv`](mehrotra-opf-2026-11-07-dc-mips.csv), [`…-dc-mehrotra.csv`](mehrotra-opf-2026-11-07-dc-mehrotra.csv), [`…-ac-mips.csv`](mehrotra-opf-2026-11-07-ac-mips.csv), [`…-ac-mehrotra.csv`](mehrotra-opf-2026-11-07-ac-mehrotra.csv) |

Every case ran to an outcome within the limits. Times are the OPF's setup and solve together,
each the median of repeated runs (one run for the large cases).

**The bar** where MATPOWER's MIPS converges:

- converge there;
- the objective within 1e-6;
- `LAM_P` and `LAM_Q` within 1e-3, relative to each price's magnitude, at least 1;
- each KKT residual within 100× MIPS's own on the same case, above floors that are the
  solvers' tolerances in unscaled units.

A miss on a converged case counts against the solver only if it can't be explained.

## Results

| | DC | AC |
|---|--:|--:|
| Cases where MATPOWER's MIPS converges | 58 | 46 |
| … where the new solver converges | **58** | **46** |
| Converged, all 73 cases: MIPS / new solver | 59 / 64 | 46 / 49 |
| Objective, worst relative difference from MATPOWER | 6.6e-8 | 2.0e-5 (explained) |
| `LAM_P`, worst | 2.6e-8 | 8.5e-4 |
| `LAM_Q`, worst | | 9.4e-3 (explained) |
| Bar misses on converged cases, unexplained | 0 | 0 |
| Iterations where both converge, median (most): MIPS / new | 10 (148) / 4 (18) | 13.5 (82) / 11 (49) |
| Iterations, new / MIPS, median (range) | 0.43 (0.10–0.73) | 0.80 (0.44–1.21) |
| Time, new / MIPS, median (range) | 0.47 (0.11–0.78) | 0.79 (0.45–1.22) |
| Time per iteration, new / MIPS, median | 1.08 | 0.99 |
| Total time where both converge: MIPS / new | 14.0 s / 4.5 s | 58.5 s / 44.4 s |
| Peak memory, any case: MIPS / new | 413 / 373 MiB | 1,091 / 1,086 MiB |

### DC

- **Iterations:** the new solver takes fewer than MIPS on every one of the 59 cases where both
  converge, a median 0.43×. The ratio falls with size: 0.43× under 1,000 buses, 0.37× up to
  10,000, and 0.26× above.
- **No misses:** prices and objectives agree with MATPOWER's to 1e-7.
- **Five cases where only the new solver converges,** MIPS failing:
  - **case2383wp, case3012wp and case3375wp** (33, 17 and 17 iterations). MATPOWER's default
    DC solver, GLPK, converges there, and the new solver's objectives match GLPK's to 1e-12,
    2e-13 and 2e-14: the same optimum.
  - **case_ACTIVSg25k and case_ACTIVSg70k** (29 and 35 iterations, against MIPS's limit of 150,
    in 5.6 s against 29.6 s and 29.9 s against 130 s). MATPOWER has no default-solver objective
    for them. Both are verified KKT points, with the largest residual 4.7e-6.
- **case_SyntheticUSA** fails on both solvers (MIPS at its iteration limit, the new solver
  numerically after 57).
- **The second solve costs little:** time per iteration is 1.08× MIPS's.

### AC

- **Iterations:** fewer on 42 of the 46 cases where both converge, equal on 2, more on 2. The
  median ratio is 0.67× under 1,000 buses and 0.91× from 1,000 to 10,000.
- **Time per iteration** equals MIPS's (0.99×): the derivatives and the factorization dominate.
- **Three cases where only the new solver converges:**
  - **case_ACTIVSg25k:** 79 iterations, largest residual 0.066 $/hr unscaled.
  - **case_ACTIVSg70k:** 84 iterations, residual 0.098 $/hr.
  - **case60nordic:** 22 iterations, by the solver's own measures (all below 1e-6). Its
    absolute feasibility is 6.7e-4 p.u., large slacks normalizing it away, and MATPOWER fails
    there too, so it's converged but unverified.
- **Where both fail:** case_ACTIVSg10k, case13659pegase and case_SyntheticUSA. The new solver
  reaches its iteration limit on case13659pegase and case_SyntheticUSA, and fails
  numerically on case_ACTIVSg10k after 74 iterations.

## The safeguard: a damped corrector on nonlinear problems

**The first AC survey** used the full corrector, `r_z = −Z μ − ΔZ_aff Δμ_aff + σ μ e`, as
Mehrotra gives it. It failed on five cases where MIPS converges: case_ACTIVSg2000 (iteration
limit), case3012wp, case3120sp, case3375wp and case8387pegase.

**What their histories showed:** they share a cause.

- **From the first iteration,** the steps were blocked at lengths of about 1e-4, with σ near 1,
  since the affine predictor made no progress.
- **The affine step** far from the solution was enormous, and its full second-order term made
  the corrector's step worse.
- **case_ACTIVSg10k and case_SyntheticUSA** reached a step of 2e12 at the first iteration, past
  `max_stepsize`.

**Experiments,** each a temporary patch on the five failures and on cases that converged:

| Change | The five fixed | Notes |
|---|--:|---|
| Separate primal and dual step lengths | 1 | breaks case9241pegase |
| Fixed σ = 0.1, as MIPS | 2 | |
| Bounds as slack-carrying rows, as MIPS | 2 | breaks case300 |
| All three together | 1 | |
| No second-order term | 5 | more iterations |
| **The term damped by `α_p,aff α_d,aff`** | **5** | adopted |
| The better of the full and damped correctors | 4 | an extra solve; case8387pegase at its limit |

**Why not merit backtracking:** it can't help when the steps are already short, so it wasn't
added.

**What was adopted,** decided with the maintainer:

- **The rule:** with nonlinear constraints, the second-order term is damped by the affine step
  lengths.
- **Unchanged elsewhere:** linear and quadratic programs keep the full term, and the DC
  survey's outcomes, iterations and residuals are unchanged by it.
- **Its origin:** this form is the crate's own safeguard, chosen from these experiments, not a
  method from the literature.
- **Its effect on the AC survey:** before it, 41 of the 46 converged, at a median 0.67× MIPS's
  iterations. After it, all 46 converge, at 0.80×.

**A cost on one nonconvex test.** From (2, 2, 2, 2), Hock and Schittkowski's problem 71 now
converges to a constrained local maximum, where the full corrector reached the optimum and
MIPS reaches its iteration limit:

- **The point:** f = 131.26, with x₁ at its upper bound.
- **The evidence:** the Lagrangian's Hessian, reduced to the active constraints' tangent space,
  has eigenvalues −18.2 and −11.9.
- **Why:** with no inertia correction, an interior-point method can't tell the two apart.
- **Where it stands:** the test records it, and inertia is left to phase O4's LDLᵀ.

## The misses, explained

24 checks miss the bar on 16 converged AC cases. Each is explained, and none counts against the
solver. A rerun of both solvers at tolerances of 1e-9 separates two kinds:

**Tolerance-limited reference prices** (11 cases): case12da, case15da, case33bw, case33mg,
case38si, case51ga, case57, case300, case_ACTIVSg2000, case3012wp and case3375wp.

- **At 1e-9,** the two solvers' results come together: `LAM_Q` within 1e-6 to 4e-4, and
  objectives within 3e-9.
- **MIPS's own prices move by as much as the miss** between its default and tight runs on 10 of
  the 11 (6.5e-3 on case51ga, 1.1e-2 on case300). The new solver's move by at most 3.9e-3
  (case300).
- **On case3375wp it's the new solver's price** that moves (2.9e-3, against MIPS's 2.5e-4).
- **On the small radial feeders,** MIPS stops with complementarity at 4.8–8.6e-7, just under
  its 1e-6 tolerance, and the new solver at 1.6e-11 to 7.8e-10.
- **So these are tolerance limits,** mostly in the reference prices (MIPS's, and so
  MATPOWER's). That includes case33mg's objective miss of 2.0e-5: MIPS's default objective is
  itself 2.0e-5 from the tight optimum, and the new solver's 1.7e-7.
- **Where the tight run itself fails:** at 1e-9 both solvers end in numerical failure on
  case3012wp, and the new solver does on case3375wp. There the agreement is of their last
  iterates.

**Different local optima** (2 cases): the differences don't shrink at 1e-9.

- **case8387pegase:** the objectives differ by 1.5e-7 (1.7e-6 from MATPOWER's), and `LAM_Q` by
  5.8e-3. A generator MIPS holds at its reactive limit of −68.2 MVAr runs at 28.7 MVAr. As O1
  found, this optimum isn't unique: MATPOWER's own two implementations differ by 9.3e-7.
- **case9241pegase:** the new solver's optimum is worse, by 2.3e-6, with `LAM_Q` 9.4e-3 away,
  and a generator at 22.0 MVAr where MIPS holds it at its limit of 177.3.
- **Verified KKT points:** the new solver's tight runs converge on both. MIPS's converges on
  case9241pegase, and ends in numerical failure at 1e-9 on case8387pegase, keeping its
  difference.
- **The rule that covers them:** distinct KKT points of a nonconvex problem, PR 2's rule for
  local optima.

**Residuals over 100× MIPS's** (case2869pegase, case8387pegase and case9241pegase):

- **Stationarity (case2869pegase, case9241pegase):** 2.8e3 and 4.5e3 $/hr on a voltage angle,
  against MIPS's 4.6 and 6.8. It comes from the constraint Jacobians times the multipliers.
- **Feasibility (case8387pegase):** 8.1e-4 p.u.² in a flow limit, against 2.1e-6.
- **Why it's allowed:** the termination measures normalize the gradient by
  `1 + ‖multipliers‖` and feasibility by `1 + max(‖x‖, ‖z‖)`, both large on the pegase cases.
  So the new solver, stopping in fewer iterations, meets the same 1e-6 measures with larger
  absolute residuals than MIPS's tighter final iterate.
- **The floor check:** the two other cases flagged by ratio, case_ACTIVSg500 and
  case1354pegase, are within the floors.

## Separate bounds

The new solver keeps x strictly inside its bounds, with their own slacks and multipliers,
where MIPS gives every bound row its own slack.

- **The KKT matrix:** in both, the bounds add the same diagonal, so its size and fill are the
  same.
- **What separate bounds save** (forming the matrix from a shorter `Jh`, and vectors over the
  bounds) wasn't measured on its own. The time per iteration (0.99× on AC, 1.08× on DC)
  includes it, the second solve and everything else.
- **A cost they may carry:** a strictly feasible x caps the primal steps when the bounds are
  narrow, as voltage bands are. The experiments above show slack-carrying bounds alone weren't
  the cause of the AC failures, and weren't better overall.

## What's left for phase O4

- **The KKT solve:** each iteration factors once and solves twice, so a cheaper factorization
  (reusing the symbolic analysis, LDLᵀ) pays off in both solvers.
- **Inertia:** an LDLᵀ factorization's inertia would let the method correct an indefinite
  Hessian, and stop converging to maxima such as hs71's from (2, 2, 2, 2).
- **The failures both solvers share:** case_ACTIVSg10k, case13659pegase and case_SyntheticUSA,
  on the AC OPF, where MATPOWER fails too.
- **Warm starts** across related problems.
