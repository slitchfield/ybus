# Line search checklist

The roadmap for long-term goal 5 in [`LongTermGoals.md`](LongTermGoals.md): a globalization for
`ybus_ipm`'s predictor-corrector method (`ybus_ipm::mehrotra`), so that it converges from poor
starting points and on nonconvex problems. Each unchecked item is sized to be one commit. The
decisions were made with the maintainer on 2026-11-22.

## Goals

1. **A filter line search for the predictor-corrector method,** as an option off by default:
   `MehrotraOptions::line_search`. `LineSearch::Off` stays today's method, bit for bit.
2. **Mehrotra's adaptive barrier kept where it works.** The method follows Nocedal, Wächter and
   Waltz's adaptive strategy:
   - **Free mode:** Mehrotra's μ, while a progress test holds.
   - **Fixed mode:** a monotone barrier globalized by Wächter and Biegler's filter line search,
     when it doesn't.
3. **Full feasibility restoration,** for when no step is acceptable, before anything is
   surveyed.
4. **A default only by a rule** written before its survey, under O3's bar.
5. **Out of scope:** MIPS. It keeps its own `step_control`, and its default stays MATPOWER's,
   bit for bit.

**Conventions carried over from the earlier roadmaps:**
- one commit per checklist item, on `dev/line-search-*` branches, with the PR's description as
  its report
- `YBUS_REQUIRE_CASES=1 ./check.sh` at each commit
- planted bugs, each tried and reverted
- `Off` and MIPS bit for bit: the existing tests, and the 72-case OPF snapshot in the default
  and `rayon` builds
- timings the same day

## Why

- **Every step is taken whole** today, at the fraction-to-the-boundary length, with no merit
  function or line search (the module docs).
- **The inertia correction (#58) can't pay off without one.** It turns a step towards descent,
  but nothing checks the step makes progress. On hs71 the method refuses the local maximum from
  2, then fails, as MIPS does. On the OPF the correction improved no outcome, so it's off by
  default.
- **On the OPF** (`benchmarks/ldlt-opf-2026-11-08-ac-mehrotra-ldlt-uncorrected.csv`, 73 AC
  cases):
  - **No outright failures:** none of the method's 24 failures is a case where MATPOWER
    converges. case13659pegase fails in MATPOWER too.
  - **Different solutions:** it converges to a different objective, where MIPS-LU matches
    MATPOWER's, on:

    | Case | Objective against MATPOWER's |
    |---|---|
    | case60nordic | 32% off |
    | case_ACTIVSg25k | 11% off |
    | case_ACTIVSg70k | 9.7% off |
    | case33mg | 2e-5 off |

    These are what a line search could change; the survey measures whether it does.

## Decided (2026-11-22)

- **The method:** the adaptive-barrier filter of Nocedal, Wächter and Waltz (2009). The
  alternatives weighed:
  - **Wächter and Biegler's filter at the iteration's target σμ:** its theory assumes μ fixed
    while it acts.
  - **A merit function with backtracking:** simpler, but less robust on nonconvex problems.
- **Feasibility restoration:** full, from the start.
- **Scope:** the predictor-corrector method only.
- **Either filter needs LDLᵀ,** whose inertia checks descent, and turns the inertia correction
  on. Asking for one with LU is a setup error.
- **The record grows:** `MehrotraIteration` gains the line search's fields and becomes
  `#[non_exhaustive]`. That breaks API, so the next release is 0.2.0.
- **No new `Status` variant,** since `Status` isn't `#[non_exhaustive]`. A failed restoration
  ends as `NumericalFailure`, with its cause in the record.

## Sources

Each was checked against its Crossref record on 2026-11-22. The method is implemented from
these papers. IPOPT's code (EPL-2.0) isn't consulted, though its options documentation may be,
for the defaults its papers name.

- A. Wächter and L. T. Biegler, "On the implementation of an interior-point filter line-search
  algorithm for large-scale nonlinear programming," *Math. Program.* 106(1):25–57, 2006.
  [doi:10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y)
  - **It gives:** the filter line search (§2.3), the second-order correction (§2.4), the
    monotone barrier update (§3.1), feasibility restoration (§3.3), and the inertia correction
    already in place (§3.1, Algorithm IC).
- A. Wächter and L. T. Biegler, "Line search filter methods for nonlinear programming:
  motivation and global convergence," *SIAM J. Optim.* 16(1):1–31, 2005.
  [doi:10.1137/S1052623403426556](https://doi.org/10.1137/S1052623403426556)
  - **It gives:** the convergence theory.
- J. Nocedal, A. Wächter and R. A. Waltz, "Adaptive barrier update strategies for nonlinear
  interior methods," *SIAM J. Optim.* 19(4):1674–1693, 2009.
  [doi:10.1137/060649513](https://doi.org/10.1137/060649513)
  - **It gives:** the free and fixed modes, and the probing oracle, which is Mehrotra's rule.

## The shape of the code

- **The option,** `MehrotraOptions::line_search: LineSearch`:
  - `Off`, the default and today's method
  - `MonotoneFilter`, Wächter and Biegler's Fiacco–McCormick barrier with the filter
  - `AdaptiveFilter`, the free and fixed modes
- **The new module,** `mehrotra/filter.rs` (with `mehrotra.rs` moving to `mehrotra/mod.rs`):
  - **The barrier problem's measures:** `θ = ‖[g; h + z]‖₁` over the nonlinear and linear rows,
    and `φ_μ = f − μ (Σ ln z + Σ ln s_l + Σ ln s_u)`.
  - **The filter:** its `(θ, φ)` pairs, with dominance and augmentation.
  - **The acceptance tests,** with the paper's defaults.
- **The direction:** today's `solve` closure. Its right-hand sides take a fixed μ in fixed mode,
  and Mehrotra's corrector in free mode.
- **Trial points:** each is evaluated with today's `evaluate`, and the accepted trial's
  evaluation is kept for the next iteration. So an iteration costs one evaluation more only when
  it backtracks.
- **The dual step:** its own fraction-to-the-boundary length (§2.2), with the κ_Σ safeguard on
  the bound multipliers. `Off` keeps today's rule.

## PR 1: the monotone filter line search

- [x] **This checklist,** and goal 5's entry in `LongTermGoals.md` pointing to it.
- [x] **The option, the record and the setup error:**
  - `LineSearch`, with `Off` as the default
  - `MehrotraIteration`'s new fields, and `#[non_exhaustive]`
  - the LU-with-filter error
  - `mehrotra.rs` moved to `mehrotra/mod.rs`, unchanged
  
  `Off` must be bit for bit.
  - **Done:** `θ` is recorded on every iteration, since it costs a norm and changes no iterate.
    The filter's other fields stay at their defaults until the filter fills them.
  - **Until the next two items,** `MonotoneFilter` only turns the inertia correction on.
  - **`Off` checked bit for bit:** the 72-case snapshot, in the default and `rayon` builds.
- [x] **The monotone barrier** (§3.1):
  - **The barrier:** a fixed `μ_j`, with `τ_j = max(τ_min, 1 − μ_j)`.
  - **Its update:** `μ_{j+1} = max(ε/10, min(κ_μ μ_j, μ_j^{θ_μ}))` once the barrier problem's
    error is at most `κ_ε μ_j`.
  - **The step:** the centered Newton step at `μ_j`.
  - **Done** (`mehrotra/filter.rs`), with the paper's defaults: `μ₀ = 0.1`, `κ_ε = 10`,
    `κ_μ = 0.2`, `θ_μ = 1.5`, `τ_min = 0.99`, and the multipliers' safeguard at `κ_Σ = 1e10`.
    - **The step lengths:** the primal and dual lengths are kept apart, with `λ` on the primal
      one, as the paper's iteration has them.
    - **An adaptation:** `μ`'s floor is `comptol / (10 · pairs)`, not the paper's `ε/10`.
      The method converges on MIPS's complementarity measure, which sums over every pair, so a
      floor per pair of `ε/10` could leave it unreachable on a large problem.
    - **A warm start** begins at its own average complementarity rather than `μ₀`.
  - **Found:** hs71 from 2, in fixed mode without a line search yet, ends at a fourth KKT
    point, `f = 27.146` at `x = (1, 5, 1.449, 3.449)`, with stationarity 2e-8. The tests
    named three KKT points until now. So PR 1's hs71 test asserts a verified KKT point that
    isn't the local maximum, not "the minimum": a line search guarantees progress to a KKT
    point, not to the global one.
- [x] **The filter line search** (§2.3–2.4):
  - the switching condition, the Armijo test, sufficient decrease, the filter and its
    augmentation
  - backtracking down to `α_min`
  - the second-order correction, up to `p_max` times
  
  Until PR 2, a solve with no acceptable step fails there.
  - **Done,** with the paper's constants (`filter.rs`):
    - **`θ_max` and `θ_min`:** `10⁴` and `10⁻⁴` times `max(1, θ₀)`.
    - **The filter resets** at each barrier update.
    - **The directional derivative** `∇φᵀd` is over `x` and the inequalities' slacks.
    - **The second-order correction** solves the same factorization with the constraints'
      residuals accumulated, and is judged at the full step's length, as §2.4 has it.
    - **The trial's evaluation** becomes the next iterate's, so an accepted first trial costs
      nothing extra.
  - **It also runs on a problem with no inequalities or bounds,** where `φ` is just `f`.
  - **The record gains the filter's size,** so that a skipped augmentation is observable.
  - **Measured on the test problems:**
    - **hs71 from 1:** where `Off` stops at the iteration limit (KKT residual 0.4), the filter
      reaches the optimum in 25 iterations, with one backtrack.
    - **Elsewhere,** the outcomes match `Off`'s.
- [ ] **Tests** (`ybus_ipm/tests/line_search_test.rs`):
  - **hs71 from its local-maximum start,** which must reach the minimum.
  - **The nonconvex problems:** `Nonlinear2` and `Nonlinear3`.
  - **A Maratos-effect problem,** on which full steps must be accepted near the solution.
  - **The filter's invariants.**
  - **Planted bugs:**
    - the Armijo test's sign
    - augmentation skipped
    - the second-order correction disabled
    - the switching condition inverted

## PR 2: feasibility restoration

- [ ] **The restoration problem** (§3.3), as a `Problem` wrapping the original:
  - **Its variables:** `x`, plus `p, n ≥ 0` for every row.
  - **Its objective:** `ρ Σ (p + n) + (ζ/2) ‖D_R (x − x_R)‖²`.
  - **Its constraints:** each row's, with `p` and `n` added.
  - **Its bounds:** the original ones.
  - **Its starting point:** the paper's closed form for `p` and `n`.
- [ ] **Solving it** with the same method in fixed mode, without nested restoration.
  - **Its stopping rule:** `θ ≤ κ_resto θ_R`, at a point acceptable to the original filter.
  - **After it returns:** the original iteration resumes, with its multipliers reset.
  - **If it fails:** the solve ends as `NumericalFailure`, recorded as local infeasibility.
- [ ] **Tests:**
  - a problem that must restore, then converge
  - a locally infeasible problem
  - **Planted bugs:** the original filter's exit test skipped, and `p` or `n` dropped from a row

## PR 3: the adaptive barrier

- [ ] **The rules, from the paper, written here before coding:**
  - the free mode's progress test
  - when to switch to fixed mode, and its μ
  - when to return
- [ ] **`LineSearch::AdaptiveFilter`:**
  - **The free mode:** today's predictor-corrector step, accepted by the free mode's test.
  - **The fixed mode:** PR 1's filter, with PR 2's restoration.
- [ ] **Tests:**
  - PR 1's problems
  - the small OPF cases (case9 to case300), within O3's bar of MIPS-LU
  - **Planted bugs:** the switch never taken, and the return never taken

## PR 4: the survey and the default

- [ ] **The survey's tooling:** `YBUS_BENCH_LINE_SEARCH` in `scripts/benchmark-opf.sh` and the
  `opf_survey` example.
- [ ] **The rule, before the survey:** a default changes only if every outcome is the same or
  better, O3's bar holds wherever today's default meets it, and the iteration and time cost
  stays within a limit stated with the rule. The objective-gap cases, the backtracks, the
  second-order corrections, the restorations and the mode switches are recorded either way.
- [ ] **The survey and the record:**
  - the AC and DC OPF on the 73 cases with costs
  - MIPS's test problems
  - a warm-started outage batch
  
  The record is `benchmarks/line-search-<date>.md`.
- [ ] **The docs:**
  - the module docs, the changelog (0.2.0) and the README's References
  - goal 5 closed in `LongTermGoals.md`
