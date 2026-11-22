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
  
  Until PR 2, a solve with no acceptable step failed there.
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
- [x] **Tests** (`ybus_ipm/tests/line_search_test.rs`):
  - **hs71 from its local-maximum start,** which must reach the minimum.
  - **The nonconvex problems:** `Nonlinear2` and `Nonlinear3`.
  - **A Maratos-effect problem,** on which full steps must be accepted near the solution.
  - **The filter's invariants.**
  - **Planted bugs:**
    - the Armijo test's sign
    - augmentation skipped
    - the second-order correction disabled
    - the switching condition inverted
  - **Done.**
    - **hs71 from 1:** the filter reaches the optimum where `Off` stalls, and its non-Armijo
      steps grow the filter.
    - **hs71 from 2:** a verified KKT point that isn't the local maximum (the found fourth
      point).
    - **`Nonlinear2` and `Nonlinear3`:** their optima, each with KKT residuals below 1e-6.
    - **The Maratos problem** (Nocedal and Wright's Example 15.4), from the circle at angle
      0.1: a second-order correction is taken and no step is shortened.
  - **Each planted bug caught,** then reverted:
    - **the Armijo test's sign:** by two unit tests of the acceptance tests
    - **augmentation skipped:** by hs71 from 1, through the record's filter size
    - **the second-order correction disabled:** by the Maratos test
    - **the switching condition inverted:** by two unit tests, and hs71 from 1

## PR 2: feasibility restoration

- [x] **The restoration problem** (§3.3), as a `Problem` wrapping the original:
  - **Its variables:** `x`, plus `p, n ≥ 0` for every row.
  - **Its objective:** `ρ Σ (p + n) + (ζ/2) ‖D_R (x − x_R)‖²`.
  - **Its constraints:** each row's, with `p` and `n` added.
  - **Its bounds:** the original ones.
  - **Its starting point:** the paper's closed form for `p` and `n`.
  - **Done** (`mehrotra/restoration.rs`), in this solver's form.
    - **The variables:** `x`, then `p_g, n_g` per nonlinear equality, `p_h` per nonlinear
      inequality (`h − p_h ≤ 0`), and `p_a, n_a` per linear row (`l ≤ A x − p_a + n_a ≤ u`).
    - **Why one variable per inequality:** a one-sided row needs only its excess, since the
      method's own slack takes the rest.
    - **The Hessian:** the original's constraint curvature, through its `lagrangian_hessian`
      with `cost_mult = 0`, plus the proximity term's.
    - **The start:** each equality's and linear row's `p` and `n` by eq. 30, so that its row
      holds exactly, and each inequality's `p` as its excess.
    - **Unit tests:** the Jacobians and gradient against finite differences, the exact start,
      and a linear row's start on its limit.
- [x] **Solving it** with the same method in fixed mode, without nested restoration.
  - **Its stopping rule:** `θ ≤ κ_resto θ_R`, at a point acceptable to the original filter.
  - **After it returns:** the original iteration resumes, with its multipliers reset.
  - **If it fails:** the solve ends as `NumericalFailure`, recorded as local infeasibility.
  - **Done** (`mehrotra/mod.rs`), through an internal entry point with hooks: the starting
    barrier parameter, no nested restoration, and the stopping predicate. Restoration's solve
    runs as `dyn Problem`, so it is one instantiation.
    - **Entering:** the current point's margin goes into the filter first. The restoration
      problem starts at `μ_R = max(μ̄, ‖c‖∞)`, `ζ = √μ̄`, with what's left of `max_it`.
    - **Also entered when no step can be computed,** or the step is longer than
      `max_stepsize`, which ended the solve before (found 2026-11-22). On Wächter and
      Biegler's counterexample the step's multipliers blow up before the line search runs,
      where the linearized constraints meet the bounds; IPOPT likewise restores when its step
      computation fails.
    - **The predicate:** the restoration point's `x`, with each slack `max(−h, μ̄)`, must have
      `θ ≤ 0.9 θ_R` and be acceptable to the original filter.
    - **Resuming:** each slack's multiplier `μ̄ / slack`, then the safeguard; `λ` by least
      squares (§3.6) from `[I Jgᵀ; Jg 0]`, or zero if larger than `λ_max = 10³`.
    - **The record:** one entry per restoration phase, `mode: Restoration`, with its
      iterations in `restoration_iterations`. `Solution.iterations` counts them too.
- [x] **Tests:**
  - a problem that must restore, then converge
  - a locally infeasible problem
  - **Planted bugs:** the original filter's exit test skipped, and `p` or `n` dropped from a row
  - **Done** (`tests/line_search_test.rs`), on Wächter and Biegler's counterexample ("Failure of
    global convergence for a class of interior point methods for nonlinear programming,"
    *Math. Program.* 88(3):565–574, 2000, doi:10.1007/pl00011386, checked against Crossref
    2026-11-22): `min x₀` subject to `x₀² − x₁ − 1 = 0`, `x₀ − x₂ = 0.5` and `x₁, x₂ ≥ 0`.
    - **From `(−0.5, 1, 1)`:** `Off` fails, and the filter restores once, then converges to the
      solution `(1, 0, 0.5)`.
    - **From `(−2, 1, 1)`, the paper's start:** the filter restores three times. The last
      restoration ends at `(−1, 0, 0)`, so the solve ends as locally infeasible (found
      2026-11-22).
      - **Why that's right:** with `x₁, x₂ ≥ 0`, the ℓ1 violation along `x₁ = x₂ = 0` is
        `1.5 − x₀ − x₀²` on `(−1, 0.5)`, and the equality's violation is zero below `−1`, so
        `(−1, 0, 0)` is a strict local minimizer of it. The barrier between it and the solution
        peaks at `x₀ = −0.5`.
      - **Starts tried:** the method converges only from starts with `x₀ ≥ −0.5`, past that
        peak.
    - **A locally infeasible problem:** `x₀² + x₁² + 1 = 0` with `x₀ ≥ 0`, which ends in
      restoration mode near the origin.
    - **hs71, from the three starts PR 1 tests:** never restores, and its count is its records'.
    - **Planted bugs, each caught, then reverted:**
      - the margin not added on entering, caught by the counterexample's filter count
      - `n` dropped from the equality rows, caught by the unit tests of the restoration problem
      - the predicate's filter test skipped, and the least-squares `λ` set to zero or solved with
        the wrong sign, each caught only once the predicate (`filter::restored`) and the
        estimate (`restoration::least_squares_multipliers`) became functions with unit tests
    - **Why the end-to-end tests didn't catch them:**
      - **The filter test:** on these problems, every point with 10% less violation is also
        acceptable to the filter.
      - **The estimate:** after restoration, the bound multipliers `μ̄ / s`, at slacks near
        `10⁻⁶`, put entries of about `2·10⁴` into the gradient. The fitted `λ` then exceeds
        `λ_max`, so it's zero either way.
    - **Not tried:** the plan's "restoration's own multipliers in place of the estimate."
      Restoration's solution isn't kept past the predicate, so there was nothing to swap in;
      the zero and the wrong sign were planted instead.

## PR 3: the adaptive barrier

- [x] **The rules, from the paper, written here before coding:**
  - the free mode's progress test
  - when to switch to fixed mode, and its μ
  - when to return
  - **Written 2026-11-22,** from §5 ("Globalization Framework") and §6 of the paper's 2005
    report: J. Nocedal, A. Wächter and R. A. Waltz, "Adaptive barrier strategies for nonlinear
    interior methods," IBM Research Report RC23563, February 2005, read from the copy on
    Optimization Online (`optimization-online.org/wp-content/uploads/2005/03/1089.pdf`). The
    2009 *SIAM J. Optim.* version wasn't available to read, so these rules are the report's.
  - **The defaults the report leaves open** come from IPOPT's options documentation
    (`coin-or.org/Ipopt/OPTIONS.html`), each read verbatim on 2026-11-22:
    `adaptive_mu_kkterror_red_fact = 0.9999`, `adaptive_mu_kkterror_red_iters = 4`,
    `adaptive_mu_kkt_norm_type = 2-norm-squared`, `fixed_mu_oracle = average_compl`,
    `adaptive_mu_monotone_init_factor = 0.8`, and `adaptive_mu_restore_previous_iterate = no`.
  - **The measure,** the original problem's KKT error (5.1):
    `Φ = ‖∇L‖² + ‖c‖² + ‖ZXe‖²`. As in IPOPT (§6), each term is divided by its vector's
    length. Here that's
    `‖r_d‖²/nx + ‖[g; h + z]‖²/m + ‖[z∘μ; s_l∘π_l; s_u∘π_u]‖²/pairs`, with an empty term
    left out.
  - **Free mode's test:** a free iteration's new point is accepted if `Φ̃ ≤ κ M_k`, with
    `M_k = max{Φ_{k−l}, …, Φ_k}` and `l = min(k, l_max)`: the KKT error falls by `κ = 0.9999`
    against the largest of the last `l_max + 1 = 5` accepted values. The values are those of
    accepted free points and of the points where fixed mode returns.
  - **Free mode's step:** the report allows any rule for μ. Its own implementations also
    line-searched each free step, with a filter reset every free iteration (§5), and dropped
    harmful correctors (§7). Here free mode is today's predictor-corrector step, taken whole
    (decided 2026-11-22), so that wherever `Φ` keeps falling the method is `Off`'s with LDLᵀ and
    the inertia correction, bit for bit.
  - **The switch:** if the test fails, fixed mode starts **from the trial point** (the report's
    framework, and IPOPT's `restore_previous_iterate = no`).
    - **Its barrier parameter:** `μ̄ = 0.8 ×` the average complementarity there, no lower than
      fixed mode's floor.
    - **The filter** is reset, and `M_k` is frozen: nothing joins the window until fixed mode
      returns.
    - **A free step that can't be computed, or is longer than `max_stepsize`,** has no trial
      point, so fixed mode starts from the current point.
  - **Fixed mode:** PR 1's monotone barrier and filter, with PR 2's restoration.
  - **The return:** "as soon as `Φ(x, y, z) ≤ κ M_k`", checked at every fixed-mode iterate,
    including after a restoration phase, without first solving its barrier problem (§5). The
    returning point's `Φ` joins the reference values.
  - **Why it converges:** `Φ` falls by `κ` at least every `l_max` iterations, so every limit
    point is a first-order point (§5). Fixed mode reaches each such point by PR 1's and PR 2's
    globalization.
- [x] **`LineSearch::AdaptiveFilter`:**
  - **The free mode:** today's predictor-corrector step, accepted by the free mode's test.
  - **The fixed mode:** PR 1's filter, with PR 2's restoration.
  - **Done** (`mehrotra/mod.rs`, `mehrotra/filter.rs`): the rules above.
    - **One flag, set per iteration:** `fixed`, which replaced PR 1's `filter` in the
      iteration. It is always on for `MonotoneFilter`, and starts off for `AdaptiveFilter`.
    - **`filter::kkt_error`** computes `Φ`, and **`filter::Progress`** keeps the reference
      window. Both have unit tests.
    - **No freeze is needed:** only free points and returning points join the window, so `M_k`
      can't change while fixed mode runs. An explicit freeze was written and then removed, once
      its planted bug (the freeze skipped) proved it changed nothing.
    - **The record:** the new field `kkt_error`, which is `Φ` on an adaptive iteration.
    - **A stalled free step,** below `alpha_min` or with complementarity out of range, starts
      fixed mode instead of failing the solve. A NaN point still fails it.
    - **A free step that can't be computed** starts fixed mode from the current point, with a
      record of its own.
- [x] **Tests:**
  - PR 1's problems
  - the small OPF cases (case9 to case300), within O3's bar of MIPS-LU
  - **Planted bugs:** the switch never taken, and the return never taken
  - **Done,** 2026-11-22.
    - **The record replays the rules:** `replay_the_modes` (`tests/line_search_test.rs`) re-applies
      them from the recorded KKT errors. It checks:
      - every switch, with fixed mode's `μ̄` at 0.8 of the average complementarity
      - every return
      - that nothing else changes mode
      
      The starting record carries `Φ₀` for this.
    - **Free mode is `Off`, bit for bit,** with LDLᵀ and the inertia correction, wherever `Φ` keeps
      falling: hs71 from `(1, 5, 5, 1)`, Nonlinear2 and Nonlinear3 never switch, and they end at
      the same `x` in the same iterations.
    - **hs71 from 2:** `Off` fails, even with the inertia correction. The adaptive barrier switches
      and returns at least twice each, and converges to PR 1's local KKT point, `f = 27.146`.
    - **hs71 from 1:** 9 iterations (free, then three fixed, then five free), against `Off`'s 13
      and `MonotoneFilter`'s 25.
    - **The Maratos problem and the counterexample** from `(−0.5, 1, 1)`: both converge. On the
      counterexample, fixed mode's first iteration is already a restoration.
    - **The small AC OPF cases,** against MIPS-LU run in the same test
      (`mehrotra_opf_test.rs`): O3's bar on the objective and the KKT residuals, wherever MIPS
      converges.

      | Case | Iterations | MIPS's | Modes |
      |---|---|---|---|
      | case9 | 9 | 11 | one fixed iteration |
      | case14 | 9 | 11 | free throughout |
      | case24_ieee_rts | 13 | 13 | two fixed |
      | case30 | 11 | 15 | free throughout |
      | case39 | 14 | 16 | two fixed |
      | case57 | 9 | 10 | free throughout |
      | case89pegase | 16 | 25 | six fixed |
      | case118 | 11 | 14 | free throughout |
      | case300 | 15 | 19 | one fixed; `Off` takes 23 |

      - **Objectives:** within 1.6e-7 of MIPS's.
      - **case145:** MIPS doesn't converge, and neither does this solver with any line search.
      - **The prices are reported, not barred.** `Off` itself misses the 1e-3 price bar against
        MIPS on case57 (LAM_Q 1.7e-3, the same here) and case300 (7.3e-3, against 4.3e-3 here).
        That bar was measured on the committed cases only.
    - **`ybus_powerflow::opf` now re-exports `LineSearch` and `IterationMode`.** Before,
      `MehrotraOptions::line_search` couldn't be set without depending on `ybus_ipm` directly
      (found 2026-11-22; that gap shipped in 0.2.0).
    - **Planted bugs, each caught, then reverted:**
      - the switch never taken
      - the return never taken
      - fixed mode's `μ̄` from `μ₀`
      
      The fourth, the reference not frozen, changed nothing, so the freeze was removed (above).

## PR 4: the survey and the default

- [x] **The survey's tooling:** `YBUS_BENCH_LINE_SEARCH` in `scripts/benchmark-opf.sh` and the
  `opf_survey` example.
  - **Done** (2026-11-22).
    - **`opf_survey --line-search off|monotone|adaptive`** adds eight columns: the line search,
      then counts from the method's record (the backtracks, the steps with a second-order
      correction, the restoration phases and their iterations, the switches to fixed mode, the
      returns, and the fixed-mode iterations).
    - **`opf_batch --line-search`,** with a `line_search` column, and
      `YBUS_BATCH_LINE_SEARCH` in `scripts/benchmark-opf-batch.sh`.
    - **`mehrotra_test.rs`'s `every_line_search_meets_the_bar_on_mipss_problems`** runs MIPS's
      test problems with each line search on LDLᵀ, under `check`'s bar without its pins
      (`meets_the_bar`). Every one passes, and the infeasible LP fails under each.
    - **`--line-search off` matches the unchanged example** on every result column: 36 runs,
      on six cases, both models, with MIPS-LU and the predictor-corrector method with
      `ldlt-uncorrected` and `ldlt`.
- [x] **The rule, before the survey:** a default changes only if every outcome is the same or
  better, O3's bar holds wherever today's default meets it, and the iteration and time cost
  stays within a limit stated with the rule. The objective-gap cases, the backtracks, the
  second-order corrections, the restorations and the mode switches are recorded either way.
  - **Written 2026-11-22, before any survey run.**
  - **What it decides:** the OPF presets' line search, `AcOpfOptions::mehrotra()` and
    `DcOpfOptions::mehrotra()`, each judged on its own model's 73 cases with costs.
    - **`MehrotraOptions::default()` keeps `Off`:** it factors by LU, which a filter can't use.
    - **MIPS,** the OPF's default solver, isn't in question.
  - **The runs compared:**
    - **Today's preset, the baseline:** the predictor-corrector method, LDLᵀ without the
      inertia correction, `Off`.
    - **Each candidate:** `monotone` or `adaptive`, on LDLᵀ with the inertia correction,
      which a filter turns on.
    - **MIPS-LU:** run beside them for O3's residual bar.
    - **An error, a time-out or a memory limit** counts as not converged.
  - **A candidate replaces `Off` in a preset only if all of these hold:**
    1. **Outcomes no worse:**
       - every case the baseline converges on, the candidate converges on
       - where both converge and MATPOWER does, the candidate's objective is no further from
         MATPOWER's than the baseline's: its relative difference is at most the baseline's, or
         at most 1e-6, O3's bar, so that a difference under the bar doesn't count
    2. **O3's bar kept:** wherever the baseline meets it, the candidate meets it.
       - The objective within 1e-6 of MATPOWER's.
       - The prices within 1e-3.
       - Each KKT residual within 100× MIPS-LU's, or below 1e-6.
    3. **The cost within 1.10× (decided 2026-11-22),** over the cases both converge on:
       - the candidate's total iterations at most 1.10× the baseline's
       - the median, over those cases, of the candidate's solve time over the baseline's, at
         most 1.10, with the runs made the same day, one after another
  - **If both candidates qualify,** the one with the fewer total iterations wins.
  - **If neither qualifies,** the presets keep `Off`, and the record says which conditions
    failed, on which cases.
  - **Recorded either way:**
    - **The objective-gap cases,** where the baseline converges to a different objective than
      MATPOWER's: case60nordic, case_ACTIVSg25k, case_ACTIVSg70k and case33mg.
    - **The cases a candidate converges on that the baseline doesn't.**
    - **Each candidate's line-search counts:** the backtracks, second-order corrections,
      restorations and mode switches.
- [ ] **The survey and the record:**
  - the AC and DC OPF on the 73 cases with costs
  - MIPS's test problems
  - a warm-started outage batch
  
  The record is `benchmarks/line-search-<date>.md`.
- [ ] **The docs:**
  - the module docs, the changelog (0.2.0) and the README's References
  - goal 5 closed in `LongTermGoals.md`
