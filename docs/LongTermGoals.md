# Long-term goals

Goals beyond the current roadmaps ([`PowerFlowChecklist.md`](PowerFlowChecklist.md), complete,
and [`OptimizationChecklist.md`](OptimizationChecklist.md), whose phases O1 to O4 are complete),
recorded 2026-11-04 with coarse plans, and goal 5 on 2026-11-18. Each becomes a checklist of its own, planned with the maintainer,
when its turn comes.

## 1. Stable foreign function interfaces: C, Python and Julia

Make `ybus_powerflow` (and the solvers beneath it) usable from other languages through stable
interfaces built with Rust's established FFI tooling. Planned with the maintainer on 2026-11-10
as a checklist of its own: [`FfiChecklist.md`](FfiChecklist.md).

- **C first,** since it's the common ground: a C ABI crate (for example `powerflow-ffi`) exposing
  opaque handles for a case, a model and a solution, functions that return status codes rather
  than panicking across the boundary, and caller-freed buffers. Headers generated with
  `cbindgen`; versioned symbols and an ABI policy written down before the first release.
- **Python** through PyO3, built and packaged with maturin: cases from files or arrays, the power
  flows and OPFs, results as NumPy arrays.
- **Julia** through the C ABI (`ccall`), wrapped in a small package, so it shares the C layer's
  stability rather than adding a second binding.
- **Open questions:** which API is stable enough to freeze (the case model and solution types
  first, the solver options later); error and panic handling at the boundary; thread safety of
  the handles; and how the bindings are tested in CI (each language's own test runner against
  the committed reference cases).

## 2. Benchmarks against published power flow libraries

Compare `ybus_powerflow` with established libraries in the languages goal 1 reaches, on the same
MATPOWER cases, by agreement and time, as the existing surveys compare it with MATPOWER.

- **Candidates to survey first** (to confirm when planned): pandapower and lightsim2grid
  (Python), PyPSA (Python), GridCal (Python), PowerModels.jl (Julia), and MATPOWER itself, already
  the reference.
- **Method:** the same cases and options where the libraries allow (Newton's method, tolerances,
  reactive limits, OPF formulation), one process per case under time and memory limits, as in
  `scripts/benchmark-*.sh`; agreement against MATPOWER's references; setup and solve timed
  separately, called both natively and through goal 1's bindings, so the FFI's own cost is
  visible.
- **Depends on** goal 1 for the in-language comparisons; native-to-native timings can start
  earlier. Licenses and citations of each library recorded, as for MATPOWER.

## 3. Security-constrained optimal power flow

Extend the OPF with contingency constraints: dispatch that stays within limits after each
credible outage.

- **Builds on** phase O1's OPF, the DC outage machinery (rank-k updates of B′'s factorization,
  PTDF and LODF, bridges), and the interior-point solvers of phases O3 and O4.
- **Coarse plan:**
  - **DC preventive SCOPF** first: post-contingency flow limits through LODF, as linear
    constraints on the base-case dispatch, with contingencies screened so that only binding
    ones enter.
  - **AC SCOPF** next, preventive then corrective (post-contingency recourse such as generator
    participation and voltage control), with contingency screening and decomposition so the
    problem doesn't grow with every outage at once.
  - Validated against published formulations and results; goal 4's data gives realistic
    instances.
- **Open questions:** the formulation to follow first (MATPOWER has no SCOPF; the GO
  Competition's Challenge 1 formulation is a natural reference), and how the KKT systems scale,
  which phase O4's work on the KKT solve feeds.

## 4. Parsing ARPA-E's SC-ACOPF case and contingency files

Read the input files of ARPA-E's Grid Optimization (GO) Competition Challenge 1, the SC-ACOPF
challenge, into the case model, so its public datasets can drive goal 3.

- **The format:** four ASCII files per problem instance
  ([GO Competition, Challenge 1 input files and format](https://gocompetition.energy.gov/challenges/challenge-1/input-data-format)):
  `case.raw` (the network, in a PSS/E-style raw format: buses, loads, fixed and switched shunts,
  generators, lines, transformers and areas), `case.rop` (generator cost curves), `case.inl`
  (generator participation factors) and `case.con` (contingencies). Only some fields of each are
  used; Appendix A of the Challenge 1 problem formulation specifies which.
- **Coarse plan:**
  - the `.raw` reader first, mapped onto the existing case model (with transformer and switched
    shunt fields the model may need to grow), checked by power flows on the datasets' base cases
  - then `.rop` into generator costs, `.con` into a contingency list (a new type), and `.inl`
    into participation factors
  - errors in the style of the MATPOWER reader, naming file and line
- **Out of scope at first:** the later challenges' formats (Challenge 3 uses a single JSON file,
  with multi-period and unit-commitment data), and the full PSS/E format beyond what Challenge 1
  uses.
- **Data licensing** to check before committing anything derived from the datasets, as for
  MATPOWER's cases.

## 5. A line search for the predictor-corrector method

Give `ybus_ipm::mehrotra` a globalization: a merit function or filter line search that accepts
or shortens each step by whether it makes progress, so that the method converges from poor
starting points and on nonconvex problems. Recorded 2026-11-18, from phase O4's records, where it
was the last open item of #59's "left for later".

- **Why it's needed:**
  - **Today the method takes every step** it computes, at the fraction-to-the-boundary length.
    It has no merit function or line search, as MIPS by default has none (MIPS's optional
    `step_control` is MATPOWER's own, and only for MIPS).
  - **The inertia correction (#58) can't pay off without one.** It makes a step a descent
    direction, but nothing then checks that the step makes progress. On hs71 it refuses the
    local maximum from 2, then fails, as MIPS does. It fails near the minima of a nonconvex
    problem with an equality. On the OPF it improved no outcome, so it's off by default.
  - **The remaining AC failures** include case13659pegase. It fails under both KKT orderings, at
    the iteration limit with AMD and by a numerical failure with matched (#82), which kept AC's
    default on AMD.
- **Builds on:**
  - the inertia from the pivoted LDLᵀ (#57, #58)
  - the predictor-corrector method's damped corrector (phase O3), whose own notes record that
    merit-function backtracking couldn't help where steps were already short
  - the KKT residuals the solvers report (`kkt_residuals`)
- **Coarse plan:**
  - **The method from the literature:** Wächter and Biegler's filter line search (IPOPT, 2006).
    It comes with its second-order correction and feasibility restoration, and it's the same
    paper whose inertia correction (Algorithm IC) #58 follows. A merit-function line search
    (ℓ1 or augmented Lagrangian) is the alternative to weigh when planned.
  - **An option first,** off by default, with the inertia correction turned on with it, since
    they're designed together.
  - **Measured** on MIPS's test problems (hs71 and the nonconvex QPs among them), then the
    full OPF survey.
  - **As a default** only by a rule written before its survey, as #59's and #82's were: O3's
    bar wherever today's default meets it, and the outcomes the same or better.
- **Open questions:**
  - **How a line search fits Mehrotra's predictor-corrector.** The literature's filter methods
    globalize a single Newton step per iteration; the corrector's combined step may need its own
    acceptance test.
  - **Whether MIPS should get it too,** beside its `step_control`. Its default must stay
    MATPOWER's, bit for bit.
  - **What it costs per iteration:** each trial step evaluates the problem again, which #79's
    plans made cheaper.

## Order

Goal 4's parser and goal 3's DC SCOPF can follow the optimization roadmap's phases O3 and O4,
which the AC SCOPF needs. Goal 1 is independent of them and can start with the C layer whenever
the API settles; goal 2 follows it. Goal 5 is independent of the others; the AC SCOPF would
benefit from it, since its problems are larger and more nonconvex than the OPF's.
