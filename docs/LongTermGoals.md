# Long-term goals

Goals beyond the current roadmaps ([`PowerFlowChecklist.md`](PowerFlowChecklist.md), complete,
and [`OptimizationChecklist.md`](OptimizationChecklist.md), whose phases O1 to O5 are complete),
recorded 2026-11-04 with coarse plans, goal 5 on 2026-11-18, and goal 6 on 2026-11-22. Each becomes a checklist of its own, planned with the maintainer,
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

**Done (2026-11-22),** in four PRs: #101, #104 (which relanded #103), #105 and this survey's.
The option, `MehrotraOptions::line_search`, has `MonotoneFilter` and `AdaptiveFilter`, with
feasibility restoration. It stays `Off` by default: under the rule written before
[the survey](benchmarks/line-search-2026-11-22.md), neither filter qualified for the OPF's
presets.
- **What's left for later:**
  - **A warm-start default.** With `Off`, warm-started outages fail where cold ones converge
    (case9241pegase AC: none of 28), and either filter converges on them. That needs a rule of
    its own, written before a survey of warm starts.
  - **The free mode's own line search** and Nocedal, Wächter and Waltz's conditional
    corrector, which their implementations used and PR 3 left out by choice (2026-11-22).
  - **MIPS:** left as it is.

**Planned (2026-11-22)** in [`LineSearchChecklist.md`](LineSearchChecklist.md): the
adaptive-barrier filter of Nocedal, Wächter and Waltz (2009), with Wächter and Biegler's filter
line search and full feasibility restoration, for the predictor-corrector method only. That
answers the open questions below: how a line search fits Mehrotra's predictor-corrector (its
free and fixed modes), and MIPS (left as it is).

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

## 6. GPU-friendly KKT formulations for the interior-point solvers

Give `ybus_ipm` KKT formulations that a GPU's sparse direct solver can factor accurately, so
that the KKT factorization can run on the GPU. Recorded 2026-11-22, from phase O5's target 4
(#95), which found cuDSS fast but wrong on today's formulation.

- **Why it's needed:**
  - **Today's KKT matrix doesn't suit a GPU solver.** The solvers factor the indefinite,
    unregularized `[W Jᵀ; J 0]`. Its zero constraint block, near half the matrix, and its
    exact zeros on W's diagonal need pivots delayed across fronts.
    - **The CPU can:** `ybus_sparse`'s multifrontal LDLᵀ delays them.
    - **GPU solvers can't:** they fix the factorization's structure at analysis for speed, and
      replace small pivots instead. cuDSS replaced a fifth to a third of all pivots, and was
      wrong on every dumped KKT matrix
      ([`o5-cudss-spike-2026-11-22.md`](benchmarks/o5-cudss-spike-2026-11-22.md)).
  - **The speed is there.** On its own ordering, cuDSS's iteration was ×2.1–5.7 faster than the
    CPU LDLᵀ, nominally. None of our orderings kept that speed while fixing the accuracy.
  - **The formulation is what fixes it.** With the zero block regularized (`[W + δI, Jᵀ; J,
    −δI]`, δ ≥ 1e-4), cuDSS perturbed no pivot and got the inertia right. The published GPU
    OPF work (MadNLP with cuDSS) condenses the KKT system into a positive definite one, which
    needs no pivoting at all. It reports up to about 10× on large OPF cases, with robustness
    limits on harder problems.
- **Builds on:**
  - **The cuDSS bindings** (`ybus_cuda::GpuLdlt`, #95): `with_ordering`, `perturbed_pivots` and
    `factor_entries`. A positive definite (Cholesky) mode, cuDSS's `CUDSS_MTYPE_SPD`, would be
    added.
  - **The KKT benchmark** (`scripts/benchmark-cudss-kkt.sh`, the `cudss_kkt` example), and
    dumping KKT matrices from real solves.
  - **Phase O4's static-pivot tools,** kept in `ybus_sparse`: the static LDLᵀ with pivot
    replacement, Ruiz scaling, normwise refinement, and their record of how regularization
    fails when δ is too large or too small.
  - **The predictor-corrector method's inertia correction (#58),** whose regularization
    parameters a quasi-definite formulation would reuse.
- **Coarse plan:**
  - **Measure before building,** on dumped KKT matrices, on the CPU and on cuDSS:
    - **Regularization:** how large δ must be for a static factorization, and how much accuracy
      each δ costs the step.
    - **A condensed system:** its accuracy and speed.
  - **The formulations to weigh, from the literature:**
    - **A regularized, quasi-definite KKT matrix** (Vanderbei, 1995; Altman and Gondzio, 1999),
      factored by a static LDLᵀ: the smallest change to the solvers.
    - **LiftedKKT** (Shin, Anitescu and Pacaud, *Electr. Power Syst. Res.*, 2024): equality
      constraints relaxed into inequalities, so the KKT system condenses to a positive definite
      one.
    - **HyKKT** (Regev et al., *Optim. Methods Softw.*, 2023): an augmented-Lagrangian Schur
      complement solved by conjugate gradient, with an inner Cholesky.
    - **The two compared on GPUs** (Pacaud et al., *Math. Program. Comput.*, 2026). The README's
      References give all three in full.
  - **An option first,** off by default, on the CPU before the GPU, so that it's judged apart
    from the device.
  - **Surveyed under O3's bar:** speedups recorded against the CPU at 1 and 14 threads, the
    same day.
  - **As a default** only by a rule written before its survey. MIPS's default stays MATPOWER's,
    bit for bit.
- **Open questions:**
  - **The inertia correction in condensed space:** a Cholesky's failure takes the inertia's
    place, and what that means for the predictor-corrector method's correction isn't
    settled.
  - **What a relaxation costs:** whether LiftedKKT's relaxation, or a regularization's δ, keeps
    the agreement with MATPOWER that today's solvers reach.
  - **Conditioning near convergence:** condensed matrices grow ill-conditioned as the barrier
    parameter falls, the robustness limit the papers report.
  - **Where the GPU still loses:** batches of solves on one device, as `opf::batch` runs them,
    and small cases, whose transfers and launches may outweigh the factorization.
  - **WSL:** GPU runs can crash the machine, so surveys start small and commit first.

## Order

Goal 4's parser and goal 3's DC SCOPF can follow the optimization roadmap's phases O3 and O4,
which the AC SCOPF needs. Goal 1 is independent of them and can start with the C layer whenever
the API settles; goal 2 follows it. Goal 5 is independent of the others; the AC SCOPF would
benefit from it, since its problems are larger and more nonconvex than the OPF's. Goal 6 is
independent too, but it touches the same KKT systems as goal 5's line search and the inertia
correction, so the two are best planned with each other in view. The AC SCOPF's larger KKT
systems are where a GPU would gain the most.
