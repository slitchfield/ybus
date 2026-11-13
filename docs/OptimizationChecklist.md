# Optimization checklist

The roadmap after [`PowerFlowChecklist.md`](PowerFlowChecklist.md), which is complete: optimal power
flow, interior-point solvers, and acceleration of the power flow and of those solvers. Each
unchecked item in phase O1 is sized to be one commit; the later phases are coarser and are broken
into commit-sized items when each is planned, as Phase 5 was. Items marked **Decide** need an
answer before the work that depends on them.

## Goals and decisions

Set with the maintainer on 2026-11-03:

1. **OPF following MIPS** (phase O1). MATPOWER's optimal power flow, solved by a Rust
   reimplementation of MIPS, MATPOWER's primal-dual interior-point solver, followed step for step
   from its BSD-licensed `mips.m`, as `newtonpf` and `fdpf` were for the power flow. It aims at
   MATPOWER's iteration counts, with the OPF formulation (costs, power balance, limits) around it.
2. **Accelerating the power flow** (phase O2). Profile first, then take the CPU's gains
   (parallelism within one solve, SIMD, memory layout), then a CUDA spike behind a feature,
   measured against the CPU path. The CPU path stays the default unless the GPU wins.
3. **A naive interior-point solver** (phase O3). A separate, simpler primal-dual method of our own
   design, written from the literature rather than from MIPS, benchmarked and verified against
   phase O1's MIPS and against MATPOWER.
4. **Accelerating the interior-point solvers** (phase O4). From profiles, with the KKT solve the
   expected hot spot.

**The order** is O1, then O2 (independent of O1, so it can overlap), then O3, then O4.

**Conventions carried over from the power flow work:**
- one commit per checklist item, on `dev/` branches, with an uncommitted report as each PR's
  description and figures taken from git and the test runner
- MATPOWER's BSD-licensed source may be consulted for the OPF formulation and MIPS, and the
  README names the files consulted; `ybus_sparse` and the phase O3 solver are written from
  published descriptions only
- MATPOWER's results committed only for case14, case30 and case118 (public-source cases), and
  generated locally for the rest, one case per process under an address-space limit
  (`ulimit -v`), since an unlimited MATPOWER job once exhausted the machine's memory
- every citation checked online before it's used, and a DOI left out if it can't be confirmed

## What there is to build on

Facts checked when this roadmap was written (2026-11-03):

- **Costs:** `Case` already parses `gencost` (`powerflow/src/case.rs`), polynomial and
  piecewise-linear, and 73 of the 78 MATPOWER cases carry it.
- **MATPOWER's OPF runs locally:** under Octave 8.4.0, `runopf` and `rundcopf` solve case9 with
  MIPS (AC cost 5296.69, DC cost 5216.03). MIPS and GLPK are the only OPF solvers installed;
  IPOPT, Knitro, Gurobi, CPLEX, MOSEK, HiGHS and the MATLAB optimizers aren't.
- **MIPS** (`mips/lib/mips.m`, 778 lines) is a general primal-dual interior-point solver for
  nonlinear programs:
  - It minimizes `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.
  - Its defaults: `feastol`, `gradtol`, `comptol` and `costtol` 1e-6, `max_it` 150, `sigma` 0.1,
    `z0` 1, `xi` 0.99995, and step control off.
  - Each iteration solves the reduced KKT system `[M dg; dgᵀ 0]` with `mplinsolve`, whose default
    is the backslash operator.
  - Its reference is Wang, Murillo-Sánchez, Zimmerman and Thomas (2007), below.
- **The OPF formulation** is about 2,400 lines across MATPOWER's `lib/opf_*.m`: power balance,
  branch flow limits (on apparent power, `opf.flow_lim = 'S'`, by default), voltage limits, costs,
  and their first and second derivatives.
- **MATPOWER's own OPF tests:** `lib/t/` holds 26 `t_opf_*` tests, with stored case9 solutions
  (`soln9_opf*.mat` and `soln9_dcopf.mat`), and `mips/lib/t` tests MIPS on its own problems.
- **External benchmarks:** PGLib-OPF publishes baseline objective values for its networks; their
  licenses need checking before anything derived from them is committed.
- **Hardware:** an NVIDIA RTX 4070 Ti SUPER (16 GB, compute capability 8.9), visible through WSL's
  `libcuda`. The CUDA toolkit isn't installed, so cuSPARSE and cuSOLVER aren't available until it
  is.
- **In `ybus_sparse`:** LU with refactoring and a given ordering, sparse right-hand sides,
  transposed solves, iterative refinement, condition estimates and parallel batch solves, all of
  which the KKT systems of an interior-point method can use.

## Phase O1: OPF following MIPS

- [x] **Decide:** where the solver lives, which MATPOWER code is the reference, which constraints
      come first, and the acceptance bar.
      Decided (2026-11-03):
      - **The solver** goes in a new workspace crate, `interior_point`: MIPS now, the phase O3
        solver later, behind one problem interface. It depends on `ybus_sparse`; `powerflow`,
        which holds the OPF formulation, depends on it.
      - **The reference is MATPOWER's default code.** Measured on case9, case14, case30, case118,
        case300 and case2383wp, AC and DC: MP-Core and the legacy code take the same iteration
        counts, with objectives within 6e-14 relative and voltages, dispatch and prices within
        3.5e-9. The legacy rule applies only where they differ, so the survey rechecks every
        case.
      - **DC references force MIPS** (`opf.dc.solver = 'MIPS'`): MATPOWER's default hands a
        linear program (costs whose quadratic terms are zero) to GLPK. Parity includes MIPS's
        failures: forced onto MIPS, case2383wp's DC OPF doesn't converge in 150 iterations,
        though GLPK solves it. The default solver's objective is kept as a cross-check.
      - **The constraints:** everything the cases use. Of the 78, 73 have costs, 42 flow limits,
        16 angle-difference limits, 8 dispatchable loads (so `makeAvl`'s constant power factor
        is in, changed from the roadmap, which deferred dispatchable loads), 2 piecewise-linear
        costs and 2 reactive-power costs. None has reactive capability curves, costs above
        quadratic or user constraints, so those are refused with an error until needed.
      - **The acceptance bar:** MATPOWER's iteration counts on the committed cases, and every
        result within 1e-6, MIPS's stopping tolerance. Then the survey records how close the
        results actually come, and the bar is tightened to what's measured, as the power flow's
        1e-9 was.
      - **Three PRs,** below.

### PR 1: the `interior_point` crate and generator costs

- [x] The `interior_point` crate and its problem interface: an objective with gradient and
      Hessian, nonlinear constraints with Jacobians, the Lagrangian's Hessian, linear constraints
      and bounds; `MipsOptions` with MIPS's defaults; and a result with the multipliers.
      `interior_point::problem` has the `Problem` trait (with nonlinear constraints and the
      Lagrangian's Hessian optional), `Constraints` (linear constraints and bounds, each
      optional) and `QuadraticProblem` (`½ xᵀ H x + cᵀ x`, as `qps_mips` builds). Constraint
      Jacobians take one row per constraint, transposed from MIPS's convention.
      `interior_point::mips` has `MipsOptions` (MIPS's defaults, including step control's), the
      exit status, the multipliers as MIPS reports them, and each iteration's record.
- [x] MIPS, ported from `mips.m`, and `qps_mips`'s quadratic-program wrapper, tested on the
      problems of `t_mips.m` and `t_qps_mips.m` with MATPOWER's tolerances and iteration counts.
      - **`interior_point::mips::mips` and `qps_mips`,** following `mips.m` step by step: the
        bounds and linear constraints split into equalities and one- or two-sided inequalities,
        the starting slacks and multipliers, the reduced KKT system solved with `ybus_sparse`'s
        LU, ordered by AMD and by COLAMD since PR 2 (where MATLAB uses its backslash; a singular
        system counts as MIPS's numerical failure, as MATLAB's NaN does), step control, the step
        lengths, the barrier
        update, the termination and failure tests, and the multipliers' filtering and mapping
        back to the bounds and linear constraints.
      - **References:** `scripts/generate-mips-reference.sh` records MATPOWER's MIPS 1.5.2 on the
        problems of `t_mips.m` and `t_qps_mips.m`, plus two starts of Hock and Schittkowski's
        problem 71 where step control changes the run and where MIPS reaches its iteration
        limit, each with step control off and on (`interior_point/tests/data/mips_reference.txt`).
      - **Agreement:** every run takes MATPOWER's exit status and iteration count (including the
        infeasible LP's numerical failure after 6, and 150 at the limit), and its solution,
        objective and multipliers agree within 6.5e-14 relative (1.0e-13 since the switch to
        COLAMD in PR 2; the bar is 1e-12), except a
        failed run's diverged multipliers (about 3e24), within 3.3e-7. Each problem also meets
        `t_mips`'s known answers to its precision.
- [x] Generator costs in `powerflow`: polynomial (with `opf_setup`'s handling) and
      piecewise-linear (with `makeAy`'s constraints), reactive-power cost rows, and the total
      cost.
      - **`powerflow::opf::cost`:** `GenCost::evaluate` and `total_cost`, as `totcost` and
        `polycost` compute them (piecewise linear costs extended along their end segments);
        `OpfCosts`, for the in-service generators: one-segment piecewise linear costs made
        linear and zero leading coefficients dropped, as `opf_setup` does, then the per-unit
        quadratic terms (`Qpg`, `cpg`, `kpg`, and the reactive `Qqg`, `cqg`, `kqg`) and
        `makeAy`'s segment constraints, real power costs first. A DC OPF leaves out the
        reactive costs. Polynomials above quadratic, and piecewise linear costs with fewer
        than two points or two at the same output, are refused (`UnsupportedCost`).
      - **Tests:** hand-computed fixtures, and the OPF form against the direct values at each
        case's dispatch, within 1e-12, on case14, case30, case118, case30pwl, case_RTS_GMLC,
        case9Q and case30Q. One cost is non-convex: case_RTS_GMLC's `gencost` row 73
        (0-based), from rounded data. Its helper variable (the largest of its segment lines)
        comes out 9.1e-5 $/hr above the interpolated cost, as it would in MATPOWER's
        formulation.

### PR 2: OPF references and the DC OPF

- [x] OPF reference results: `dcopf` (MIPS forced) and `opf` runs in `generate_reference.m`,
      committed for case14, case30 and case118 and generated locally for every case with costs.
      - **Runs:** `rundcopf` with `opf.dc.solver = 'MIPS'`, and `runopf` with its defaults
        (MIPS). Each records bus, generator and branch results with their prices and
        multipliers, a summary (objective, MIPS's iterations, success, whether Octave's
        least-squares fallback was taken, and for DC the default solver's objective and name), and MIPS's per-iteration history. A run where MIPS fails
        is still recorded.
      - **Committed:** case14, case30 and case118, 30 files, 152 KB (33 KB compressed).
        Regenerating them gives byte-for-byte the same files. All six runs converge: DC in 10, 9
        and 11 iterations, AC in 11, 15 and 14.
      - **Provenance:** case14's and case118's costs were added by MATPOWER's authors, with no
        stated source (the IEEE format has none). case30's come from Ferrero et al. The
        reference README says so.
      - **Running under `ulimit -v`:** Octave hung at no CPU in `rundcopf` under an
        address-space limit until BLAS was made single-threaded. `generate-reference.sh` now
        sets `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` unless they're already set.
- [x] The DC OPF: `[Va, Pg, y]` with the power balance, flow limits, angle limits and the
      piecewise-linear rows, solved through the quadratic-program wrapper, against `rundcopf`.
      - **`powerflow::opf::dc`:** `dc_opf` follows the legacy `opf_setup`'s DC model, after
        `ext2int`:
        - every bus of type REF has its angle fixed; there's no `bustypes` and no island check,
          as in MATPOWER's OPF
        - the constraint blocks `Pmis`, `Pf`, `ang` (`makeAang`, including its literal 0° limits)
          and `ycon`, in that order
        - `dcopf_solver`'s starting point, and MATPOWER's MIPS options (`feastol` 5e-6)
        - the constant cost terms added after the solve, as `opt_model.solve` adds them
        - the results mapped back as `dcopf_solver` and `opf` map them, with out-of-service
          elements zeroed
      - **Shared with the DC power flow:** the `ext2int` half of the topology (`Network`) and
        `makeBdc` (`Bdc`), now crate-private pieces of their own.
      - **Fixtures**, worked out by hand: a binding flow limit, a binding angle limit,
        piecewise-linear costs, a dispatchable load, and two islands with a reference each.
      - **Committed references:** MATPOWER's iteration counts on case14, case30 and case118,
        with every quantity within 6.9e-12 (the bar is 1e-9) and objectives within 1.7e-15.
      - **Every local case with costs (73):** `local_cases_match_rundcopf`, run with
        `YBUS_DCOPF_CASES`. Its criteria follow what was found:
        - **The outcome** always matches: MATPOWER converges on 58 cases and fails on 15, and so
          does `powerflow`, except on case3120sp (below).
        - **Iteration counts** match on 56 of the 57 runs MATPOWER finishes without Octave's
          fallback (below). The exception is case2736sp, at 44 to MATPOWER's 45, since MIPS
          switched to COLAMD (`KKT_ORDERING`).
        - **Objectives and prices** agree wherever both converge: objectives within 1e-10
          relative, prices and multipliers within 1e-8, and the optimum within 1e-6 of the
          default solver's (GLPK on 52 cases).
        - **Dispatch, angles and flows aren't compared locally.** On degenerate problems
          (linear costs, several generators at the same marginal cost) the optimum is a set,
          and the point found differs, by up to 61.5 MW on case8387pegase, while objectives
          agree within 5.5e-14 and prices within 3e-11. MATPOWER's own MP-Core and legacy
          implementations differ from each other the same way: 3.84 MW on case1354pegase, with
          objectives equal within 4e-16. The committed cases' optima are unique, and are
          compared in full.
      - **Octave's least-squares fallback:** when a KKT system is singular to machine precision,
        Octave's backslash warns and returns a least-squares solution by QR, not the LU
        solution. MIPS then continues where an exact solve would stop it on too long a step.
        UMFPACK's own LU, called directly, gives our step (case16am's third iteration: a norm of
        1.7e14 from both, against the fallback's 0.33). The references now record when this
        happened (a `singular` column): on 6 cases, 5 of which MATPOWER fails anyway. There the
        iteration counts aren't compared. On case3120sp `powerflow` converges at 148 iterations
        where MATPOWER fails at 150.

### PR 3: the AC OPF

Split from the survey (now PR 4) with the maintainer on 2026-11-04. The references stay MP-Core's,
with legacy results generated locally to explain differences.

- [x] The AC OPF: `[Va, Vm, Pg, Qg, y]` with the power balance and flow limits (and their first
      and second derivatives), angle limits, dispatchable loads and the piecewise-linear rows,
      from `nlpopf_solver`'s starting point, against `runopf`.
      - **Derivatives** (`opf::derivatives`, public, as MATPOWER's are), in Technical Note 2's
        complex matrix notation. Every quantity is `[C V] conj(Y V)`, so one formula gives the
        first derivatives, and one Hessian of the bilinear form `Vᵀ M conj(V)` gives the second
        derivatives of the bus injections and both ends' flows. The squared apparent power adds
        `2 Re(dSᴴ [μ] dS)`. Tested (`tests/opf_derivatives_test.rs`) against central finite
        differences at 1e-8 relative (a planted sign error fails them), for symmetry, and against
        the Newton Jacobian at 1e-12.
      - **`opf::ac`:** `ac_opf` follows `opf_setup`'s AC model for MIPS (`legacy_formulation`
        0):
        - REF angles fixed; `Pmis` and `Qmis`; `Sf` and `St` as squared apparent power on the
          limited branches; `makeAvl`, `makeAang` and `makeAy` rows, in that order
        - the cost's constant terms inside MIPS's objective, as `nlp_costfcn` puts them
        - MIPS's options for the AC OPF (`feastol` 5e-6, `cost_mult` 1e-4) and
          `nlpopf_solver`'s start
        - the results and multipliers mapped back as `nlpopf_solver` and `opf` map them,
          including `MU_SF = 2 μ RATE_A/baseMVA`
        - PQ capability curves, and dispatchable loads `makeAvl` rejects, are refused
          (`UnsupportedGenerator`)
      - **Shared with the DC OPF and the power flow:** `makeYbus`'s stamping (`Admittances`),
        and the OPF's reference buses, starting point and linear rows.
      - **Fixtures:** every solution must balance power at each bus from its reported flows,
        respect every limit, satisfy complementarity (multiplier × slack ≈ 0), price power at
        each generator's marginal cost, and cost what `totcost` says. Cases cover a binding
        apparent-power limit, binding voltage limits, reactive costs, a dispatchable load at its
        power factor, and the refusals.
      - **Committed references:** MATPOWER's iteration counts (11, 15 and 14), every quantity
        within 5.4e-9 (case30's voltage-limit multipliers; the bar is 1e-8), and objectives
        within 1.1e-15.
      - **Every local case with costs (73)** (`YBUS_ACOPF_CASES`). Every case finishes, the
        largest in 678 s (`case_SyntheticUSA`):
        - **Outcome:** matches on all 73 (MATPOWER converges on 46, fails on 27).
        - **Iteration counts:** match on all 46 converged runs but case8387pegase, at 82 to 80.
          Its optimum isn't unique: MATPOWER's legacy code differs from MP-Core there by as much
          (246 MVAr of reactive dispatch, objective 9.3e-7).
        - **Objectives and bus prices:** objectives within 1.1e-11 relative, `LAM_P` within 2.4e-7
          and `LAM_Q` within 9.4e-5 (the bars: 1e-10, 1e-6 and 2e-4).
        - **Individual bound multipliers aren't compared locally.** Where near-zero-impedance
          branches lock neighbouring buses at a voltage limit, the split of the limit's price
          between their multipliers is ill-determined. On case3012wp `mu_vmax` differs by up to
          0.43, at neighbouring buses whose differences cancel (99% of their magnitude), and
          `LAM_Q` shifts by 1.16e-4 over the region, with the objective within 1.3e-12. MATPOWER's
          two implementations differ less there (5e-7 in `LAM_Q`), since they share Octave's
          backslash for the KKT solve.

### PR 4: the OPF survey

- [x] The OPF survey over every case with costs, DC and AC: an `opf_survey` example and
      `scripts/benchmark-opf.sh`, rechecking MP-Core against the legacy code, recorded in
      `docs/benchmarks/`, with the bars set from it.
      - **Recorded** in [`docs/benchmarks/opf-2026-11-04.md`](benchmarks/opf-2026-11-04.md). The
        survey reproduces the local tests:
        - **Outcomes:** the same as MATPOWER's on 72 of 73 DC cases and all 73 AC.
        - **Iteration counts:** MATPOWER's on 56 of 57 and 45 of 46 converged runs.
        - **Agreement:** objectives within 7.5e-12 and 1.1e-11, and `LAM_P` within 4.5e-11 and
          2.4e-7, outside the explained cases. The bars stay as PR 2 and PR 3 set them from the
          same measurements.
      - **MP-Core against legacy:** the two agree on every outcome and on every converged run's
        iteration count, with objectives within 3.9e-13 (DC) and 6.8e-12 (AC, outside
        case8387pegase). They differ on the degenerate dispatch (86 MW on case8387pegase, DC) as
        `powerflow` does, and fail after different numbers of iterations on 7 AC runs that fail
        in both. So the MP-Core references stand, and no case needs the legacy rule.
      - **Time against MATPOWER under Octave:** far faster on small cases (Octave's fixed costs),
        and close on large ones (0.64–1.27 times MATPOWER's time for AC). On the three largest
        cases, which fail in both, 1.4–2.9 times, since each iteration analyzes and factors the KKT
        system afresh: a baseline for phase O4.
      - **Generation:** `generate_reference.m` gains `dcopf_legacy` and `opf_legacy` runs and
        local timing files (not committed, since the time varies), so the recheck is
        reproducible.

**Phase O1 is complete** (PRs #38–#40 and this one, 2026-11-04).

## Phase O2: accelerating the power flow

- [x] **Decide** what to accelerate first, from profiles of the existing surveys: Newton,
      fast-decoupled, DC, N-1 outages and reactive limits. Record where the time goes.
      - **Recorded** in [`docs/benchmarks/profile-2026-11-04.md`](benchmarks/profile-2026-11-04.md):
        samply profiles (`scripts/profile-pf.sh`, the `profile_pf` example, the `profiling`
        build profile) of seven workloads on case2383wp, case9241pegase, case_ACTIVSg70k and
        case_SyntheticUSA, folded into stages.
      - **Decided with the maintainer (2026-11-04):** results stay bit-identical on the CPU path;
        samply for the CPU, Nsight Systems for the CUDA spike.
      - **The targets, in order:**
        1. **N-1 outages without per-outage allocation:** allocating, zeroing and copying each
           outage's result vectors take 35–55% of the N-1 time, against 7–9% for the sparse
           solves (bound 1.5–2.2×).
        2. **The LU refactorization kernel:** 41–64% of one-at-a-time reactive limits and 25% of
           a large Newton solve (halving it bounds the gain at 1.3–1.5× for reactive limits).
        3. **Parsing:** 25–48% of a single DC or fast-decoupled run.
      - **Smaller follow-ups:** the sorting in `solve_round`, allocation in the Newton loop, and
        reusing AMD's ordering.
- [ ] **The CPU.** Parallelism within one solve (refills, mismatches and flows with rayon), SIMD
      in the hot kernels, cache-friendlier layouts, and level-scheduled or supernodal triangular
      solves if the profiles point there. Each measured against the current records, with the
      results unchanged or within an agreed bar.
      - [x] **N-1 outages without per-outage allocation** (target 1;
        [`dc-contingency-allocation-2026-11-04.md`](benchmarks/dc-contingency-allocation-2026-11-04.md)):
        - **Changes:** reachability from a bridge forest built once per model, instead of an
          adjacency rebuilt and searched per outage; `ybus_sparse`'s `solve_sparse_into` with a
          reusable `SparseSolveWorkspace`; per-thread outage workspaces; and `finish` computing
          flows into the result and `Bbus θ` only at the reference buses.
        - **Results:** bit for bit the same. All outages of all cases take 4.5× less time
          serially and 4.6× in parallel; on case_SyntheticUSA, 157 s serially, from 680 s. The
          reachability check gave most of it: its breadth-first search, inlined into the outage,
          hid from the profile's allocation share.
        - **What's left:** each outage's dense updates and back substitutions, which a streaming
          API returning only the flows a caller wants could cut (deferred).
      - [x] **The LU refactorization kernel** (target 2;
        [`lu-refactor-2026-11-05.md`](benchmarks/lu-refactor-2026-11-05.md)):
        - **The kernel is latency-bound:** 0.7–1.4 ns per multiply-subtract, on indirect
          gathers and scatters.
        - **Fusing nested column pairs** (85–98% of the work on large Jacobians: each bus's angle
          and magnitude columns) halves the scatter traffic with the factors bit for bit
          unchanged: refactorization 1.27× faster on case_ACTIVSg70k and 1.19× on
          case_SyntheticUSA; one-at-a-time reactive limits there 1.10× (565 → 514 s), Newton
          1.04–1.06×. Its buffers are also kept between calls.
        - **Set aside with the maintainer (2026-11-05):** level scheduling (at most 1.4× even
          ideally on case_ACTIVSg70k, since most work is a serial chain of separator columns),
          unchecked indexing (about 10%, needs `unsafe`), parallelism within large columns, and
          low-rank updates instead of refactoring for reactive limits (not bit-identical).
      - [x] **Parsing** (target 3; [`parsing-2026-11-05.md`](benchmarks/parsing-2026-11-05.md)):
        - **Byte-level scans:** the reader scans lines and tokenizes matrices by bytes (every
          byte it stops on is ASCII), with a vectorizable precheck that rules out the lines
          holding nothing to stop on, and sizes rows like the first. Float parsing is still
          `str::parse`.
        - **Results:** the same results for 18,749 inputs, malformed ones included; reading a
          case 1.54× faster over all cases (median 1.51× from 20 KB up), DC runs 1.17–1.30×.
        - **Left:** float parsing is now its largest part (about 28%). A fused single pass or
          our own float parser were not tried.
- [x] **A CUDA spike** behind a `cuda` feature (through `cudarc`). Batched DC and N-1 solves and
      the Jacobian refill on the GPU, with transfers counted, against the CPU path. The CUDA
      toolkit gets installed only when cuSPARSE or cuSOLVER is needed. Record the outcome.
      - **Recorded** in [`cuda-spike-2026-11-07.md`](benchmarks/cuda-spike-2026-11-07.md).
      - **Built:**
        - **#46:** the `ybus_cuda` crate, the workspace's only `unsafe` code, and SpSM
          solves against CPU factorizations.
        - **#47:** `DcStudy::gpu` (batched solves, single-branch outages solved whole or
          screened with `outage_loadings`) and `DcSensitivity::gpu` (PTDF columns and
          transfers).
      - **Results:**
        - **N-1 screening:** 1.3–2.1× faster than 28 threads from 25,000 buses up, and even
          around 9,000–14,000.
        - **Whole results:** batched DC, outage solutions and PTDF columns run at 0.1–0.3× of
          28 threads from 1,354 buses up, since their per-result host work runs on one thread.
        - **The solves alone:** they beat 28 threads only on the two largest Jacobians
          (1.8–1.9×).
        - **Agreement:** to within 9.7e-13, with identical errors.
      - **Decided with the maintainer (2026-11-07):**
        - The CPU path stays the default, and `cuda` stays opt-in.
        - **The Jacobian refill is descoped,** with Newton on the GPU (cuDSS): one right-hand
          side per factorization, a latency-bound refactorization, and a new dependency. A GPU
          factorization comes back in phase O4's GPU KKT item.
        - The whole-result workloads are recorded as they are.
      - **Left:**
        - parallel host work around the device's solves;
        - rank-k updates for sets of branches, and LODF columns on the device;
        - alternatives to SpSM.

**Phase O2 is complete** (PRs #42–#47 and this one, 2026-11-07).

## Phase O3: a naive interior-point solver

- [x] **Decide** the method, for example a textbook primal-dual interior-point method, possibly
      with Mehrotra's predictor-corrector, and its references (Nocedal and Wright; Mehrotra
      1992). It takes the same problem interface as the MIPS crate.
      - **Decided with the maintainer (2026-11-07):**
        - **Method:** a primal-dual method with Mehrotra's predictor-corrector
          (`interior_point::mehrotra`).
          - **Predictor:** an affine step (σ = 0) sets the centering σ = (μ_aff/μ)³.
          - **Corrector:** adds the predictor's second-order term.
          - **Per iteration:** one KKT factorization serves both solves.
          - **Steps:** at most `xi` of the way to the boundary. A linear program takes
            separate primal and dual lengths. With curvature in the objective, both take the
            shorter one, as in Nocedal and Wright's §16.6: there the dual residual depends on
            x, and separate lengths left a small QP with two-sided bounds oscillating. With
            nonlinear constraints, they always share one length (decided 2026-11-07, for PR 2),
            since the Jacobians couple x to the multipliers too. MIPS takes separate lengths
            there.
          - **What MIPS does instead:** it fixes σ at 0.1 and solves once per iteration.
          - **References:** Mehrotra (1992); Nocedal and Wright, ch. 14 (Algorithm 14.3), 16.6
            and 19; S. J. Wright (1997). Torres and Quintana (1998) used a predictor-corrector
            method on the AC OPF.
        - **Interface:** MIPS's `Problem` and `Constraints`, unchanged, so the comparison is
          fair. Any extension (warm starts, fixed sparsity patterns) is additive, with
          defaults, and left for phase O4.
        - **Bounds:** `xmin ≤ x ≤ xmax` are kept apart from the linear rows, with their own
          multipliers. x stays strictly inside them, and they add only diagonal terms to the
          KKT matrix. MIPS stacks them as identity rows of its inequality Jacobian, which comes
          to the same diagonal. So the KKT matrix's size and fill don't change; what's saved is
          forming it, and the vectors over the bounds. A variable with `xmin = xmax` stays an
          equality row, as in MIPS.
        - **Globalization:** none at first, like MIPS. A merit function with backtracking is
          added only if AC cases fail where MIPS converges, decided from those failures.
        - **Termination:** MIPS's four measures, defined the same way (with the bounds' slacks
          and multipliers counted as MIPS counts its bound rows), so iteration counts compare
          fairly.
        - **The bar:**
          - **Must:** converge wherever MIPS converges, with the objective within 1e-6
            (relative) of MIPS's and MATPOWER's; KKT residuals below the tolerances; prices
            matching where they're unique, as in phase O1.
          - **Reported, not required:** iterations and time.
          - **Later:** the bar tightens after measuring.
        - **PRs:**
          1. this decision, and the solver for problems without nonlinear constraints;
          2. nonlinear constraints;
          3. the OPF's choice of solver, with MIPS staying the default;
          4. the survey and the record.
- [x] **The solver,** written from the literature, tested on the MIPS test problems.
      - **`interior_point::mehrotra`** (`mehrotra`, `qps_mehrotra`), for MIPS's whole problem
        class. It reports MIPS's multipliers and is checked by the solver-independent
        `interior_point::solution::kkt_residuals`.
      - **On MIPS's test problems:** every problem converges where MATPOWER's MIPS does.
        - **The bar:** the objective within 1e-6 of MIPS's, x and the multipliers within 1e-4,
          t_mips's known answers, and KKT residuals below 1e-5.
        - **Iterations against MIPS:** fewer on qp2 (8 to 11), lp3 (10 to 11) and nonlinear3
          (6 to 8); more on qp4 (13 to 10), nonlinear2 (12 to 8) and hs71 (11 to 8); equal on
          the rest.
        - **Pinned:** the counts are fixed in the tests, so any change to the method shows.
      - **On nonconvex starts** (decided 2026-11-07), the bar is a verified KKT point at a known
        optimum wherever MIPS converges.
        - **From (4.5, 4.5, 4.5, 4.5):** hs71 reaches the optimum (f = 17.014, 13 iterations),
          where MIPS stops at a local point (35.28).
        - **From (2, 2, 2, 2):** it reaches the optimum in 13, where MIPS reaches its iteration
          limit.
- [x] **Verified and benchmarked** against phase O1's MIPS and against MATPOWER, on DC and AC OPF
      over every case with costs: objective, KKT residuals, iterations and time. Recorded.
      - **Recorded** in [`mehrotra-opf-2026-11-07.md`](benchmarks/mehrotra-opf-2026-11-07.md):
        both solvers, all 73 cases with costs, DC and AC.
      - **The bar:** the new solver converges wherever MATPOWER's MIPS does (58 DC, 46 AC
        cases). It has no unexplained misses: the 24 AC checks over the bar are
        tolerance-limited reference prices, two non-unique optima (case8387pegase and
        case9241pegase), and residuals that the termination measures' normalization allows.
      - **Beyond MATPOWER:** it also converges on 5 DC cases where MIPS fails (matching GLPK's
        optimum to 1e-12 where GLPK has one), and on 3 AC cases.
      - **Iterations against MIPS,** median: 0.43× on DC, 0.80× on AC.
      - **Time against MIPS,** median: 0.47× on DC, 0.79× on AC. The second solve per
        factorization costs little (time per iteration 1.08× on DC, 0.99× on AC).
      - **Decided with the maintainer (2026-11-07):** on problems with nonlinear constraints,
        the corrector's second-order term is damped by the affine step lengths.
        - **Why:** the full term failed 5 AC cases where MIPS converges, with steps blocked
          from the first iteration. Merit backtracking couldn't help, since the steps were
          already short.
        - **Its origin:** this form is the crate's own safeguard, chosen from the experiments
          in the record.
        - **Its cost:** hs71 from (2, 2, 2, 2) now converges to a constrained local maximum.
          With no inertia correction, a maximum can't be told from a minimum; inertia is
          phase O4's.
      - **Left for O4:** the KKT solve (factored once, solved twice), inertia through LDLᵀ, the
        AC cases both solvers fail (case_ACTIVSg10k, case13659pegase, case_SyntheticUSA), and
        warm starts.

**Phase O3 is complete** (PRs #49–#51 and this one, 2026-11-07).

## Phase O4: accelerating the interior-point solvers

- [x] **Decide** what to accelerate, from profiles: the KKT solve is the expected hot spot.
      - **Recorded** in [`ipm-profile-2026-11-07.md`](benchmarks/ipm-profile-2026-11-07.md): samply
        profiles of both solvers on the DC and AC OPF, three DC and five AC cases, folded into
        stages, with the KKT matrices' size and fill.
      - **What it found:**
        - **The numeric LU factorization of the KKT matrix dominates,** at 38–80%, growing with
          size, with fill of 5.5–9.1×.
        - **The analysis (COLAMD and the symbolic LU) comes second:** 27–35% on mid-size DC and
          5–15% on AC. It's redone every iteration, though the pattern is fixed.
        - **Then allocation (5–18%) and the OPF's evaluation** (7–19% on AC). The second solve
          costs 0.6–2.2%.
        - **Both solvers share these hot spots.**
      - **Decided with the maintainer (2026-11-07):**
        - **The rule on results:** MIPS's default path stays bit for bit the same, keeping its
          MATPOWER parity. Arithmetic changes come as options, or as the predictor-corrector
          solver's default, judged by O3's bar.
        - **The targets, in order:**
          1. analyze the KKT matrix once per solve, bit for bit the same (up to 1.5× on
             mid-size DC);
          2. a symmetric indefinite LDLᵀ with Bunch–Kaufman pivoting, with its inertia (up to
             1.7× on the largest cases if it halves the factorization);
          3. allocation in the KKT assembly and the loops;
          4. the OPF's evaluation.
      - **Re-profiled after targets 1 and 2** in
        [`ipm-profile-2026-11-08.md`](benchmarks/ipm-profile-2026-11-08.md). The new cases are a
        small one, DC case300, and DC on LDLᵀ.
        - **On AC, on the predictor-corrector method with LDLᵀ:**
          - the OPF's evaluation is now the largest stage, at 28–38%;
          - the factorization takes 30–42%, and its arithmetic, the dense kernel, only 8–19%;
          - its bookkeeping takes 18–31%: front assembly and extend-add, a symmetry check and
            supernode partition redone every factorization, and allocation.
        - **On DC with LDLᵀ:** allocation is 18–27% and the bookkeeping 24–31%. The fronts are
          tiny (a median of 4 rows, 0.9 per variable). That is the cost per front behind its
          small-case loss.
        - **MIPS on LU** is 49–89% factorization.
        - **The order, decided with the maintainer (2026-11-08):**
          1. the LDLᵀ's bookkeeping (bound DC 1.3–1.5×, AC 1.2–1.3× from the bookkeeping
             alone), in two PRs:
             - [x] **The per-pattern work kept with the analysis, and the buffers reused,**
               bit for bit the same
               ([`ldlt-bookkeeping-2026-11-08.md`](benchmarks/ldlt-bookkeeping-2026-11-08.md)).
               - **The change:** `SymbolicIndefinite` holds the supernodes, the gather map and
                 the mirror map.
               - **Speed:** 1.16–1.32× on DC with LDLᵀ, 1.20–1.25× on mid-size AC, and
                 1.06–1.10× on the largest AC cases.
             - [x] **Relaxed amalgamation of the tiny fronts**
               ([`ldlt-amalgamation-2026-11-08.md`](benchmarks/ldlt-amalgamation-2026-11-08.md)).
               - **The change:** postorder the ordering, then merge by CHOLMOD's rule.
               - **The bar:** O3's, kept on DC and AC wherever LU meets it, with the same
                 iterations on every converged case.
               - **Against LU:** DC on LDLᵀ is now 1.15× faster in the median (0.79× before #61)
                 and 3.14× in total; AC 1.16× and 3.07×.
               - **Amalgamation's own gain:** 1.02–1.23×. The dense kernel pays per entry for
                 the larger fronts' zeros.
               - **Decided with the maintainer (2026-11-08):** the DC OPF's
                 predictor-corrector default switches to LDLᵀ, which it now qualifies for under
                 #59's rule. MIPS and `MehrotraOptions::default()` keep LU.
          2. the OPF's evaluation (bound 1.4–1.6× on AC);
          3. allocation and KKT assembly;
          4. the kernel and the KKT ordering.
- [ ] **The KKT solve.** Reusing the symbolic analysis across iterations while the pattern stays
      fixed; a symmetric indefinite LDLᵀ factorization in `ybus_sparse`, with Bunch–Kaufman-style
      pivoting, written from the literature and compared with the LU; and reduced against
      augmented KKT systems.
      - [x] **The analysis kept while the pattern holds** (target 1;
        [`kkt-analysis-2026-11-07.md`](benchmarks/kkt-analysis-2026-11-07.md)):
        - **The change:** both solvers share `KktAnalysis`, which analyzes the KKT matrix once
          and reanalyzes only when `SymbolicLu::factor` reports a different pattern.
        - **Bit for bit the same:** a snapshot of 72 cases, and every result column of the full
          OPF survey, match.
        - **Reuse held everywhere:** all 292 survey solves analyzed once.
        - **Faster, median:** DC 2.15× (MIPS) and 1.66× (the predictor-corrector method), AC
          1.38× and 1.33×. The largest AC cases gain 1.05–1.16×, where the factorization
          dominates.
        - **More than the profile's bound on DC,** whose innermost-frame rule had counted the
          analysis's own allocation as allocation.
      - [x] **A symmetric indefinite LDLᵀ** (target 2), with its inertia
        ([`ldlt-opf-2026-11-08.md`](benchmarks/ldlt-opf-2026-11-08.md)).
        - **Decided with the maintainer (2026-11-07):**
          - **The method:** static pivots (1×1, no pivoting during the factorization), on the
            Cholesky's symbolic analysis, rather than Bunch–Kaufman's dynamic 1×1 and 2×2 pivots.
          - **Why static pivots work:** the solvers regularize the KKT matrix into a
            quasi-definite one, `[W + δ_w I, Jgᵀ; Jg, −δ_c I]`, which factors under any
            symmetric ordering (Vanderbei, 1995).
          - **The inertia** comes from D's signs.
          - **Bunch–Kaufman** stays a possible follow-up if the numbers call for it.
        - **Three PRs:**
          - [x] **The factorization in `ybus_sparse`:** `SymbolicCholesky::factor_ldlt` and
            `CsrMatrix::ldlt` give `Ldlt`, with its solves and `inertia()`. It's Davis's
            up-looking LDLᵀ on the Cholesky's exact pattern.
            - **Tested:** against the Cholesky on SPD matrices, against a known D (L₀ D₀ L₀ᵀ),
              on real and complex quasi-definite matrices under every ordering, and on a matrix
              that needs pivoting.
            - **Accuracy:** it degrades as 1/δ with the regularization (about 3e-10 at
              δ = 1e-8).
          - [x] **The tools from the static-pivot attempt in the solvers,** kept in `ybus_sparse`:
            `CsrMatrix::symmetric_scaling` (Ruiz), `SymbolicCholesky::factor_ldlt_replacing`
            (tiny pivots replaced with their block's sign), and `Ldlt::solve_refined` and
            `Ldlt::solve_refined_normwise`.
            - **Why the static pivots were set aside (2026-11-07):** in the solvers they failed
              where LU converges.
              - **With δ = 1e-8 on the scaled matrix,** a test problem (nonlinear2) failed: the
                scaling's factors, about 1e-7, made δ about 1e6 on the original matrix, far
                above W's small eigenvalues.
              - **With no static δ** (scaling and pivot replacement only), the test suites
                passed. But AC case9241pegase failed on both solvers, and DC case_ACTIVSg25k
                on the predictor-corrector method. About 9,700 of W's zero-diagonal pivots were
                replaced per factorization, D's entries reached 7e12, and refinement couldn't
                recover.
              - **Where it worked,** it was fast: the same result 2.2× faster on AC
                case_ACTIVSg25k, on the predictor-corrector method.
            - **Refinement is normwise:** the componentwise backward error is 1.0 on every OPF
              KKT system, even LU's. A row with a zero right-hand side and almost nothing in
              |K||x| sets it.
        - **Decided with the maintainer (2026-11-07): a pivoted factorization instead.** A
          multifrontal LDLᵀ on the Cholesky's elimination tree, with threshold Bunch–Kaufman
          1×1 and 2×2 pivots in each dense front. A pivot too small to take is delayed to the
          parent front (Duff and Reid, 1983; MA57). It's written from the literature, for real
          and Hermitian complex matrices.
          - [x] **The dense kernel:** threshold Bunch–Kaufman on a front, with a public dense
            factorization to test it (`DenseLdlt`).
            - **No fallback to Bunch and Kaufman's own rule:** with u ≤ 1/2, a fully summed
              nonsingular matrix always has a threshold pivot (Duff and Reid's argument, in
              `DenseLdlt`'s docs).
          - [x] **The multifrontal driver:** `SymbolicCholesky::factor_indefinite` and
            `CsrMatrix::ldlt_indefinite` give `IndefiniteLdlt`.
            - **What it has:** fundamental-supernode fronts, extend-add, delayed pivots, the
              inertia, solves, refinement, and the analysis reused.
            - **On the OPF's KKT matrices** (DC and AC case9241pegase and case_ACTIVSg25k,
              iterations 1 and 10, u = 1e-4):
              - **Fill and speed:** 2.4–4.9× fewer entries than LU, and 2–5.3× faster on all
                but DC case9241pegase at iteration 10 (0.83×), with no regularization.
              - **Accuracy:** refined to LU's residual.
              - **Inertia:** correct, except where W isn't positive definite on the null space
                (AC case_ACTIVSg25k, off by 1–2).
            - **What u = 0.01 (MA57's default) costs:** more delays, 53k on DC
              case9241pegase at iteration 10, and so 1.6× LU's fill. u = 1e-8 loses accuracy.
            - **The rest of the delays:** about 8,400–48,000 at any threshold, from
              zero-diagonal constraint rows that AMD makes leaves without a partner. An ordering
              that keeps each with a variable it touches would remove them (the KKT ordering
              study below).
          - [x] **The solvers' option** (`KktFactorization::Ldlt` on `MipsOptions` and
            `MehrotraOptions`, `OpfSolver::with_kkt_factorization`, `opf_survey --kkt`), on the
            pivoted factorization at u = 1e-4, with up to 3 steps of normwise refinement.
            - **The predictor-corrector method corrects the inertia** (Wächter and Biegler,
              2006, Algorithm IC, with their defaults).
            - **The default is unchanged:** LU everywhere, bit for bit (a 72-case snapshot).
            - **On the test problems:** both solvers end as they do with LU.
            - **On the large cases:** both solvers end as with LU on DC and AC case9241pegase
              and case_ACTIVSg25k, mostly in the same iterations, 1.2–2.3× faster per solve.
            - **The correction:**
              - **Where it helps:** on a nonconvex QP it reaches the minimum where LU converges
                to a saddle.
              - **Where it does nothing:** AC case_ACTIVSg25k on the predictor-corrector method
                needs none, its W carrying the bounds' barrier terms.
              - **Where it costs:** AC case9241pegase takes 37 corrections, 58 iterations
                against 38, which costs LDLᵀ's speed there.
              - **What it can't do:** without a line search, it refuses hs71 from 2's local
                maximum but then fails, as MIPS does, and it fails near the minima of a
                nonconvex problem with an equality. A globalization is a separate question.
          - [x] **The survey, the record and the defaults**
            ([`ldlt-opf-2026-11-08.md`](benchmarks/ldlt-opf-2026-11-08.md)): ten runs on all 73
            cases.
            - **The default path is unchanged:** the LU runs match #54's survey on every
              result column.
            - **The bar holds:** no LDLᵀ setting loses O3's bar where LU meets it.
            - **Speed:** the large AC cases are 2.4–3.0× faster on the predictor-corrector
              method, and the small DC cases slower (0.55× in the median below 1,000 buses), from
              the factorization's cost per front.
            - **Decided with the maintainer (2026-11-08), by the rule set before the survey,
              per model:**
              - **AC:** `AcOpfOptions::mehrotra()` factors by LDLᵀ, without the correction:
                1.12× faster in the median, 2.70× in total, with the same outcomes and
                iterations.
              - **DC keeps LU:** 0.79× in the median, though 2.05× in total.
              - **MIPS keeps LU,** and so does `MehrotraOptions::default()`.
              - **The inertia correction** is now off by default. It improved no OPF outcome,
                and cost iterations on 16 AC cases.
            - **Left for later:**
              - the kernel's cost per front, which would let DC qualify;
              - a KKT-aware ordering for the structural delays;
              - a globalization for the predictor-corrector method.
- [ ] **A fuller study of the KKT ordering and factorization,** once phase O3's solver exists to
      compare against. MIPS's LU switched from AMD to COLAMD in O1's PR 2: with zero diagonals
      in the equality block, AMD's factors filled 14 million entries on `case_ACTIVSg10k`
      (24 s per factorization) against COLAMD's 640,000 (0.036 s). COLAMD lost one iteration-count
      match (case2736sp's DC OPF, 44 to MATPOWER's 45). The study would cover AMD, COLAMD,
      orderings that keep the KKT structure (constraints after the variables they touch), LDLᵀ
      and UMFPACK-like strategies, by fill, time and agreement with MATPOWER.
- [ ] **Parallel and GPU.** Batched OPF solves with rayon, and a GPU KKT factorization if phase
      O2's spike supports one; warm starts across related problems, such as outage sequences.
      Recorded.

## References

Each was checked online (Crossref) on 2026-11-03.

- H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman and R. J. Thomas, "On computational issues of
  market-based optimal power flow," *IEEE Trans. Power Syst.* 22(3):1185–1193, 2007.
  [doi:10.1109/TPWRS.2007.901301](https://doi.org/10.1109/TPWRS.2007.901301). MIPS's own
  reference.
- J. Nocedal and S. J. Wright, *Numerical Optimization*, Springer, 2006.
  [doi:10.1007/978-0-387-40065-5](https://doi.org/10.1007/978-0-387-40065-5)
- S. Mehrotra, "On the implementation of a primal-dual interior point method," *SIAM J. Optim.*
  2(4):575–601, 1992. [doi:10.1137/0802028](https://doi.org/10.1137/0802028)
- T. A. Davis, "Algorithm 849: a concise sparse Cholesky factorization package," *ACM Trans.
  Math. Softw.* 31(4):587–591, 2005.
  [doi:10.1145/1114268.1114277](https://doi.org/10.1145/1114268.1114277) (checked on Crossref on
  2026-11-07)
- R. J. Vanderbei, "Symmetric quasidefinite matrices," *SIAM J. Optim.* 5(1):100–113, 1995.
  [doi:10.1137/0805005](https://doi.org/10.1137/0805005) (checked on Crossref on 2026-11-07)
- P. A. Knight, D. Ruiz and B. Uçar, "A symmetry preserving algorithm for matrix scaling,"
  *SIAM J. Matrix Anal. Appl.* 35(3):931–955, 2014.
  [doi:10.1137/110825753](https://doi.org/10.1137/110825753) (checked on Crossref on 2026-11-07)
- A. Wächter and L. T. Biegler, "On the implementation of an interior-point filter line-search
  algorithm for large-scale nonlinear programming," *Math. Program.* 106(1):25–57, 2006.
  [doi:10.1007/s10107-004-0559-y](https://doi.org/10.1007/s10107-004-0559-y) (checked on
  Crossref on 2026-11-07)
- A. Altman and J. Gondzio, "Regularized symmetric indefinite systems in interior point methods
  for linear and quadratic optimization," *Optim. Methods Softw.* 11(1–4):275–302, 1999.
  [doi:10.1080/10556789908805754](https://doi.org/10.1080/10556789908805754) (checked on
  Crossref on 2026-11-07)
- S. J. Wright, *Primal-Dual Interior-Point Methods*, SIAM, 1997.
  [doi:10.1137/1.9781611971453](https://doi.org/10.1137/1.9781611971453) (checked on Crossref
  on 2026-11-07)
- G. L. Torres and V. H. Quintana, "An interior-point method for nonlinear optimal power flow
  using voltage rectangular coordinates," *IEEE Trans. Power Syst.* 13(4):1211–1218, 1998.
  [doi:10.1109/59.736231](https://doi.org/10.1109/59.736231) (checked on Crossref on
  2026-11-07)
