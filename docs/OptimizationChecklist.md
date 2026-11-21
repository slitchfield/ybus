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

- **Costs:** `Case` already parses `gencost` (`ybus_powerflow/src/case.rs`), polynomial and
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
      - **The solver** goes in a new workspace crate, `ybus_ipm`: MIPS now, the phase O3
        solver later, behind one problem interface. It depends on `ybus_sparse`; `ybus_powerflow`,
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

### PR 1: the `ybus_ipm` crate and generator costs

- [x] The `ybus_ipm` crate and its problem interface: an objective with gradient and
      Hessian, nonlinear constraints with Jacobians, the Lagrangian's Hessian, linear constraints
      and bounds; `MipsOptions` with MIPS's defaults; and a result with the multipliers.
      `ybus_ipm::problem` has the `Problem` trait (with nonlinear constraints and the
      Lagrangian's Hessian optional), `Constraints` (linear constraints and bounds, each
      optional) and `QuadraticProblem` (`½ xᵀ H x + cᵀ x`, as `qps_mips` builds). Constraint
      Jacobians take one row per constraint, transposed from MIPS's convention.
      `ybus_ipm::mips` has `MipsOptions` (MIPS's defaults, including step control's), the
      exit status, the multipliers as MIPS reports them, and each iteration's record.
- [x] MIPS, ported from `mips.m`, and `qps_mips`'s quadratic-program wrapper, tested on the
      problems of `t_mips.m` and `t_qps_mips.m` with MATPOWER's tolerances and iteration counts.
      - **`ybus_ipm::mips::mips` and `qps_mips`,** following `mips.m` step by step: the
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
        limit, each with step control off and on (`ybus_ipm/tests/data/mips_reference.txt`).
      - **Agreement:** every run takes MATPOWER's exit status and iteration count (including the
        infeasible LP's numerical failure after 6, and 150 at the limit), and its solution,
        objective and multipliers agree within 6.5e-14 relative (1.0e-13 since the switch to
        COLAMD in PR 2; the bar is 1e-12), except a
        failed run's diverged multipliers (about 3e24), within 3.3e-7. Each problem also meets
        `t_mips`'s known answers to its precision.
- [x] Generator costs in `ybus_powerflow`: polynomial (with `opf_setup`'s handling) and
      piecewise-linear (with `makeAy`'s constraints), reactive-power cost rows, and the total
      cost.
      - **`ybus_powerflow::opf::cost`:** `GenCost::evaluate` and `total_cost`, as `totcost` and
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
      - **`ybus_powerflow::opf::dc`:** `dc_opf` follows the legacy `opf_setup`'s DC model, after
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
          does `ybus_powerflow`, except on case3120sp (below).
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
        iteration counts aren't compared. On case3120sp `ybus_powerflow` converges at 148 iterations
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
        `ybus_powerflow` does, and fail after different numbers of iterations on 7 AC runs that fail
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
      results unchanged or within an agreed bar. **The untried parts moved to phase O5**
      (2026-11-17): rayon within one solve, SIMD, layouts and the triangular solves.
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
          (`ybus_ipm::mehrotra`).
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
      - **`ybus_ipm::mehrotra`** (`mehrotra`, `qps_mehrotra`), for MIPS's whole problem
        class. It reports MIPS's multipliers and is checked by the solver-independent
        `ybus_ipm::solution::kkt_residuals`.
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
          2. [x] the OPF's evaluation (bound 1.4–1.6× on AC), through plans kept for the solve
             that replay today's arithmetic bit for bit. Decided with the maintainer
             (2026-11-16): bit for bit; the plans public in `ybus_sparse`; the Hessian, then the
             constraints' Jacobians. In four commits:
             - [x] **Product and sum plans in `ybus_sparse`:** `ProductPlan`
               (`plan_product`) records each product term in `matmul`'s order, and `SumPlan`
               (`plan_sum`) each entry's sources as `try_add` merges them.
               - **Bit for bit** with `matmul` and `try_add`, pattern and values, on random real
                 and complex matrices in both orientations, with `±0.0`, cancellations and
                 magnitudes from 1e-8 to 1e16.
               - **Planted bugs caught:** terms in reverse order, a copied entry turned into an
                 add of zero (which turns `-0.0` into `+0.0`), and an output not zeroed first.
                 Swapping `a + b` to `b + a` is not a bug: IEEE addition commutes bit for bit.
             - [x] **The planned Lagrangian Hessian:** `WeightedHessianPlan` and
               `SquaredMagnitudeHessianPlan` in `opf::derivatives`, kept by the AC OPF's
               problem with two sum plans and `widen`'s assembly as an `AssemblyMap`.
               - **The one pattern that depends on values:** `bilinear_hessian`'s diagonal
                 terms go through `from_diagonal`, which drops zeros. Where `M + Mᵀ` stores the
                 position, the plan copies its value for a zero term, as `try_add` does. Where
                 only the term can, a zero term made the plan decline, and the Hessian was
                 computed by the functions, as before. A bus no limited branch touches has no
                 term at all. ~~Declines would be rare.~~ They weren't: see the -0.0 item
                 below, which replaced the decline in the OPF.
               - **Bit for bit:** against the functions on case14, case30, case118 and a 4-bus
                 case, at random voltages with random, zero and signed-zero weights. The
                 72-case OPF snapshot hashes identically to `main`'s on all four solver and
                 model pairs.
               - **Planted bugs caught:** the scalings reassociated, and a zero term added
                 rather than copied, which only a unit test of the diagonal step reaches. A
                 search of 200,000 small random Hessians found no input where it shows, since
                 the products' sums start at `+0.0`.
             - [x] **The planned constraint Jacobians:** `FlowsPlan` (`flows`' scalings and
               merge) and `SquaredMagnitudesPlan` (`squared_magnitudes`' scaling, real part and
               `hstack`). The constraints' two `widen`s are kept `AssemblyMap`s, and
               `branch_flows` goes through the plans, so the Hessian's flows do too. There's
               no zero-dropping step here, so these plans never decline.
               - **Bit for bit:** against `flows` and `squared_magnitudes` on the same cases.
                 The 72-case snapshot still hashes as `main`'s.
               - **Planted bugs caught:** `Y`'s scaled term reassociated, by the tests, and two
                 of `dg`'s blocks swapped, by the reference tests. Reassociating `C`'s term
                 changes nothing: `C`'s entries are all 1, which multiplies exactly.
             - [x] **Dropped terms held as `-0.0`, after the surveys found the fallback
               running on every iteration of the large cases.**
               - **Why:** near-zero-impedance ties (case2383wp's branch row 82 has `r = 0`,
                 `x = 1e-4`) keep their two buses at exactly equal voltages, so the tie's flow
                 is exactly `0` at every iteration. The flow Hessian's weights `2μ Re S` and
                 `2μ Im S` are then `0`, and its diagonal term at the far bus is dropped. All
                 32 evaluations of case2383wp's solve declined. Each computed the Hessian
                 twice, so the cases above 2,000 buses ran 0.95–1.00×, while the small ones
                 gained 1.13–1.28×.
               - **Decided with the maintainer (2026-11-16):** hold a dropped term as `-0.0`.
                 `x + (-0.0)` is `x` bit for bit for every `x`, `-0.0` included, so every later
                 sum gives `try_add`'s copy of the other operand. That replaces the copy step
                 for a covered position too.
               - **The pattern:** the plans report the positions really absent, where both the
                 real and the imaginary part dropped the term and nothing else stores it. The
                 bus Hessian stores every diagonal position, since Ybus stores its diagonal,
                 which the plans check when they're made. So the Lagrangian's pattern never
                 changes, and it falls back only if that check fails, or for a NaN term.
                 `WeightedHessianPlan::evaluate` and `SquaredMagnitudeHessianPlan::evaluate`
                 still return `None` for a really absent position, and only then.
               - **Bit for bit:** the 72-case snapshot still hashes as `main`'s. case2383wp on
                 the predictor-corrector method now takes 512 ms, from 765 ms (1.49×).
               - **Planted bugs caught:** a `+0.0` sentinel, a covered zero term added as is,
                 and the absent positions taken as the union of the two parts' drops rather
                 than their intersection.
             - [x] **The record**
               ([`opf-evaluation-2026-11-16.md`](benchmarks/opf-evaluation-2026-11-16.md)):
               - **Bit for bit:** the 72-case snapshot matches `main`'s, and none of the 1,825
                 result columns of either solver's AC survey differs.
               - **Speed, interleaved on one machine:** the predictor-corrector method on LDLᵀ
                 is 1.38–1.48× faster on case2383wp, case9241pegase and case_ACTIVSg25k, and
                 1.28× on case300. MIPS on LU is 1.08–1.22× faster, its LU dominating.
               - **The surveys' medians:** 1.26× (predictor-corrector, 1.28× in total) and 1.22×
                 (MIPS).
               - **The evaluation itself** is 2.6–2.8× faster on every profiled workload: 19.3%
                 of case9241pegase's predictor-corrector solve, from 39.3%.
               - **Next:** the KKT assembly, now 18.1% there, which `from_triplets` rebuilds at
                 every factorization. That's target 3.
          3. [x] allocation and KKT assembly, bit for bit. Decided with the maintainer
             (2026-11-17): kept plans and assembly maps, and buffer-taking `Problem` methods. In
             four commits:
             - [x] **The Newton system planned in both solvers** (`ybus_ipm`'s `newton` and
               `kkt` modules).
               - **`WeightedGram`** keeps `W = H + Jhᵀ (D Jh)` planned: `Jhᵀ` as a traced
                 gather, the product and the sum through #79's plans.
               - **`KktSolver`** keeps the KKT matrix's assembly as `AssemblyMap`s: one in CSC
                 for LU; for LDLᵀ, the lower triangle's, the mirror's gather and the full
                 matrix's.
               - **Both are keyed by their inputs' patterns,** and planned afresh when one
                 differs, as the kept analysis is.
               - **`factor` now takes `W` as a matrix** with the bounds' diagonal terms beside
                 it, rather than triplets.
               - **Bit for bit:** the 72-case snapshot matches `main`'s. Tests against the
                 triplet assembly and against `try_add` and `matmul` cover reused plans and
                 patterns that change, `H`'s alone and `Jh`'s alone.
               - **Planted bugs caught:** an assembly that's never replanned, and plans that
                 ignore `Jh`'s pattern or `H`'s. The first version of the `W` test drew new
                 patterns every round, so it never reused a plan, and missed the second. It now
                 runs in phases.
             - [x] **Buffer-taking `Problem` methods:** `objective_hessian_into`,
               `constraints_into` and `lagrangian_hessian_into`.
               - **Each defaults to the allocating method,** so implementors are unchanged.
               - **The solvers keep their buffers:** the Hessian's across iterations; the
                 constraints' in a `RefCell` beside each solver's evaluation, which takes
                 `&self` in MIPS.
               - **`AcOpfProblem` refills in place** when the buffer holds the planned pattern,
                 compared as slices (`ybus_ipm::problem::same_pattern`). Its allocating methods
                 now call these.
               - **`QuadraticProblem`** copies `H`'s values rather than cloning it.
               - **Not done here:** #79's plans still allocate their intermediate vectors. Those
                 move into scratch in the next commit only if the allocation profile names them.
               - **Bit for bit:** the 72-case snapshot still matches `main`'s.
             - [x] **The rest of the loops' allocation, by profile.**
               - **What the profile found:** samply on case9241pegase after the first two
                 commits, attributing each "allocation and copying" sample to its nearest frame
                 from our crates.
                 - **The KKT assembly** was already down to 2.2% of the predictor-corrector
                   solve, from 18.1%.
                 - **#79's plans' intermediate vectors** came to about 0.1%, so they stay as they
                   are.
                 - **MIPS's allocation** was almost all `lu::gilbert_peierls` growing `L`'s and
                   `U`'s arrays one push at a time: 5.9% of the solve.
                 - **The LDLᵀ's** was spread over sites of 0.4% or less, left for the kernel's
                   PR.
               - **The change:** `SymbolicLu` keeps the last factorization's sizes of `L` and
                 `U`, in atomics, so the analysis stays `Send`, `Sync` and `Clone`. Each
                 factorization reserves them, and a repivoting refactorization reserves its own
                 old sizes. Only the capacity changes, never the arithmetic.
               - **Bit for bit:** the 72-case snapshot still matches `main`'s. A kept analysis
                 whose fill changes between factorizations factors as a fresh one does.
             - [x] **The record**
               ([`kkt-assembly-2026-11-17.md`](benchmarks/kkt-assembly-2026-11-17.md)):
               - **Bit for bit:** the 72-case snapshot matches `main`'s, and none of the 1,825
                 result columns differs in any of the four surveys (AC and DC, both solvers).
               - **Speed, interleaved:**
                 - the predictor-corrector method 1.20–1.34× on AC and 1.09–1.26× on DC
                 - MIPS 1.03–1.10× on AC and 1.12–1.16× on DC
               - **The KKT assembly** is 10–20× smaller everywhere, at 0.2–3.2% of each solve.
               - **Next:** target 4. The LDLᵀ's dense kernel is now 27.6% of case9241pegase's
                 AC predictor-corrector solve.
          4. [x] the kernel and the KKT ordering. Planned (2026-11-17) as two PRs: the dense kernel,
             then the ordering.
             - **The kernel, bit for bit:** found while planning, both of its changes keep the
               same arithmetic. Each trailing entry's updates are a fixed sequence of
               subtractions in pivot order, which a blocked kernel can keep while changing only
               the loop order and the memory traffic.
               - [x] **The lower triangle only, mirrored once.**
                 - **Today's waste:** each pivot copied the whole trailing lower triangle into
                   the upper, a strided write per entry, though only the threshold test and the
                   swap read the upper.
                 - **The change:** they now read `conj(a_ji)` from the lower triangle, which has
                   the same magnitude. The swap is a Hermitian one on the lower triangle, as
                   LAPACK's `xHETRF` does. `factor_front` writes the upper once, at the end, as
                   it was left before, so the multifrontal driver and `DenseLdlt` are
                   unchanged.
                 - **Bit for bit:** against a test-only copy of the old kernel, on random
                   Hermitian fronts, real and complex, of 1 to 150 rows. The comparison covers
                   the whole array, the order and D's blocks. The cases include 2×2 pivots,
                   delays and singular fronts, which the test requires it to reach. The
                   72-case snapshot matches `main`'s.
                 - **Planted bugs caught:** a swap that doesn't conjugate the joining entry, or
                   the stretch between, and a mirror that loses a 2×2 block's upper entry.
               - [x] ~~A blocked, left-looking panel.~~ Tried and dropped (2026-11-17), bit for
                 bit but slower.
                 - **Only the contribution block can be deferred exactly:** in the fully summed
                   part, swaps conjugate entries across the diagonal, and
                   `conj(w_j)·l_i ≠ conj(w_j·conj(l_i))` bit for bit.
                 - **Deferring it was 3–15% slower** than the lower-triangle kernel alone,
                   whether swept by column, tiled by 8 columns, or cut off by block size. The
                   fronts are small, and the eager update already streams contiguous columns.
               - [x] **The record**
                 ([`ldlt-kernel-2026-11-17.md`](benchmarks/ldlt-kernel-2026-11-17.md)):
                 - **Bit for bit:** the fronts, the 72-case snapshot, and the predictor-corrector
                   surveys, DC and AC.
                 - **The kernel stage** is 1.23–2.00× faster.
                 - **The solves:** AC 1.11–1.26× faster, and 1.24× in the survey's total (case_SyntheticUSA
                   94.9 s → 75.2 s). DC 1.13× on case_ACTIVSg25k, and 1.19× in total.
                 - **The temporary reference kernel is removed.** It's at `9734167`, to
                   restore for any later change to the kernel.
               - **Next:** the KKT-aware ordering, the second PR.
             - **The KKT-aware ordering** (planned 2026-11-17): Duff and Pralet's matched
               ordering as an option, made a default only by a rule written before its survey.
               The survey also closes the ordering study. In six commits:
               - [x] **Analyses with a given permutation:** `analyze_cholesky_with_permutation`
                 and `analyze_ldlt_with_permutation` in `ybus_sparse`, as LU's
                 `analyze_lu_with_permutation` already was. Given AMD's own permutation, each
                 gives AMD's analysis.
               - [x] **The KKT factorizations' statistics:** `KktStatistics` on `Solution`
                 (`solution.kkt`), with the factorizations, delayed pivots, 2×2 pivots and
                 factor entries totalled over a solve, the inertia correction's
                 refactorizations included. `opf_survey` gains the four columns. These are
                 committed rather than a temporary instrument, decided with the maintainer, so
                 any later ordering change surveys the same way. The 72-case snapshot is
                 unchanged.
               - [x] **The matched ordering:** `KktOrdering { Amd, Matched }` on both solvers'
                 options (`kkt_ordering`, AMD by default; LU ignores it),
                 `OpfSolver::with_kkt_ordering`, `opf_survey --kkt-ordering` and
                 `YBUS_BENCH_KKT_ORDERING`.
                 - **The method:**
                   - each constraint row is matched to a variable, greedily by `|Jg|` and then by
                     augmenting paths, so the matching is maximum
                   - each pair is merged into one node, and the merged graph ordered by AMD
                   - each pair is placed as its variable, then its constraint
                 - **On case118, against AMD,** with the same iterations and objective:
                   - **DC:** 44 → 20 delays, 13,516 → 13,255 factor entries
                   - **AC:** 303 → 115 delays, 787 → 112 2×2 pivots, 90,629 → 83,894 factor
                     entries
                 - **With AMD the default,** the 72-case snapshot is unchanged.
                 - **Planted bugs caught:**
                   - a pair placed constraint first
                   - a greedy pass that reuses a variable
                   - the ordering given to LU, caught by the MATPOWER reference test
               - [x] **The rule, written before the survey** (2026-11-17, decided with the
                 maintainer as #59's). For each model, DC and AC, the matched ordering becomes
                 the predictor-corrector default (`DcOpfOptions::mehrotra()`,
                 `AcOpfOptions::mehrotra()`) if, on all 73 cases with costs, against the same
                 method with AMD:
                 1. **O3's bar holds wherever it holds with AMD:**
                    - convergence wherever MIPS converges
                    - the objective within 1e-6 (relative) of MIPS's and MATPOWER's
                    - KKT residuals within the tolerances
                    - prices matching where they're unique
                 2. **The outcomes are the same,** case by case.
                 3. **The median `opf_ms` is faster.**

                 Iterations are reported, not required. A model whose survey fails any condition
                 keeps AMD, and the record says which condition failed.
               - [x] **The survey, and the defaults the rule gives** (all 73 cases, DC and AC,
                 the predictor-corrector method on LDLᵀ with AMD and with matched, and on LU,
                 and MIPS on LU).
                 - **DC: the rule holds, so `DcOpfOptions::mehrotra()` now orders by the matched
                   ordering.**
                   1. O3's bar is met on 58 of 58 where AMD meets it.
                   2. The outcomes are the same on all 73.
                   3. It's 1.14× faster in the median, and 1.27× in total.
                 - **AC: the rule fails on its second condition, so `AcOpfOptions::mehrotra()`
                   keeps AMD.** case13659pegase ends at the iteration limit with AMD (150) and in
                   a numerical failure with matched (147). Neither converges, but the outcomes
                   differ, and the rule asks them to be the same.
                   - **The rule's other conditions held:** O3's bar on 30 of 30, and 1.09×
                     faster in the median and 1.17× in total.
                   - **Not revised after the fact.** A later rule, written before another survey,
                     could treat two failures as the same outcome.
                 - **The snapshot changes only in its DC predictor-corrector column,** on 67 of 72
                   cases. MIPS and AC are untouched.
               - [x] **The record**
                 ([`kkt-ordering-2026-11-17.md`](benchmarks/kkt-ordering-2026-11-17.md)), which
                 closes O4:
                 - **Delays:** 7.4× fewer on DC and 22× fewer on AC, with 19–21% less fill.
                 - **The DC default:** 1.14–1.30× faster.
                 - **The study:** AMD and matched on LDLᵀ, and COLAMD on LU, over 64 DC and 49 AC
                   cases.
- [x] **The KKT solve.** Reusing the symbolic analysis across iterations while the pattern stays
      fixed; a symmetric indefinite LDLᵀ factorization in `ybus_sparse`, with Bunch–Kaufman-style
      pivoting, written from the literature and compared with the LU; and reduced against
      augmented KKT systems.
      - [x] **Reduced against augmented: closed by decision (2026-11-17).** The solvers keep
        MIPS's reduced `[W Jgᵀ; Jg 0]`, with the inequalities folded into `W` and the right-hand
        side, because:
        - MIPS's parity with MATPOWER needs it
        - the augmented system is larger, by the inequalities' rows
        - it needs the same 2×2 pivoting

        Not implemented.
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
              - a globalization for the predictor-corrector method, now goal 5 of
                [`LongTermGoals.md`](LongTermGoals.md) (2026-11-18). The other two are done:
                the kernel's cost per front (#61, #62, #81) and the KKT-aware ordering (#82).
- [x] **A fuller study of the KKT ordering and factorization,** once phase O3's solver exists to
      compare against. MIPS's LU switched from AMD to COLAMD in O1's PR 2: with zero diagonals
      in the equality block, AMD's factors filled 14 million entries on `case_ACTIVSg10k`
      (24 s per factorization) against COLAMD's 640,000 (0.036 s). COLAMD lost one iteration-count
      match (case2736sp's DC OPF, 44 to MATPOWER's 45). The study would cover AMD, COLAMD,
      orderings that keep the KKT structure (constraints after the variables they touch), LDLᵀ
      and UMFPACK-like strategies, by fill, time and agreement with MATPOWER.
      - **Done in target 4's ordering survey**
        ([`kkt-ordering-2026-11-17.md`](benchmarks/kkt-ordering-2026-11-17.md)):
        - **The matched ordering** beats AMD on LDLᵀ on every measure but AC's one differing
          outcome.
        - **COLAMD on LU** fills 2–3× more than either, and stays MIPS's, for MATPOWER parity.
        - **UMFPACK-like strategies weren't pursued:** the multifrontal LDLᵀ already plays their
          part.
- [ ] **Parallel and GPU.** Batched OPF solves with rayon, and a GPU KKT factorization if phase
      O2's spike supports one; warm starts across related problems, such as outage sequences.
      Recorded. **Moved to phase O5** (rayon, SIMD and GPU), with O2's untried "The CPU" ideas:
      see "Phase O5" below.

**Phase O4 is complete** (PRs #53–#62 and #79–#81, and this one, 2026-11-17). Its open item,
"Parallel and GPU", continues in O5.

## Phase O5: rayon, SIMD and GPU

Planned 2026-11-19, after the warm starts, as the last optimization phase before the publish
(P5). It gathers O2's untried "The CPU" ideas and O4's "Parallel and GPU" item.

- **Decided with the maintainer (2026-11-19):**
  - **The targets:** all four, in this order unless the bounds say otherwise:
    1. batched OPF solves
    2. parallelism within one solve
    3. SIMD and layouts
    4. a GPU KKT factorization
  - **Bit for bit, independent of the thread count:**
    - **Parallel paths:** they give the serial path's results bit for bit, whatever the
      thread count, through fixed partitions (never by thread count) and ordered reductions.
    - **Reordered arithmetic:** a change that reorders it (SIMD reductions, FMA) must meet O3's
      bar and be opt-in, or be the predictor-corrector default within the bar, by the standing
      rule.
    - **MIPS's default** stays bit for bit.
  - **Profile and targets first,** as O2 and O4 began.
- **Standing constraints:**
  - **rayon stays an optional feature,** and the serial path the default. `ybus_ipm` gains a
    `rayon` feature, forwarded by `ybus_powerflow/rayon`. The C build keeps default features
    only (F0).
  - **The GPU comes last.** GPU runs have crashed WSL here, so state stays in the repo, commits
    come before a GPU run, and runs start small.
- **How each target is checked:**
  - the 72-case snapshot
  - results at 1, 2, 7 and 28 threads, bit for bit against the serial path
  - planted bugs, such as a reduction in thread order or an unordered extend-add
  - a survey whose rule is written before it runs
- [x] **The profile and the targets** (O5 PR 1;
  [`o5-profile-2026-11-19.md`](benchmarks/o5-profile-2026-11-19.md)):
  - **What it measured:**
    - re-timed baselines, with outcomes and iterations as archived, and times 1–16% slower
      than 2026-11-17's
    - samply profiles, all at 1.0 threads busy
    - each solve's front tree (`FrontTree`, `KktStatistics::tree`, and the `o5_parallelism`
      example)
    - the machine's ceiling for independent work: about ×10–11 on 14 threads, and ×5 on
      case_ACTIVSg70k
  - **The bounds:**
    - **Batched solves:** about ×10, for both solvers.
    - **Within one solve:** ×1.5–2.1 (DC) and ×2.6–3.0 (AC) for the predictor-corrector method.
      For MIPS it's at most ×1.12, since its LU isn't tree-parallel.
    - **SIMD in the dense kernel:** at most ×1.14–1.23.
  - **The order stands.**
- [x] **Batched OPF solves** (O5 target 1): many solves sharing a base (N-1 OPF warm from the
  base with `map_start`, or a day's hours), each on its own thread with its own model, results
  in input order. Bound: about ×10.
  - **The survey's rule, written 2026-11-20 before it runs:**
    - **Bit for bit at every thread count.** `opf_batch` compares every parallel batch with
      the serial one and fails on a difference, so a difference is a bug, not a result.
    - **The speedups are recorded, with nothing chosen.** Target 1 adds APIs and changes no
      default.
    - **The record compares them with the ×10 bound,** and explains any batch far below it.
  - **The survey:**
    - **The cases:** case2383wp, case3120sp and case9241pegase, DC and AC, both solvers, both
      batches, cold and warm at each solver's own `warm_shift`. Also case_ACTIVSg25k on the
      predictor-corrector method, with 7 outages of each kind.
    - **The threads:** 1, 2, 4, 8, 14 and 28.
  - **Changed from the plan:**
    - **14 branch and 14 generator outages per case,** not 50 each, to keep the survey near two
      hours. 28 fill 28 threads in one wave.
    - **Each process times the serial batch once,** then the `par_` form on each pool, rather
      than one process per thread count, which would repeat the serial batch six times.
  - **The result** ([`o5-batched-solves-2026-11-20.md`](benchmarks/o5-batched-solves-2026-11-20.md)):
    - **Bit for bit:** every parallel batch matches its serial one, all 321 rows at 1–28
      threads.
    - **MIPS:** a median ×9.35 (DC) and ×7.69 (AC) at 28 threads, up to ×11.6, reaching the
      bound.
    - **The predictor-corrector method:** ×5.56 (DC) and ×4.13 (AC) at 14 threads, and less at
      28. The balance doesn't explain it, since its cold batches could reach ×22–27. Memory
      bandwidth or allocation are the likely causes, for target 3.
    - **Warm outage batches:** held to ×2–4 by their slowest outage, which was measured.
    - **Memory:** case_ACTIVSg25k's AC batches went over 12 GiB beyond 8 and 14 threads.
- [x] **Within one solve** (O5 target 2): the LDLᵀ's independent subtrees with `rayon::join`,
  with each child's update extended into its parent in the serial order; and the AC
  evaluation's plans split by fixed element ranges. Bound: ×1.5–3.0 on the predictor-corrector
  method.
  - **Decided (2026-11-20):**
    - **Two PRs:** PR 2a, the LDLᵀ, then PR 2b, the AC evaluation, each with its own survey.
    - **Automatic with `rayon`,** by fixed work cutoffs, with no new option.
    - **Node-level parallelism** (the large fronts' trailing update, split by columns) in PR 2a.
  - **PR 2a's cutoffs, from a calibration** (decided 2026-11-20, before the survey):
    - **The plan** was to calibrate on case2383wp alone. There the parallel factorization gains
      nothing on DC and ×1.14 on AC, so that rule would have turned it off.
    - **The large cases do gain:** at 14 threads, case_ACTIVSg25k ×1.44 (AC) and ×1.42 (DC),
      and case9241pegase ×1.24 and ×1.11.
    - **The choice:** over all six workloads, subtree work 2²⁰ with no column split. It has a
      geometric-mean speedup of ×1.15, and the best worst case, ×0.97 on case2383wp's DC.
    - **The column split** stays in the code, tested bit for bit, but off: it never beat the
      noise, and from trailing order 256 it made solves 10–20× slower.
    - **Calibrated on the survey's cases,** so the record's speedups on them are at tuned
      cutoffs. The record says so.
  - **PR 2a's survey rule, written before it runs:**
    - **Bit for bit:** every OPF survey result column matches the serial build's, at every
      thread count.
    - **The speedups are recorded, with nothing chosen.**
    - **The record compares them with the bounds:** the trees' ×4–9 for the factorization, and
      the whole-solve ×1.5–3.0, part of which is PR 2b's.
  - [x] **PR 2a, the parallel LDLᵀ**
    ([`o5-parallel-ldlt-2026-11-20.md`](benchmarks/o5-parallel-ldlt-2026-11-20.md)):
    - **Bit for bit:** 1,740 survey columns at 1–28 threads, and the 72-case snapshot in both
      builds.
    - **The predictor-corrector method's single solves:** AC a median ×1.22 at 4–8 threads (up
      to ×1.34), DC about ×1.0 (×1.36 on case_ACTIVSg25k), with only 1.3–1.6 threads busy at 14.
      The root fronts, on the critical path, hold most of the work, and the column split
      didn't pay. A blocked kernel is target 3's.
    - **Off inside batches and on a one-thread pool** (decided 2026-11-20): `opf::batch`'s
      `par_` forms run each solve's factorizations serially (`ybus_sparse`'s `serially`). The
      batch then matches `main`'s same-day speed, and the serial batch loop gains ×1.13.
    - **A correction:** the batch slowdown first reported (10%) compared two days. Same-day, it
      was 0–4%.
  - [x] **PR 2b, the AC evaluation's plans**
    ([`o5-parallel-evaluation-2026-11-21.md`](benchmarks/o5-parallel-evaluation-2026-11-21.md)):
    - **Decided (2026-11-21):** coarse tasks only. The evaluation's independent parts run
      concurrently with `rayon::join`, and the kernels aren't split.
    - **The cutoff,** calibrated before the survey at 14 threads, each case with the
      evaluation on and off:
      - case300 is slower with it on (×0.75–0.87)
      - case1354pegase and larger all gain: ×1.05–1.18 on the predictor-corrector method, and
        ×1.01–1.02 on MIPS
      
      By PR 2a's rule, it goes concurrent from 1,000 buses (`EVALUATION_BUSES`). The
      calibration's timings are beside the record.
    - **The survey's rule, written before it runs:**
      - **Bit for bit:** every result column matches the serial build's at every thread count.
      - **The speedups are recorded, with nothing chosen.**
      - **The record compares them with the bound:** the evaluation's share, about ×1.1–1.25 on
        the AC solve. It separates this PR's gain from PR 2a's with a same-day run of PR 2a's
        branch.
    - **The result:**
      - **Bit for bit:** 2,320 survey columns, and the 72-case snapshot in both builds.
      - **With PR 2a,** the predictor-corrector method's AC solves reach a median ×1.50 at
        8 threads, up to ×1.67, from PR 2a's ×1.15 alone the same day, so this PR adds
        about ×1.3. MIPS gains ×1.03.
      - **Batches:** unchanged, and their serial loop gains ×1.27–1.37.
      - **The rest of the bound** is in the evaluation's largest part, the root fronts, the
        serial triangular solves and rayon's runtime. The solves and kernel are target 3's.
  - **Target 2 is complete:** within one solve, the predictor-corrector method's AC solves
    reach ×1.50 at 8 threads, about half of the ×2.6–3.0 bound.
- [ ] **SIMD and layouts** (O5 target 3):
  - **The dense kernel's rank-k update,** vectorized without reordering any element's
    operations, if it isn't auto-vectorized already.
  - **Supernodal or blocked triangular solves.**
  - **Less allocation and copying.**
  
  Reductions that reorder are opt-in under O3's bar. Bound: ×1.14–1.23 from the kernel.
- [ ] **A GPU KKT factorization** (O5 target 4): a cuDSS spike behind `cuda`. It's gated: only
  if the KKT factorization is still the largest stage of a large AC solve after targets 1–3,
  and O2's transfer measurements leave it room.

## Warm starts for OPF sequences

Planned 2026-11-18, before O5, from O3's and O4's open items ("warm starts across related
problems, such as outage sequences"). Many uses of an OPF solve a sequence of related problems.
Each solve now starts cold, rebuilding its problem and starting from the same point.

- **Decided with the maintainer (2026-11-18):**
  - **Both kinds of sequence, in stages:** the same network with new data first, then outages.
  - **What carries over:** the primal point and the duals, shifted off the boundary (Yildirim
    and Wright, 2002; Gondzio and Grothey, 2008). Measured against a primal-only warm start and
    a cold start.
  - **Both solvers, opt-in:** a cold start stays today's path, MIPS's bit for bit.
  - **A kept OPF model first,** as `AcModel` is for the power flow.
- [x] **Kept OPF models** (PR A): `AcOpfModel` and `DcOpfModel`. `new(case)` builds what the
  network fixes, and `solve(case, options)` reads only values. `ac_opf_with_options` and
  `dc_opf_with_options` are now `new(case)?.solve(case, options)`.
  - **What AC keeps:** the network, its admittances, and the problem with #79's plans and #80's
    assemblies. Each solve copies the problem with its case's values.
  - **What AC refuses,** with `TopologyMismatch` naming the field: a changed topology, and any
    value that fixes a pattern:
    - which branches are limited (`rate_a_mva`)
    - which have angle limits (`angmin_deg`)
    - which generators are dispatchable loads (`pmin_mw`)
    - a piecewise cost's points (`ncost`)
    - which quadratic coefficients are 0 (`quadratic`)
  - **What DC keeps:** the network and its B matrices. Nothing kept depends on the limits or
    costs, so only the topology is refused.
  - **Not kept: the solver's KKT analysis,** a change from the plan (2026-11-18). Keeping it
    needs a public workspace type in `ybus_ipm`, and the analysis is about 1% of a solve in the
    profiles.
  - **Bit for bit:**
    - the 72-case snapshot is unchanged
    - a model solving a copy with loads scaled ×0.95 or ×1.04 matches the one-shot solve of the
      copy, on case9, case30, case118 and case2383wp, both solvers, AC and DC
  - **Planted bug caught:** a solve keeping the model's own loads.
  - **What it saves on its own: little**
    ([`opf-models-2026-11-18.md`](benchmarks/opf-models-2026-11-18.md)). On an 11-solve load
    sweep the speedup is 0.98–1.04×, since a solve's own work dwarfs its setup. The models are
    for what they enable: warm starts, and sequences through one value with pattern changes
    refused.
- [x] **Warm starts in the solvers** (PR B, #85): `ybus_ipm::start` with `Start::Cold` and
  `Start::Warm(&WarmStart)`. The solvers take them through `mips_from`, `mehrotra_from` and their
  `qps_` forms, and the original functions delegate with a cold start, bit for bit.
  - **The warm start:**
    - It inverts each solver's multiplier report, a function whose round trip a unit test checks.
    - It recomputes the slacks from `h(x)`.
    - It raises every inequality slack and multiplier to at least `warm_shift`.
    - The predictor-corrector method pushes its point `warm_shift` inside its bounds.
  - **Three defects, found by PR C's OPF sequences and fixed there:**
    - the multipliers weren't scaled back by `cost_mult`
    - the predictor-corrector method's bound slacks were raised, letting a step cross its bound
    - it read the thresholded report instead of the final iterate's multipliers
- [x] **Warm OPF sequences on one network** (PR C): `solve_from` with `OpfStart` on both
  models, the `opf_sequence` survey, and `warm_shift`'s default from it.
  - **The result** ([`opf-warm-start-2026-11-18.md`](benchmarks/opf-warm-start-2026-11-18.md)):
    no shift qualifies under the corrected rule, for either solver, so `warm_shift` stays at
    the provisional 1e-4.
    - **MIPS at 1e-6** comes closest. It has no regressions and 38 rescues, with median
      speedups of ×2.11 (DC) and ×3.40 (AC), and misses on 3 AC steps of 616.
    - **The predictor-corrector method** fails most on case9241pegase's AC, where an absolute
      shift moves the start far from stationary.
  - **Left for later:** a shift relative to the problem (for example, to the previous
    solution's complementarity, as Gondzio and Grothey's is), surveyed by the same tool.
  - **The rule, written 2026-11-18 before the survey,** as #59's and #82's were:
    - **The sequences:** each case's loads swept ×0.90 to ×1.10 (21 steps), and a seeded random
      walk of 25 steps, each bus ±2% a step. Every step is solved cold, primal-only and warm
      at each shift, each chain from its own previous solution.
    - **The shifts:** 1, 1e-1, 1e-2, 1e-4 and 1e-6, widened from the planned 1e-2, 1e-4 and
      1e-6 (decided 2026-11-18). On case2383wp's DC OPF, the predictor-corrector method warm
      from ×1.00 to ×1.02 converges from 1, 1e-1 and 1e-2 in 17, 15 and 22 iterations (31 cold),
      and fails from 1e-4 and 1e-6 on the step-size limit. Its multipliers reach 1e5.
    - **A warm step is acceptable** if it has the cold step's outcome and, where the cold step
      converged, meets O3's bar against it without `LAM_Q`:
      - the objective within 1e-6
      - `LAM_P` within 1e-3
      - each KKT residual within 100× the cold step's or below 1e-6
    - **Why `LAM_Q` is reported but not barred** (decided 2026-11-18): two cold solvers already
      differ in it. MIPS and the predictor-corrector method, each against MATPOWER, differ by
      3.5e-4 on case118, 7.3e-3 on case300 and 9.4e-3 on case9241pegase. On case118's AC,
      the warm and primal steps' only misses were `LAM_Q`'s, at 1.0e-3 to 1.8e-3.
    - **On AC, a warm step reaching a different local optimum** is counted and listed per case,
      but doesn't disqualify a shift. It means both steps converged with every residual within
      the bar, but the objective beyond 1e-6 (its prices differ with it). The problem is
      nonconvex, the warm start is opt-in, and neither optimum is wrong.
    - **`warm_shift`'s default** is the shift with the fewest median iterations over all warm
      steps, among those with no unacceptable step on DC. On AC, the steps where the cold step
      converged count too, apart from the different local optima. Both solvers get the same
      default. Ties go to the larger shift, which is safer. If no shift qualifies, the default
      stays provisional and the record says why.
    - **Corrected 2026-11-18, after the first tier-1 run and before the rerun** (decided with
      the maintainer). The rule as written above chose no shift, for reasons outside the warm
      start:
      - **The outcome:** "the cold step's outcome" counted a warm step that converged where the
        cold one failed as unacceptable: 38 such rescues on MIPS's DC alone. Now only a
        **regression** is unacceptable, meaning the cold step converged and the warm one didn't.
        Rescues are counted as gains.
      - **The residuals' reference:** "100× the cold step's" flagged warm steps that stopped at
        the solver's tolerance where the cold solve over-converged: on case9, a stationarity of
        5e-6 against 1e-11, with objectives agreeing to 1e-11. They're now compared with
        **MIPS-LU's cold solve of the step**, as O3's bar compared them. They aren't compared
        where that solve didn't converge.
      - **Per solver:** each solver gets the default this rule chooses from its own steps, since
        the best shifts differ by solver (1e-6 for MIPS on DC and AC; 1e-2 on DC and 1e-6 on AC
        for the predictor-corrector method). A solver with no qualifying shift keeps the
        provisional 1e-4.
      - **The first run also found a defect, fixed before the rerun:** the warm start read the
        thresholded multipliers (`mu_threshold`). It now reads the final iterate's,
        `Solution::iterate_multipliers`.
    - **The cases:** tier 1 runs every shift on case9, case30, case118, case300, case2383wp,
      case3120sp and case9241pegase, DC and AC, both solvers, both sequences. Tier 2 runs the
      chosen shift on case13659pegase (AC), case_ACTIVSg10k (DC, AC) and case_ACTIVSg25k (DC,
      and AC on the predictor-corrector method only).
- [x] **Warm starts across outages** (PR D): `map_start` on both models, with
  `OpfStart::Mapped` and `MappedPrimal`.
  - **The mapping:** each model keys its variables and constraints by case row, built from the
    same lists that build its problem. A solution of one model carries onto another model of
    the same case's tables, such as with a branch or a generator out of service.
  - **The survey** ([`opf-outage-warm-start-2026-11-19.md`](benchmarks/opf-outage-warm-start-2026-11-19.md)):
    up to 50 sampled branch outages (not bridges) and 50 generator outages per tier-1 case,
    each warm from the base solution, as an N-1 study runs, at all five shifts.
    - **MIPS at 1e-6:** ×2.05–2.24 (DC) and ×2.99–3.38 (AC), with 49 rescues, and 4
      regressions and 6 misses in 1,106 steps.
    - **The predictor-corrector method:** it fails where it failed in PR C, on case9241pegase's
      AC at the large shifts.
    - **A model per outage** costs about 2% of a cold solve.
  - **Decided (2026-11-19):** `warm_shift`'s default wasn't revisited, so it stays at 1e-4.
- **The warm-start section is complete.** Left for later is a shift relative to the problem
  (PR C), which both surveys' tools can measure.

## References

Each was checked online (Crossref) on 2026-11-03, and the last on 2026-11-17.

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
- I. S. Duff and S. Pralet, "Strategies for scaling and pivoting for sparse symmetric indefinite
  problems," *SIAM J. Matrix Anal. Appl.* 27(2):313–340, 2005.
  [doi:10.1137/04061043X](https://doi.org/10.1137/04061043X). The matched ordering of the KKT
  matrix (`ybus_ipm::kkt::KktOrdering::Matched`).
