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
  README names the files consulted; `sparse_mat` and the phase O3 solver are written from
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
- **In `sparse_mat`:** LU with refactoring and a given ordering, sparse right-hand sides,
  transposed solves, iterative refinement, condition estimates and parallel batch solves, all of
  which the KKT systems of an interior-point method can use.

## Phase O1: OPF following MIPS

- [x] **Decide:** where the solver lives, which MATPOWER code is the reference, which constraints
      come first, and the acceptance bar.
      Decided (2026-11-03):
      - **The solver** goes in a new workspace crate, `interior_point`: MIPS now, the phase O3
        solver later, behind one problem interface. It depends on `sparse_mat`; `powerflow`,
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
        the starting slacks and multipliers, the reduced KKT system solved with `sparse_mat`'s
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
          adjacency rebuilt and searched per outage; `sparse_mat`'s `solve_sparse_into` with a
          reusable `SparseSolveWorkspace`; per-thread outage workspaces; and `finish` computing
          flows into the result and `Bbus θ` only at the reference buses.
        - **Results:** bit for bit the same. All outages of all cases take 4.5× less time
          serially and 4.6× in parallel; on case_SyntheticUSA, 157 s serially, from 680 s. The
          reachability check gave most of it: its breadth-first search, inlined into the outage,
          hid from the profile's allocation share.
        - **What's left:** each outage's dense updates and back substitutions, which a streaming
          API returning only the flows a caller wants could cut (deferred).
      - [ ] **The LU refactorization kernel** (target 2).
      - [ ] **Parsing** (target 3).
- [ ] **A CUDA spike** behind a `cuda` feature (through `cudarc`). Batched DC and N-1 solves and
      the Jacobian refill on the GPU, with transfers counted, against the CPU path. The CUDA
      toolkit gets installed only when cuSPARSE or cuSOLVER is needed. Record the outcome.

## Phase O3: a naive interior-point solver

- [ ] **Decide** the method, for example a textbook primal-dual interior-point method, possibly
      with Mehrotra's predictor-corrector, and its references (Nocedal and Wright; Mehrotra
      1992). It takes the same problem interface as the MIPS crate.
- [ ] **The solver,** written from the literature, tested on the MIPS test problems.
- [ ] **Verified and benchmarked** against phase O1's MIPS and against MATPOWER, on DC and AC OPF
      over every case with costs: objective, KKT residuals, iterations and time. Recorded.

## Phase O4: accelerating the interior-point solvers

- [ ] **Decide** what to accelerate, from profiles: the KKT solve is the expected hot spot.
- [ ] **The KKT solve.** Reusing the symbolic analysis across iterations while the pattern stays
      fixed; a symmetric indefinite LDLᵀ factorization in `sparse_mat`, with Bunch–Kaufman-style
      pivoting, written from the literature and compared with the LU; and reduced against
      augmented KKT systems.
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
