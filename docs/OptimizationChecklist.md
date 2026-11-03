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
- [ ] MIPS, ported from `mips.m`, and `qps_mips`'s quadratic-program wrapper, tested on the
      problems of `t_mips.m` and `t_qps_mips.m` with MATPOWER's tolerances and iteration counts.
- [ ] Generator costs in `powerflow`: polynomial (with `opf_setup`'s handling) and
      piecewise-linear (with `makeAy`'s constraints), reactive-power cost rows, and the total
      cost.

### PR 2: OPF references and the DC OPF

- [ ] OPF reference results: `dcopf` (MIPS forced) and `opf` runs in `generate_reference.m`,
      committed for case14, case30 and case118 and generated locally for every case with costs.
- [ ] The DC OPF: `[Va, Pg, y]` with the power balance, flow limits, angle limits and the
      piecewise-linear rows, solved through the quadratic-program wrapper, against `rundcopf`.

### PR 3: the AC OPF and the survey

- [ ] The AC OPF: `[Va, Vm, Pg, Qg, y]` with the power balance and flow limits (and their first
      and second derivatives), angle limits, dispatchable loads and the piecewise-linear rows,
      from `nlpopf_solver`'s starting point, against `runopf`.
- [ ] The OPF survey over every case with costs, DC and AC, rechecking MP-Core against the
      legacy code, recorded in `docs/benchmarks/`, with the bar tightened to the measured
      agreement. Then the docs.

## Phase O2: accelerating the power flow

- [ ] **Decide** what to accelerate first, from profiles of the existing surveys: Newton,
      fast-decoupled, DC, N-1 outages and reactive limits. Record where the time goes.
- [ ] **The CPU.** Parallelism within one solve (refills, mismatches and flows with rayon), SIMD
      in the hot kernels, cache-friendlier layouts, and level-scheduled or supernodal triangular
      solves if the profiles point there. Each measured against the current records, with the
      results unchanged or within an agreed bar.
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
