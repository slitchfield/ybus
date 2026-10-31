# Implementation report: `dev/dc-power-flow`

This branch adds the DC power flow and finishes Phase 1.5: the model, the solve, the comparison
with MATPOWER's results, and benchmarks. It has nine commits of work on top of `main` at
`b7b2f36`, and this report with its updates. The first two are from an earlier session, which
ended before writing a report. This session reviewed them against MATPOWER's source, fixed the
one problem the review found, and added the rest. The last two, on backward error, answer a
question the maintainer asked about the first version of this report. This report covers the
whole branch.

| | `main` | `dev/dc-power-flow` |
|---|---|---|
| `powerflow` tests | 62 | 97 (adds 13 DC model, 15 DC solve, 4 reference comparison and 3 backward error) |
| MATPOWER 8.1 cases the DC power flow solves | none | all 78, within 1.1e-10 of MATPOWER's results relative to each quantity's largest value, and backward stable (at most 10.7 ε) |
| Benchmarks | none | criterion on case2383wp, and a recorded survey of all 78 cases |
| Diff, without this report | | 22 files, +2,245 / −24: 481 lines of library code, 1,043 of tests, 368 of benchmark code, scripts and manifest, and 353 of documentation |

## Commits

| Commit | Change |
|---|---|
| `6484fd3` | The DC model: Bbus, Bf and phase-shifter injections (earlier session) |
| `35c3cc3` | The DC solve: angles, branch flows and slack generation (earlier session) |
| `7036a35` | Validates the case before building the model, a fix from this session's review |
| `eec0722` | Compares case14, case30 and case118 with MATPOWER's committed results |
| `0e31586` | Criterion benchmarks on case2383wp |
| `da2b907` | A survey of every MATPOWER case, and the recorded results |
| `97a1648` | The README |
| `2659525` | Tests the solve's backward error on six MATPOWER cases |
| `0b9eb18` | Records the backward error in the survey, and tracks its link to fill-in |
| `4fd2918` and later | This report, with corrections and updates |

## How the work was driven

1. **Review before building on it.** The earlier session's two commits were checked against
   MATPOWER 8.1's `makeBdc.m`, `dcpf.m` and `runpf.m`, then probed with invalid cases, which
   found the validation gap.
2. **Measured tolerances.** The reference tests' tolerance was set after measuring the actual
   differences, and states them.
3. **Measurements over inferences in the benchmarks.** A first pass at 20 s per case checked the
   tools and informed the choice of benchmark case. Claims about fill-in rest on counted factor
   entries, and the survey was rerun with the final tool so that the committed tool reproduces
   the record.
4. **Separating the solver from the problem.** Asked whether differences from MATPOWER that grow
   with the case were acceptable, I measured backward errors and condition numbers before
   answering, and the guess that came out of it is recorded with a baseline, to be tested
   rather than assumed.

## Details

### The DC power flow (`6484fd3`, `35c3cc3`, earlier session)

- **The model,** `DcModel::new`: MATPOWER's `makeBdc` after `ext2int`'s eliminations (isolated
  buses, out-of-service generators and branches, and those touching an isolated bus). Branch
  stamps go through `from_triplets`, so parallel branches sum, and each branch's two
  off-diagonal stamps are pushed together, so Bbus is exactly symmetric, as `cholesky` requires.
  Bus roles follow `bustypes`, with several reference buses allowed.
- **The solve,** `DcModel::solve` and `dc_power_flow`: the injections and solve of `runpf` and
  `dcpf`, with the slack taken by the first in-service generator at each reference bus. Results
  are in case-table order, with 0.0 for out-of-service elements, as `runpf` reports them. B′ is
  solved by Cholesky, falling back to LU when it isn't positive definite, and
  `DcSolution::solver` records which.
- **Errors:** `ZeroReactance`, `NoReferenceBus`, `IslandWithoutReference` (found by a
  breadth-first search, not left to a factorization that rounding can fool), `SingularNetwork`
  and `CaseMismatch`.
- **Tests:** 27 on hand-built fixtures with hand-computed results, covering what the reference
  cases can't: phase shifters, shunt conductance, isolated buses, out-of-service elements,
  several reference buses, the LU fallback and the errors. Each solution's physics is checked
  too: branch equations, and power balance at every bus.

### The review, and its fix (`7036a35`)

- **Against MATPOWER's source:** the susceptance, tap and phase-shifter formulas (`makeBdc.m`
  lines 64–86), the 1e5 failure threshold (`dcpf.m` line 24), the slack (`runpf.m` lines
  224–229) and the zeroed results (lines 474–480) all match.
- **One problem:** `DcModel::new` trusted the case to be valid, but only `Case::parse`
  validates, and `Case`'s fields are public. A probe on case14:

| Edit | Before | After |
|---|---|---|
| A generator moved to bus 999 | `Ok`, with its 40 MW silently moved onto the slack | `UnknownBus` |
| Bus 14 relabelled 1 | `IslandWithoutReference` at bus 1, the reference bus | `DuplicateBus` |
| `baseMVA` of 0 | `SingularNetwork` | `InvalidBaseMva` |

- **The fix:** `DcModel::new` runs `Case::validate` first, as its doc comment now states. The
  new test fails without the fix.

### Comparison with MATPOWER (`eec0722`)

- **The tests:** `tests/dc_reference_test.rs` compares angles, flows and generator outputs for
  case14, case30 and case118 with the committed `rundcpf` results. The largest differences
  measured are 2.3e-13° and 3.6e-12 MW (case118), and the stated tolerance is 1e-9.
- **Shared code:** the reader and the comparison are in `tests/common/reference.rs`, so the
  survey uses them too. The comparison first checks that the results describe the case, and
  treats a NaN difference as a failure, which `f64::max` would skip. A test shows it rejects
  another case's results, a 1e-6 MW error and a NaN.

### Benchmarks (`0e31586`, `da2b907`)

- **Criterion,** `benches/dc_power_flow.rs`, times parsing, the model and the solve separately
  on case2383wp: 2.3 ms, 0.50 ms and 317–330 ms. Chosen because the reference data's README
  already names it as the benchmark case, it parses straight from MATPOWER's files, it has 6
  phase shifters, and its numbering suits an unordered factorization badly. Its Cholesky factor
  holds 40 times the entries of B′'s lower triangle, against 13 times for the similarly sized
  case_ACTIVSg2000, so Phase 3's gains will show.
- **The survey:** `examples/dc_survey.rs` and `scripts/benchmark-dc.sh` run every case in its
  own process, under a time and a memory limit. They compare its results with MATPOWER's,
  generated locally by `scripts/generate-reference.sh <dir> <case>...`, which are never
  committed. The record, `docs/benchmarks/dc-power-flow-2026-10-31.md`, found:
  - all 78 cases solve, in 18 minutes, within 1.1e-10 of MATPOWER relative to each quantity's
    largest value
  - 59 by Cholesky and 19 by LU, which all have negative series susceptances
  - run time follows fill-in, not size: unordered factors hold up to 109 times the entries
    they would without fill-in, and case_ACTIVSg10k solves in 1.2 s where the similarly sized
    case9241pegase takes 69 s
  - for Phase 3: a fill-reducing ordering is the main gain; LU costs about four times Cholesky
    at similar fill, and its pivot search visits every remaining row at every step; the failed
    Cholesky attempt before LU costs nothing
- **Octave versions:** with no arguments, `generate-reference.sh` regenerates the committed
  results byte for byte under Octave 8.4. That settles the handoff's open question about the
  6.4 they were made with.

### Backward error (`2659525`, `0b9eb18`)

- **Why:** the differences from MATPOWER grow with the case, from 8e-16 relative on case14 to
  1e-10 on case_SyntheticUSA, and a comparison alone can't tell a solver that loses accuracy
  from an ill-conditioned problem. The normwise backward error can: a backward-stable solve
  keeps it within a small multiple of ε whatever the conditioning. On nine cases measured with
  condition estimates, both solvers were backward stable, and every difference was 30 to 7,000
  times below what conditioning allows, so the growth comes from the problems.
- **The test,** `tests/dc_backward_error_test.rs`, requires at most 100 ε on four Cholesky and
  two LU cases, which measure 0.07–0.38 ε, and shows that a 1e-9° error in one angle
  (3,700 ε) and a NaN fail it. The measure itself, in `tests/common/accuracy.rs`, derives P′
  independently of `DcModel::solve`.
- **The survey** records our backward error and MATPOWER's for every case. Ours is at most
  10.7 ε (case2383wp) and MATPOWER's 1.4 ε, but ours is up to 30 times MATPOWER's on the cases
  whose factors fill in most, and it correlates with fill (Spearman 0.46) where MATPOWER's
  doesn't (−0.03).
- **Tracked:** the record states the guess that fill-in drives our backward error, with this as
  its baseline. Size rises with fill when nothing is ordered, so only an ordering can test it,
  and a new Phase 3 item asks for that comparison.

## Tests

- **`powerflow`:** 97 tests, from 62. `dc_model_test.rs` has 13 (the earlier session's 12, plus
  the validation test), `dc_solve_test.rs` 15, `dc_reference_test.rs` 4 and
  `dc_backward_error_test.rs` 3. The other 62 are unchanged: 23 data model, 18 locator, 4 fuzz,
  8 syntax and 9 normalized cases.
- **The benchmark** runs once as a test under `cargo test --all-targets`. In debug that takes
  about 6 s per pass locally and 10.5 s on CI's runners, until Phase 3 makes the solve fast.
- **`sparse_mat`:** unchanged.

## Verification

- **Locally:** `./check.sh` with `YBUS_REQUIRE_CASES=1`, the Rust 1.87 build and tests, and
  `cargo doc` with warnings as errors all pass. `cargo-deny` isn't installed locally; instead,
  `cargo metadata` shows that all 60 packages in the resolved graph, criterion's included, have
  licenses the allow list accepts.
- **The survey's limits:** the time limit was exercised by a first pass at 20 s per case, which
  9 cases hit; the error path on a case that doesn't exist; and the memory limit by hand, at
  256 MiB.
- **CI:** Check, MSRV and Dependencies pass on this pull request, on both its push and pull
  request runs, so `cargo-deny` accepts criterion's dependencies. The Check logs show the
  reference tests, the backward-error tests and the benchmark's test runs in all four passes,
  with nothing skipped. The MSRV job doesn't fetch MATPOWER, so there those tests return early,
  and its log shows the benchmark skipping.

## Mistakes caught

- **A misreported status:** in this session's first review, my wrapper around `check.sh`
  reported `tail`'s exit status rather than the script's. The log showed it had passed; later
  runs record the script's own status.
- **A check that tested nothing:** I ran the reference tests with `MATPOWER_DIR` empty to see
  them skip, but the locator still found the fetched cases. I didn't count it as verification;
  CI's MSRV job, which has no MATPOWER, exercises the skip.
- **Stale figures from the earlier session:** the checklist said 14 cases need LU, and gave 36 s
  and 166 s for case8387pegase and case13659pegase. The survey finds 19 of 78 cases, 15 s and
  83 s. The earlier count covered 65 cases, and its machine wasn't recorded. The checklist is
  corrected.
- **An overstated summary:** in conversation I summarized the agreement with MATPOWER as about
  ten digits everywhere. Small cases agree to about 15 digits and the largest to about 10,
  because the condition number grows with the case. The committed text, a bound, was right;
  checking the claim led to the backward-error work.
- **A local figure given as CI's:** the report first said the benchmark's test runs cost CI
  about 6 s per pass, which I'd measured locally. CI's logs show about 10.5 s, so I corrected it
  after the first CI run.
- **A survey rerun:** I first surveyed without factor sizes. Once the timings showed that
  fill-in drives the cost, I added them and reran the survey, so that the committed tool
  reproduces the record.

## Decide

- **Case facts in the record.** `docs/benchmarks/` names every case with its size, timings,
  factor size and agreement with MATPOWER, case2383wp included. I read the license decision
  ("nothing derived from case2383wp is committed") as covering case data and results, not facts
  like these, which reveal nothing about the grid. Say if you'd rather drop those rows.
- **CI time.** The benchmark's test-mode run takes about 10.5 s in each of the four test passes
  in CI's `check.sh` step, about 42 s per Check run, until Phase 3 makes the solve fast. With a
  warm cache, the `check.sh` step took 185–187 s on this pull request, against 81 s for this
  branch before this session's commits. I can't attribute the rest from these runs: there are
  new targets to compile in each pass, and runners vary (PR #13's two Check runs took 1m35s and
  2m30s). The obvious
  alternative case, case_ACTIVSg2000, solves in 0.44 s in a local debug build, against
  case2383wp's 5.0 s, but it would show less of Phase 3's gains.
- **`DcSolver::None`** hides `Option::None` wherever someone writes `use DcSolver::*`. A name
  such as `NoSolve` would avoid that. Not changed, since it's an API decision.

## Next

Phase 2 (a fixed-pattern `CsrMatrix`) and Phase 3 (a Gilbert–Peierls LU with orderings and
refactorization). The survey shows where Phase 3's gains are: a fill-reducing ordering for the
PEGASE and RTE cases, and a pivot search that only looks at nonzeros for the largest cases.
Once the orderings land, rerunning the survey tests whether fill-in drives the backward error.
