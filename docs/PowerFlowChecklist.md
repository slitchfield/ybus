# Power flow checklist

The roadmap for using `sparse_mat` as the numerical base of a DC and AC power flow library. Each
unchecked item is sized to be one commit. Items marked **Decide** need an answer before the work
that depends on them.

## Review notes: changes from the original plan

The plan was first written in conversation on 2026-10-28. Reviewing it against the code after
Phase 0 merged led to these changes:

1. **`num-complex` becomes a dependency in Phase 1, not Phase 4.** The complex helpers (`conj`,
   `real`, `imag`) need it, and so does the solver generalization, because `ComplexFloat` lives
   in `num_complex`, not in `num-traits`. That location was a correction from the original plan.
2. **The new LU is generic over `ComplexFloat` from the start.** Writing it for `Float` in
   Phase 3 and generalizing it in Phase 4 would mean writing the numeric core twice. Phase 4 now
   only migrates the existing solvers.
3. **Complex Cholesky and CG need Hermitian, not symmetric, semantics.** Complex Cholesky factors
   `A = L Lᴴ`, so the symmetry check must become a Hermitian check. CG's inner products need a
   conjugate. The original plan said only "rerun the tests for `Complex<f64>`".
4. **A DC power flow prototype moves up to right after Phase 1.** It only needs `select` (to drop
   the slack bus) and the existing `cholesky`. Building it early creates the case reader, the
   MATPOWER reference results and the benchmark harness before the LU rewrite, so Phase 3 is
   measured against real cases instead of synthetic ones.
5. **Open items from the PR #1 pre-push review are tracked here:** parser fuzzing, `cargo audit`,
   package metadata and the leftover template `add` function in `lib.rs`.

## Phase 0: Hardening

Merged in PR #2.

- [x] `from_csr` rejects `rows == usize::MAX` instead of overflowing.
- [x] Compressed-storage allocation guarded against huge dimensions (clear panic, and half the
      pointer memory).
- [x] The Matrix Market reader caps pre-allocation from the untrusted size line.
- [x] `check.sh` lints and tests each optional feature on its own.
- [x] MSRV declared as 1.87 (measured), with the MSRV-aware resolver and a CI job.

## Phase 1: Cheap primitives

- [x] **Decide:** whether `num-complex` is always a dependency or behind a `complex` feature.
      Decided: always a dependency, re-exported as `sparse_mat::num_complex` and
      `sparse_mat::Complex`.
- [x] `scale_rows(&[T])` and `scale_cols(&[T])` (diag(v)·A and A·diag(v)), with
      `ShapeMismatch` on a wrong length. Test against `from_diagonal(v) * A`.
- [x] `select(rows: &[usize], cols: &[usize])` for non-contiguous index sets. Rows and columns
      come out in the order given; `OutOfRange` on a bad index; repeated indices are allowed and
      repeat the row or column (decided and documented). Test that `select` with ranges matches
      `submatrix`.
- [x] `permute(row_perm, col_perm)` and `permute_symmetric(perm)`, validating that each is a
      permutation. Test that a permutation followed by its inverse is the identity.
- [x] Complex helpers: `conj()` (identity for real types via `ComplexFloat`), `adjoint()`,
      `real()`, `imag()`, and `from_parts(re, im)`. Test that `from_parts(a.real(), a.imag())`
      gives back `a`.

## Phase 1.5: DC power flow prototype (power flow crate)

- [x] **Decide:** a Cargo workspace in this repo (`sparse_mat/` and `powerflow/`) or a separate
      repo. Recommended: a workspace, so one CI run covers both. Converting means moving the
      crate into a subdirectory; `check.sh` already uses `--workspace`.
      Decided: a workspace. The repository is renamed `ybus`; the crates are `sparse_mat/` and
      `powerflow/`.
- [x] **Decide:** confirm the MATPOWER license terms before vendoring its case files as test data.
      Decided: don't vendor the case files. MATPOWER's code is BSD 3-clause, but its LICENSE
      says the case files are not covered by it, so `scripts/fetch-matpower.sh` downloads a
      pinned release (8.1, SHA-256 checked) into the gitignored `testdata/` instead. Reference
      results computed from case14, case30 and case118 are committed, since those cases come
      from the public UW test case archive and published papers. Nothing derived from
      case2383wp is committed: its data is distributed with the provider's permission, whose
      scope isn't stated.
- [x] Test-data locator (`powerflow/tests/common`): finds cases in `MATPOWER_DIR` or the
      fetched `testdata/`, skips when they're missing unless `YBUS_REQUIRE_CASES=1`, and fails
      on a broken installation. CI fetches MATPOWER with a cache keyed on the fetch script
      (version and checksum) and requires the cases.
- [x] Case model: buses (type, demand, shunts, voltage limits), branches (r, x, b, tap, shift,
      status) and generators, in per unit. Done in `powerflow/src/case.rs`, with generator costs
      and bus names too. Changed from the plan: the case stores the file's units (MW, MVAr,
      degrees, with units in the field names) rather than per unit, and the solvers convert, as
      MATPOWER does, so the case stays a faithful copy of the file. `Case::validate` checks bus
      labels and references.
- [x] MATPOWER case reader for the `mpc.bus`, `mpc.branch`, `mpc.gen` and `mpc.baseMVA` matrices
      in `.m` files. Done in `powerflow/src/case/matpower.rs` (version 2 format, plus `gencost`
      and `bus_name`). It reads only `mpc.<field> = <literal>;` statements and rejects any other
      code, because 24 of MATPOWER 8.1's data files rescale their data with code after the
      literals; ignoring that would load wrong values. 52 of its 84 data files parse.
- [x] Cases that compute their data with code: `scripts/normalize-cases.sh` evaluates every
      MATPOWER case with MATPOWER's own `loadcase` under Octave and writes it as literals at
      full precision (`%.17g`; MATPOWER's `savecase` would round to `%.9g`) into the gitignored
      `testdata/`. All 78 of MATPOWER 8.1's cases then parse; the other 6 data files aren't
      cases. Chosen over interpreting MATLAB in the reader, which would be a partial language
      implementation. Tests use `common::normalized_case`, and CI caches the output, installing
      Octave only on a cache miss. Along the way the reader learned to accept version 2 rows
      without the optional columns, filled in as MATPOWER's own conversion does.
- [x] Case reader hardening for user-written files (`dev/parser-hardening`):
      - **Encoding:** `Case::from_bytes` (used by `from_path`) ignores a UTF-8 byte order mark
        and falls back to Latin-1 for files that aren't UTF-8.
      - **Fuzzing:** deterministic fuzz tests (`matpower_fuzz_test.rs`) found a `gencost` length
        overflow panic, now fixed. Whole numbers are limited to 2^53 instead of saturating.
      - **Syntax:** several statements per line, double-quoted strings, `%{ %}`/`#{ #}` block
        comments and `#` comments, the struct named by the `function` line, and skipping any
        assignment to a field the model doesn't use (`matpower_syntax_test.rs`).
      - **Clearer errors:** version 1 files are reported as such. NaN in a read column, and a
        case with no buses, are errors.
      - **No regressions:** the same 52 MATPOWER files parse directly as before.
- [x] Bbus assembly with `from_triplets` (summing branch stamps), honoring branch status and
      phase shifters. `powerflow::dc::DcModel::new`, following MATPOWER's `ext2int`,
      `bustypes` and `makeBdc`. Fields are private, with accessors (decided). Beyond the plan:
      a case may have several reference buses, as in `bustypes`; and an island that can't
      reach a reference bus is found by a breadth-first search and reported by bus
      (`IslandWithoutReference`), rather than left to a singular factorization, which rounding
      can hide. Zero reactance is an error (`ZeroReactance`), where MATPOWER gives `Inf`.
      The case is validated first, since only `Case::parse` validates: in a review, a case
      edited after parsing with an unknown generator bus label solved with that generator
      silently dropped.
- [x] DC power flow: drop the slack row and column with `select`, solve B′θ = P with
      `cholesky`, then compute branch flows. `DcModel::solve` and `dc_power_flow`, following
      `runpf` and `dcpf`; results indexed in case-table order (decided). Out-of-service and
      disconnected generators and branches report 0.0, as `runpf.m` (lines 474–480) does,
      confirmed with `rundcpf`; isolated buses keep their angle and their load is not served.
      Cholesky falls back to LU when B′ isn't positive definite (negative reactances), which
      19 of MATPOWER's 78 cases need. **Found:** without a fill-reducing ordering the solve is
      slow on large cases (case13659pegase 83 s, case8387pegase 15 s in release; see
      `docs/benchmarks/`), which Phase 3's orderings address. First recorded as 14 cases (of
      65 surveyed), 166 s and 36 s, on a machine that wasn't recorded.
- [x] Reference results: MATPOWER `rundcpf` outputs for case14, case30 and case118 checked into
      test data. Tests compare angles and flows within a stated tolerance. The results are
      generated and committed (`powerflow/tests/data/matpower-reference/`, by
      `scripts/generate-reference.sh`); the comparison tests come with the DC power flow.
      Done in `powerflow/tests/dc_reference_test.rs`, which also compares generator outputs.
      Angles, flows and outputs agree within 1e-9 (degrees, MW); the largest differences
      measured are 2.3e-13° and 3.6e-12 MW, on case118. The comparison first checks that the
      results describe the case (buses, branch ends and statuses, generators), and a test shows
      it fails on another case's results, a 1e-6 MW error or a NaN.
- [x] Benchmark harness (criterion) running DC power flow on the largest available case, so
      later phases have a baseline. Changed from the plan (decided): criterion benchmarks a
      medium case, since sampling the largest repeatedly is impractical without a
      fill-reducing ordering (nine cases took over 20 s each in a first pass), and a survey
      times every case once.
      - [x] Criterion harness, `powerflow/benches/dc_power_flow.rs`, on case2383wp (2,383
            buses), the benchmark case named in the reference data's README. It times parsing,
            the model (which Phase 2 targets) and the solve (which Phase 3 targets) separately.
            Baseline on the maintainer's machine: parse 2.3 ms, model 0.50 ms, solve 317–330 ms
            (two runs; the machine varies by about 4% between runs). Its numbering suits an
            unordered factorization badly, so Phase 3's orderings will show: the Cholesky factor
            has 40 times the entries of B′'s lower triangle.
      - [x] Survey of every MATPOWER case, with the results recorded in
            `docs/benchmarks/dc-power-flow-2026-10-31.md`. `scripts/benchmark-dc.sh` runs the
            `dc_survey` example on each case in its own process, under a time and memory limit,
            and compares the results with MATPOWER's, generated locally by
            `scripts/generate-reference.sh`, which now takes an output directory and cases. All
            78 cases solve in 18 minutes and agree with MATPOWER to within 1.1e-10 of each
            quantity's largest value in its case, and every solve is backward stable (at most
            10.7 ε, checked on six cases by `dc_backward_error_test.rs`). Run time follows
            fill-in, not size: unordered factors hold up to 109 times the entries they would
            without fill-in, and LU, which 19 cases need, costs about four times Cholesky at
            similar fill.

## Phase 2: Fixed-pattern `CsrMatrix`

- [ ] **Decide:** whether `CsrMatrix` becomes the main type for numerical work, with
      `SparseMatrix` kept for editing and assembly. Recommended: yes.
- [ ] Promote the internal `Compressed<T>` to a public `CsrMatrix<T>` with an immutable pattern
      and `values()`/`values_mut()`. Build it from `SparseMatrix` and from validated raw arrays
      (reusing `from_csr`'s checks), and convert back.
- [ ] Zero-copy row access (`row(i)` returning column and value slices), `position(i, j)` by
      binary search, and cached diagonal positions.
- [ ] `mul_vec` on `CsrMatrix`, reusing `mul_vec_csr`, plus `par_mul_vec` under `rayon`.
- [ ] An assembly map built once from a triplet list's positions that records each triplet's
      slot in `values`, so a refill with new values is an O(nnz) write with no hashing. Test that
      a refill matches a fresh `from_triplets`.
- [ ] Benchmark: refilling Ybus and Jacobian-sized patterns through the assembly map versus
      `from_triplets`.

## Phase 3: Factorization rewrite (the largest piece)

- [ ] Reverse Cuthill–McKee ordering, with a test that bandwidth doesn't increase on banded and
      grid matrices.
- [ ] Symbolic analysis, `analyze(pattern)`: ordering, plus the elimination tree for Cholesky.
- [ ] Gilbert–Peierls left-looking LU with threshold partial pivoting, generic over
      `ComplexFloat` from the start (review note 2). Keep the current `Lu` API (`solve`,
      `determinant`, `l`, `u`, `permutation`) as the front end.
- [ ] `refactor(values)`, which reuses the analysis and pivot sequence when only values change.
      If a pivot becomes too small, fall back to a full factorization and report that.
- [ ] Keep the current `BTreeMap` LU as a test-only reference implementation, and compare the
      two on pseudo-random and power flow matrices.
- [ ] `solve_in_place` and `solve_many` (several right-hand sides).
- [ ] AMD ordering for Cholesky, and COLAMD (or AMD on the pattern of AᵀA) for LU. Benchmark fill
      and time against RCM on the large case.
- [ ] Test the hypothesis from the Phase 1.5 survey that fill-in drives the DC solve's backward
      error: unordered, ours was up to 30 times MATPOWER's on high-fill cases, and correlated
      with fill where MATPOWER's, from an ordered solve, didn't. Rerun `scripts/benchmark-dc.sh`
      with the orderings, and compare `backward_error` with `factor_entries` and
      `matpower_backward_error` against the baseline in
      `docs/benchmarks/dc-power-flow-2026-10-31.md`. Threshold partial pivoting trades some
      stability for sparsity, so `dc_backward_error_test.rs` guards its limit meanwhile.
- [ ] Cholesky after ordering, reusing the existing up-looking algorithm.

## Phase 4: Complex-capable solvers (migrating what Phase 3 didn't rewrite)

- [ ] Change the bounds of the triangular solves, `solve` and the norms from `Float` to
      `Scalar + ComplexFloat`. Pivot choice compares `abs()` (a real value), and norms return
      `Real`.
- [ ] Hermitian Cholesky (`A = L Lᴴ`): replace the symmetry check with a Hermitian check, and
      take real square roots of the diagonal (review note 3).
- [ ] Hermitian conjugate gradient: use conjugated inner products.
- [ ] Run every solver test for both `f64` and `Complex<f64>`.
- [ ] Note in the docs that user-defined `Float` types lose solver support, since `ComplexFloat`
      is sealed.

## Phase 4.5: AC power flow (power flow crate)

- [ ] Ybus assembly with `from_triplets` (complex branch π-models, taps, shifts, shunts), then
      converted to `CsrMatrix` with an assembly map for re-stamping.
- [ ] Newton–Raphson in polar form: mismatch, then the dS/dV Jacobian built with
      `scale_rows`/`scale_cols`, `conj` and `select` over the PV/PQ index sets, then a solve
      using `refactor` from the second iteration on.
- [ ] Convergence control: mismatch tolerance, iteration limit, and reported iteration count and
      final mismatch.
- [ ] Fast-decoupled power flow (XB and BX variants), which factors B′ and B″ once.
- [ ] Validation against MATPOWER `runpf` for case14, case30, case118 and a large case such as
      case2383wp: voltage magnitudes, angles and branch flows within a stated tolerance.
- [ ] Benchmarks against the Phase 1.5 baseline (`docs/benchmarks/`; rerun its tools on the same
      machine to compare).

## Phase 5: DC studies

- [ ] Sparse right-hand-side solves, using the Gilbert–Peierls reach to touch only the needed
      part of the factor.
- [ ] PTDF and LODF computation using `solve_many` or sparse right-hand sides.
- [ ] Rank-k updates (Woodbury) for N-1 contingencies without refactoring. Test against a full
      refactorization.
- [ ] Island detection: connected components on the pattern, with an error naming the islands
      instead of `SingularMatrix`.
- [ ] Iterative refinement and a condition estimate.
- [ ] Parallel batch solves over a shared factorization (`rayon`). Check that `Lu` is
      `Send + Sync` and that results are deterministic.

## Housekeeping (carried over from the PR #1 review)

- [x] Fuzz the Matrix Market parser with truncated, malformed and oversized input, and confirm
      it always errors instead of panicking. Done as deterministic tests
      (`tests/matrix_market_fuzz_test.rs`); fixed a skew-symmetric negation overflow.
- [x] **Decide:** integer overflow when summing duplicate entries (`from_triplets`, and so the
      Matrix Market reader) panics in debug builds. Keep integer `+` semantics (documented), or
      add a checked-add bound so it can return an error. Decided: keep integer `+` semantics.
- [x] Run `cargo audit` or `cargo deny` over the dependency tree, and consider adding it to CI.
      `cargo deny` (advisories, licenses, bans, sources) is configured in `deny.toml` and runs
      in a CI job.
- [x] License: dual MIT OR Apache-2.0 (`LICENSE-MIT`, `LICENSE-APACHE`, and `license` in
      `Cargo.toml`).
- [x] Add the remaining package metadata (`description`, `repository`) and check
      `cargo package --list`, which currently warns that they are missing.
- [x] Remove the leftover `add(left, right)` template function from `src/lib.rs`.
- [x] Decide whether to doctest the README by including it as crate docs
      (`#![doc = include_str!("../README.md")]`), so its example can't drift. Decided: doctest
      it through a `#[cfg(doctest)]` item in `lib.rs` instead, which keeps the crate docs
      separate from the README.
