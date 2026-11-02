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

## Phase 2: Fixed-pattern compressed matrices

- [x] **Decide:** whether `CsrMatrix` becomes the main type for numerical work, with
      `SparseMatrix` kept for editing and assembly. Recommended: yes.
      Decided (2026-10-31):
      - **Two types:** `SparseMatrix` builds and edits; a compressed type does the numerical
        work.
      - **The numerical methods move** to the compressed type, rather than staying on
        `SparseMatrix` as forwarding methods. Breaking changes are acceptable at 0.1.
      - **One compressed type, oriented either way, with the orientation in the type:**
        `CompressedMatrix<T, O>`, with aliases `CsrMatrix<T>` and `CscMatrix<T>`. Each algorithm
        states the orientation it needs, so a mismatch fails to compile and conversions are
        always explicit. Chosen over a runtime orientation flag, which would hide conversions
        inside algorithms.
      - **Stored zeros are kept:** the pattern is fixed when the matrix is built, so refills,
        refactorization and Newton's flat start (where some Jacobian entries are exactly zero)
        keep one pattern. Equality is mathematical, ignoring stored zeros and orientation.
- [x] `CompressedMatrix<T, O>`, promoted from the internal `Compressed<T>`: a fixed pattern with
      `values()`/`values_mut()` and stored zeros, built from `SparseMatrix`, from triplets and
      from validated raw arrays (reusing `from_csr`'s checks), and converted back. Transposing
      swaps the orientation in O(1); reorienting is an explicit O(nnz) conversion.
      Done in `sparse_mat/src/compressed_matrix.rs`, with mathematical equality (across
      orientations and with `SparseMatrix`), formatting and `serde`. `SparseError::InvalidCsr`
      became `InvalidCompressed`, since raw arrays can now be either orientation, and `from_csr`
      validates through `CsrMatrix::try_from_raw_parts`.
- [x] Move the numerical methods from `SparseMatrix` to the compressed type, each on the
      orientation it needs, one group at a time:
      - [x] products, including `mul_vec` in both orientations and the `rayon` versions. A dot
            kernel and a scatter kernel serve both orientations, and every result entry adds its
            terms in the same order either way, so CSR and CSC results are bit-for-bit equal.
            `matmul` keeps cancellations as stored zeros, so its pattern doesn't depend on the
            values. `par_mul_vec` is on CSR and `par_vec_mul` on CSC, where each splits into
            independent dot products.
      - [x] solvers, with factors returned as compressed matrices. `lu`, `solve`, `cholesky`
            and `solve_spd` are on `CsrMatrix`, since they work by rows, and build their factors
            straight into CSR arrays instead of through `HashMap`s. The triangular solves and
            conjugate gradient work in either orientation. Changed from the plan: `solve_spd`
            back-substitutes with `Lᵀ` reoriented to CSR, an O(nnz) step, rather than the O(1)
            transpose with a column-oriented solve, which adds terms in another order. That keeps
            every result bit-for-bit what it was: the DC survey's factor sizes, backward errors
            and differences from MATPOWER match the recorded ones exactly on case14, case300,
            case2383wp and case2848rte.
      - [x] access, selection, permutation and properties. `row` and `row_iter` are on CSR,
            `col` and `col_iter` on CSC, where each is a slice; `diagonal`, `submatrix`,
            `select`, the permutations, `density`, `sum`, `trace`, `is_symmetric` and the norms
            work in either orientation. Stored zeros are carried along, count in `density`,
            and don't break symmetry.
      - [x] arithmetic, block construction and complex parts. Patterns follow from the
            operation: sums and differences store the union of their operands' patterns,
            `hadamard` their intersection, and scaling, negation, `map`, `apply`, row and column
            scaling and complex parts keep the pattern. `adjoint`, like `transpose`, returns the
            other orientation. Stacking along the compressed axis concatenates the arrays.
- [x] `powerflow` builds and stores its DC model's matrices as `CsrMatrix`. Every stamped
      position stays stored, so Bbus's pattern depends only on the topology. On the 70 cases
      that solve in under 30 s, every survey column but time and memory (sizes, solver, B′ and
      factor entries, backward errors and differences from MATPOWER) matches the recorded
      survey exactly.
- [x] `SparseMatrix` becomes the editing type only: its compressed caches, the cache API and the
      moved methods are removed, and the README describes both types. It keeps construction,
      element edits and reads, `iter` (which now sorts on demand, so it works for any shape),
      transposes, `to_dense`, Matrix Market I/O, `serde`, and `to_csr`/`to_csc`. `from_csr` is
      gone in favour of `CsrMatrix::try_from_raw_parts`.
- [x] Zero-copy major-axis slices (a CSR row or a CSC column), `position(i, j)` by binary search,
      and cached diagonal positions. `outer(i)` borrows a row (CSR) or column (CSC) as its
      index and value slices, and `outer_mut(i)` makes the values mutable for refills.
      `position(row, col)` returns the entry's slot in `values()`, or `None` if it isn't
      stored. `diagonal_positions()` computes the diagonal's slots on first use and caches
      them in a `OnceLock`; a fixed pattern means the cache never goes stale, so values can
      change freely, and `transpose` keeps it, since a diagonal entry keeps its slot.
- [x] An assembly map built once from a triplet list's positions that records each triplet's
      slot in `values`, so a refill with new values is an O(nnz) write with no hashing. Test that
      a refill matches a fresh `from_triplets`.
      `CompressedMatrix::from_triplets_mapped` returns the matrix and an `AssemblyMap<O>`, and
      `refill(&map, &values)` takes one value per triplet. Each position's first triplet is
      written and the rest added in input order, so a refill equals `from_triplets` bit for bit,
      down to the sign of a zero (tested on 200 inexact triplets with repeats, in both
      orientations). The map carries the orientation in its type, so it can't refill a matrix of
      the other orientation, and checks the shape and stored count (`PatternMismatch`) and the
      number of values (`LengthMismatch`, both new errors) before writing anything.
- [x] Benchmark: refilling Ybus and Jacobian-sized patterns through the assembly map versus
      `from_triplets`.
      `powerflow/benches/assembly.rs`, on case2383wp, on the maintainer's machine (medians):

      | Pattern | Size, triplets | `from_triplets` | `refill` | Speedup |
      |---|---|---|---|---|
      | B stamps (Ybus-shaped) | 2,383², 11,584 | 69.0 µs | 5.39 µs | 12.8× |
      | 2×2 blocks (Jacobian-shaped) | 4,766², 46,336 | 309.0 µs | 34.9 µs | 8.9× |

      The block pattern expands each B stamp into the four quadrants, as the polar Jacobian's
      blocks share Ybus's pattern. Its refill costs more per triplet: its map and values
      (about 1.1 MB) outgrow the cache, so the scattered writes miss. Routing `from_triplets`
      through the shared core left the DC model step unchanged (382.6 µs, against 387.7 µs
      before).

## Phase 3: Factorization rewrite (the largest piece)

Done in PRs #19 to #22 (2026-10-31 to 2026-11-01). Ordered Cholesky and the Gilbert–Peierls LU
took the 78 MATPOWER cases' DC solves from 496 s to 0.35 s together, and brought every case's
backward error below 0.8 ε. COLAMD stays deferred, below.

- [x] **Decide:** the shape of the ordered factorizations, where AMD comes from, whether to
      implement COLAMD, and the LU's pivot threshold.
      Decided (2026-10-31):
      - **Symbolic/numeric split:** `analyze_cholesky`/`analyze_lu` return a symbolic object
        (ordering, elimination tree, counts), whose `factor(&A)` returns `Cholesky`/`Lu` factor
        types that carry their permutations. `Lu::refactor` reuses a factorization's pivots and
        patterns. `cholesky()`, `lu()`, `solve` and `solve_spd` stay as one-call front ends,
        ordered by default, and return the new types; breaking changes are acceptable.
      - **AMD is written from the papers** (Amestoy, Davis and Duff, 1996; Davis, 2006, ch. 7),
        not ported, so `sparse_mat` stays written from published descriptions only.
      - **COLAMD is deferred.** LU orders with AMD on the pattern of A+Aᵀ, which suits every
        matrix power flow factors (B′, B″, Ybus and the polar Jacobian are structurally
        symmetric), as KLU does by default.
      - **The pivot threshold defaults to 0.1** and is configurable: a pivot on the diagonal is
        kept if its magnitude is at least 0.1 times the column's largest.
- [x] Reverse Cuthill–McKee ordering, with a test that bandwidth doesn't increase on banded and
      grid matrices.
      `CompressedMatrix::ordering(Ordering::Rcm)` in `compressed_matrix/ordering.rs`, on the
      pattern of A+Aᵀ, so it's the same in either orientation and for any values. Each component
      starts from a pseudo-peripheral node (George and Liu). On shuffled k×k grids it recovers
      bandwidth k exactly, with a smaller envelope than the natural order (10,700 against
      15,024 at k = 25). `bandwidth()` and `profile()` measure it.
- [x] Symbolic analysis, `analyze_cholesky`: ordering, the elimination tree and L's column
      counts.
      `CsrMatrix::analyze_cholesky(Ordering)` returns a `SymbolicCholesky` holding the
      permutation, the elimination tree (Liu's algorithm with path compression), a postorder,
      and the exact pattern of L, found as each row's subtree of the tree, from which come the
      column counts and `predicted_entries()`. It reads only the pattern of A+Aᵀ, stored zeros
      included. Tested against dense symbolic elimination on random patterns, under both
      orderings.
- [x] Cholesky after ordering, reusing the existing up-looking algorithm, with the elimination
      tree's reach in place of its `BTreeMap` work row. Moved up from the end of the phase: it
      costs little once an ordering exists, and gives the first measurement.
      `SymbolicCholesky::factor` returns a `Cholesky` (`l`, `permutation`, `solve`) of P A Pᵀ,
      filling the analyzed pattern with a dense work vector, so the factor stores its whole
      symbolic pattern (numerical zeros included) and its memory is known in advance. Each row
      still takes its entries in increasing column order, so under `Ordering::Natural` the
      factor and solve match the old algorithm bit for bit, apart from the sign of a zero
      (tested against it, kept in `tests/common`). `cholesky()` and `solve_spd` order by RCM
      until AMD arrives, and `NotPositiveDefinite` names the failing row of A in its own
      numbering. The factor keeps Lᵀ in CSR form as well, for a back substitution in the same
      order as before, at twice the factor's memory.
- [x] The DC power flow solves with ordered Cholesky, and the survey can choose its ordering.
      Rerun the survey with RCM and make a first check of the hypothesis below on the Cholesky
      cases.
      `DcModel::solve` orders B′ by RCM, and `solve_with_ordering` takes any ordering; the
      survey's `--ordering` (and the script's `YBUS_BENCH_ORDERING`) chooses it, and a new
      `ordering` column records it. Recorded in `docs/benchmarks/dc-power-flow-rcm-2026-10-31.md`:
      - **Check:** ordered naturally, the 70 cases that solve in under 30 s match the baseline
        in every column but time and memory, backward errors included.
      - **Speed:** the 59 Cholesky cases solve in 19 ms together instead of 17.5 s;
        case2383wp's solve takes 1.56 ms in the criterion benchmark instead of 292 ms. No case
        fills more, and the 19 LU cases are unchanged.
      - **Hypothesis, first check:** on the Cholesky cases, the median ratio of our backward
        error to MATPOWER's falls from 1.23 to 1.00 and its rank correlation with fill from 0.39
        to 0.07; case2383wp's falls from 10.7 ε to 0.72 ε (MATPOWER's: 0.65 ε). But
        case8387pegase's barely moves (2.44 ε to 2.35 ε, MATPOWER's 0.94 ε) though its fill
        falls 15.5 times, so fill explains the large excesses, not the whole gap.
- [x] AMD ordering, for Cholesky and (on the pattern of A+Aᵀ) for LU. Benchmark fill and time
      against RCM on the large case.
      `Ordering::Amd` (`compressed_matrix/ordering.rs`), written from the paper and the book:
      quotient graph, approximate degrees, aggressive absorption, supervariables found by
      hashing, and rows with more than max(16, 10√n) neighbours ordered last. Tested for no
      fill on random trees and arrowheads, dense rows last, and less fill than RCM on shuffled
      grids: 184,412 fill entries against 651,750 on a 100×100 grid, ordered in 7.2 ms
      (119 ms for 300×300).
      `cholesky()`, `solve_spd` and `DcModel::solve` now order by AMD. Recorded in
      `docs/benchmarks/dc-power-flow-amd-2026-11-01.md`, on the 59 Cholesky cases against RCM:
      - **Fill:** never more than RCM's; the large cases hold 2.9–4.7 times fewer entries, and
        the worst fill falls from 11.2× to 2.4×.
      - **Time:** a wash for one solve (19.0 ms against 18.6 ms over the 59 cases;
        case2383wp 1.36 ms against 1.52 ms in criterion). AMD's ordering costs about three times
        RCM's, but its factor is about twice as fast and its solve two to four times, so it pays
        off wherever an analysis is reused.
      - **Backward error:** see the hypothesis item below.
- [x] Keep the current `BTreeMap` LU as a test-only reference implementation, and compare the
      two on pseudo-random and power flow matrices.
      `reference_lu` in `sparse_mat/tests/common`, rewritten against the public API, so no test
      code ships in the library. It matches the library's LU bit for bit (factors, permutation
      and determinant) before the new LU replaces it; the comparisons with the new LU come with
      it.
- [x] Gilbert–Peierls left-looking LU with threshold partial pivoting, generic over
      `ComplexFloat` from the start (review note 2). Keep the current `Lu` API (`solve`,
      `determinant`, `l`, `u`, `permutation`) as the front end.
      `CscMatrix::analyze_lu(ordering)` and `SymbolicLu::factor(&A, LuOptions)` in
      `compressed_matrix/lu.rs`; `lu()` and `solve` (AMD, threshold 0.1) moved from `CsrMatrix`
      to `CscMatrix`, since the algorithm works by columns. `Lu` keeps `solve`, `determinant`,
      `l`, `u` (now CSC) and `permutation`, and adds `column_permutation`: P A Q = L U.
      - **Pivoting:** the diagonal of the ordered matrix is kept while its magnitude is at least
        the threshold times the column's largest; otherwise the largest is taken. Checked: 0.1
        is UMFPACK's default partial pivot tolerance (its separate 0.001 applies only to its
        symmetric strategy's diagonal pivots).
      - **Tested:** P A Q = L U under every ordering and thresholds 0.001, 0.1 and 1.0; solves
        matching the reference LU; at threshold 1.0 with no ordering, the same pivots as the
        reference; complex matrices (a Ybus-like grid, and an unsymmetric one needing
        pivoting, with its determinant); the determinant's sign from both permutations.
      - **Errors:** `SingularMatrix` names the column of A with no pivot, which the DC solve
        maps to a bus, as before.
- [x] `refactor(values)`, which reuses the analysis and pivot sequence when only values change.
      If a pivot becomes too small, fall back to a full factorization and report that.
      `Lu::refactor(&A)` reuses the column ordering, the pivot sequence and the patterns of L
      and U: each column is a sparse triangular solve over U's known pattern, with no search.
      Each reused pivot is checked against the factorization's threshold; if one fails (or is
      zero), it factors afresh under the same column ordering and returns
      `RefactorOutcome::Repivoted`, otherwise `RefactorOutcome::Reused` (named `Refactor` until
      the PR #21 review asked for a name that reads as a call's result). A different pattern is
      `PatternMismatch`, and on any error the factorization is unchanged. Tested on new values
      (same pivots and patterns, reconstructing P A Q, solving and the determinant), on pivots
      that shrink or vanish, on complex values, and on errors.
- [x] The DC power flow's LU fallback uses the new LU.
      With the same ordering as Cholesky (AMD by default). Recorded in
      `docs/benchmarks/dc-power-flow-lu-2026-11-01.md`, on the 19 LU cases:
      - **Speed:** their solves take 0.33 s together instead of 478.5 s (case13659pegase 15 ms
        instead of 83 s). Unordered, the new LU fills about as much as the old one but runs
        26 times faster at the median; AMD gives another factor of 40, cutting the worst fill
        from 109× to 2.5×. All 78 cases now solve in 0.35 s together, against 496 s.
      - **Agreement with MATPOWER:** every case within 2.5e-12 of each quantity's largest value,
        against 1.1e-10 before.
- [x] `solve_in_place` and `solve_many` (several right-hand sides).
      On both `Cholesky` and `Lu`. `solve_in_place(&mut rhs)` overwrites its argument with the
      solution (permuting still takes one scratch vector), and `solve_many(&[b])` returns one
      solution per right-hand side, sharing one scratch vector. `solve` now runs through the same
      code, so all three agree bit for bit (tested, with complex right-hand sides for LU). A
      wrong length is `ShapeMismatch` and leaves an in-place argument untouched.
- [x] Test the hypothesis from the Phase 1.5 survey that fill-in drives the DC solve's backward
      error: unordered, ours was up to 30 times MATPOWER's on high-fill cases, and correlated
      with fill where MATPOWER's, from an ordered solve, didn't. Rerun `scripts/benchmark-dc.sh`
      with the orderings, and compare `backward_error` with `factor_entries` and
      `matpower_backward_error` against the baseline in
      `docs/benchmarks/dc-power-flow-2026-10-31.md`. Threshold partial pivoting trades some
      stability for sparsity, so `dc_backward_error_test.rs` guards its limit meanwhile.
      **Cholesky verdict (2026-11-01): supported.** Ordered by AMD, the 59 Cholesky cases' median
      ratio of our backward error to MATPOWER's is 1.00, the largest backward error is 0.78 ε
      (10.7 ε unordered), and its rank correlation with fill is −0.03 (0.45 unordered).
      case8387pegase, which RCM left at 2.35 ε, falls to 0.71 ε, below MATPOWER's 0.94 ε.
      **LU verdict (2026-11-01): supported.** The same new LU, unordered, gives the old excess
      back (median ratio 2.57, largest 28); ordered by AMD, the median ratio is 1.23, the
      largest 2.7, and the largest backward error 0.49 ε (6.5 ε before). So the excess followed
      the fill, not the algorithm, and threshold pivoting at 0.1 cost no stability here. Every
      case's backward error is now below 0.8 ε.
- [ ] Deferred: COLAMD, for LU of a matrix whose pattern is far from symmetric. Nothing in power
      flow needs it yet.

## Phase 4: Complex-capable solvers (migrating what Phase 3 didn't rewrite)

Done in PR #23 (2026-11-01). Every solver and norm now works for real and complex element
types through one code path, and every real result is bit for bit what it was.

- [x] **Decide:** what Cholesky reports for a matrix that isn't Hermitian.
      Decided (2026-11-01): keep `SparseError::NotSymmetric`, documented as "not symmetric, or
      for complex element types, not Hermitian (A ≠ Aᴴ)", so nothing breaks for real callers
      such as the DC solve's LU fallback. Every real result stays bit for bit what it was: for
      `f64`, `conj` is the identity and re² + im² is x·x + 0, so each change keeps the same
      operations in the same order.
- [x] Change the bounds of the triangular solves, conjugate gradient and the norms from `Float`
      to `Scalar + ComplexFloat`, with norms returning `Real`. (`solve`, the general one, is the
      LU's, which Phase 3 already made generic.)
      Done for the triangular solves and the norms; conjugate gradient changes with its
      Hermitian item below. The norms return `T::Real`: Frobenius sums re² + im² (exactly x·x for
      a real x), and the 1- and infinity norms sum `abs()`. New `is_hermitian()` (A equals Aᴴ
      mathematically; `is_symmetric` for real types). `check_system` needs no numeric bound, so it
      moved to an impl bound only by `Scalar`. Tested with complex triangular solves in both
      orientations, complex norms against hand-computed values, and `is_hermitian` on Hermitian,
      symmetric-but-not-Hermitian and real matrices.
- [x] Hermitian Cholesky (`A = L Lᴴ`): replace the symmetry check with a Hermitian check, and
      take real square roots of the diagonal (review note 3).
      `SymbolicCholesky::factor`, `Cholesky`, `cholesky()` and `solve_spd` now take any
      `ComplexFloat` element type. The check is `is_hermitian()` (still `NotSymmetric`); each
      update subtracts conj(l_mj)·l_ij; the remaining pivot is real (`T::Real`), losing
      re² + im² of each entry; and the stored back-substitution factor is Lᴴ. For real types
      these are the same operations as before, and the natural-order check against the old
      Cholesky still passes bit for bit. Tested on a Hermitian (not symmetric) grid matrix under
      every ordering (P A Pᵀ = L Lᴴ, a real positive diagonal, solves agreeing with LU), and on
      a complex symmetric matrix (`NotSymmetric`) and a Hermitian indefinite one.
- [x] Hermitian conjugate gradient: use conjugated inner products.
      `conjugate_gradient` takes any `ComplexFloat` element type, with ⟨u, v⟩ = Σ conj(uᵢ) vᵢ.
      The residual's squared norm, the curvature and the step lengths are real (`T::Real`), and
      `tolerance` is `T::Real` (still `f64` for `f64` callers). A real solution's bits, recorded
      before the change, are unchanged (a permanent test). Tested on a Hermitian grid matrix,
      matching the Cholesky solution, with CSR and CSC agreeing bit for bit.
- [x] Run every solver test for both `f64` and `Complex<f64>`.
      Done as one generic suite in `factorization_test.rs`, run for each type: triangular solves
      (both orientations), Cholesky under every ordering with `solve_many` and `solve_in_place`,
      LU under every ordering with a refactor, conjugate gradient against Cholesky, and the
      norms against sums taken from the dense matrix. A small `Element` trait lifts the real
      fixtures (off-diagonals `v ↦ v·(1 + 0.5i)`, mirrored conjugated where Hermitian is
      needed). Changed from the plan: the existing real-only tests stay as they are rather than
      being rewritten generically, since many check exact values or bit-for-bit agreement with
      the old algorithms, which only make sense for `f64`.
- [x] Note in the docs which float types keep solver support. Corrected while planning:
      `ComplexFloat` is sealed, but its seal is implemented for every `Float + FloatConst` type
      (`num-complex` 0.4.6, `complex_float.rs`), so a user-defined float type keeps solver support
      if it also implements `FloatConst`; only a type implementing `Float` alone loses it.
      Done in the Element types section of `sparse_mat/README.md`, which also lists the
      complex-capable solvers and norms.

## Phase 4.5: AC power flow (power flow crate)

- [x] **Decide:** how to split the phase, whether to enforce reactive power limits, how strictly
      to validate against MATPOWER, and how to build the Jacobian.
      Decided (2026-11-01):
      - **Four PRs:** the AC model and Newton–Raphson with validation; fast-decoupled XB and BX;
        an AC survey of every case and benchmarks; then a review of the references, attribution
        and `AI-DISCLAIMER.md`.
      - **Reactive power limits are off,** as in `runpf`'s default (`pf.enforce_q_lims = 0`);
        enforcing them is a deferred item below.
      - **Validation is strict:** the same iteration count as MATPOWER, and agreement within 1e-9
        in every result, as for the DC power flow. A case that can't meet it is investigated,
        then its measured gap and tolerance are recorded.
      - **The Jacobian is built by matrix operations** each iteration (`scale_rows`/`scale_cols`,
        `conj`, `select`, stacking), as MATPOWER's `dSbus_dV` does. Its pattern stays fixed, so
        the LU refactors. An assembly-map Jacobian is a follow-up only if the benchmark shows the
        build matters.
      Like the DC power flow, the AC code reproduces MATPOWER 8.1's `runpf` with its default
      options, following `ext2int`, `bustypes`, `makeYbus`, `makeSbus`, `newtonpf`, `dSbus_dV`,
      `makeB`, `fdpf` and `pfsoln`.
- [x] Extract the topology `DcModel` derives (connected buses, in-service elements, bus types,
      the island check) into a module both power flows share, with no change to the DC results.
      `powerflow/src/topology.rs` (crate-private). Each model passes its own branch check
      (`x = 0` for DC), run as branches are found, so errors keep their precedence. The DC tests
      pass unchanged, and a survey rerun matches the AMD and LU records in every column but
      time and memory, on all 78 cases.
- [x] Ybus assembly with `from_triplets` (complex branch π-models, taps, shifts, shunts), then
      converted to `CsrMatrix` with an assembly map for re-stamping.
      `AcModel::new` in `powerflow/src/ac.rs` builds Ybus, Yf and Yt as
      `CsrMatrix<Complex<f64>>` with `from_triplets_mapped`, keeping each assembly map, and
      `AcModel::restamp` re-stamps changed impedances, charging, taps, phase shifts and shunts
      through them; the result equals a fresh build bit for bit (tested). `sbus_pu` gives
      `makeSbus`'s injections. A zero-impedance in-service branch is the new
      `PowerflowError::ZeroImpedance`. Tested against the π-model computed by hand, rows summing
      to the shunts, symmetry, and dropped elements. Corrected while testing: the plan expected
      Ybus to be Hermitian without phase shifts, but it's symmetric (Yft = Ytf for a real tap),
      and never Hermitian.
- [x] MATPOWER's AC reference results (`runpf`, Newton–Raphson) for case14, case30 and case118,
      committed beside the DC ones.
      `scripts/generate_reference.m` also runs `runpf` with its defaults and writes
      `<case>_acpf_{bus,branch,gen,summary}.csv` (voltages, all four branch flows, generator P
      and Q, and the iteration count: 2, 3 and 3). Regenerated under Octave 8.4.0, the DC files
      came out byte-for-byte unchanged; `VERSIONS` now names 8.4.0. `tests/common/reference.rs`
      reads them (`read_ac`), and the reference README documents the files.
- [x] Newton–Raphson in polar form: mismatch, then the dS/dV Jacobian built with
      `scale_rows`/`scale_cols`, `conj` and `select` over the PV/PQ index sets, then a solve
      using `refactor` from the second iteration on. Convergence control: mismatch tolerance,
      iteration limit, and reported iteration count and final mismatch. `pfsoln`'s
      post-processing (generator reactive power, the slack, branch flows). Validated against
      the reference results within 1e-9, with MATPOWER's iteration counts.
      `AcModel::solve(&case, AcOptions)` and `ac_power_flow` in `powerflow/src/ac.rs`, following
      `runpf`'s start, `newtonpf`, `dSbus_dV` and `pfsoln`. The Jacobian is factored by
      AMD-ordered LU once and refactored afterwards; `AcSolution::full_factorizations` counts
      fresh factorizations (1 on the reference cases). Non-convergence is the new
      `PowerflowError::NoConvergence { iterations, mismatch_pu }`, and `SingularNetwork` and
      `CaseMismatch` now describe either power flow.
      - **Validation:** case14, case30 and case118 take MATPOWER's 2, 3 and 3 iterations and
        agree within 6.3e-12 in every result (case118's flows the largest).
      - **Found while validating:** MATPOWER 8.1's `runpf` uses its MP-Core implementation by
        default under Octave too (`have_feature('mp_core')` is true), so the reference results
        come from MP-Core, while this code follows the legacy functions (`newtonpf`, `pfsoln`).
        On case118 the two agree within 6.4e-12, with the same iteration count, measured with
        `exp.use_legacy_core`.
      - **Tests:** power balance at every bus and each branch's π-model flows; dropped elements;
        non-convergence; zero iterations from a solved start; a flat start (agreeing to the
        tolerance, 1.3e-8 MW measured); and `pfsoln`'s reactive sharing between generators at a
        bus, with finite and infinite limits.
- [x] Fast-decoupled power flow (XB and BX variants), which factors B′ and B″ once, validated
      the same way.
      `AcOptions::method` (`AcMethod::Newton`, `FastDecoupledXb`, `FastDecoupledBx`), with
      `AcOptions::for_method` giving MATPOWER's defaults (30 iterations for fast-decoupled).
      `makeB`'s B′ and B″ are factored once each, then `fdpf`'s alternating P and Q
      half-iterations run, with mismatches divided by `Vm`. All six runs (case14, case30,
      case118, XB and BX) take MATPOWER's iteration counts and agree within 7.8e-12. A branch
      with `x = 0` is `ZeroReactance` here, where MATPOWER would divide by zero.
      Changed from the plan:
      - **LU, not Cholesky:** B′ keeps the phase shifts (only taps, charging and shunts are
        removed), so it's unsymmetric whenever a phase shifter is in service. Both matrices are
        factored by AMD-ordered LU, as MATPOWER does.
      - **References from the legacy code:** against MATPOWER's default MP-Core results, case14
        FDXB differed by 2.3e-7 at the same iteration count. MATPOWER's own MP-Core and legacy
        fast-decoupled results differ by exactly that (and by up to 9.7e-9 on the other runs),
        while this code follows the legacy `fdpf`. So the fast-decoupled references are generated
        with `exp.use_legacy_core`, and the tolerance stays 1e-9.
- [x] An AC survey of every MATPOWER case against MATPOWER's `runpf`, recorded in
      `docs/benchmarks/`, and criterion benchmarks on case2383wp, including the Jacobian build's
      share of the Newton–Raphson time.
      Recorded in `docs/benchmarks/ac-power-flow-2026-11-01.md` (`ac_survey`,
      `scripts/benchmark-ac.sh`, and the `ac_power_flow` criterion bench):
      - **Convergence:** NR converges on 77 cases, FDXB on 67 and FDBX on 66: exactly the cases
        MATPOWER converges on, each in MATPOWER's number of iterations.
      - **Agreement:** within 2.3e-9 (NR), 4.0e-9 (FDXB) and 1.9e-10 (FDBX) of each quantity's
        largest value. In 10 cases MATPOWER's MP-Core and legacy post-processing share
        generator reactive output differently (MP-Core even reports NaN on case6470rte); there
        this code matches the legacy `pfsoln`, which it follows, within 2.9e-9.
      - **Time:** NR solves all 77 in 2.17 s; fast-decoupled solves its cases about three times
        faster than NR. On case2383wp, NR takes 12.9 ms, FDXB 4.15 ms and FDBX 3.87 ms.
      - **The Jacobian build is about half of Newton's time** (6 builds of 1.12 ms in 12.9 ms),
        so an assembly-map Jacobian becomes the follow-up item below. `AcModel::jacobian` is
        public, so the benchmark measures the solver's own build.
- [x] Review the README's references, the attribution of MATPOWER and of committed data, and
      `AI-DISCLAIMER.md` against the repository as it stands.
      - **Provenance corrected:** the session transcripts show which MATPOWER files were read.
        `makeSbus` and `dSbus_dV` were not, though the README (from PR #24) said the AC code was
        written with reference to them; the injections follow the User's Manual's definition and
        the Jacobian Technical Note 2's formulas. The README, `ac.rs`'s docs and the references
        now say so (Technical Note 2's DOI as MATPOWER's own source files cite it). No
        SuiteSparse source was read: only its licence file, for the licence notes.
      - **References:** every code name in the implemented table exists; nothing done is still
        listed as planned. Stott and Alsaç and KLU, newly implemented in this phase, were
        rechecked online (Stott and Alsaç has no DOI in the table, and none was confirmed).
      - **Committed data:** the reference results' README now also cites MATPOWER 8.1 by its
        version DOI, as MATPOWER's `CITATION` file recommends.
      - **`AI-DISCLAIMER.md`:** brought up to date on the work done with AI (through Phase 4.5),
        how it was checked (references kept in tests, MATPOWER validation and surveys), the CI
        jobs, the test data, and the dependencies' licences, which it listed slightly wrong
        (the Unlicense only ever appears as an alternative to MIT).
- [x] Moved to the next section: an assembly-map Jacobian, and enforcing generator reactive
      power limits.

## Before Phase 5: AC follow-ups

Moved ahead of Phase 5 by the maintainer (2026-11-01): the Phase 4.5 follow-up and deferred
item.

- [x] **Decide:** which reactive-limit modes, how to split the work, and which MATPOWER results
      to validate against.
      Decided (2026-11-01):
      - **Both of MATPOWER's modes:** all violators converted at once (`enforce_q_lims = 1`)
        and only the largest per round (`= 2`), off by default.
      - **Two PRs:** the Jacobian first, then reactive limits stacked on it.
      - **References:** if MATPOWER's MP-Core and legacy results agree within 1e-9, use
        `runpf`'s defaults; otherwise the legacy code this implementation follows, recording the
        gap, as for fast-decoupled.
- [x] Build the Newton Jacobian through an assembly map from Ybus's entries, refilling its values
      in place each iteration, bit for bit equal to the matrix-operation build it replaces.
      `JacobianAssembly` in `powerflow/src/ac.rs` lays the Jacobian out once per set of bus
      roles with `from_triplets_mapped` (each Ybus entry feeds up to one position per block)
      and refills it each iteration, computing each entry with the operations the old build
      performed on it. Tested bit for bit against the old build (kept in
      `tests/ac_jacobian_test.rs` as `reference_jacobian`) on the AC fixtures, case14, case30,
      case118 and case2383wp, at four voltage states each, and a refill against a fresh layout.
- [x] Measure it: criterion before and after, and an AC survey rerun that matches the record.
      - **Results unchanged:** the AC survey rerun (NR, FDXB, FDBX) matches
        `docs/benchmarks/ac-power-flow-2026-11-01.csv` in every column but time and memory, on
        all 234 case-and-method rows.
      - **Faster:** case2383wp's Newton solve takes 7.15 ms instead of 11.8 ms in criterion
        (−40%), and the survey's Newton solves 1.67 s together instead of 2.17 s. Fast-decoupled
        is unchanged, as it should be.
      - **Inferred refill cost:** the solve saved 4.7 ms over 6 iterations, so a refill costs
        about 0.18 ms where the old build took 0.96 ms. The `jacobian` benchmark measures a
        layout and one fill (0.86 ms), since `AcModel::jacobian` makes a fresh layout each call.
- [x] Let the AC solve take changing bus roles (REF, PV, PQ), with a `bustypes` shared with the
      topology, and no change in results.
      `topology.rs` has a `bustypes` function returning `Roles` (the ref, pv and pq lists) from a
      case's bus types, which `Topology::new` uses. `AcModel::solve` is split into `start`,
      `solve_round` (one solve under given roles and a case's injections) and `results` (which
      takes the roles for its slack). The DC survey (78 cases) and the AC survey (234
      case-and-method rows) match their records in every column but time and memory.
- [x] Enforce generator reactive power limits as `runpf`'s `pf.enforce_q_lims` does (both
      modes), with fixture tests.
      `AcOptions::reactive_limits` (`ReactiveLimits::Ignore` by default, `AllAtOnce`,
      `OneAtATime`), following legacy `runpf`'s loop: solve, `pfsoln`, find in-service generators
      past `Qmax + 5e-6` or `Qmin − 5e-6`, hold them (or the largest) at the limit and make their
      buses PQ, rerun `bustypes`, and solve again from the last voltages, with `pfsoln`'s outputs
      carried into the next round's injections. Violations are listed in MATPOWER's internal
      generator order (by bus, then row). `AcSolution::limited_generators` lists the held
      generators, and `iterations` totals every round. New errors: `ReactiveLimitsInfeasible`,
      `ReferenceAtReactiveLimit` (MATPOWER's multiple-slack error) and `InvertedReactiveLimits`.
      Changed from MATPOWER: a generator with `Qmin > Qmax` is refused, since holding it at
      either limit violates the other and `runpf`'s loop alternates forever (found when a test
      fixture did exactly that). With proper limits the loop always ends: a bus made PQ stays
      PQ, and a generator held at a limit can't violate it again. Fixture tests: upper and
      lower limits, two generators sharing a bus, a reference handing over to a PV bus (with
      the original angle restored), several references, infeasibility, inverted limits, and
      one at a time against all at once, each checked for power balance.
- [x] MATPOWER's reactive-limit reference results, and the tests against them.
      - **Which code:** MP-Core and legacy agree within 5.1e-12 on case14, case30 and case118,
        but over all cases they differ (up to 71 MW on case_SyntheticUSA; case6470rte converges
        only in legacy; case_ACTIVSg10k takes 17 iterations in MP-Core against 16). So, by the
        agreed rule, the `qlim1`/`qlim2` references come from the legacy `runpf`.
      - **Files:** `generate_reference.m` writes them, plus `_limited.csv`, the generators held at
        a limit (in service at a bus made PQ, at exactly one limit), since `runpf` doesn't report
        them. Its `write_csv` no longer leaves a blank line in a table with no rows.
      - **Validation:** all six runs (three cases, two modes) take MATPOWER's iteration counts
        (case14 5 and 5, case30 3 and 3, case118 6 and 16), hold the same generators, and agree
        within 8.1e-12. case14 holds its reference generator, so the reference moves to a PV bus;
        case30's limits don't bind; case118 holds six generators.
- [x] A survey with reactive limits enforced, recorded in `docs/benchmarks/`, and the docs.
      - **Survey:** `ac_survey --qlim 1|2` and `YBUS_BENCH_QLIM` in `benchmark-ac.sh`, with
        columns comparing the held generators. Over all 78 cases, in both modes, we converge on
        exactly MATPOWER's 74, with its iteration counts (up to 14,225 one at a time on
        case_SyntheticUSA), the same held generators, and within 1.7e-10 relative; the same four
        cases fail (`docs/benchmarks/ac-reactive-limits-2026-11-01.md`). The survey without
        limits is unchanged from its record.
      - **The reference's held generators,** inferred since `runpf` doesn't report them, wrongly
        included generators with equal limits beside a held one on three cases. The rule now
        leaves them out; nothing else changed, and the committed files are identical.
        `generate_reference.m` also writes each run's summary last, so that it marks a complete
        run.
      - **Not done:** comparing MP-Core with legacy one at a time over all cases. A machine
        restart interrupted it, and it needs many hours on the largest cases.
      - **Docs:** the README's status and provenance (`runpf`'s limit loop, and the legacy
        references), and the benchmarks README.
- [ ] **Follow-up:** reuse the ordering, LU analysis and Jacobian assembly across the rounds of
      a limit loop. One at a time spends 2,450 s on case_SyntheticUSA over 6,337 rounds of about
      2.2 iterations, each round repeating that setup, though each moves only one bus from PV to
      PQ.

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
