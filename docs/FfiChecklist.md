# FFI checklist

The roadmap for long-term goal 1 in [`LongTermGoals.md`](LongTermGoals.md): a stable C interface
to `powerflow`, then Python and Julia bindings. Each unchecked item in phases F1 to F6 is sized to
be one commit; phase F7's are coarser and are broken into commit-sized items when each is planned.
Phase F0's decisions were made with the maintainer on 2026-11-10.

## Goals

1. **A C ABI first.** A crate of its own that wraps `powerflow` in `extern "C"` functions: opaque
   handles for cases, options, models and solutions, status codes instead of panics, and buffers
   the caller owns. It's the common ground the other languages can build on.
2. **Python and Julia after it.** Python through PyO3 and maturin, binding `powerflow` directly;
   Julia through the C library with `ccall`, so it shares the C layer's stability.
3. **Out of scope at first:** `sparse_mat`'s matrices, the interior-point solvers on their own,
   the DC outage studies and sensitivities, and the `rayon` and `cuda` features. Phase F7 lists
   them.

**Conventions carried over from the earlier roadmaps:**
- one commit per checklist item, on `dev/` branches, with the PR's description as its report
- unsafe code only in the new crate, as `ybus_cuda` holds the GPU path's; `powerflow` keeps
  `#![forbid(unsafe_code)]`
- every exported function tested from C as well as from Rust, against the committed reference
  cases

## What there is to build on

- **Entry points.** `Case::from_path`, `from_bytes` and `parse`; `dc_power_flow`,
  `ac_power_flow`, `dc_opf` and `ac_opf`, each with a `_with_options` partner; and `DcModel` and
  `AcModel`, built once and solved for many scenarios. Since PR #65 they all take options by
  reference and return one solution type per solve.
- **Data shapes.** `Case` holds vectors of structs (`Bus`, `Gen`, `Branch`, `GenCost`) and
  optional bus names. Each solution is a struct of columns, mostly `Vec<f64>`, indexed like the
  case's tables, plus a few scalars.
- **Errors.** `PowerflowError` has 26 variants, many carrying data, and is `#[non_exhaustive]`.
  Its `Display` messages are written for people.
- **Options.** `DcOptions`, `AcOptions`, `DcOpfOptions` and `AcOpfOptions` are all
  `#[non_exhaustive]`. The OPF ones choose a solver through `OpfSolver`, which wraps
  `MipsOptions` or `MehrotraOptions`.
- **Panics.** Invariants are checked with `unwrap`, `expect` and `panic!`, so a bug behind any of
  them panics. The fuzz tests (`powerflow/tests/matpower_fuzz_test.rs`) already use
  `catch_unwind` to check that the MATPOWER reader returns errors instead of panicking. The
  workspace's profiles keep the default `panic = "unwind"`.
- **Two hazards for a C caller.**
  - **Mismatched cases.** A model is solved against a case passed again, and `check_case` in
    `powerflow/src/topology.rs` compares only the tables' lengths, so a different case of the
    same size is solved silently with the wrong network.
  - **A borrowed model.** `DcStudy<'a>` borrows its `DcModel`, which a C handle can't express.

## Phase F0: Decisions

- [x] **Decide** the first API's scope and the shape of its handles. Decided (2026-11-10):
      - loading a case from a path or from bytes
      - the DC and AC power flows and both OPFs, with default or custom options
      - the DC and AC models, for building once and solving many scenarios
      - reading every column of each solution

      Cases, options, models and solutions are separate handles, as in the Rust API. A solution
      then outlives later solves, read-only handles can be shared between threads, and the
      Python and Julia bindings stay thin. One session handle holding all of them was the
      alternative: fewer objects to free, but one thread per session, and each run replacing
      the last one's results. The outage studies, sensitivities, matrices and features wait for
      phase F7.
- [x] **Decide** the crate's name and form. Decided (2026-11-10): `ybus_ffi` in the workspace,
      with `[lib] name = "ybus"` so C links `-lybus` and includes `ybus.h`, matching the `ybus_`
      symbol prefix; `crate-type = ["cdylib", "staticlib", "rlib"]`, the `rlib` so Rust tests can
      call the exports; `publish = false` until a release; `powerflow`'s default features only.
      The long-term goal suggested `powerflow-ffi`, but the other crates use underscores, and
      the library is the whole project's.
- [x] **Decide** how a case crosses the boundary. Decided (2026-11-10): an opaque handle with
      whole-column getters and setters, for example every bus's `pd_mw` read into or written
      from one array. A C caller can run scenarios on one case, the ABI doesn't depend on the
      structs' layouts, and whole arrays suit NumPy and Julia. Rejected:
      - **loaders only:** C couldn't change a case
      - **`#[repr(C)]` mirror structs:** every field of `Bus`, `Gen` and `Branch` would be part of
        the ABI
      - **one element at a time:** slow across the boundary on the largest cases

      Changing part of a column, by range or by index set, waits for phase F7.
- [x] **Decide** how results cross the boundary. Decided (2026-11-10):
      - each solution is an opaque handle
      - each column has an accessor that copies it into a buffer the caller allocates, with a
        length check (too short is an error, never a partial write)
      - the scalars (`objective`, `iterations`, `mismatch_pu`) have plain getters
      - the OPF's `status` becomes a code
      - the OPF's `solver` and `residuals` fields stay out, since they expose `interior_point`'s
        internals

      Pointers into the solution's own memory were the alternative: no copies, but a C caller
      could use them after freeing the solution.
- [x] **Decide** how errors cross the boundary. Decided (2026-11-10):
      - every fallible function returns a `ybus_status`, an `int32_t` enum
      - `YBUS_OK` is 0, with one code per `PowerflowError` variant and codes for the boundary's
        own failures: a null argument, invalid UTF-8, a buffer too short, and a caught panic
      - `YBUS_ERROR_OTHER` covers variants added to `PowerflowError` later, until each gets its
        own code
      - codes are never renumbered
      - the message is kept per thread and read by `ybus_last_error_message`, which keeps every
        signature short and works when a constructor fails before any handle exists

      Rejected: an error handle returned through an out-parameter (an extra argument on every
      call, and one more object to free), and a message kept on the failing handle (awkward for
      constructors and for functions that take several handles).
- [x] **Decide** how options cross the boundary. Decided (2026-11-10): opaque handles created
      at the defaults (`ybus_ac_options_new`), with one typed setter per field (for example
      `ybus_ac_options_set_tolerance`), and a null handle meaning the defaults. The compiler
      checks each name and type, and a new option is a new function, which doesn't break the ABI.
      String-keyed setters were the alternative: no new symbols for new options, but a
      misspelled name or wrong type would only fail at run time.
- [x] **Decide** the ABI and versioning policy, to be written in the crate's README. Decided
      (2026-11-10):
      - every symbol is prefixed `ybus_`
      - `ybus_version` returns the library's version, and the header defines it as macros
      - before 1.0, any minor version may break the ABI
      - from 1.0, minor versions only add functions and codes
- [x] **Decide** the thread-safety contract, to be written in the crate's README. Decided
      (2026-11-10):
      - a handle can move between threads
      - a handle that isn't being changed (a case, a model, a solution) can be read from many
        threads at once
      - a case being changed, or an options handle being set, needs the caller's own lock

      Phase F1 checks that the Rust types underneath are `Send` and `Sync`.
- [x] **Decide** how the header is made. Decided (2026-11-10): generated with the `cbindgen`
      command-line tool and committed as `ybus_ffi/include/ybus.h`, with CI regenerating it and
      failing when it differs. `cbindgen` is MPL-2.0, which `deny.toml`'s license list doesn't
      include, so as a build dependency it would fail the Dependencies job; installed in CI as a
      tool, it stays out of the dependency tree. Rejected: generating it in `build.rs` (a build
      dependency and a license exception), and writing it by hand (signatures could drift in
      ways the C compiler doesn't catch, such as a changed integer width).

## Phase F1: Preparing `powerflow`

- [x] `check_case` compares the case's topology as well as the tables' lengths, so a model
      solved against a different case of the same size returns an error. Decided with the
      maintainer (2026-11-10):
      - **What's compared:** everything topological, which decides the elements a model covers
        and each bus's role: each bus's label and type, each generator's bus and status, and
        each branch's ends and status. The injections, setpoints, voltages and (AC) impedances
        and shunts may still differ, as each model's `solve` and `restamp` document.
      - **The error:** a new variant, `TopologyMismatch { table, row, field }`, naming the first
        difference (1-based row, MATPOWER's table names), while `CaseMismatch` stays for tables
        of different lengths. Phase F2's status codes give each its own code.
- [x] A test asserts at compile time that `Case`, the models, the options and the solutions are
      `Send` and `Sync`, so the thread-safety contract is checked rather than assumed
      (`powerflow/tests/thread_safety_test.rs`, which also covers `DcStudy` and
      `PowerflowError`).

## Phase F2: The crate and its boundary

- [x] The crate: a workspace member named and built as decided in F0, with `ybus_version`, and a
      README holding the ABI policy and thread-safety contract. `check.sh` and CI build and lint
      it, through their `--workspace` commands. `ybus_version` returns the version encoded as
      SQLite encodes its own, major × 1 000 000 + minor × 1 000 + patch, as does the header's
      `YBUS_VERSION`.
- [x] The header: generated by the `cbindgen` command-line tool, committed, and regenerated in
      CI, which fails when it differs. `ybus_ffi/generate-header.sh` runs cbindgen 0.29.4 (CI
      installs that version); `check.sh` runs it with `--verify`, and skips the check, saying
      so, on a machine without cbindgen, unless `YBUS_REQUIRE_CBINDGEN=1`, as CI sets.
- [x] A C test harness: a small C program per area, compiled against the header and linked
      against the static library, run by `check.sh` and CI. Later items each add their C tests.
      `ybus_ffi/tests/c/run.sh` builds each `test_*.c` with `-std=c11 -Wall -Wextra -Wpedantic
      -Werror`, using one `CHECK` macro (`check.h`) rather than a framework; `check.sh` runs it
      once, after the Rust tests.
- [x] One guard that every export goes through:
      - `catch_unwind` turns a panic into `YBUS_ERROR_PANIC` and records its message
      - null and alignment checks on pointer arguments
      - the per-thread last error recorded on failure

      A panic reaching an `extern "C"` function aborts the process (since Rust 1.81), so the
      guard turns that into a status code. A test forces a panic and checks the code. Done in
      `ybus_ffi/src/boundary.rs`, with the boundary's own codes (`YBUS_OK` 0, then
      `YBUS_ERROR_OTHER`, `PANIC`, `NULL_ARGUMENT`, `MISALIGNED_ARGUMENT`, `INVALID_UTF8` and
      `BUFFER_TOO_SMALL`, 1 to 6) as a `typedef int32_t ybus_status` and `#define`s, as
      decided with the maintainer (2026-11-10). The last error is set on failure only, so a
      later success doesn't clear it.
- [x] `ybus_status` and the mapping from `PowerflowError`, with `YBUS_ERROR_OTHER` as the
      mapping's catch-all arm, since `#[non_exhaustive]` requires one outside `powerflow`.
      `ybus_last_error_message` copies the message into a caller's buffer and returns the length
      it needs. Each variant's code is from 100, in declaration order (`YBUS_ERROR_GPU` is
      reserved, since the C build has no `cuda`). No export can fail with a `PowerflowError`
      yet, so as decided with the maintainer (2026-11-10), a Rust test maps every current
      variant to its own code, and C tests trigger real errors from phase F3 on, starting with
      a missing case file.

## Phase F3: Cases

- [x] Loading and freeing: `ybus_case_from_path` (a NUL-terminated UTF-8 path),
      `ybus_case_from_bytes`, `ybus_case_clone` and `ybus_case_free` (a no-op on null). Each
      constructor takes an out-parameter, set to `NULL` first and left `NULL` on an error.
      Decided with the maintainer for phase F3 (2026-11-10):
      - **every getter** returns a `ybus_status` and writes through an out-parameter
      - **columns** are read and written by field code, one getter and one setter per table and
        value type, rather than about 94 functions written by hand (cbindgen can't see
        functions generated by macros)
      - **integer columns** are `int64_t`, and setters reject values that can't be represented
      - **costs** are flattened, with offsets into one array of values

      The C tests' first real errors are here: a missing file, malformed text, a path that
      isn't UTF-8 and a null argument. Their case, `ybus_ffi/tests/c/three_bus.m`, is a file
      the Rust tests read too.
- [x] Sizes and scalars: the bus, generator and branch counts, `base_mva`, and whether the case
      has costs and bus names. `ybus_case_set_base_mva` sets the base, unchecked until a solve.
- [x] Column getters for every numeric field of the bus, generator and branch tables, with bus
      types and statuses as integer columns using MATPOWER's codes. Bus names, if any, are read
      one at a time by index. A field code is the field's 1-based column number in MATPOWER's
      matrix for the table, so the format fixes the codes; a test reads every code against the
      case file's column of that number. New codes: `YBUS_ERROR_INVALID_VALUE` 7,
      `LENGTH_MISMATCH` 8, `INVALID_FIELD` 9 and `MISSING_DATA` 10.
- [x] Column setters for the same fields. A set isn't validated; invalid data is reported when
      the case is solved, as in Rust, where building a model runs `Case::validate`. Only what an
      `int64_t` column can't hold is refused when set, the whole column checked first; a test
      sets each column to its values rotated and checks that no other column changes.
- [x] Generator costs, read and set in a form for both cost models. The form (flattened
      coefficients with offsets, or one generator at a time) is decided with the item. Decided:
      flattened, with offsets. `ybus_case_cost_sizes` gives the rows and values;
      `ybus_case_get_costs` and `ybus_case_set_costs` take each row's model (MATPOWER's code),
      startup and shutdown, offsets, and one array of values, coefficients or interleaved
      points. Setting checks the structure before changing anything, and 0 rows removes the
      costs.

## Phase F4: Power flows

- [x] Options handles for the DC and AC power flows: create at the defaults, free, and one typed
      setter per field: the ordering; the method, tolerance, iteration limit and reactive
      limits. Decided with the maintainer for phase F4 (2026-11-10):
      - **The AC method is chosen when the options are made,** as `AcOptions::for_method`
        does: `ybus_ac_options_new(method, &options)` starts at that method's defaults (30
        fast-decoupled iterations, 10 Newton), and `set_method` changes only the method.
      - **Each setting also has a getter,** an addition that lets a binding show the settings.
      - **The codes:** the orderings are 1 to 4, and the methods 1 to 3 (MATPOWER's `NR`,
        `FDXB` and `FDBX`). The reactive limits take MATPOWER's `enforce_q_lims` codes, 0 to 2.
- [x] `ybus_dc_power_flow`: a case and options (null for the defaults) in, a solution handle out,
      with its column accessors and free. Tests against the committed MATPOWER references.
      Decided with the maintainer (2026-11-10):
      - **A solution's columns are read by the case's field codes,** one getter per table, at
        the column MATPOWER writes each result to. The branch table gains the result codes 14
        to 17 (`PF`, `QF`, `PT` and `QT`), which a case doesn't have.
      - **The tests are fidelity, plus a C reference.** A Rust test checks that the FFI's
        result is `powerflow`'s, bit for bit, under every ordering, on `three_bus.m` and on
        case14, case30 and case118. A C test compares case14 with the committed CSVs, at
        `powerflow`'s 1e-9.

      That test found that `Ordering::Colamd` panicked in `DcModel::study`. COLAMD orders
      columns only, so Cholesky refuses it; B′ now goes straight to LU under it, a fix in
      `powerflow` in its own commit.
- [x] `ybus_ac_power_flow`, in the same form, with the iteration count, mismatch, limited
      generators and factorization count. Tested as the DC power flow is: in Rust, against
      `powerflow`'s own results under each method and both kinds of reactive limits; in C,
      case14 against MATPOWER's `runpf`, without and with reactive limits.
- [x] The DC and AC models: build from a case, solve against a case (with F1's stricter
      `check_case`), restamp (AC), and free. A C test solves several load scenarios against one
      model. Each solve on a model is, bit for bit, the power flow of the same case, and a
      branch switched out is `YBUS_ERROR_TOPOLOGY_MISMATCH`. The DC model refactors B′ on
      each solve, as `DcModel::solve` does; factoring once waits for `DcStudy` in phase F7.

## Phase F5: Optimal power flows

- [x] OPF options handles: the solver (MIPS or the predictor-corrector method), and its main
      tolerances, iteration limit and KKT factorization, each by a typed setter. Which of the
      solvers' other fields to expose is decided with the item; the rest keep their defaults.
      Decided with the maintainer for phase F5 (2026-11-10):
      - **Two handles,** `ybus_dc_opf_options` and `ybus_ac_opf_options`, each made for a
        solver at its own problem's defaults, so options can't reach the wrong OPF.
      - **Exposed:** the four termination tolerances, the iteration limit, the KKT
        factorization, and `cost_mult`, which the AC OPF already sets apart. All seven are
        common to both solvers, so every setter applies to either.

      Each setting has a getter too, as for the power flows' options.
- [x] `ybus_dc_opf`: a solution handle with dispatch, flows, prices, multipliers, objective,
      status code and iteration count. Tests against the committed references.
      - **Columns:** read by the case's field codes, at MATPOWER's OPF result columns. These
        add `LAM_P`, `LAM_Q`, `MU_VMAX` and `MU_VMIN` (bus, 14 to 17), `MU_PMAX`, `MU_PMIN`,
        `MU_QMAX` and `MU_QMIN` (gen, 22 to 25), and `MU_SF`, `MU_ST`, `MU_ANGMIN` and
        `MU_ANGMAX` (branch, 18 to 21).
      - **The status** takes MIPS's `exitflag` codes. A solve that doesn't converge still
        returns its solution.
      - **Tests:** in Rust, `powerflow`'s results bit for bit under each solver and KKT
        factorization; in C, case14 against MATPOWER's `rundcopf`.
- [x] `ybus_ac_opf`, in the same form, with the voltage, reactive power and branch multipliers.
      Tested as the DC OPF is. In C, case14 is compared with MATPOWER's `runopf`, every column
      at `powerflow`'s 1e-8, and the generators' PQ capability curves in `three_bus.m` are
      refused with `YBUS_ERROR_UNSUPPORTED_GENERATOR` until zeroed.

## Phase F6: Hardening and a first release

- [x] CI runs the C tests under AddressSanitizer and UndefinedBehaviorSanitizer, so leaks,
      double frees and out-of-bounds writes at the boundary fail the build.
      - **How:** `check.sh` runs `tests/c/run.sh --sanitize` after the plain run.
      - **Coverage:** the Rust library isn't instrumented on stable Rust, but its allocations go
        through ASan's `malloc`. A leaked case, a double free and a C-side overrun, each
        planted in a scratch program, fail the run.
      - **No suppressions were needed:** Rust's runtime reported nothing.
- [x] A C example in the crate (load a case, run each solve, print a few results), and its use
      in the crate's README. `ybus_ffi/examples/c/solve_case.c` runs the four solves with
      default options, reporting each failure and carrying on. `tests/c/run.sh` builds it and
      runs it on case14 in both modes, so the sanitizers check it too.
- [ ] A record of the exported symbols, compared with the built library in CI, so a symbol
      removed or renamed is a visible decision. Whether this starts before 1.0 is decided with
      the item.
- [ ] The boundary's own cost, measured on mid-size cases: each solve through C against the same
      solve in Rust. Recorded in `docs/benchmarks/`; it feeds long-term goal 2.

## Phase F7: Later

Coarse, and planned when its turn comes:

- **Partial case edits:** setting part of a column by range or by index set, so a scenario that
  changes a few loads on a large case doesn't copy every column.
- **DC outage studies** through C, after `DcStudy` owns or shares its model (an `Arc`, say)
  instead of borrowing it.
- **Sensitivities** (PTDF and LODF), and builds with the `rayon` and `cuda` features.
- **Python** through PyO3 and maturin, binding `powerflow` directly: cases from files or arrays,
  the power flows and OPFs, and results as NumPy arrays. Tested with pytest in CI.
- **Julia** through the C library with `ccall`, wrapped in a small package, tested with Julia's
  own test runner in CI.
