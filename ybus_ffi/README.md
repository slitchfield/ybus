# ybus_ffi

The C interface to [`powerflow`](../powerflow/): the `ybus` library, with the header
[`include/ybus.h`](include/). It's being built phase by phase from
[`docs/FfiChecklist.md`](../docs/FfiChecklist.md), and isn't released yet.

## Building and linking

```sh
cargo build -p ybus_ffi --release
```

This builds `target/release/libybus.so` and `target/release/libybus.a`. Compile against the
header and link with `-lybus`. A program linking the static library on Linux also needs the
system libraries Rust's standard library uses, which
`cargo rustc -p ybus_ffi --lib --crate-type staticlib -- --print native-static-libs` lists:

```sh
cc -I ybus_ffi/include program.c target/release/libybus.a \
    -lgcc_s -lutil -lrt -lpthread -lm -ldl -lc
```

## Tests

The Rust tests run with `cargo test -p ybus_ffi`. The C tests in [`tests/c/`](tests/c/) are
small programs, compiled against the header and linked with the static library as a C program
would be, then run, by `tests/c/run.sh`.

`check.sh` runs the C tests twice: plainly, then with `tests/c/run.sh --sanitize`, under
AddressSanitizer (with LeakSanitizer) and UndefinedBehaviorSanitizer.
- **What it catches:** the Rust library isn't instrumented, since stable Rust can't build with
  `-Zsanitizer`, but its allocations go through the `malloc` that ASan replaces. So a handle
  that's leaked or freed twice fails the run, as do a C buffer overrun and undefined behaviour
  in the C code.
- **What it doesn't:** Rust writing past a C buffer. The exports' length checks guard against
  that, and the Rust tests cover them.

## Cases

A `ybus_case` is an opaque handle to a case, made by one of three functions and freed by
`ybus_case_free`:

- **`ybus_case_from_path(path, &case)`** reads a MATPOWER file. The path is a NUL-terminated
  UTF-8 string.
- **`ybus_case_from_bytes(bytes, len, &case)`** parses the file's text from memory.
- **`ybus_case_clone(case, &copy)`** copies a case, so the copy can change without changing the
  original.

Each sets `*out` to `NULL` first, and leaves it `NULL` on an error. `ybus_case_free(NULL)` does
nothing. Rows are numbered from 0, as in C.

Every getter returns a status and writes its value through an out-parameter, such as
`ybus_case_bus_count(case, &buses)`. The sizes and scalars are:
- `ybus_case_bus_count`, `ybus_case_gen_count` and `ybus_case_branch_count`
- `ybus_case_base_mva`, and `ybus_case_set_base_mva`
- `ybus_case_has_costs` and `ybus_case_has_bus_names`

### Columns

Each column of the bus, generator and branch tables is read whole, one entry per row in case
order, by a field code. The code is the field's 1-based column number in MATPOWER's matrix for
that table. For example, `YBUS_BUS_PD_MW` is 3, `YBUS_GEN_STATUS` is 8 and `YBUS_BRANCH_X` is
4; the header lists them all.

There's one getter per table and value type:

```c
double pd[3];
ybus_case_get_bus_f64(case, YBUS_BUS_PD_MW, pd, 3);
int64_t status[3];
ybus_case_get_branch_i64(case, YBUS_BRANCH_STATUS, status, 3);
```

- **`int64_t` columns:** labels, areas, zones, bus types (MATPOWER's codes 1 to 4) and statuses
  (1 in service, 0 not). Every other field is a `double` column.
- **A code that's wrong for the function** gives `YBUS_ERROR_INVALID_FIELD`, whether the table
  doesn't have it or it's of the other type.
- **Too short a buffer** gives `YBUS_ERROR_BUFFER_TOO_SMALL` and is left untouched. A longer
  one keeps its tail.

Each getter has a setter, such as `ybus_case_set_bus_f64(case, YBUS_BUS_PD_MW, values, 3)`.
- **Length:** the array must have exactly one entry per row, or it's
  `YBUS_ERROR_LENGTH_MISMATCH`.
- **`int64_t` columns** are checked whole before anything is written. A negative label, area or
  zone, a bus type outside 1 to 4, or a status other than 0 or 1 is
  `YBUS_ERROR_INVALID_VALUE`, and the column is left as it was.
- **`double` columns** aren't checked. Data a solve can't use, such as a branch to a bus that
  doesn't exist, is reported when the case is solved.

`ybus_case_bus_name(case, row, buffer, len, &needed)` reads one bus name, as a NUL-terminated
string. Call it with `NULL, 0` to learn the size. A case without names gives
`YBUS_ERROR_MISSING_DATA`.

### Costs

The generator costs are read and set whole, as flat arrays. Each row has:
- a model: `YBUS_COST_POLYNOMIAL` (2) or `YBUS_COST_PIECEWISE_LINEAR` (1), MATPOWER's codes
- a startup cost and a shutdown cost
- a run of values, `values[offsets[i]]` up to `values[offsets[i + 1]]`. A polynomial's values
  are its coefficients, highest order first. A piecewise linear cost's values are its points,
  interleaved as `p0, cost0, p1, cost1, …`.

The rows are one per generator for real power, optionally followed by one per generator for
reactive power, as in MATPOWER.

```c
size_t rows, count;
ybus_case_cost_sizes(case, &rows, &count);
/* model, startup and shutdown have rows entries, offsets rows + 1, values count */
ybus_case_get_costs(case, model, startup, shutdown, offsets, rows, values, count);
ybus_case_set_costs(case, model, startup, shutdown, offsets, rows, values, count);
```

Setting replaces the whole table, and setting 0 rows removes it. The structure is checked first:
the models, an even number of values for each piecewise linear row, and offsets that start at
0, never decrease, and end at the number of values. What a solve can't use is reported when the
case is solved.

## Power flows

### Options

Options are handles made at MATPOWER's defaults and freed with their own `_free`. Wherever a
solve takes options, `NULL` means the defaults.

- **DC:** `ybus_dc_options_new(&options)` makes the options. The one setting is
  `ybus_dc_options_set_ordering`, B′'s fill-reducing ordering: `YBUS_ORDERING_AMD` by default,
  or `NATURAL`, `RCM` or `COLAMD`.
- **AC:** `ybus_ac_options_new(method, &options)` takes the method: `YBUS_AC_NEWTON`,
  `YBUS_AC_FAST_DECOUPLED_XB` or `YBUS_AC_FAST_DECOUPLED_BX`. It starts at that method's
  defaults: tolerance 1e-8, at most 10 Newton or 30 fast-decoupled iterations, and reactive
  limits ignored.
  - The setters are `set_method` (which changes only the method), `set_tolerance` (positive and
    finite), `set_max_iterations` and `set_reactive_limits`.
  - The reactive limits use MATPOWER's `enforce_q_lims` codes: `YBUS_REACTIVE_LIMITS_IGNORE` 0,
    `ALL_AT_ONCE` 1 and `ONE_AT_A_TIME` 2.

Every setter has a getter, and a code that isn't one of the setting's is
`YBUS_ERROR_INVALID_VALUE`.

### Solutions

A solve makes a solution handle, freed with its own `_free`. Its columns are read whole, one
getter per table, by the case's field codes. The codes are the columns MATPOWER writes each
result to, so a result reads where MATPOWER would put it:

```c
ybus_dc_solution *solution;
ybus_dc_power_flow(case, NULL, &solution);
ybus_dc_solution_get_bus(solution, YBUS_BUS_VA_DEG, angles, buses);
ybus_dc_solution_get_branch(solution, YBUS_BRANCH_PF_MW, flows, branches);
ybus_dc_solution_free(solution);
```

- **DC:** bus `YBUS_BUS_VA_DEG`, gen `YBUS_GEN_PG_MW`, and branch `YBUS_BRANCH_PF_MW` (14).
  `ybus_dc_solution_solver` says which factorization solved B′.
- **AC:** bus `VM_PU` and `VA_DEG`; gen `PG_MW` and `QG_MVAR`; and branch `PF_MW`, `QF_MVAR`,
  `PT_MW` and `QT_MVAR` (14 to 17).
  - The scalars are `ybus_ac_solution_iterations`, `_mismatch_pu` and `_full_factorizations`.
  - `ybus_ac_solution_limited_generators` gives the generator rows held at a reactive limit,
    numbered from 0. Call it with `NULL, 0` to ask the count.
- **A code the solution doesn't have** in that table gives `YBUS_ERROR_INVALID_FIELD`.

### Models

A model is built once from a case's network and solves many scenarios of it, giving the same
result as the power flow of each:

- **`ybus_dc_model_new(case, &model)` and `ybus_dc_model_solve(model, case, options,
  &solution)`.** The DC model refactors B′ on each solve, as `DcModel::solve` does.
- **`ybus_ac_model_new(case, &model)` and `ybus_ac_model_solve(model, case, options,
  &solution)`.** `ybus_ac_model_restamp(model, case)` takes in changed impedances, charging,
  taps, shifts and shunts.

The case solved can change its injections, setpoints and starting voltages, but not its
network. A case whose tables' lengths differ is `YBUS_ERROR_CASE_MISMATCH`, and one whose
topology differs (a bus's label or type, a generator's bus or status, or a branch's ends or
status) is `YBUS_ERROR_TOPOLOGY_MISMATCH`.

## Optimal power flows

### Options

Each OPF has its own options handle, made for a solver and freed with its own `_free`:
- `ybus_dc_opf_options_new(solver, &options)` and `ybus_ac_opf_options_new(solver, &options)`
- the solver is `YBUS_OPF_MIPS`, MATPOWER's MIPS, or `YBUS_OPF_MEHROTRA`, the predictor-corrector
  method

**The defaults:**
- MIPS starts at MATPOWER's settings: `feastol` 5e-6, with the AC OPF's cost scaled by
  `cost_mult` 1e-4.
- The predictor-corrector starts at its tuned settings, which factor the KKT system by the
  pivoted LDLᵀ.
- The solver can't change after the options are made.

**The settings** apply to either solver, and each has a setter and a getter:
- `feastol`, `gradtol`, `comptol` and `costtol`, the termination tolerances, positive and
  finite
- `cost_mult`, the objective's scale, positive and finite
- `max_iterations`
- `kkt_factorization`, `YBUS_KKT_LU` or `YBUS_KKT_LDLT`

The solvers' other settings keep their defaults.

### Solutions

`ybus_dc_opf(case, options, &solution)` and `ybus_ac_opf(case, options, &solution)` solve the
OPFs, with `NULL` options for MATPOWER's defaults. Their solutions' columns are read as the
power flows' are, by the case's field codes, at the columns MATPOWER writes the OPF's results
to:

| Table | DC OPF columns | AC OPF columns |
|---|---|---|
| bus | `VA_DEG`, and `LAM_P` (14), the price | `VM_PU`, `VA_DEG`, `LAM_P` (14), `LAM_Q` (15), `MU_VMAX` (16) and `MU_VMIN` (17) |
| gen | `PG_MW`, and `MU_PMAX` (22) and `MU_PMIN` (23) | `PG_MW`, `QG_MVAR`, `VG_PU`, and `MU_PMAX` to `MU_QMIN` (22 to 25) |
| branch | `PF_MW`, `MU_SF` (18), `MU_ST` (19), `MU_ANGMIN` (20) and `MU_ANGMAX` (21) | `PF_MW` to `QT_MVAR` (14 to 17), and `MU_SF` to `MU_ANGMAX` (18 to 21) |

The scalars are `_objective`, `_iterations` and `_status`. The status uses MIPS's `exitflag`
codes: `YBUS_OPF_CONVERGED` (1), `YBUS_OPF_MAX_ITERATIONS` (0) and
`YBUS_OPF_NUMERICAL_FAILURE` (-1). A solve that doesn't converge still returns `YBUS_OK` with
its solution, as in MATPOWER, so check the status.

## ABI policy

- Every symbol is prefixed `ybus_`, and every macro `YBUS_`.
- `ybus_version()` returns the version of the library linked at run time, encoded as the
  header's `YBUS_VERSION` is (major × 1 000 000 + minor × 1 000 + patch). The header also
  defines `YBUS_VERSION_MAJOR`, `YBUS_VERSION_MINOR` and `YBUS_VERSION_PATCH`.
- Before 1.0, any minor version may break the ABI.
- From 1.0, minor versions only add functions and status codes.

## Errors

- **Status codes.** Every fallible function returns a `ybus_status` (an `int32_t`): `YBUS_OK`,
  which is 0, or the code of what failed. The header lists the codes. A code is never
  renumbered, and a failure without a code of its own, such as one added in a later version,
  is `YBUS_ERROR_OTHER`.
  - 1 to 6 are the boundary's own failures, such as a null argument or a buffer too short.
  - From 100, each `PowerflowError` variant has its own code, in declaration order, named
    after it: `YBUS_ERROR_NO_CASE_FILE` is 100, and `YBUS_ERROR_TOPOLOGY_MISMATCH` is 116.
- **The last error.** Each thread keeps the message of its last failed call, which
  `ybus_last_error_message(buffer, len)` copies out. It returns the bytes the message needs,
  counting the NUL, or 0 if no call on the thread has failed. It writes nothing unless the
  buffer is big enough, so call it with `NULL, 0` first to learn the size. The message is set
  only on failure, so it describes the most recent failure, not the most recent call.
- **Panics.** A panic in the library is a bug. It's caught at the boundary and returned as
  `YBUS_ERROR_PANIC`, rather than aborting the process. Afterwards, only free the handles that
  were passed to that call. Rust's default panic hook prints the panic to stderr first; the
  library doesn't replace the process-wide hook.

## Thread safety

- Any handle can move between threads.
- A handle that isn't being changed (a case, a model, a solution) can be read from many threads
  at once.
- A case being changed, or an options handle being set, needs the caller's own lock.

The Rust types behind the handles are checked to be `Send` and `Sync` by
`powerflow/tests/thread_safety_test.rs`.
