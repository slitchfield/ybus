# The C interface's own cost, 2026-11-10

Phase F6 of [`FfiChecklist.md`](../FfiChecklist.md): what calling a solve through the C
interface costs, against the same solve in Rust. This feeds long-term goal 2, whose
comparisons with other libraries will run through the bindings.

**Result:** the boundary's own work is microseconds. Reading every column of a solution, which
is most of what the boundary adds, takes 0.8 to 8.7 µs on cases of 1,354 to 3,120 buses. That's
at most 0.03% of a power flow and 0.007% of an OPF. Timed side by side, the whole solves differ
by −3.8% to +4.5%, in both directions, which is this machine's run-to-run noise rather than the
boundary.

## Setup

| | |
|---|---|
| Machine | Intel Core i7-6600U (2 cores, 4 threads), 11 GiB of memory, Ubuntu 22.04 (Linux 5.15) |
| Code | `dev/ffi-f6-hardening` at `a9415e1`, with this record's example |
| Build | `--release`, Rust 1.99.0 |
| Cases | MATPOWER 8.1's case1354pegase, case2383wp and case3120sp (the power flows), and case1354pegase (the OPFs) |
| Options | Each solve's defaults: Newton for the AC power flow, and MIPS at MATPOWER's settings for the OPFs |

This laptop isn't the machine of the earlier records. The comparison is within one process on
one machine, so it should carry over, but the absolute times don't.

## Method

The example [`ybus_ffi/examples/boundary_cost.rs`](../../ybus_ffi/examples/boundary_cost.rs)
times each solve two ways, alternating, in one process.

- **Rust:** `powerflow`'s function, such as `dc_power_flow_with_options`, on a case it loaded
  itself. The result is kept, through `black_box`.
- **Through C:** the export, such as `ybus_dc_power_flow`, on a case loaded with
  `ybus_case_from_path`.
  - Every column the solution has is read into buffers allocated beforehand, as a C caller's
    would be, and then the solution is freed.
  - That's the code a C program runs: the guard, the solution's handle, the copies and the
    free. The exports are called through the rlib, and a real C caller adds only a plain
    function call.

**Before timing,** each pair of results is checked to agree bit for bit. That check doubles as
each side's warm-up.

**Runs:**
- 21 runs of each power flow and 5 of each OPF, reporting the medians.
- "Reading the columns" is the boundary's work alone: every column of one kept solution, read
  1001 times, reporting the median.

To rerun, on another machine or other cases:

```sh
cargo run --release -p ybus_ffi --example boundary_cost -- \
    testdata/matpower/matpower8.1/data/case1354pegase.m \
    testdata/matpower/matpower8.1/data/case2383wp.m \
    testdata/matpower/matpower8.1/data/case3120sp.m
cargo run --release -p ybus_ffi --example boundary_cost -- --opf \
    testdata/matpower/matpower8.1/data/case1354pegase.m
```

## Results

| Case | Solve | Columns | Rust (median) | Through C (median) | Difference | Relative | Reading the columns |
|---|---|---:|---:|---:|---:|---:|---:|
| case1354pegase | DC power flow | 3 | 3.818 ms | 3.915 ms | +97.3 µs | +2.55% | 0.8 µs |
| case1354pegase | AC power flow | 8 | 11.730 ms | 11.873 ms | +142.9 µs | +1.22% | 2.9 µs |
| case2383wp | DC power flow | 3 | 5.655 ms | 5.477 ms | −178.2 µs | −3.15% | 1.6 µs |
| case2383wp | AC power flow | 8 | 25.381 ms | 25.892 ms | +511.4 µs | +2.01% | 4.4 µs |
| case3120sp | DC power flow | 3 | 12.162 ms | 11.702 ms | −459.6 µs | −3.78% | 1.9 µs |
| case3120sp | AC power flow | 8 | 37.095 ms | 36.988 ms | −107.6 µs | −0.29% | 6.0 µs |
| case1354pegase | DC OPF | 10 | 73.864 ms | 77.200 ms | +3336.0 µs | +4.52% | 4.8 µs |
| case1354pegase | AC OPF | 21 | 1887.389 ms | 1907.614 ms | +20224.5 µs | +1.07% | 8.7 µs |

## Reading the numbers

**The side-by-side differences are noise.**
- They change sign from case to case: the DC power flow is faster through C on two cases and
  slower on one.
- An earlier run of the same example, before the column timing was added, gave −4.1% and
  +11.7% for case1354pegase's DC power flow in two runs.
- They're hundreds of times larger than the boundary's measured work.

The laptop's clock scaling and the solves' own variation dominate at these sizes.

**What the boundary adds**, per solve:
- a guard, which is a `catch_unwind` that costs nothing unless a panic is caught
- a check on each pointer
- one box for the solution, and one free
- one copy of each column into the caller's buffer

The copies are the measurable part: about 0.3 ns per value here. case3120sp's AC solution is
22,022 values (2 columns over 3,120 buses, 2 over 505 generators and 4 over 3,693 branches),
read in 6.0 µs. Next to a solve's milliseconds, that's at most 0.03%, and a caller that reads
fewer columns pays less.

**Not measured:**
- loading a case through C, which is `Case::from_path` behind the same guard
- the case's own column getters and setters, which copy the same way
- the models, which add nothing over the power flows' boundary

For goal 2, the bindings' cost can be taken as nothing next to the solves. A comparison with
another library through these bindings measures the solver, not the boundary.
