# Benchmarks

Recorded performance of the `powerflow` crate, so that later phases can measure what they
change. Each record names the code, the machine and the method it used.

| Record | What |
|---|---|
| [`dc-power-flow-2026-10-31.md`](dc-power-flow-2026-10-31.md) | The DC power flow on every MATPOWER 8.1 case, before Phase 2's fixed-pattern matrices and Phase 3's ordered LU: the Phase 1.5 baseline. The survey's raw output is in the `.csv` beside it. |
| [`dc-power-flow-rcm-2026-10-31.md`](dc-power-flow-rcm-2026-10-31.md) | The same survey with B′'s Cholesky factorization ordered by reverse Cuthill–McKee (Phase 3), with the first test of the hypothesis that fill-in drives the backward error. |
| [`dc-power-flow-amd-2026-11-01.md`](dc-power-flow-amd-2026-11-01.md) | The Cholesky cases ordered by approximate minimum degree against RCM (Phase 3): fill, the cost of each stage, and the hypothesis's verdict for Cholesky. |
| [`dc-power-flow-lu-2026-11-01.md`](dc-power-flow-lu-2026-11-01.md) | The LU cases with the Gilbert–Peierls LU (Phase 3), unordered and AMD-ordered, against the old LU, and the hypothesis's verdict for LU. |
| [`ac-power-flow-2026-11-01.md`](ac-power-flow-2026-11-01.md) | The AC power flow on every case by Newton's method and fast-decoupled XB and BX (Phase 4.5), against MATPOWER's `runpf`, with criterion benchmarks on case2383wp. |
| [`ac-reactive-limits-2026-11-01.md`](ac-reactive-limits-2026-11-01.md) | Newton's method on every case with generator reactive limits enforced, all at once and one at a time, against MATPOWER's legacy `runpf`, including the held generators. |

## Measuring

- **One case, repeatably:** `cargo bench -p powerflow --bench dc_power_flow` runs the criterion
  benchmarks in `powerflow/benches/dc_power_flow.rs`, which time parsing, the model and the
  solve on case2383wp. Criterion compares each run with the previous one on the same machine.
- **Reassembly:** `cargo bench -p powerflow --bench assembly` times refilling case2383wp's
  B stamps, and a Jacobian-shaped pattern of 2×2 blocks built from them, through an assembly
  map against building them afresh with `from_triplets` (`powerflow/benches/assembly.rs`).
- **The AC power flow:** `cargo bench -p powerflow --bench ac_power_flow` times the AC model,
  the three solution methods and one Newton Jacobian on case2383wp, and
  `scripts/benchmark-ac.sh [<output.csv>]` surveys every case with the `ac_survey` example, the
  method chosen by `YBUS_BENCH_METHOD` (`nr`, `fdxb` or `fdbx`) and reactive limits by
  `YBUS_BENCH_QLIM` (0, 1 or 2, as `pf.enforce_q_lims`).
- **Every case, once:** `scripts/benchmark-dc.sh [<output.csv>]` runs the `dc_survey` example
  on each MATPOWER case in its own process, under a time and memory limit, and compares the
  results with MATPOWER's own. `YBUS_BENCH_ORDERING` chooses the Cholesky ordering. Its header comment lists what it needs and its settings.

Timings depend on the machine, so compare runs made on the same one, and record the machine
with any new record.

## What a record may hold

Case names and sizes, the time and memory this crate's code takes on them, the backward error
of its results, and how far they are from MATPOWER's. Never case data or power flow results: MATPOWER's license doesn't
cover its case files, and only results from public-source cases are committed
(`docs/PowerFlowChecklist.md`, Phase 1.5).
