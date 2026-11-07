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
| [`dc-power-flow-condition-2026-11-02.md`](dc-power-flow-condition-2026-11-02.md) | The DC survey with B′'s estimated condition number and the backward error after iterative refinement (Phase 5): condition doesn't explain the backward error, and refinement has little left to do. |
| [`dc-contingency-parallel-2026-11-02.md`](dc-contingency-parallel-2026-11-02.md) | The N-1 survey with the outages solved in parallel (Phase 5): the same results, and 11× faster on cases of 1,000 buses or more. |
| [`dc-contingency-2026-11-02.md`](dc-contingency-2026-11-02.md) | Every single-branch outage of every case by rank-1 updates of B′'s kept factorization (Phase 5), against refactoring. |
| [`colamd-2026-11-03.md`](colamd-2026-11-03.md) | COLAMD against AMD on every case's Jacobian and B′: AMD fills less on every case, so it stays the default. |
| [`ac-reactive-limits-reuse-2026-11-03.md`](ac-reactive-limits-reuse-2026-11-03.md) | The reactive-limit survey with a limit loop's later rounds solving against one superset Jacobian: the same results as MATPOWER, one at a time 2.9× faster in all. |
| [`lu-refactor-2026-11-05.md`](lu-refactor-2026-11-05.md) | The LU refactorization kernel (phase O2): latency-bound; fused nested column pairs make it 1.19–1.27× faster on the largest Jacobians, bit for bit, and why level scheduling and other options were set aside. |
| [`parsing-2026-11-05.md`](parsing-2026-11-05.md) | Reading MATPOWER cases (phase O2): byte-level scans and a vectorizable precheck, with the same results for every input; 1.54× faster over all cases, and what's left. |
| [`cuda-spike-2026-11-07.md`](cuda-spike-2026-11-07.md) | The CUDA spike (phase O2): triangular solves with many right-hand sides and the DC workloads on the GPU against 28 CPU threads, transfers and host work counted. N-1 screening wins from about 25,000 buses (up to 2.1×); whole-result workloads are host-bound. The solve and workload CSVs are beside it. |
| [`dc-contingency-allocation-2026-11-04.md`](dc-contingency-allocation-2026-11-04.md) | The N-1 survey after phase O2's first CPU target: reachability from a bridge forest instead of a search per outage, sparse solves into reused buffers, and per-thread outage workspaces. The same results bit for bit, 4.5× faster serially and 4.6× in parallel over all cases. |
| [`profile-2026-11-04.md`](profile-2026-11-04.md) | Where the power flows spend their time (phase O2): samply profiles of DC, Newton, fast-decoupled, reactive limits and N-1 outages on four large cases, folded into stages, and the targets chosen from them. |
| [`opf-2026-11-04.md`](opf-2026-11-04.md) | The DC and AC optimal power flows on every case with costs, solved by MIPS (phase O1), against MATPOWER's MP-Core and legacy implementations: outcomes, iterations, agreement and time against MATPOWER. |
| [`ac-reactive-limits-2026-11-01.md`](ac-reactive-limits-2026-11-01.md) | Newton's method on every case with generator reactive limits enforced, all at once and one at a time, against MATPOWER's legacy `runpf`, including the held generators. |

## Measuring

- **Profiling:** `--profile profiling` builds the release profile with debug line tables
  (`Cargo.toml`), so a sampling profiler names functions and lines in the same optimized code
  the surveys time, for example `cargo build --profile profiling -p powerflow --example
  ac_survey`, then a profiler on `target/profiling/examples/ac_survey`.
  `scripts/profile-pf.sh [<output.md>]` profiles the power flows with samply through the
  `profile_pf` example, which repeats each workload's whole pipeline so every stage weighs as in
  one run, and `scripts/profile_summary.py` folds the samples into stages. Its `parse` workload
  reads a case alone, for timing the reader.

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
- **Outages:** `scripts/benchmark-n1.sh [<output.csv>]` takes out each in-service branch of
  each case alone with the `n1_survey` example, by rank-k updates of B′'s kept factorization,
  and compares a sample with refactoring (every outage up to 3,000 buses, about 300 above).
  `YBUS_BENCH_PARALLEL=1` solves them in parallel (`rayon`).
- **The optimal power flow:** `scripts/benchmark-opf.sh [<output.csv>]` runs the `opf_survey`
  example on each case with costs, the model chosen by `YBUS_BENCH_MODEL` (`dc` or `ac`), and
  compares it with MATPOWER's MP-Core and legacy results, generating any that are missing.
- **Orderings:** `scripts/benchmark-ordering.sh [<output.csv>]` factors each case's Jacobian and
  B′ under AMD and COLAMD with the `ordering_survey` example.
- **Every case, once:** `scripts/benchmark-dc.sh [<output.csv>]` runs the `dc_survey` example
  on each MATPOWER case in its own process, under a time and memory limit, and compares the
  results with MATPOWER's own. `YBUS_BENCH_ORDERING` chooses the Cholesky ordering. Its header comment lists what it needs and its settings.
- **The GPU (feature `cuda`):** `scripts/benchmark-cuda-dc.sh [<output.csv>]` runs the
  `cuda_dc_survey` example's DC workloads on each case, each in its own process under a time
  limit and a cgroup memory limit through `systemd-run` (`ulimit -v` can't limit CUDA).
  `cuda_survey <case.m>` times B′ and Jacobian solves for many right-hand sides by block width.
  Both need the `cuda` and `rayon` features and a CUDA device.

Timings depend on the machine, so compare runs made on the same one, and record the machine
with any new record.

## What a record may hold

Case names and sizes, the time and memory this crate's code takes on them, the backward error
of its results, and how far they are from MATPOWER's. Never case data or power flow results: MATPOWER's license doesn't
cover its case files, and only results from public-source cases are committed
(`docs/PowerFlowChecklist.md`, Phase 1.5).
