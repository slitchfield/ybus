# The CUDA spike, 2026-11-07

Phase O2's last item: a spike behind a `cuda` feature, measuring the GPU against the CPU path
with transfers counted. Two PRs built it:

- **#46:** the `ybus_cuda` crate (the workspace's only `unsafe` code) and triangular solves with
  many right-hand sides against a CPU factorization.
- **#47:** the DC power flow's workloads on the GPU: batched solves, PTDF columns and
  single-branch outages, solved whole or screened for their most loaded branch.

As decided when the spike started (2026-11-05), the GPU path is judged against the reference
tolerances, not bit for bit, and the CPU path stays the default unless the GPU wins.

## Setup

| | |
|---|---|
| Code | `main` at `1afb08a` (#47 merged) |
| Machine and system | Intel Core i7-14700K (28 hardware threads), 15 GB for WSL2, Ubuntu 24.04.3 under WSL2 (kernel 6.6.87.2), Rust 1.99.1, release profile |
| GPU | NVIDIA GeForce RTX 4070 Ti SUPER, 16 GB; driver 595.95 (CUDA 13.2); cuSPARSE 12.8.2 and NVRTC from the CUDA 13.3 toolkit, loaded at run time by cudarc 0.19 |
| Solves | The `cuda_survey` example (#46): many right-hand sides against B′ (Cholesky, or LU where it's indefinite) and the Newton Jacobian at the case voltages (LU). Each matrix gets up to 4,096 right-hand sides (at most 1 GB of them), solved once serially (`solve_many`), once on rayon's 28 threads (`par_solve_many`), and on the GPU at block widths 64–4,096; per-matrix results in [`cuda-spike-2026-11-07-solves.csv`](cuda-spike-2026-11-07-solves.csv) |
| Workloads | The `cuda_dc_survey` example through `scripts/benchmark-cuda-dc.sh` (#47), described below; per-workload results in [`cuda-spike-2026-11-07-dc.csv`](cuda-spike-2026-11-07-dc.csv) |
| Cases | Ten, from case118 to case_SyntheticUSA (82,000 buses), each in its own process |

Each case ran under `systemd-run` with a 10 GB memory limit and no swap, and a time limit.
`ulimit -v` can't limit CUDA, which reserves far more address space than it uses; a cgroup
limits what's resident, page-locked memory included. CUDA runs under WSL2 can crash the whole
instance, so every result went to `target/` as each case finished, not to `/tmp`.

**Every time includes the transfers and the host's work,** except where a column says it's the
solves alone. Each pass is timed once, after a warm-up. Times on the GPU varied between runs: N-1
screening's speed against 28 threads was 1.00× on case9241pegase here, against 1.16× in #47's
run of the same code, and 1.31× against 1.59× on case_ACTIVSg25k. Ratios near 1 can go either
way.

## The solves alone

#46's question was whether the GPU's sparse triangular solves (cuSPARSE's SpSM, the factors
uploaded once) are competitive with the CPU's. At each matrix's fastest block width:

| Case | Matrix | Order | Right-hand sides | CPU, 1 thread | CPU, 28 threads | GPU solves | GPU in total | Solves vs 28 threads | Total vs 28 threads |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| case2383wp | B′ (Cholesky) | 2,382 | 4,096 | 86.2 ms | 9.28 ms | 19.8 ms | 61.8 ms | 0.47× | 0.15× |
| | Jacobian | 4,438 | 4,096 | 227 ms | 22.7 ms | 41.3 ms | 120 ms | 0.55× | 0.19× |
| case9241pegase | B′ (LU) | 9,240 | 4,096 | 430 ms | 49.3 ms | 62.5 ms | 231 ms | 0.79× | 0.21× |
| | Jacobian | 17,036 | 4,096 | 1,157 ms | 110 ms | 122 ms | 440 ms | 0.90× | 0.25× |
| case_ACTIVSg25k | B′ (LU) | 24,999 | 4,096 | 1,273 ms | 135 ms | 143 ms | 589 ms | 0.95× | 0.23× |
| | Jacobian | 47,246 | 2,840 | 2,682 ms | 317 ms | 291 ms | 932 ms | 1.09× | 0.34× |
| case_ACTIVSg70k | B′ (LU) | 69,999 | 1,917 | 2,073 ms | 238 ms | 257 ms | 855 ms | 0.92× | 0.28× |
| | Jacobian | 134,104 | 1,000 | 3,372 ms | 544 ms | 303 ms | 1,029 ms | 1.80× | 0.53× |
| case_SyntheticUSA | B′ (LU) | 81,997 | 1,636 | 2,259 ms | 241 ms | 257 ms | 875 ms | 0.94× | 0.28× |
| | Jacobian | 156,255 | 858 | 3,347 ms | 620 ms | 325 ms | 993 ms | 1.90× | 0.62× |

(The CSV has every case and width. "GPU in total" is `solve_many` from host vectors to host
vectors.)

- **The solves are competitive,** but they win only on the largest Jacobians. From about 9,000
  unknowns up, the GPU's solves run at 0.8–1.1× of 28 threads on B′ and the smaller Jacobians,
  and 1.8–1.9× on the two largest Jacobians, whose factors hold 2.4–2.8 million entries. On
  smaller matrices they're slower: 0.4–0.6× at 2,000–4,000 unknowns.
- **In total, the GPU loses every comparison from case300 up** (0.15–0.71×). Uploading and
  downloading take about 50 ms each on the largest cases, and the host's packing into the
  factorization's order and unpacking out of it cost more than the solves themselves.
- **Block width matters little past 256.** From case1354pegase up, it changed the GPU's solve
  time by at most 1.8× (case_ACTIVSg2000's Jacobian, 64 against 4,096), most of it between 64
  and 256. Widths of 256 and more differed by at most 1.19×. On case118 and case300, narrow
  blocks cost up to 5.5×.
- **Agreement:** to within 5.4e-14 relative, on every matrix and width.

## The DC workloads

#46's conclusion was to keep the work on the device. #47 builds right-hand sides `e_plus −
e_minus` there, applies each outage's rank-one update there, and computes flows and per-column
maxima there, so only what's asked for comes back. Four workloads, each against the CPU's own
API, serially and on rayon's 28 threads:

- **N-1 screening:** each single-branch outage's most loaded branch (`|Pf| / rate_a`).
  - CPU: the outages' solutions from `solve_outages` and `par_solve_outages`, reduced on the
    host.
  - GPU: `DcGpuStudy::outage_loadings`. The flows and their maximum are computed on the device,
    so one value per outage comes back.
- **Batched DC:** copies of the case with scaled loads.
  - CPU: `DcStudy::solve` for each.
  - GPU: `DcGpuStudy::solve_many`. Only the solves run on the GPU; the right-hand sides, flows
    and slack are the CPU path's code, on the host.
- **PTDF columns:** columns for evenly spaced buses, at MATPOWER's slack.
  - CPU: `ptdf_columns` and `par_ptdf_columns`.
  - GPU: `DcGpuSensitivity::ptdf_columns`, with the flows computed on the device and
    downloaded.
- **N-1 solutions:** the same outages, solved whole.
  - CPU: `solve_outages` and `par_solve_outages`.
  - GPU: `DcGpuStudy::solve_outages`. The angles are downloaded, and the flows and slack are
    computed on the host.

Each workload has up to 8,192 items (fewer on the largest cases, to keep right-hand sides times
unknowns under 2²⁸). They run in chunks of the GPU's block width (1,024 on every case) on both
paths, and each chunk's results are dropped as they come. A first and a middle chunk of the
GPU's results are compared with the serial CPU's.

### N-1 screening: the GPU's win

| Case | Outages | CPU, 1 thread | CPU, 28 threads | GPU | GPU vs 1 thread | GPU vs 28 threads |
|---|--:|--:|--:|--:|--:|--:|
| case1354pegase | 1,991 | 29.9 ms | 7.53 ms | 10.0 ms | 2.98× | 0.75× |
| case2383wp | 2,896 | 82.1 ms | 16.2 ms | 20.7 ms | 3.97× | 0.78× |
| case_ACTIVSg2000 | 3,206 | 123 ms | 22.5 ms | 38.9 ms | 3.15× | 0.58× |
| case9241pegase | 8,192 | 1,742 ms | 324 ms | 323 ms | 5.40× | 1.00× |
| case13659pegase | 8,192 | 1,661 ms | 242 ms | 289 ms | 5.74× | 0.84× |
| case_ACTIVSg25k | 8,192 | 2,936 ms | 714 ms | 546 ms | 5.38× | 1.31× |
| case_ACTIVSg70k | 3,834 | 5,817 ms | 1,515 ms | 792 ms | 7.34× | 1.91× |
| case_SyntheticUSA | 3,273 | 5,862 ms | 1,566 ms | 761 ms | 7.70× | 2.06× |

- **The GPU's lead over 28 threads grows with size:** it's even around 9,000–14,000 buses, and
  1.3–2.1× faster from 25,000 up. In #47's run of the same code, the figures were 1.16×, 0.97×,
  1.59×, 1.90× and 1.87× on the five largest cases.
- **Over one thread, it's 3–8× faster from 1,354 buses up.** Below that, the GPU loses to both
  (case118 and case300, in the CSV): a launch costs more than the whole CPU pass.
- **Where the time goes:** on case9241pegase, a block of 1,024 outages took about 20 ms of
  triangular solves, 3 ms of flows and maxima, and 40 µs of the host's checks (#47's stage
  timings). SpSM is the bottleneck: #46's Nsight Systems trace put 91% of the kernel time in
  its solve kernel.

### Batched DC, PTDF columns and N-1 solutions: host-bound

GPU against 28 threads (and against one thread):

| Case | Batched DC | PTDF columns | N-1 solutions |
|---|--:|--:|--:|
| case1354pegase | 0.10× (0.99×) | 0.15× (1.07×) | 0.11× (0.98×) |
| case2383wp | 0.10× (1.22×) | 0.22× (1.33×) | 0.12× (0.92×) |
| case_ACTIVSg2000 | 0.11× (1.16×) | 0.14× (1.40×) | 0.10× (0.92×) |
| case9241pegase | 0.13× (1.66×) | 0.15× (1.79×) | 0.10× (0.98×) |
| case13659pegase | 0.12× (1.43×) | 0.14× (1.78×) | 0.12× (0.98×) |
| case_ACTIVSg25k | 0.12× (1.18×) | 0.18× (1.92×) | 0.12× (1.20×) |
| case_ACTIVSg70k | 0.25× (1.38×) | 0.33× (2.42×) | 0.20× (1.33×) |
| case_SyntheticUSA | 0.27× (1.42×) | — | 0.22× (1.30×) |

(case_SyntheticUSA has several islands, each with its own reference, so no single slack reaches
every bus and it has no PTDF, on either path.)

- **The GPU runs at 0.9–2.4× the speed of one CPU thread,** because each result costs O(buses +
  branches) on the host, and on the GPU path that host work runs on one thread. Batched DC and
  N-1 solutions build each case's right-hand side, place its angles, and compute its flows and
  slack. PTDF downloads and places a column per bus. The CPU path spreads the same work over 28
  threads.
- **PTDF shows where:** on case9241pegase, a block of 1,024 columns took 18 ms to solve, 49 ms
  to compute, download and copy the flows (131 MB), and 14 ms to place them by case branch row
  (#47's stage timings). The download itself is about 11 ms at the page-locked rate; the rest
  is copying into newly allocated vectors.

### Agreement

On every workload and case, the GPU's results agree with the serial CPU's to within 9.7e-13,
relative to each result's largest magnitude. No item fails on one path and not the other, and
every error is the same.

An outage whose update the GPU flags as singular, or whose angles fail the 1e5-radian check, is
solved again on the CPU, which decides its result. That keeps near-singular outages from
differing between the paths. The tests (`powerflow/tests/cuda_dc_test.rs`) exercise both
fallbacks.

## Engineering notes

- **cuSPARSE 12.8's SpSM** gives wrong results for a 1×1 matrix with more than 64 right-hand
  sides per block, with no error (found by #46's tests). A single-row factor is solved on the
  host, and its device blocks are at most 64 wide.
- **Kernels** are compiled once per device by NVRTC, into a CUBIN for the device's
  architecture. This driver (CUDA 13.2) rejects the 13.3 toolkit's PTX, while a CUBIN loads on
  any CUDA 13 driver. They're compiled without fused multiply-adds, so the outage update rounds
  as the CPU's does.
- **Page-locked staging:** every copy goes through page-locked host memory. On
  case_ACTIVSg70k's Jacobian, uploads took 55 ms against 87 ms from pageable memory, and
  downloads 62 ms against 93 ms (1.07 GB each way).
- **A block always solves its full width.** The survey's first run passed chunks of 512 to
  1,024-wide blocks, and screening on case9241pegase took 455 ms instead of 291 ms. Callers pass
  multiples of `block_width()`.
- **Errors:** `PowerflowError::Gpu` carries the GPU's error as text, and is converted
  explicitly, not by a `From` impl. A second `From` into `PowerflowError` left `?` in the
  reader's closures without a type to infer, so enabling the feature would have broken the
  build.
- **Nsight Systems under WSL2:** CUDA tracing and per-process CPU sampling work, and
  system-wide sampling doesn't (#46).

## Decisions

Decided with the maintainer on 2026-11-07:

- **The CPU path stays the default,** and `cuda` stays an opt-in feature. The GPU path is
  worth using for N-1 screening on large cases, from about 25,000 buses. Whole results (batched
  solves, outage solutions, PTDF columns) are faster on the CPU's parallel methods.
- **The Jacobian refill on the GPU is descoped,** with Newton on the GPU (cuDSS). Newton solves
  one right-hand side per factorization, where the GPU's solves are latency-bound and win only
  in wide blocks on the largest Jacobians. The refactorization kernel is latency-bound on the
  CPU too ([`lu-refactor-2026-11-05.md`](lu-refactor-2026-11-05.md)), since most of its work is
  a serial chain of separator columns. And cuDSS would be a new run-time dependency. A GPU
  factorization comes back as a question in phase O4's GPU KKT item, where one factorization
  serves an interior-point iteration's larger systems.
- **Record only:** the whole-result workloads stay as they are, documented as host-bound.

## What's left

- **Parallel host work on the GPU path:** with `rayon`, the whole-result workloads could
  prepare and finish their results on many threads around the device's solves. Overlapping
  that with the next block's solve might then beat the CPU path, but it's unmeasured.
- **Sets of branches:** the device's update is rank-one, so sets of branches stay on the CPU.
  A rank-k update, and LODF columns on the device, would follow the same pattern.
- **The solves themselves:** SpSM dominates screening. Level-scheduled or supernodal kernels of
  our own, or another library, weren't tried.
