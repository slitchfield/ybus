# Fused arithmetic in the LDLᵀ: the spike, 2026-11-22

This is phase O5's target 3, PR 3c (`OptimizationChecklist.md`, "Phase O5"). Fused
multiply-add changes rounding, so it would be opt-in under O3's bar. It was measured on a
throwaway patch before any option was built (decided 2026-11-22), against a gate written
beforehand.

**Outcome: no go.** Every fused build was slower than the default build, on both of the gate's
conditions.
- **The dense kernel:** ×1.24 as long with `+fma`, against the required 10% cut.
- **The AC solves:** a median ×0.94 at one thread and ×0.96 at eight, against the required
  ×1.05.

No option was built, and PR 3c closes with this record.

## Why measure first

- **No runtime dispatch:** `forbid(unsafe_code)` rules it out, since calling a
  `#[target_feature]` function needs `unsafe`. So fused instructions need a build flag.
- **PR 3a's slowdown:** its `x86-64-v3` (AVX2) build made the dense kernel take 1.83× as long,
  on short columns.
- **A memory-bound kernel:** the kernel updates the whole trailing triangle once per pivot, so
  fusing might save little.

## The spike

- **The patch:** every multiply-subtract in the LDLᵀ's hot loops became a fused `mul_add`, for
  real element types:
  - the dense kernel's rank-1 and rank-2 updates (`eliminate_one`, `eliminate_two`)
  - the forward and backward sweeps, serial and parallel
  
  It wasn't committed.
- **A shortcut in the patch:** the parallel forward sweep's recorded terms stayed unfused. So
  case_ACTIVSg25k at 8 threads, which takes the parallel path, gave different fused results
  from its serial run. Every other run of each build agreed with itself. A real option would
  fuse at the replay.
- **The builds,** all with `rayon`, each in its own worktree:

  | Build | Flags | Fused instructions |
  |---|---|---|
  | default | none, unpatched | none |
  | no flag | none | `mul_add` calls `compiler_builtins`' `fma`, which picks a hardware or software path at run time: a call per element |
  | `+fma` | `-C target-feature=+fma` | packed, all 256-bit (`ymm`) |
  | `+fma`, 128-bit | `+fma` and `-C llvm-args=-force-vector-width=2` | packed, all 128-bit (`xmm`) |
  | `x86-64-v3` | `-C target-cpu=x86-64-v3` | packed, 256-bit |

- **A premise that didn't hold:** the plan assumed `+fma` kept the default build's SSE2 width.
  It doesn't. FMA3 instructions need AVX, so `+fma` widens the vectors to 256 bits as AVX2
  does. The 128-bit build was added to measure what the plan meant.
- **The runs:**
  - `scripts/benchmark-opf.sh`, the predictor-corrector method's AC solves (`ldlt-uncorrected`,
    AMD), on case9241pegase, case_ACTIVSg10k and case_ACTIVSg25k, at 1 and 8 threads
  - two rounds, alternating the builds, the same day
  - one-thread profiles of case_ACTIVSg10k's AC
- **The rows:** [`o5-fused-spike-2026-11-22.csv`](o5-fused-spike-2026-11-22.csv).
- **The machine:** an Intel Core i7-14700K (28 threads, 15 GB under WSL2), Linux 6.6.

## The solves

The mean of the two rounds, against the default build. Each round's ratio was within 0.05 of
the other's.

| Build | Threads | case9241pegase | case_ACTIVSg10k | case_ACTIVSg25k | Median |
|---|---|---|---|---|---|
| default | 1 | 1.86 s | 3.77 s | 11.88 s | |
| no flag | 1 | ×0.61 | ×0.53 | ×0.49 | ×0.53 |
| `+fma` | 1 | ×0.97 | ×0.94 | ×0.92 | **×0.94** |
| `+fma`, 128-bit | 1 | ×0.94 | ×0.88 | ×0.91 | ×0.91 |
| `x86-64-v3` | 1 | ×0.92 | ×0.82 | ×0.84 | ×0.84 |
| default | 8 | 1.36 s | 2.43 s | 6.78 s | |
| no flag | 8 | ×0.70 | ×0.62 | ×0.62 | ×0.62 |
| `+fma` | 8 | ×0.98 | ×0.96 | ×0.96 | **×0.96** |
| `+fma`, 128-bit | 8 | ×0.96 | ×0.94 | ×0.97 | ×0.96 |
| `x86-64-v3` | 8 | ×0.96 | ×0.91 | ×0.91 | ×0.91 |

## Where the time went

case_ACTIVSg10k's AC at one thread, each stage's share of the profile times its wall time:

| Build | Solve | Dense kernel | Triangular solves | Front assembly | Evaluation |
|---|---|---|---|---|---|
| default | 3,768 ms | 991 ms | 411 ms | 384 ms | 1,044 ms |
| `+fma` | 4,177 ms | 1,232 ms (×1.24) | 459 ms | 430 ms | 1,065 ms |
| `+fma`, 128-bit | 4,620 ms | 1,529 ms (×1.54) | 517 ms | 420 ms | 1,123 ms |
| `x86-64-v3` | 4,840 ms | 1,815 ms (×1.83) | 479 ms | 421 ms | 1,094 ms |

- **The loss is the kernel's, not only the wider vectors'.** The 128-bit build, at the default
  build's width, lost more than `+fma`.
- **The cause wasn't measured** instruction by instruction. One reading: `eliminate_two` chains
  two fused operations per entry, and fusing replaces an independent multiply with an
  operation on the critical path. On short columns, that latency isn't hidden. The other
  stages barely moved.

## Accuracy, not gating

- **Unchanged:** every fused build gave the default's outcomes and iteration counts on all
  three cases.
- **The changes are in the last digits:** case9241pegase's objective error went from
  2.3327593789e-6 to 2.3327593773e-6, and its stationarity residual from 4545.23301967 to
  4545.23301974. So O3's bar would have held. (case_ACTIVSg10k fails, and case_ACTIVSg25k
  misses MATPOWER's objective, in the default build too: they're the uncorrected LDLᵀ's known
  AC failures.)
- **The same results in every flagged build:** the no-flag, `+fma` and `x86-64-v3` builds gave
  identical results, as the design expected, since `mul_add` is correctly rounded however
  it's computed.
- **One exception, unexplained:** the 128-bit build's results differ from the other fused
  builds' in the last digits. `-force-vector-width` is an internal LLVM option, so this wasn't
  pursued: it's one more reason not to rely on it.

## The gate

Written into the checklist before measuring:

| Condition | Required | Measured (`+fma`) | |
|---|---|---|---|
| The dense kernel's profile time, case_ACTIVSg10k's AC, one thread | at least 10% less | 24% more | not met |
| The AC solves' median, at 1 and at 8 threads | at least ×1.05 | ×0.94 and ×0.96 | not met |

**No go.** PR 3c closes with no option built, and target 3 closes with it.

**Outside the gate:** the 128-bit variant fails both conditions too (×1.54, ×0.91 and ×0.96).
A faster fused kernel would need a different structure, such as a blocked update that keeps
more of each column in registers. That's not planned.
