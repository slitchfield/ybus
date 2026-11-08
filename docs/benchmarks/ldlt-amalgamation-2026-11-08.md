# Relaxed amalgamation of the LDLᵀ's fronts, 2026-11-08

This is the second half of the pivoted LDLᵀ's bookkeeping, the target
[`ipm-profile-2026-11-08.md`](ipm-profile-2026-11-08.md) proposed first. The first half,
[`ldlt-bookkeeping-2026-11-08.md`](ldlt-bookkeeping-2026-11-08.md) (#61), kept the per-pattern work
with the analysis.

**What was left:** the fronts are tiny. On DC they had a median of 4 rows, at 0.9 fronts per
variable, so each factorization paid its cost per front about once per variable.

**This half:** it merges small supernodes into larger fronts, which changes the pivots and so the
results. It's judged by O3's bar, with #59's survey rerun for the predictor-corrector method.

## The change

- **`SymbolicIndefinite` postorders its ordering.** It's the same elimination tree, renumbered so
  that each subtree's columns are contiguous and a parent follows its last child.
- **It merges a child into its parent** when two things hold:
  - **the child's columns come just before the parent's;**
  - **the merged front passes a relaxed-amalgamation rule** (Ashcraft and Grimes, 1989). The
    default is CHOLMOD's: under 4 columns always, under 16 with under 80% explicit zeros, under
    48 with under 10%, and any size with under 5%.
- **`Amalgamation::NONE`** keeps the fundamental supernodes on the analysis's own ordering, as
  in #61.

**The fronts at the first factorization,** DC case1197 (order 2,396) and case9241pegase (order
19,928):

| | Fronts | Front rows: median / max | Delays | L entries |
|---|--:|---|--:|--:|
| case1197, none | 2,394 | 3 / 4 | 1,193 | 5,981 |
| case1197, merged | 613 | 4 / 18 | 65 | 14,897 |
| case9241pegase, none | 17,999 | 4 / 96 | 8,982 | 131,809 |
| case9241pegase, merged | 3,727 | 6 / 103 | 1,012 | 268,477 |

**Merging also absorbs most of the delays.** A delayed pivot is often a zero-diagonal constraint
row whose partner sits in its parent. Once the two share a front, they make a 2×2 pivot at once.

## A bug found by the survey's speed

The survey's first attempt ran up to 1,500× slower than #59: case1197 took 5.5 s per solve
against 3.6 ms.

- **The cause:** after a merge, the merged front's entries were recorded as its dense
  trapezoid. That counted its explicit zeros as true entries, so every later merge into it
  underestimated its zero fraction. The merges cascaded into dense fronts: case1197 became one
  front of 2,396 rows, with 2.9 million L entries.
- **The fix:** it now counts the true entries, as CHOLMOD does. A test bounds the fill at under 4×
  the unmerged factorization's (2.4× measured), and the bug fails it.
- **The runs** were stopped, and the survey below was rerun from the start on the fixed code.

## The same answers, by O3's bar

The rerun survey: the predictor-corrector method, DC and AC, LU against LDLᵀ (uncorrected, the AC
default), all 73 cases with costs, same machine and morning. The CSVs are beside this file, and
MIPS-LU's residuals for the bar are #59's, since MIPS is unchanged.

| | Outcomes differing from LU | LU meets the bar | LDLᵀ meets it | Met by LU, not by LDLᵀ | Iterations |
|---|--:|--:|--:|--:|---|
| DC | 0 | 58 | 58 | 0 | equal on all 64 converged |
| AC | 0 | 30 | 30 | 0 | equal on all 49 converged |

The LU runs match #54's survey on every result column, as before.

## Speed

**Against LU** (LU's time over LDLᵀ's, where both converged), with #59's figures for the same
comparison before #61 and this PR:

| | Median | Under 1,000 buses | 1,000–10,000 | Over 10,000 | Total |
|---|--:|--:|--:|--:|--:|
| DC, now | **1.15×** | 0.82× (39) | 1.79× (21) | 2.46× (4) | 28.8 → 9.2 s (3.14×) |
| DC, #59 | 0.79× | 0.57× | 1.32× | 1.72× | 2.02× |
| AC, now | 1.16× | 1.09× (34) | 1.63× (13) | 3.08× (2) | 402.7 → 131.2 s (3.07×) |
| AC, #59 | 1.12× | 1.04× | 1.53× | 2.70× | 2.70× |

- **The largest cases:**
  - **DC case_ACTIVSg70k:** 21.8 → 5.4 s, against #59's 8.9 s.
  - **AC case_ACTIVSg70k:** 301 → 87 s, against #59's 96.5 s.
- **One DC outlier is left:** case3375wp at 0.27×, where delays still enlarge its fronts.
- **Peak memory:** the median ratio to LU is 1.00. On the largest cases LDLᵀ now uses less than
  LU: 280 against 328 MiB on DC case_ACTIVSg70k, and 889 against 944 MiB on AC.

**Amalgamation's own share,** timed on one binary with it on and off (`profile_pf`, a mean over
5 s of repeated runs, at least one):

| Case | Off | On | Gain |
|---|--:|--:|--:|
| DC case300 | 2.23 ms | 2.19 ms | 1.02× |
| DC case1197 | 3.86 ms | 3.67 ms | 1.05× |
| DC case2383wp | 120.4 ms | 115.4 ms | 1.04× |
| DC case9241pegase | 99.8 ms | 94.0 ms | 1.06× |
| DC case_ACTIVSg10k | 361 ms | 342 ms | 1.06× |
| DC case_ACTIVSg70k | 7.32 s | 5.96 s | 1.23× |
| AC case2383wp | 792 ms | 766 ms | 1.03× |
| AC case9241pegase | 3.99 s | 3.76 s | 1.06× |
| AC case_ACTIVSg25k | 26.4 s | 25.1 s | 1.05× |

- **What the gain is:** modest, 2–23%. Most of the turnaround against LU on DC comes from #61 and
  this PR together.
- **Why it isn't more:** the profiles show front assembly down 43–67% and allocation down 12–66%,
  but the dense kernel spending 27–101% more, on fronts that are larger and carry explicit zeros.
  The kernel isn't blocked, so a larger front costs it per entry about what a smaller one did.

## Decisions for review

1. **The DC OPF's predictor-corrector default.** It now meets #59's rule:
   - **The bar:** kept wherever LU meets it.
   - **The median:** 1.15× faster.
   - **The total:** 3.14× faster.

   Switching `DcOpfOptions::mehrotra()` to LDLᵀ is left to the maintainer. It's a change of
   results, though the same outcomes and iterations on every case.
2. **The amalgamation rule.** CHOLMOD's defaults were taken as they are. Its zero allowances were
   set for a blocked supernodal kernel. Ours pays for zeros per entry, so a tighter rule may do
   better. That's a tuning question for a later PR, with the kernel's blocking.

## What's left

- **The dense kernel** is now the factorization's largest part. Blocked updates would let the
  larger fronts pay off.
- **The delays,** which still enlarge some fronts (case3375wp): a KKT-aware ordering.
- **On AC, the OPF's evaluation** (#60's second target).
