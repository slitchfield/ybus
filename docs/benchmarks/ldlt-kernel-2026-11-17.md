# The LDLᵀ's dense kernel on the lower triangle only, 2026-11-17

This is the kernel half of O4's fourth target. After #80, the LDLᵀ's dense kernel, the threshold
Bunch–Kaufman factorization the multifrontal LDLᵀ runs on each front (`factor_front`), was the
largest stage of the AC predictor-corrector solve: 27.6% on case9241pegase.

**One change was kept: the kernel updates only the lower triangle**, and writes the upper once at
the end. A blocked update was tried and dropped, below. Both stayed **bit for bit the same.**

## The change

- **Before,** each pivot's elimination copied the whole trailing lower triangle into the upper:
  `a[(j, i)] = conj(a[(i, j)])`. That's a strided write per entry, as many writes as the update
  computed. Yet during the factorization only the threshold test and the swap read the upper
  triangle.
- **The readers now read the lower triangle:** `conj(a_ji)` for `a_ij` above the diagonal, which
  has the same magnitude. That gives the same maximum, and the same first row on ties.
- **The swap is a Hermitian swap of the lower triangle,** as LAPACK's `xHETRF` does it, with moves
  and conjugations only.
- **The upper triangle is written once, at the end,** as before, so the multifrontal driver and
  `DenseLdlt` are unchanged:
  - zero above the diagonal in the eliminated rows, except each 2×2 block's `conj(d21)`
  - the conjugate of the lower triangle in the trailing square

## The same results, checked

- **Fronts:** a test-only copy of the old kernel compared the new one on random Hermitian fronts,
  real and complex, of 1 to 150 rows.
  - **The variations:** partial (`k < m`) and complete, thresholds 1e-8, 1e-4 and 0.5.
  - **Required of the cases:** they reach 2×2 pivots, delays and singular fronts.
  - **Compared:** the whole front array, the order, the eliminated count and D's blocks, bit for
    bit, or the same error.
  - **Removed with this record,** as planned. The copy and its tests are at `9734167`
    (`bunch_kaufman/reference.rs`, and `bunch_kaufman.rs`'s test module), to restore for any
    later change to the kernel.
- **Planted bugs,** each caught:
  - a swap that doesn't conjugate the joining entry, or the stretch between (seen only in
    complex values, as conjugation does nothing to reals)
  - a mirror that loses a 2×2 block's upper entry
- **The 72-case OPF snapshot** matches `67d936b`'s, #80's tip. MIPS factors by LU, so its hashes
  are unchanged too.
- **The predictor-corrector surveys** on LDLᵀ, all 73 cases on DC and AC: none of the 1,241
  result columns differs in either.

## Speed

**Mean time per run** of `profile_pf`'s whole pipeline, for the predictor-corrector method on
LDLᵀ:
- **Both binaries:** `67d936b` and `9734167`.
- **Runs:** interleaved before, after, two rounds, on an otherwise idle machine.
- **Repeats:** 15 s per case, and 40 s on case_ACTIVSg25k.

| Model | case300 | case2383wp | case9241pegase | case_ACTIVSg25k |
|---|--:|--:|--:|--:|
| AC | 1.11× | 1.14× | 1.12× | 1.26× |
| DC | 0.97× | 1.05× | 0.99× | 1.13× |

AC case_ACTIVSg25k: 13.7 s → 10.9 s. DC's small cases are within noise; their fronts are tiny,
so the kernel is little of the time.

**The surveys' times,** from `opf_ms`, run back to back:

| Model | Median | Under 1,000 buses | 1,000–10,000 | Over 10,000 | Total |
|---|--:|--:|--:|--:|--:|
| AC | 1.05× | 1.04× | 1.11× | 1.24× | 1.24× |
| DC | 1.05× | 1.04× | 1.07× | 1.17× | 1.19× |

- **AC case_SyntheticUSA:** 94.9 s → 75.2 s.
- **AC case_ACTIVSg70k:** 49.3 s → 39.3 s.

**The kernel stage itself:** samply at 1 kHz, 10 s per pipeline, folded by `profile_summary.py`.
The "LDLᵀ: dense kernel" stage's time per run:

| Model | case2383wp | case9241pegase | case_ACTIVSg25k |
|---|--:|--:|--:|
| AC | 110 → 70 ms (1.56×) | 542 → 329 ms (1.65×) | 5,191 → 2,598 ms (2.00×) |
| DC | 32 → 26 ms (1.23×) | 14 → 11 ms (1.30×) | 472 → 318 ms (1.48×) |

**The stages after,** on AC case9241pegase:

| Stage | Share |
|---|--:|
| OPF evaluation | 25.7% |
| Dense kernel | 18.2% |
| Triangular solves | 11.9% |
| Front assembly and extend-add | 10.0% |

## Tried and dropped: a blocked update

The plan's second commit was a blocked update, also bit for bit. It was measured three ways, all
bit for bit against the reference kernel, and every one was slower than the lower-triangle
kernel alone. Mean per run, 10 s each, idle machine:

| Variant | AC case2383wp | AC case9241pegase | DC case_ACTIVSg25k |
|---|--:|--:|--:|
| Lower triangle only (kept) | 338–342 ms | 1,730–1,770 ms | 1,060–1,082 ms |
| Contribution block deferred, swept by column | 367 ms | 1,859 ms | 1,192 ms |
| The same, tiled 8 columns per row | 367–373 ms | 1,913–1,921 ms | 1,253–1,260 ms |
| Tiled, only for blocks of at least 32, 64 or 128 rows | | 1,778–1,816 ms | 1,107–1,122 ms |

What was learned:

- **Only the contribution block can be deferred bit for bit.** That's the rows and columns that
  aren't fully summed, which no threshold test reads and no swap moves.
  - **Inside the fully summed part,** a swap moves entries across the diagonal, conjugating them.
    A lazily updated entry would then take `conj(w_j)·l_i` where the old kernel had
    `conj(w_j·conj(l_i))`.
  - **Those are equal mathematically,** but not always bit for bit, so a full LAPACK-style lazy
    panel there would change results.
- **Deferring didn't pay at these front sizes.** The median AC front has 7 rows, so the
  contribution blocks are small. Once the strided mirror was gone, the eager update already ran
  down contiguous columns that stay in cache.
  - **The deferred sweep** only moved the reads from `a` to the stored `w` columns.
  - **Tiling's transposed copies** cost more than the reuse saved.
  - **No size cutoff** beat the kernel without deferral.

So the code was not kept. Further kernel gains need reordered or fused arithmetic, such as SIMD
reductions and FMA, which change bits. That belongs to O5's exploration, judged by O3's bar.

## What's left

- **For target 4:** the KKT-aware ordering (the next PR), which should shrink the fronts' delays
  and fill.
- **For O5:**
  - the dense kernel with reordered arithmetic
  - the triangular solves (11.9%)
  - the front assembly (10.0%)
