# The matched KKT ordering, and the ordering study, 2026-11-17

This is the ordering half of O4's fourth target, and O4's last item.

**The problem:** the pivoted LDLᵀ delayed 8,400–48,000 pivots on the OPF's KKT matrices at any
threshold.
- **The cause:** AMD makes the constraints' zero-diagonal rows leaves of the elimination tree.
  Each is eliminated before any variable it touches, with nothing in its front to pair with, so
  it's delayed to the parent.
- **The fix:** the matched ordering (Duff and Pralet, 2005). It pairs each constraint with a
  variable it touches and places the constraint right after it. By the constraint's turn, its
  diagonal holds the Schur complement's entry, and it can pivot.

## The change

- **`KktOrdering { Amd, Matched }`** on both solvers' options (`kkt_ordering`), for the pivoted
  LDLᵀ. LU ignores it.
- **The method:**
  - each constraint row is matched to a variable with `Jg(k, j) ≠ 0`, greedily by `|Jg|` and
    then by augmenting paths, so the matching is maximum
  - each pair is merged into one node, and the merged graph ordered by AMD
  - each pair is placed as its variable, then its constraint
- **Kept as the LDLᵀ's own ordering would be:** it's made once per analysis, and kept while the
  pattern holds. The new `analyze_ldlt_with_permutation` takes it, and the relaxed amalgamation's
  postorder keeps each constraint after its variable.
- **Measured with the new statistics** on `Solution` (`solution.kkt`): the factorizations, delayed
  pivots, 2×2 pivots and factor entries, totalled over a solve, and surveyed by `opf_survey`.

## The rule, and the defaults

The rule was written into the checklist before the survey ran (`f81e7d9`). For each model,
matched becomes the predictor-corrector default if, on all 73 cases against AMD:
1. O3's bar holds wherever it holds with AMD
2. the outcomes are the same
3. the median `opf_ms` is faster

O3's bar was evaluated as #59 evaluated it, by the same functions. Where MATPOWER converged:
- the objective within 1e-6
- the prices within 1e-3
- each KKT residual within 100× MIPS-LU's, or below 1e-6

| Model | 1. Bar where AMD meets it | 2. Outcomes differing | 3. Median speedup (total) | Default |
|---|--:|--:|--:|---|
| DC | 58 of 58 | 0 | 1.14× (1.27×) | **matched** (`DcOpfOptions::mehrotra()`) |
| AC | 30 of 30 | 1 | 1.09× (1.17×) | **AMD**, unchanged |

**AC's one difference:** case13659pegase ends at the iteration limit with AMD (150) and in a
numerical failure with matched (147). Neither converges, but the outcomes differ, which the rule
doesn't allow. So AC keeps AMD, and the rule wasn't revised after the fact. A later rule, written
before another survey, could count two failures as the same outcome.

**The snapshot:** the 72-case OPF snapshot changes only in its DC predictor-corrector column (67
of 72 cases), as the switch should. MIPS's columns and AC's are unchanged.

## The study

Totals over the cases that converged in all three settings: 64 on DC, 49 on AC. The three
settings are the predictor-corrector method with LDLᵀ under AMD, with LDLᵀ under matched, and with
LU under COLAMD (MIPS's KKT ordering).

| Model | Setting | Factor entries per factorization | Delays | 2×2 pivots | Time | Iterations |
|---|---|--:|--:|--:|--:|--:|
| DC | LDLᵀ, AMD | 310,090 | 891,548 | 3,040,121 | 7.0 s | 499 |
| DC | LDLᵀ, matched | 252,228 | 120,379 | 1,692,951 | 5.7 s | 499 |
| DC | LU, COLAMD | 609,994 | — | — | 27.0 s | 499 |
| AC | LDLᵀ, AMD | 1,118,891 | 5,121,799 | 10,721,643 | 61.2 s | 915 |
| AC | LDLᵀ, matched | 887,919 | 236,029 | 328,327 | 51.2 s | 915 |
| AC | LU, COLAMD | 2,868,759 | — | — | 356.5 s | 915 |

- **The matched ordering removes almost all the delays:** 7.4× fewer on DC, and 22× fewer on AC.
- **Fill falls with them:** 19–21% fewer factor entries.
- **The 2×2 pivots fall most on AC,** 33× fewer. Once a constraint follows its variable, it can
  usually be a 1×1 pivot.
- **COLAMD on LU fills 2–3× more** than either LDLᵀ ordering. It's 4–7× slower, MIPS's
  MATPOWER-faithful default.
- **The iterations are the same** in every setting, case by case, on these cases.

**The large cases:**

| Model | Case | Delays, AMD → matched | Factor entries per factorization | Time |
|---|---|--:|--:|--:|
| DC | case_ACTIVSg25k | 145,308 → 13,147 | 919,033 → 731,844 | 1,063 → 798 ms |
| DC | case_ACTIVSg70k | 538,315 → 41,053 | 2,712,876 → 2,186,146 | 4,143 → 3,309 ms |
| DC | case_SyntheticUSA | 1,057,422 → 77,011 | 3,196,639 → 2,543,801 | 8,182 → 6,327 ms |
| AC | case_ACTIVSg25k | 949,064 → 4,061 | 2,626,985 → 2,027,771 | 11.7 → 9.4 s |
| AC | case_ACTIVSg70k | 3,223,912 → 12,670 | 7,826,047 → 6,157,452 | 40.2 → 33.5 s |
| AC | case_SyntheticUSA | 7,361,008 → 11,322 | 9,403,101 → 7,150,767 | 78.9 → 67.6 s |

**The study's conclusions,** which close "a fuller study of the KKT ordering and factorization":
- **On LDLᵀ, the matched ordering is better than AMD** on every measure but AC's one outcome.
- **COLAMD stays LU's,** and MIPS's, for MATPOWER parity.
- **UMFPACK-like strategies weren't pursued:** the multifrontal LDLᵀ already plays their part, at
  less than half LU's fill.

## The DC default's speed

**Mean time per run** of `profile_pf`'s pipeline for the DC predictor-corrector default:
- **Both binaries:** `9a6600e` (AMD) and `f5da62f` (matched).
- **Runs:** interleaved, two rounds, on an otherwise idle machine.

| case300 | case2383wp | case9241pegase | case_ACTIVSg25k |
|--:|--:|--:|--:|
| 1.87 → 1.62 ms (1.16×) | 79.7 → 61.5 ms (1.30×) | 79.7 → 69.3 ms (1.15×) | 1,171 → 927 ms (1.26×) |

**The dense kernel's time per run** (samply, 10 s per pipeline):
- case2383wp: 28.7 → 16.4 ms (1.75×)
- case_ACTIVSg25k: 329 → 199 ms (1.66×)

The smaller fronts, freed of their delayed pivots, make most of the difference.

## Closing O4

- **"Reduced against augmented" is closed by decision.** The solvers keep MIPS's reduced
  `[W Jgᵀ; Jg 0]`, with the inequalities folded into `W` and the right-hand side, because:
  - MIPS's parity with MATPOWER needs it
  - the augmented system is larger, by the inequalities' rows
  - it needs the same 2×2 pivoting
- **The ordering study is closed** by the table above.
- **"Parallel and GPU" and O2's untried "The CPU" ideas** move to O5, whose plan comes next:
  rayon, SIMD and GPU acceleration.
