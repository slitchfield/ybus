# AC power flow with reactive limits on every MATPOWER case, 2026-11-01

The survey for enforcing generator reactive power limits: `powerflow`'s Newton power flow on all
78 MATPOWER 8.1 cases, with limits enforced all at once (`ReactiveLimits::AllAtOnce`, MATPOWER's
`pf.enforce_q_lims = 1`) and one at a time (`OneAtATime`, `= 2`), against MATPOWER's own `runpf`
results for each mode. The survey's output, one row per case and mode, is in
[`ac-reactive-limits-2026-11-01.csv`](ac-reactive-limits-2026-11-01.csv).

## Setup

| | |
|---|---|
| Code | The `dev/reactive-limits` branch, at the commit that adds this record |
| Machine and system | As in the earlier records: Intel Core i7-14700K, Ubuntu 24.04.3 under WSL2, Rust 1.97.1, release profile |
| Reference | MATPOWER 8.1 under Octave 8.4.0, run locally by `scripts/generate-reference.sh` (results not committed): `runpf` by Newton's method with `pf.enforce_q_lims` 1 and 2, through the legacy code (`exp.use_legacy_core`), which `powerflow` follows. MP-Core's results with limits differ from it (below) |
| Run | `scripts/benchmark-ac.sh` with `YBUS_BENCH_QLIM` set to 1 (21:22–21:23 EST) and then 2 (21:23–22:40 EST), with the default limits (1 hour and 10 GiB per case) |

Each step is timed by its median over repeated runs, for up to 1 s or 101 runs; a solve longer
than 1 s runs once. An earlier run of both modes gave the same results in every column but time,
memory and the held generators of three cases, which came from the reference (below).

## Results

| | All at once | One at a time |
|---|--:|--:|
| Cases that converge | 74 | 74 |
| Cases where MATPOWER converges | 74 | 74 |
| Converged with MATPOWER's iteration count | all 74 | all 74 |
| Converged holding exactly MATPOWER's generators | all 74 | all 74 |
| Cases holding any generator | 36 | 36 |
| Iterations, median (most) | 4 (20) | 5 (14,225) |
| Worst relative difference from MATPOWER | 1.7e-10 | 1.7e-10 |
| Worst absolute difference, any quantity | 2.7e-9 | 9.1e-9 |
| Solve time, all converged cases | 6.54 s | 4,397 s |
| Slowest case | 2.93 s (case_SyntheticUSA) | 2,450 s (case_SyntheticUSA) |
| Peak memory, any case | 370 MiB | 328 MiB |

- **Convergence matches MATPOWER exactly,** in both modes, as does the iteration count, which
  totals every solve. That includes case_ACTIVSg10k's 16 iterations all at once, where MP-Core
  takes 17, and the one-at-a-time counts in the thousands on the largest cases.
- **The same four cases fail in both,** in both modes:
  - case4gs and case118zh with `ReactiveLimitsInfeasible`: every generator left at a REF or PV
    bus violates the same limit. (case118zh has one generator.)
  - case9target and case16am don't converge. Without limits, case16am fails too and case9target
    converges.
- **The same generators are held.** The comparison counts each mode's held generators and
  compares them as sets of case rows; with limits, 36 cases hold at least one, up to 7,214 on
  case_SyntheticUSA all at once.
- **The results agree with MATPOWER's** within 1.7e-10 relative, the largest absolute
  difference over the largest magnitude of that quantity in MATPOWER's results for the case, as
  in the earlier records. The worst is case141's flows.
  - **Closer than without limits:** these references come from the legacy code `powerflow`
    follows, while the earlier Newton references came from MP-Core. case141's angles, for one,
    differ by 5.0e-11° here against 6.8e-10° there.

## Held generators: what `runpf` doesn't report

`runpf` doesn't say which generators it held, so `generate_reference.m` infers them: generators
in service at a bus `runpf` made PQ, at exactly one of their limits.

- **The first run disagreed on three cases,** case2736sp, case2746wp and case2746wop, by 2 to 8
  generators, in both modes, with every other column matching.
- **The inference was wrong there.** Each extra generator has equal limits (Qmax = Qmin) and
  shares its bus with a generator both lists hold. That one made the bus PQ, and a generator with
  equal limits then sits at them whatever happened, without ever violating one.
- **The rule now leaves such a generator out** when a generator with unequal limits at its bus is
  held. The three cases were regenerated, and the lists match ours. Applied to every other
  case's list, the rule removes nothing. Rerunning case14, case30 and case118 under it gave
  committed files identical to the bytes.

## One at a time is slow on the largest cases

`OneAtATime` makes one bus PQ per round, so a case holding thousands of generators takes
thousands of rounds, as in MATPOWER:

| Case | Iterations | Rounds (held + 1) | Solve |
|---|--:|--:|--:|
| case_SyntheticUSA | 14,225 | 6,337 | 2,450 s |
| case_ACTIVSg70k | 10,993 | 4,969 | 1,699 s |
| case_ACTIVSg25k | 3,880 | 2,171 | 185 s |
| case_ACTIVSg10k | 2,612 | 1,245 | 36.6 s |
| case9241pegase | 436 | 197 | 5.94 s |

- **The four largest cases take 99.4% of the time.** The other 70 converged cases take 26.5 s
  together.
- **Each round repeats the setup.** Rounds average about 2.2 iterations, and each one builds a
  new AMD ordering, LU analysis and Jacobian assembly, since the Jacobian's pattern changes with
  the bus roles. So case_SyntheticUSA takes 172 ms per iteration, against 143 ms for an
  iteration of its Newton solve without limits in the earlier record, which spread one setup
  over 6 iterations and predates the assembly-map Jacobian.
- **A follow-up could reuse it:** each round moves one bus from PV to PQ, adding one row and
  column to the Jacobian, so much of the previous round's ordering and assembly carries over
  (`docs/PowerFlowChecklist.md`).
- **MATPOWER is slower still.** Generating its references took about 2.5 hours for
  case_SyntheticUSA one at a time and 1.5 hours for case_ACTIVSg70k. That is incidental: the
  generation wasn't set up as a benchmark, and Octave spent most of its CPU time in the kernel,
  coordinating threads over small steps.

## MP-Core against legacy

**All at once:** MP-Core and the legacy code were compared on every case MATPOWER solves:

- They agree within 5.1e-12 on case14, case30 and case118, the committed cases.
- Over all cases they differ, by up to 71 MW on case_SyntheticUSA. case6470rte converges only in
  legacy, and case_ACTIVSg10k takes 17 iterations in MP-Core against 16.

By the rule agreed for this work, the references therefore come from the legacy code.

**One at a time** (rerun 2026-11-02, after a machine restart interrupted the first attempt):
each case ran once in MP-Core and once in legacy, under single-threaded BLAS, an address-space
limit and a 5-hour timeout per run.

| Outcome | Cases |
|---|--:|
| Both converge | 69 |
| Neither converges (case4gs, case9target, case16am, case118zh) | 4 |
| Only legacy converges (case6470rte) | 1 |
| Legacy converges, MP-Core runs out of memory (case_ACTIVSg10k, case_ACTIVSg25k, case_ACTIVSg70k, case_SyntheticUSA) | 4 |

- **Where both converge:** the final bus types agree in all 69. 53 agree within 1e-9. Of the
  other 16, two take different iteration counts (case2848rte 184 in MP-Core against 183 in
  legacy, case3012wp 482 against 481), and case2868rte differs most, by 0.040 MW/MVAr in flows
  and generator outputs and 0.0054° in angles. The rest are within 8.4e-9. case14, case30 and
  case118 agree within 5.1e-12, at the same iteration counts.
- **The four largest cases run MP-Core out of memory under Octave.** Legacy converged on each,
  with the iteration counts above (2,612 to 14,225), in 98 s to 92 minutes. MP-Core failed on
  each within about 3.5 minutes at a 12 GiB limit, and its resident memory on case_ACTIVSg10k
  rose steadily from 0.34 GB to 11.3 GB in 207 s, while legacy's stayed near 0.3 GB.
- **The cause is a leak, not the problem's size.** MP-Core's math model stores a closure that
  captures the model itself (`add_node_balance_constraints` in `math_model_pf_acps.m`), a
  reference cycle among handle objects, which Octave's reference counting never frees. One
  plain `runpf` call leaks 2.7 MB on case_ACTIVSg2000, against nothing for legacy, and one at a
  time rebuilds the models every round. Without the memory limit it would exhaust this
  machine's 16 GB, so these four can't be compared here. As of 2026-11-02 the leak isn't
  reported in MATPOWER's issues, pull requests, discussions or change log.

So in both modes MP-Core's results with limits differ from legacy's, and one at a time it can't
run the largest cases under Octave: the legacy references stand.
