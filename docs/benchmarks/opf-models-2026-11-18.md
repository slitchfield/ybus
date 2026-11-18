# Kept OPF models, 2026-11-18

This is the first stage of the warm-start plan in `OptimizationChecklist.md` ("Warm starts for
OPF sequences"). `AcOpfModel` and `DcOpfModel` build what a network fixes once, and solve many
cases on it, reading only values. `ac_opf_with_options` and `dc_opf_with_options` are now
`new(case)?.solve(case, options)`.

## The same results, checked

- **The 72-case OPF snapshot** is unchanged. It hashes both solvers on DC and AC: their
  solutions, iterations, points and multipliers.
- **A changed case:** a model solving a copy of its case with every load scaled by 0.95 or 1.04
  matches the one-shot solve of that copy, bit for bit. That covers case9, case30, case118 and
  case2383wp, on both solvers, AC and DC (`tests/opf_model_test.rs`).
- **Refusals:** changes the model can't solve are refused with the field named:
  - a branch out of service
  - for AC, a flow limit removed, or a quadratic cost coefficient set to 0

## What it saves on its own: little

**The measurement:** a load sweep of 11 solves, loads ×0.95 to ×1.05 in steps of 0.01, solved
by one-shot calls and by one model. Two rounds alternating; mean time per solve:

| Model and solver | Case | One-shot | Model | Speedup |
|---|---|--:|--:|--:|
| AC, predictor-corrector | case2383wp | 288.4 ms | 281.7 ms | 1.02× |
| AC, predictor-corrector | case9241pegase | 1,828.7 ms | 1,820.6 ms | 1.00× |
| AC, MIPS | case2383wp | 923.1 ms | 937.2 ms | 0.98× |
| DC, predictor-corrector | case2383wp | 58.6 ms | 58.2 ms | 1.01× |
| DC, predictor-corrector | case9241pegase | 71.5 ms | 68.9 ms | 1.04× |
| DC, predictor-corrector | case_ACTIVSg25k | 824.4 ms | 814.7 ms | 1.01× |
| DC, MIPS | case9241pegase | 325.1 ms | 324.1 ms | 1.00× |

**Within 4%, and within noise for MIPS.** The setup a model keeps, the network, the admittances
or B matrices, and #79's plans, is a small part of a solve. The solve itself, its factorizations
and evaluations, is almost all of it.

**So the models are for what they enable,** not for this saving:
- **warm starts** (the plan's next stages), which need a fixed problem to carry a solution across
- **sequences and batches** through one value, with the pattern-deciding changes refused
  explicitly rather than handled silently

**The solver's KKT analysis isn't kept across solves,** a change from the plan.
- **Why:** keeping it needs a public workspace type in `ybus_ipm`, and the analysis is about 1% of
  a solve in the profiles.
- **Its turn comes** with O5's batched solves, if the measurements there want it.
