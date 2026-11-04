# MATPOWER power flow and optimal power flow reference results

Expected results for the `powerflow` power flow tests, computed by MATPOWER: the DC power flow
by `rundcpf`, and the AC power flow by `runpf` with its default options (Newton's method,
tolerance 1e-8 on the largest mismatch, reactive power limits not enforced), with its
fast-decoupled XB and BX methods, and by Newton's method with reactive power limits enforced
(`pf.enforce_q_lims` 1 and 2); and the optimal power flow, DC by `rundcopf` and AC by `runopf`,
both solved by MIPS. The versions used are in [`VERSIONS`](VERSIONS).

The optimal power flow results use MATPOWER's defaults (MP-Core), which agree with the legacy
code on OPF: the same iteration counts on case9, case14, case30, case118, case300 and
case2383wp, AC and DC, with objectives within 6e-14 relative and voltages, dispatch and prices
within 3.5e-9 (measured while planning phase O1). The DC runs force MIPS
(`opf.dc.solver = 'MIPS'`), since `powerflow` reproduces MIPS: by default MATPOWER hands a DC OPF
with linear costs to GLPK. The objective with the default solver is recorded too, as a check of
the optimum. A run where MIPS fails is still recorded, with `success` 0, since `powerflow`
should fail the same way.

MATPOWER 8.1's `runpf` uses its newer MP-Core implementation by default, under Octave too. The
Newton results come from it, and agree with the legacy `newtonpf` within 6.4e-12 (case118). The
fast-decoupled results come from the legacy `fdpf` (`exp.use_legacy_core`), which `powerflow`
follows: MP-Core's fast-decoupled results differ from it by up to 2.3e-7 (case14, FDXB) at the
same iteration counts. The reactive-limit results also come from the legacy code: on these
three cases the two agree within 5.1e-12, but over all of MATPOWER's cases MP-Core's differ (all
at once by up to 71 MW on case_SyntheticUSA, and in convergence on case6470rte; one at a time by
up to 0.040 MW, and on the four largest cases MP-Core runs out of memory under Octave).

These files are **generated**. Don't edit them by hand; regenerate them with
`scripts/generate-reference.sh`, which fetches MATPOWER (`scripts/fetch-matpower.sh`) and runs
`scripts/generate_reference.m` under Octave. The output is deterministic, so a regeneration
with the same versions leaves no diff. `YBUS_REFERENCE_RUNS` chooses which results to make
(for example `islands`), so new ones don't need the slow reactive-limit runs again.

## Files

For each case:

| File | Columns |
|---|---|
| `<case>_dcpf_bus.csv` | `bus` (bus number), `va_deg` (voltage angle, degrees) |
| `<case>_dcpf_branch.csv` | `row` (1-based row in the case file's branch matrix), `from`, `to`, `status`, `pf_mw` (real power flow at the from end, MW) |
| `<case>_dcpf_gen.csv` | `row` (1-based row in the gen matrix), `bus`, `pg_mw` (real power output, MW) |
| `<case>_acpf_bus.csv` | `bus`, `vm_pu` (voltage magnitude, p.u.), `va_deg` (voltage angle, degrees) |
| `<case>_acpf_branch.csv` | `row`, `from`, `to`, `status`, `pf_mw`, `qf_mvar` (power flowing into the branch at the from end), `pt_mw`, `qt_mvar` (at the to end) |
| `<case>_acpf_gen.csv` | `row`, `bus`, `pg_mw`, `qg_mvar` (real and reactive power output) |
| `<case>_acpf_summary.csv` | `iterations` (Newton iterations), `success` (1 if `runpf` converged) |
| `<case>_fdxb_*.csv`, `<case>_fdbx_*.csv` | The same four tables for fast-decoupled XB and BX; `iterations` counts P half-iterations |
| `<case>_qlim1_*.csv`, `<case>_qlim2_*.csv` | The same four tables for Newton's method with reactive limits enforced, all at once (1) or one at a time (2); `iterations` totals every solve |
| `<case>_islands.csv` | `bus`, `island`: each bus's island from `find_islands` (1-based, in MATPOWER's order, largest first), or 0 for a bus with no in-service branch |
| `<case>_bridges.csv` | `row`: the branches `find_bridges` reports as bridges (1-based rows, increasing). Only where `find_bridges` succeeds: case14 and case30, but not case118, whose parallel branches make it fail |
| `<case>_ptdf.csv` | `makePTDF`'s PTDF at the default slack, after `ext2int`: a header `row` and the columns' bus numbers, then one line per in-service branch, its 1-based row and its entries. case14 and case30 only: case118's would be 1.6 MB |
| `<case>_ptdf_distributed.csv` | The same with the slack distributed in proportion to each bus's load |
| `<case>_lodf.csv` | `makeLODF`'s LODF, unmasked: a header `row` and the outaged branches' rows, then one line per in-service branch. A bridge's column is rounding noise there (division by 1 − h ≈ 0), which the tests check by `powerflow`'s own arithmetic instead |
| `<case>_dcopf_bus.csv` | `bus`, `va_deg`, `lam_p` (the bus's real power price, $/MWh) |
| `<case>_dcopf_gen.csv` | `row`, `bus`, `status`, `pg_mw`, `mu_pmin`, `mu_pmax` (multipliers of the real power limits, $/MWh) |
| `<case>_dcopf_branch.csv` | `row`, `from`, `to`, `status`, `pf_mw`, `mu_sf`, `mu_st` (multipliers of the flow limit at the from and to ends, $/MWh), `mu_angmin`, `mu_angmax` (of the angle difference limits, $/degree) |
| `<case>_dcopf_summary.csv` | `objective` ($/hr), `iterations` (MIPS's), `success`, `default_objective` and `default_solver` (`rundcopf` with its default solver; NaN if it failed) |
| `<case>_dcopf_history.csv` | `iteration` (0 for the starting point), then MIPS's `feascond`, `gradcond`, `compcond`, `costcond`, `gamma`, `stepsize`, `obj`, `alphap` and `alphad` after it |
| `<case>_opf_bus.csv` | `bus`, `vm_pu`, `va_deg`, `lam_p`, `lam_q` ($/MVAr-hr), `mu_vmin`, `mu_vmax` (multipliers of the voltage magnitude limits, $/p.u.) |
| `<case>_opf_gen.csv` | `row`, `bus`, `status`, `pg_mw`, `qg_mvar`, `mu_pmin`, `mu_pmax`, `mu_qmin`, `mu_qmax` |
| `<case>_opf_branch.csv` | `row`, `from`, `to`, `status`, `pf_mw`, `qf_mvar`, `pt_mw`, `qt_mvar`, `mu_sf`, `mu_st` (multipliers of the apparent power limits, $/MVA), `mu_angmin`, `mu_angmax` |
| `<case>_opf_summary.csv`, `<case>_opf_history.csv` | As for the DC runs, without the default-solver columns |
| `<case>_qlim1_limited.csv`, `<case>_qlim2_limited.csv` | `row`: the generators held at a reactive limit (1-based rows, increasing), found as those in service at a bus `runpf` made PQ, at exactly one of their limits, since `runpf` doesn't report them. A generator with equal limits is left out when one with unequal limits at its bus is held: it sits at its limits whenever its bus is PQ |

Values have 17 significant digits, so they parse back to the exact doubles MATPOWER computed.
In a DC power flow the flow at the to end is `-pf_mw`; in an AC one, `pf_mw + pt_mw` is the
branch's loss. The DC results came out byte-for-byte the same when regenerated under Octave
8.4.0, the version that added the AC results.

## Sources and attribution

The results were computed with MATPOWER:

> R. D. Zimmerman, C. E. Murillo-Sánchez and R. J. Thomas, "MATPOWER: Steady-state operations,
> planning, and analysis tools for power systems research and education," *IEEE Trans. Power
> Syst.* 26(1):12–19, 2011. [doi:10.1109/TPWRS.2010.2051168](https://doi.org/10.1109/TPWRS.2010.2051168)

> R. D. Zimmerman and C. E. Murillo-Sánchez, MATPOWER, Version 8.1 [software], 2025.
> <https://matpower.org>, [doi:10.5281/zenodo.15871662](https://doi.org/10.5281/zenodo.15871662)

MATPOWER's `CITATION` file asks results derived from it to cite the 2011 paper, and recommends
citing the specific software version for reproducibility; both are given above.

MATPOWER's code is under the BSD 3-clause license. Its case files are not covered by that
license, so they are not included in this repository; the fetch script downloads them. Only
results are committed here, and only for cases that come from public sources, as recorded in
each case file's header:

| Case | Origin |
|---|---|
| `case14` | IEEE 14-bus test case (1962), converted by MATPOWER from the IEEE Common Data Format file in the University of Washington Power Systems Test Case Archive, <https://labs.ece.uw.edu/pstca/>. |
| `case30` | Compiled by MATPOWER's authors from O. Alsaç and B. Stott, "Optimal load flow with steady-state security," *IEEE Trans. Power App. Syst.* PAS-93(3):745–751, 1974, [doi:10.1109/TPAS.1974.293972](https://doi.org/10.1109/TPAS.1974.293972), and R. W. Ferrero, S. M. Shahidehpour and V. C. Ramesh, "Transaction analysis in deregulated power systems using game theory," *IEEE Trans. Power Syst.* 12(3):1340–1347, 1997, [doi:10.1109/59.630479](https://doi.org/10.1109/59.630479), with modifications described in the case file. |
| `case118` | IEEE 118-bus test case (1961), from the same University of Washington archive, with base kV data and two transformer corrections added by MATPOWER's authors. |

The optimal power flow results also depend on each case's costs and limits. case30's come from
Ferrero, Shahidehpour and Ramesh, cited above. The IEEE Common Data Format has no cost data, so
case14's and case118's generator costs were added by MATPOWER's authors, and the case files
don't say where they came from. As with the power flow results, only values MATPOWER computed
are committed, not the case data.

The larger `case2383wp` (Polish grid), used for benchmarks, is distributed by MATPOWER with
its data provider's permission. No results derived from it are committed; anything that needs
it fetches MATPOWER.
