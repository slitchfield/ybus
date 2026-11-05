# ybus

Sparse linear algebra, power flow and interior-point optimization for Rust. A Cargo workspace
with three crates:

| Crate | Status | Contents |
|---|---|---|
| [`sparse_mat`](sparse_mat/) | Usable | Sparse matrices: a hash-map type for building and editing, and compressed (CSR/CSC) types with fixed patterns for arithmetic, products, fill-reducing orderings (AMD, RCM, COLAMD), sparse Cholesky and LU (complex too, with refactoring, sparse right-hand sides and transposed solves), iterative refinement and condition estimates, conjugate gradient, Matrix Market I/O, optional `serde` and `rayon` support. Generic over the element type, including complex numbers. |
| [`powerflow`](powerflow/) | In progress | DC and AC power flow on top of `sparse_mat`. Done: the case model, a MATPOWER case reader, DC power flow, which agrees with MATPOWER on all 78 of its cases, and AC power flow by Newton–Raphson and fast-decoupled (XB and BX), optionally enforcing generator reactive limits, which match MATPOWER's `runpf` on its reference cases. On every MATPOWER case, it converges exactly where MATPOWER does, in the same number of iterations, and with limits holds the same generators. It also lists a network's islands and bridges, computes PTDF and LODF from a kept factorization of B′, as MATPOWER does, and solves branch outages by rank-k updates of that factorization, matching a refactored solve, optionally in parallel (`rayon`). Also the optimal power flow, DC and AC, solved by `interior_point`'s MIPS, which matches MATPOWER's `rundcopf` and `runopf` on their reference cases, in the same number of iterations; over every case with costs, recorded in [`docs/benchmarks/`](docs/benchmarks/), it ends as MATPOWER does on every AC case and all but one DC case. |
| [`interior_point`](interior_point/) | In progress | Interior-point solvers on top of `sparse_mat`. Done: MIPS, MATPOWER's primal-dual interior-point solver, reimplemented step for step from `mips.m`, with its quadratic-program wrapper. On the problems of MATPOWER's MIPS tests, it ends as MATPOWER does in the same number of iterations, and its results agree within 1.0e-13. |

The name is the bus admittance matrix, Y<sub>bus</sub>: the sparse matrix every power flow
formulation is built around.

See [`sparse_mat/README.md`](sparse_mat/README.md) for the matrix library's quick start and
API overview, [`docs/PowerFlowChecklist.md`](docs/PowerFlowChecklist.md) for the completed power
flow roadmap, and [`docs/OptimizationChecklist.md`](docs/OptimizationChecklist.md) for the next
one: optimal power flow, interior-point solvers and acceleration.

## Development

- **MSRV:** Rust 1.88, set in `[workspace.package]` in the root `Cargo.toml` and checked by a CI
  job.
- **`./check.sh`:** runs everything the Check CI job runs, over the whole workspace: fmt,
  clippy with `-D warnings`, and the tests with no features, all features and each
  `sparse_mat` feature alone, plus doc tests (including the `sparse_mat` README's example).
  The benchmarks run once, as tests, in the first of those passes only.
- **Dependencies:** `cargo deny --all-features check` checks advisories, licenses, banned
  crates and sources against `deny.toml`. The Dependencies CI job runs it.
- **MATPOWER test data:** `scripts/fetch-matpower.sh` downloads a pinned MATPOWER release
  (checksum verified) into the gitignored `testdata/`. It needs `curl` and `unzip`. Set
  `MATPOWER_DIR` to use an existing installation instead.
  - **How tests use it:** `powerflow` tests find cases through `powerflow/tests/common`. They
    skip when MATPOWER hasn't been fetched. Set `YBUS_REQUIRE_CASES=1` to make missing cases
    fail instead; the Check CI job fetches MATPOWER (cached) and sets it.
  - **Cases that compute their data with code:** some MATPOWER cases rescale their data with
    MATLAB code, which the `powerflow` reader rejects rather than misread.
    `scripts/normalize-cases.sh` evaluates every case with MATPOWER's own `loadcase` and writes
    it out as plain literals, at full precision, into the gitignored `testdata/`. It needs
    Octave (for example `apt install octave`). Tests find these cases through
    `common::normalized_case`. The Check CI job installs Octave and normalizes the cases only
    when its cache of them is out of date.
  - **Regenerating the results:** `scripts/generate-reference.sh` regenerates the committed DC
    power flow reference results in
    [`powerflow/tests/data/matpower-reference/`](powerflow/tests/data/matpower-reference/). It
    needs Octave. Given an output directory and cases, it generates MATPOWER's results for any
    case instead, for local comparisons; never commit those.
- **Benchmarks:** `cargo bench -p powerflow` times, with criterion, the DC and AC power flows on
  case2383wp and reassembling matrices through an assembly map, and `scripts/benchmark-dc.sh`
  and `scripts/benchmark-ac.sh` survey every MATPOWER case once, comparing the results with
  MATPOWER's. [`docs/benchmarks/`](docs/benchmarks/) records the results.

## Documentation

- [`docs/PowerFlowChecklist.md`](docs/PowerFlowChecklist.md) is the completed roadmap toward
  DC/AC power flow.
- [`docs/OptimizationChecklist.md`](docs/OptimizationChecklist.md) is the roadmap that follows:
  optimal power flow by a reimplementation of MATPOWER's MIPS, a simpler interior-point solver of
  our own, and accelerating the power flow and both solvers.
- [`docs/LongTermGoals.md`](docs/LongTermGoals.md) records the goals beyond these roadmaps, with
  coarse plans: foreign function interfaces to C, Python and Julia, benchmarks against other
  power flow libraries, security-constrained OPF, and reading ARPA-E's SC-ACOPF data.
- [`docs/BaseFunctionalityChecklist.md`](docs/BaseFunctionalityChecklist.md) is the completed
  checklist of `sparse_mat`'s core functionality.
- [`docs/benchmarks/`](docs/benchmarks/) records performance baselines, such as the DC power
  flow on every MATPOWER case.
- [`AI-DISCLAIMER.md`](AI-DISCLAIMER.md) describes how AI was used to build this repository and
  how to report missing attribution.
- API documentation: `cargo doc --workspace --open --all-features`.

## License

Licensed under either of

- Apache License, Version 2.0 ([`LICENSE-APACHE`](LICENSE-APACHE) or
  <https://www.apache.org/licenses/LICENSE-2.0>)
- MIT license ([`LICENSE-MIT`](LICENSE-MIT) or <https://opensource.org/licenses/MIT>)

at your option. See [`AI-DISCLAIMER.md`](AI-DISCLAIMER.md) for how this licensing relates to the
use of AI in developing the crate.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in
this work by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without
any additional terms or conditions.

## References

The sparse matrix algorithms in `sparse_mat` were written from the standard published
descriptions below; no source code from other libraries was copied or ported.

`powerflow` reproduces MATPOWER's behaviour, so that its results can be checked against
MATPOWER's to within 1e-9. Its DC power flow (`powerflow::dc`) was written from MATPOWER's
documentation and with reference to its BSD-licensed source code (`ext2int`, `bustypes`,
`makeBdc`, `dcpf` and `runpf`), to match its handling of every case, including isolated
buses, out-of-service elements and phase shifters. Its AC power flow (`powerflow::ac`) was
written the same way, with reference to the source of `makeYbus`, `newtonpf`, `makeB`, `fdpf`
and `pfsoln` as well, and to `runpf`'s loop for enforcing reactive limits. Its island and
bridge listing (`powerflow::islands`) follows `find_islands`, `connected_components` and
`find_bridges`, except that parallel branches count as separate edges, which MATPOWER's
`find_bridges` mishandles. Its PTDF and LODF (`DcSensitivity`) follow `makePTDF` and
`makeLODF`, with bridges' LODF columns NaN, as `makeLODF`'s `mask_bridge` makes them. Two of its pieces come from MATPOWER's published formulas instead, since
their source wasn't consulted: the bus injections (what `makeSbus` computes, as the User's
Manual defines them) and the Newton Jacobian's derivatives (what `dSbus_dV` computes, from
Technical Note 2, cited below). Both power flows are reimplementations in Rust, not
translations, and the module documentation names the MATPOWER functions each step follows.

`interior_point`'s MIPS (`interior_point::mips`) is different: it is ported step for step from
MATPOWER's `mips.m` and `qps_mips.m` (MIPS 1.5.2), so that its iterations can be checked
against MATPOWER's one by one, and so it is derived from them and carries MIPS's BSD notice, as
do its tests and the script that records MIPS's results on the same problems
([`LICENSE-MIPS`](LICENSE-MIPS) lists the files). The optimal power flow's generator costs
(`powerflow::opf::cost`) follow `opf_setup`, `makeAy`, `pqcost`, `totcost` and `polycost`, and
its DC optimal power flow (`powerflow::opf::dc`) follows `opf_setup`'s DC model, `makeAang`,
`dcopf_solver`'s starting point and multiplier mapping, and `opf`'s reporting of out-of-service
elements. Its AC optimal power flow (`powerflow::opf::ac`) follows `opf_setup`'s AC model for MIPS,
`makeAvl`, and `nlpopf_solver`'s starting point and multiplier mapping, with the constraints'
first and second derivatives written from Technical Note 2's formulas (cited below), not from
MATPOWER's `dSbus_dV`, `d2Sbus_dV2` and related functions. All three are reimplementations, as
the power flows are.

MATPOWER 8.1 runs its newer MP-Core implementation by default, while `powerflow` follows the
legacy functions named above. The committed DC and Newton reference results come from
MP-Core; for Newton's method it agrees with the legacy code within 6.4e-12 on case118. Its
fast-decoupled results differ from the legacy `fdpf`'s by up to 2.3e-7, so the fast-decoupled
reference results come from `fdpf`. With reactive limits enforced, MP-Core differs from the
legacy `runpf` over all cases in both modes (by up to 71 MW), and one at a time it runs out of
memory on the largest cases under Octave, so those reference results come from the legacy code
too (see
[the reference results' README](powerflow/tests/data/matpower-reference/README.md)). Across all
MATPOWER cases, MP-Core's Newton post-processing also shares generator reactive output
differently from the legacy `pfsoln` in 10 cases
([AC survey](docs/benchmarks/ac-power-flow-2026-11-01.md)); there `powerflow` matches `pfsoln`.

### Implemented

| Where (in `sparse_mat` unless a crate is named) | Algorithm | Reference |
|---|---|---|
| `solve_lower_triangular`, `solve_upper_triangular`, `Lu::solve`, `Lu::determinant` | Forward and back substitution; LU factorization with partial pivoting | G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013, ch. 3. |
| `lu`, `analyze_lu`, `SymbolicLu::factor` | Left-looking sparse LU: a depth-first search for each column's reach, then a sparse triangular solve, with threshold partial pivoting | J. R. Gilbert and T. Peierls, "Sparse partial pivoting in time proportional to arithmetic operations," *SIAM J. Sci. Stat. Comput.* 9(5):862–874, 1988. [doi:10.1137/0909058](https://doi.org/10.1137/0909058) |
| | | T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006, ch. 6. |
| `Lu::refactor` | Refactoring with the previous pivot sequence and factor patterns, checking each reused pivot | T. A. Davis and E. Palamadai Natarajan, "Algorithm 907: KLU, a direct sparse solver for circuit simulation problems," *ACM Trans. Math. Softw.* 37(3), article 36, 2010. [doi:10.1145/1824801.1824814](https://doi.org/10.1145/1824801.1824814) |
| `Cholesky::solve_lower_sparse`, `Lu::solve_lower_sparse`, `solve_sparse` | Sparse right-hand sides: the forward solve computes only the reach of the right-hand side's entries in the factor's graph | J. R. Gilbert, "Predicting structure in sparse matrix computations," *SIAM J. Matrix Anal. Appl.* 15(1):62–79, 1994. [doi:10.1137/S0895479887139455](https://doi.org/10.1137/S0895479887139455) |
| `cholesky`, `analyze_cholesky` | Up-looking (row-by-row) sparse Cholesky, with symbolic analysis by the elimination tree and row subtrees | T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006, ch. 4. |
| | | J. W. H. Liu, "The role of elimination trees in sparse factorization," *SIAM J. Matrix Anal. Appl.* 11(1):134–172, 1990. |
| `Lu::solve_refined`, `Cholesky::solve_refined` | Fixed-precision iterative refinement, measured by the componentwise backward error | R. D. Skeel, "Iterative refinement implies numerical stability for Gaussian elimination," *Math. Comp.* 35(151):817–832, 1980. [doi:10.1090/S0025-5718-1980-0572859-4](https://doi.org/10.1090/S0025-5718-1980-0572859-4) |
| | | W. Oettli and W. Prager, "Compatibility of approximate solution of linear equations with given error bounds for coefficients and right-hand sides," *Numer. Math.* 6(1):405–409, 1964. [doi:10.1007/BF01386090](https://doi.org/10.1007/BF01386090) |
| | | N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, 2002, ch. 12. [doi:10.1137/1.9780898718027](https://doi.org/10.1137/1.9780898718027) |
| `inverse_norm_1_estimate`, `condition_estimate_1` (on `Lu<f64>` and `Cholesky<f64>`) | Estimating ‖A⁻¹‖₁ from a few solves with A and Aᵀ, as LAPACK's `xLACON` does | W. W. Hager, "Condition estimates," *SIAM J. Sci. Stat. Comput.* 5(2):311–316, 1984. [doi:10.1137/0905023](https://doi.org/10.1137/0905023) |
| | | N. J. Higham, "FORTRAN codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation," *ACM Trans. Math. Softw.* 14(4):381–396, 1988. [doi:10.1145/50063.214386](https://doi.org/10.1145/50063.214386) |
| `conjugate_gradient` | Conjugate gradient | M. R. Hestenes and E. Stiefel, "Methods of conjugate gradients for solving linear systems," *J. Res. Natl. Bur. Stand.* 49(6):409–436, 1952. [doi:10.6028/jres.049.044](https://doi.org/10.6028/jres.049.044) |
| | | Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, 2003. |
| | | J. R. Shewchuk, "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain," Tech. Rep. CMU-CS-94-125, Carnegie Mellon University, 1994. |
| `matmul`, `par_matmul`; `CompressedMatrix::from_triplets`, `to_csr`/`to_csc` and `SparseMatrix::to_csr`/`to_csc` | Row-by-row sparse matrix product; counting-sort transposition into compressed form, also used to group triplets and to reorient between CSR and CSC | F. G. Gustavson, "Two fast algorithms for sparse matrices: multiplication and permuted transposition," *ACM Trans. Math. Softw.* 4(3):250–269, 1978. [doi:10.1145/355791.355796](https://doi.org/10.1145/355791.355796) |
| `Ordering::Rcm` | Reverse Cuthill–McKee, from a pseudo-peripheral node found by George and Liu's method | E. Cuthill and J. McKee, "Reducing the bandwidth of sparse symmetric matrices," *Proc. 24th ACM National Conference*, pp. 157–172, 1969. |
| | | A. George and J. W. H. Liu, *Computer Solution of Large Sparse Positive Definite Systems*, Prentice-Hall, 1981. |
| `Ordering::Amd` | Approximate minimum degree on the quotient graph, with aggressive absorption and supervariables | P. R. Amestoy, T. A. Davis and I. S. Duff, "An approximate minimum degree ordering algorithm," *SIAM J. Matrix Anal. Appl.* 17(4):886–905, 1996. [doi:10.1137/S0895479894278952](https://doi.org/10.1137/S0895479894278952) |
| | | T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006, ch. 7. |
| `Ordering::Colamd` | Column approximate minimum degree: a column ordering for LU from the row–column structure, with row elements, approximate scores, aggressive absorption and supercolumns | T. A. Davis, J. R. Gilbert, S. I. Larimore and E. G. Ng, "A column approximate minimum degree ordering algorithm," *ACM Trans. Math. Softw.* 30(3):353–376, 2004. [doi:10.1145/1024074.1024079](https://doi.org/10.1145/1024074.1024079) |
| | | T. A. Davis, J. R. Gilbert, S. I. Larimore and E. G. Ng, "Algorithm 836: COLAMD, a column approximate minimum degree ordering algorithm," *ACM Trans. Math. Softw.* 30(3):377–380, 2004. [doi:10.1145/1024074.1024080](https://doi.org/10.1145/1024074.1024080) |
| Matrix Market I/O | Coordinate exchange format | R. F. Boisvert, R. Pozo and K. A. Remington, "The Matrix Market Exchange Formats: Initial Design," NISTIR 5935, NIST, 1996. [doi:10.6028/NIST.IR.5935](https://doi.org/10.6028/NIST.IR.5935) |
| `powerflow`'s tests and `dc_survey` example (`tests/common/accuracy.rs`) | Normwise backward error, η = ‖b − Ax‖ / (‖A‖ ‖x‖ + ‖b‖) | J. L. Rigal and J. Gaches, "On the compatibility of a given solution with the data of a linear system," *J. ACM* 14(3):543–548, 1967. |
| | | N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, 2002, ch. 7. [doi:10.1137/1.9780898718027](https://doi.org/10.1137/1.9780898718027) |
| `powerflow::islands` (`find_bridges`) | Bridges by depth-first search, iteratively, with parallel branches as separate edges | R. E. Tarjan, "A note on finding the bridges of a graph," *Inf. Process. Lett.* 2(6):160–161, 1974. [doi:10.1016/0020-0190(74)90003-9](https://doi.org/10.1016/0020-0190(74)90003-9) |
| `powerflow::dc` (`DcSensitivity`: `ptdf_*`, `lodf_columns`, `lodf`) | Power transfer and line outage distribution factors, from transfer solves against B′ | A. J. Wood, B. F. Wollenberg and G. B. Sheblé, *Power Generation, Operation, and Control*, 3rd ed., Wiley, 2014, ch. 7 and app. 7B. |
| | | T. Güler, G. Gross and M. Liu, "Generalized line outage distribution factors," *IEEE Trans. Power Syst.* 22(2):879–881, 2007. [doi:10.1109/TPWRS.2006.888950](https://doi.org/10.1109/TPWRS.2006.888950) |
| `powerflow::dc` (`DcStudy::solve_outage`, `solve_outages`) | Outages by rank-k updates of B′'s kept factorization, the Sherman–Morrison–Woodbury identity | J. Sherman and W. J. Morrison, "Adjustment of an inverse matrix corresponding to a change in one element of a given matrix," *Ann. Math. Statist.* 21(1):124–127, 1950. [doi:10.1214/aoms/1177729893](https://doi.org/10.1214/aoms/1177729893) |
| | | M. A. Woodbury, "Inverting modified matrices," Memorandum Report 42, Statistical Research Group, Princeton University, 1950. |
| | | W. W. Hager, "Updating the inverse of a matrix," *SIAM Rev.* 31(2):221–239, 1989. [doi:10.1137/1031049](https://doi.org/10.1137/1031049) |
| `powerflow::case` (`Case::parse`) | The MATPOWER case format, as `caseformat.m` defines it | R. D. Zimmerman and C. E. Murillo-Sánchez, *MATPOWER User's Manual*, Version 8.1, 2025. [doi:10.5281/zenodo.15871661](https://doi.org/10.5281/zenodo.15871661) |
| `powerflow::ac` (`AcModel::jacobian`, Newton's method) and `powerflow::opf::ac` | The power flow equations' derivatives in complex matrix form, `dS/dVa` and `dS/dVm`, and for the AC optimal power flow the second derivatives of the bus injections and of the branch flows' squared apparent power | R. D. Zimmerman, "AC power flows, generalized OPF costs and their derivatives using complex matrix notation," MATPOWER Technical Note 2, February 2010. [doi:10.5281/zenodo.3237866](https://doi.org/10.5281/zenodo.3237866) |
| `powerflow::ac` (`AcMethod::FastDecoupledXb`, `FastDecoupledBx`) | Fast-decoupled power flow, XB and BX, as MATPOWER's `fdpf` and `makeB` compute it | B. Stott and O. Alsaç, "Fast decoupled load flow," *IEEE Trans. Power App. Syst.* PAS-93(3):859–869, 1974. |
| `interior_point::mips` (`mips`, `qps_mips`) | MIPS, MATPOWER's primal-dual interior-point solver, as its `mips.m` and `qps_mips.m` implement it (see above) | H. Wang, C. E. Murillo-Sánchez, R. D. Zimmerman and R. J. Thomas, "On computational issues of market-based optimal power flow," *IEEE Trans. Power Syst.* 22(3):1185–1193, 2007. [doi:10.1109/TPWRS.2007.901301](https://doi.org/10.1109/TPWRS.2007.901301) |
| `powerflow::dc`, `powerflow::ac`, `powerflow::islands` and `powerflow::opf`, and the reference results their tests compare against | DC power flow as MATPOWER's `rundcpf` computes it, AC power flow as its `runpf` does, PTDF and LODF as its `makePTDF` and `makeLODF` do, islands and bridges as its `find_islands` and `find_bridges` find them, and generator costs and the DC and AC optimal power flows as its `opf_setup`, `makeAy`, `makeAang`, `makeAvl`, `dcopf_solver` and `nlpopf_solver` build and solve them (see above) | R. D. Zimmerman, C. E. Murillo-Sánchez and R. J. Thomas, "MATPOWER: Steady-state operations, planning, and analysis tools for power systems research and education," *IEEE Trans. Power Syst.* 26(1):12–19, 2011. [doi:10.1109/TPWRS.2010.2051168](https://doi.org/10.1109/TPWRS.2010.2051168) |
| | | R. D. Zimmerman and C. E. Murillo-Sánchez, MATPOWER, Version 8.1 [software], 2025. <https://matpower.org>, [doi:10.5281/zenodo.15871662](https://doi.org/10.5281/zenodo.15871662) |

### Planned (see [`docs/OptimizationChecklist.md`](docs/OptimizationChecklist.md))

| Planned work | Reference |
|---|---|
| A simpler primal-dual interior-point solver of our own | J. Nocedal and S. J. Wright, *Numerical Optimization*, Springer, 2006. [doi:10.1007/978-0-387-40065-5](https://doi.org/10.1007/978-0-387-40065-5) |
| | S. Mehrotra, "On the implementation of a primal-dual interior point method," *SIAM J. Optim.* 2(4):575–601, 1992. [doi:10.1137/0802028](https://doi.org/10.1137/0802028) |

### Licenses of related software

Algorithms are not covered by copyright, so implementing them from the publications above
carries no license obligations. Porting the source code of a reference implementation would.
The current terms of the most relevant implementations:

- **SuiteSparse AMD, CAMD, COLAMD and CCOLAMD:** BSD 3-clause. Porting is possible if their
  copyright notices are kept.
- **SuiteSparse CSparse, CXSparse, KLU and BTF:** LGPL 2.1 or later. Porting their code would
  bring LGPL obligations, so implement from the papers instead.
- **MATPOWER:** the code is BSD 3-clause. Its `LICENSE` says the case files are not covered by
  that license, so this repository doesn't include them: `scripts/fetch-matpower.sh` downloads
  them. Results computed from the public-source cases (case14, case30, case118) are committed,
  with each case's origin credited in
  [their README](powerflow/tests/data/matpower-reference/README.md). MIPS, MATPOWER's solver,
  is under the same BSD 3-clause license; `interior_point::mips` is ported from it, so it and
  the files that restate MIPS's tests keep MIPS's notice ([`LICENSE-MIPS`](LICENSE-MIPS)).

Before relying on these, check them against each project's current `LICENSE` file.
