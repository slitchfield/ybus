# ybus

Sparse linear algebra and power flow for Rust. A Cargo workspace with two crates:

| Crate | Status | Contents |
|---|---|---|
| [`sparse_mat`](sparse_mat/) | Usable | Sparse matrices: construction, arithmetic, products, LU/Cholesky/CG solvers, Matrix Market I/O, optional `serde` and `rayon` support. Generic over the element type, including complex numbers. |
| [`powerflow`](powerflow/) | In progress | DC and AC power flow on top of `sparse_mat`. Done: the case model, a MATPOWER case reader, and DC power flow, which agrees with MATPOWER on all 78 of its cases. Next: Newton–Raphson and fast-decoupled AC power flow. |

The name is the bus admittance matrix, Y<sub>bus</sub>: the sparse matrix every power flow
formulation is built around.

See [`sparse_mat/README.md`](sparse_mat/README.md) for the matrix library's quick start and
API overview, and [`docs/PowerFlowChecklist.md`](docs/PowerFlowChecklist.md) for the roadmap.

## Development

- **MSRV:** Rust 1.87, set in `[workspace.package]` in the root `Cargo.toml` and checked by a CI
  job.
- **`./check.sh`:** runs everything the Check CI job runs, over the whole workspace: fmt,
  clippy with `-D warnings`, and the tests with no features, all features and each
  `sparse_mat` feature alone, plus doc tests (including the `sparse_mat` README's example).
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
- **Benchmarks:** `cargo bench -p powerflow` times the DC power flow on case2383wp with
  criterion, and `scripts/benchmark-dc.sh` surveys every MATPOWER case once, comparing the
  results with MATPOWER's. [`docs/benchmarks/`](docs/benchmarks/) records the results.

## Documentation

- [`docs/PowerFlowChecklist.md`](docs/PowerFlowChecklist.md) is the roadmap toward DC/AC power
  flow.
- [`docs/BaseFunctionalityChecklist.md`](docs/BaseFunctionalityChecklist.md) is the completed
  checklist of `sparse_mat`'s core functionality.
- [`docs/benchmarks/`](docs/benchmarks/) records performance baselines, such as the DC power
  flow on every MATPOWER case.
- [`docs/ai-reports/`](docs/ai-reports/) holds implementation reports for past branches.
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

The algorithms in this repository were written from the standard published descriptions below. No
source code from other libraries was copied or ported.

### Implemented

| Where (in `sparse_mat` unless noted) | Algorithm | Reference |
|---|---|---|
| `solve_lower_triangular`, `solve_upper_triangular`, `lu`, `Lu::determinant` | Forward and back substitution; Gaussian elimination with partial pivoting | G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013, ch. 3. |
| `cholesky` | Up-looking (row-by-row) sparse Cholesky. This crate's version does no elimination-tree symbolic analysis. | T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006, ch. 4. |
| `conjugate_gradient` | Conjugate gradient | M. R. Hestenes and E. Stiefel, "Methods of conjugate gradients for solving linear systems," *J. Res. Natl. Bur. Stand.* 49(6):409–436, 1952. [doi:10.6028/jres.049.044](https://doi.org/10.6028/jres.049.044) |
| | | Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, 2003. |
| | | J. R. Shewchuk, "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain," Tech. Rep. CMU-CS-94-125, Carnegie Mellon University, 1994. |
| `matmul`, `par_matmul`, CSR/CSC construction | Row-by-row sparse matrix product; counting-sort transposition into compressed form | F. G. Gustavson, "Two fast algorithms for sparse matrices: multiplication and permuted transposition," *ACM Trans. Math. Softw.* 4(3):250–269, 1978. [doi:10.1145/355791.355796](https://doi.org/10.1145/355791.355796) |
| Matrix Market I/O | Coordinate exchange format | R. F. Boisvert, R. Pozo and K. A. Remington, "The Matrix Market Exchange Formats: Initial Design," NISTIR 5935, NIST, 1996. [doi:10.6028/NIST.IR.5935](https://doi.org/10.6028/NIST.IR.5935) |

### Planned (see [`docs/PowerFlowChecklist.md`](docs/PowerFlowChecklist.md))

| Planned work | Reference |
|---|---|
| Left-looking sparse LU (Phase 3) | J. R. Gilbert and T. Peierls, "Sparse partial pivoting in time proportional to arithmetic operations," *SIAM J. Sci. Stat. Comput.* 9(5):862–874, 1988. [doi:10.1137/0909058](https://doi.org/10.1137/0909058) |
| Refactoring circuit-like matrices (Phase 3) | T. A. Davis and E. Palamadai Natarajan, "Algorithm 907: KLU, a direct sparse solver for circuit simulation problems," *ACM Trans. Math. Softw.* 37(3), article 36, 2010. [doi:10.1145/1824801.1824814](https://doi.org/10.1145/1824801.1824814) |
| Reverse Cuthill–McKee ordering (Phase 3) | E. Cuthill and J. McKee, "Reducing the bandwidth of sparse symmetric matrices," *Proc. 24th ACM National Conference*, pp. 157–172, 1969. |
| | A. George and J. W. H. Liu, *Computer Solution of Large Sparse Positive Definite Systems*, Prentice-Hall, 1981. |
| AMD ordering (Phase 3) | P. R. Amestoy, T. A. Davis and I. S. Duff, "An approximate minimum degree ordering algorithm," *SIAM J. Matrix Anal. Appl.* 17(4):886–905, 1996. [doi:10.1137/S0895479894278952](https://doi.org/10.1137/S0895479894278952) |
| COLAMD ordering (Phase 3) | T. A. Davis, J. R. Gilbert, S. I. Larimore and E. G. Ng, "A column approximate minimum degree ordering algorithm," *ACM Trans. Math. Softw.* 30(3):353–376, 2004. |
| Low-rank updates for contingencies (Phase 5) | J. Sherman and W. J. Morrison, "Adjustment of an inverse matrix corresponding to a change in one element of a given matrix," *Ann. Math. Statist.* 21(1):124–127, 1950. [doi:10.1214/aoms/1177729893](https://doi.org/10.1214/aoms/1177729893) |
| | M. A. Woodbury, "Inverting modified matrices," Memorandum Report 42, Statistical Research Group, Princeton University, 1950. |
| Newton–Raphson power flow and its Jacobian; reference results (Phases 1.5 and 4.5) | R. D. Zimmerman, C. E. Murillo-Sánchez and R. J. Thomas, "MATPOWER: Steady-state operations, planning, and analysis tools for power systems research and education," *IEEE Trans. Power Syst.* 26(1):12–19, 2011. [doi:10.1109/TPWRS.2010.2051168](https://doi.org/10.1109/TPWRS.2010.2051168) |
| Fast-decoupled power flow (Phase 4.5) | B. Stott and O. Alsaç, "Fast decoupled load flow," *IEEE Trans. Power App. Syst.* PAS-93(3):859–869, 1974. |

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
  [their README](powerflow/tests/data/matpower-reference/README.md).

Before relying on these, check them against each project's current `LICENSE` file.
