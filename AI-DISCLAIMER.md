# AI disclaimer

This repository was developed with substantial help from AI coding tools. This document
describes how AI was used, how sources were attributed, how to report attribution that is
missing or incomplete, and how the licensing relates to AI-assisted development.

## How AI was used

- **Original implementation:** the first version of the library, up to commit `db09ecd`
  (February 2025 to July 2026), was written by the maintainer without AI assistance. It contains
  the dictionary-of-keys matrix type, row iteration, `Display`, transpose and addition.
- **AI-assisted development:** from October 2026 on, most new code, tests and documentation were
  written with **Claude Code**, Anthropic's agentic coding tool, using the Claude Sonnet 5.5 and
  Claude Opus 5.5 models. That covers the large majority of the current source and tests,
  including, in `sparse_mat`:
  - the generic element type, the products, Matrix Market I/O, and the optional `serde` and
    `rayon` features
  - the fixed-pattern compressed matrix type and its assembly maps (Phase 2)
  - the fill-reducing orderings (RCM and AMD), the symbolic analysis, the ordered Cholesky and
    the Gilbert–Peierls LU with refactoring (Phase 3)
  - the complex-capable solvers (Phase 4)

  and the whole of `powerflow`: the case model and MATPOWER reader, the DC power flow, and the AC
  power flow by Newton's method and fast-decoupled (Phases 1.5 and 4.5), with their tests,
  survey tools and benchmark records. [`docs/PowerFlowChecklist.md`](docs/PowerFlowChecklist.md)
  records each phase.
- **Commits:** commits produced with AI assistance carry a `Co-Authored-By: Claude …` trailer.
  Some of them also contain code the maintainer wrote and the AI committed. The trailer marks
  AI involvement, not exclusive AI authorship.
- **Maintainer's role:** the maintainer directed the work, made the design decisions, wrote parts
  of the code, and reviewed every pull request before merging it.
- **Reports:** the description of each AI-assisted pull request is an implementation report,
  recording what was changed, the decisions made, and the mistakes caught along the way. Until
  October 2026 the reports were also committed, under `docs/ai-reports/`, and they remain in the
  commit history.

### How AI-generated code was checked

AI-generated code can be plausible and still wrong, so it was not accepted on trust:
- **Tests:** numerical routines (products, orderings, LU, Cholesky, conjugate gradient,
  Kronecker products) are tested against independent dense reference implementations and
  algebraic identities, such as P A Q = L U, P A Pᵀ = L Lᴴ, (AB)ᵀ = BᵀAᵀ and solution residuals,
  for real and complex element types.
- **Old algorithms kept as references:** when an algorithm was replaced, the old one was kept in
  the tests, and the new one was required to reproduce its results bit for bit where they should
  agree (for example the ordered Cholesky without an ordering, and the new LU's pivots under true
  partial pivoting).
- **Validation against MATPOWER:** the power flows are compared with MATPOWER's own results.
  Committed reference cases must agree within 1e-9 and, for the AC power flow, take MATPOWER's
  number of iterations. Surveys of all 78 MATPOWER cases, recorded in
  [`docs/benchmarks/`](docs/benchmarks/), check every case against MATPOWER, and differences
  were investigated rather than tolerated. A backward-error test checks that the DC solve is
  backward stable.
- **Regressions:** bugs were reproduced with a failing test before they were fixed.
- **CI:** every change runs `./check.sh` in CI (formatting, clippy with warnings as errors, and
  tests with each feature combination), plus a minimum-supported-Rust-version job.
- **Review:** the maintainer reviewed each pull request before merging it.

None of this guarantees the absence of errors. Please report any you find.

## Attribution of sources

Considerable care went into attributing the ideas this code is built on:
- **No copied code:** the sparse matrix algorithms were implemented from published
  descriptions in papers and textbooks, and no source code from other libraries was knowingly
  copied or ported. The DC and AC power flows are the exception: they reproduce MATPOWER's
  results, and were written with reference to MATPOWER's BSD-licensed source code as well as
  its documentation. They are reimplementations in Rust, not translations; the README's
  [References](README.md#references) section describes how they relate to MATPOWER.
- **References:** the [References](README.md#references) section of the README lists the
  publication behind each implemented algorithm, and behind each planned one.
- **Checked citations:** each citation's details (authors, venue, pages, DOI) were checked against
  online records rather than taken from the AI model's memory. Where a record could not be
  confirmed, the citation says so or was left out.
- **Licenses of related software:** the README records the license terms of the most relevant
  reference implementations (SuiteSparse, MATPOWER). That way future work can avoid porting code
  whose license would conflict with this repository's.

### A limitation of AI-generated code

AI models learn from large amounts of public text and code. An AI-written routine was not copied
from any particular source, but it may still resemble an existing implementation, especially for
short, standard algorithms that are written the same way in many places. This repository can't
rule that out completely, which is why the next section exists.

## Reporting missing or incomplete attribution

If you believe any part of this repository should credit you or someone else, please get in
touch. That includes code that resembles yours, an idea that deserves a citation, or a reference
that is wrong.

- **How:** open an issue at <https://github.com/slitchfield/ybus/issues>. If the matter is
  sensitive, open an issue asking for a private contact and the maintainer will follow up.
- **What helps:** the file and lines concerned, the original source (a link, publication or
  repository), and what you think the correct attribution would be.
- **What happens next:** reports are taken seriously. Depending on the case, the response will be
  to add or correct the attribution, or to rewrite or remove the affected code, and to update this
  document if the process itself needs to change.

## Licensing and AI-assisted development

This repository is dual-licensed under the MIT license or the Apache License, Version 2.0, at
your option (see [`LICENSE-MIT`](LICENSE-MIT) and [`LICENSE-APACHE`](LICENSE-APACHE)). This
licensing follows current best practices for AI-assisted software development:

- **A standard permissive license:** the same dual license as Rust itself and most of its
  ecosystem, applied to the whole repository, including AI-assisted portions.
- **No third-party code:** no source code from other projects is included, so there are no
  inherited license obligations to track.
- **Third-party test data:** MATPOWER's case files are not included; a script downloads them
  for testing. The repository does include DC power flow results that MATPOWER computed from
  three public test cases, as expected values for tests. Their sources are credited in
  [`powerflow/tests/data/matpower-reference/README.md`](powerflow/tests/data/matpower-reference/README.md).
- **Compatible dependencies:** every dependency, including transitive ones, is under a compatible
  permissive license (MIT, Apache-2.0, Unlicense or Unicode-3.0).
- **Transparency:** AI involvement is disclosed in this document, in commit trailers and in the
  pull requests' implementation reports, and sources are attributed as described above.

The legal treatment of AI-generated code is still developing, for example the extent to which it
is protected by copyright, which varies by jurisdiction. The licenses grant rights to this
repository's contents to the extent the maintainer holds them. They will be reviewed and updated
as necessary to keep up with best practices and with legal and community guidance as it changes.
Any such change will be recorded in the commit history.

This document describes the maintainer's practices and intentions. It is not legal advice.
