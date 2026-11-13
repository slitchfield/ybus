# Publishing checklist

The roadmap to a first release of the workspace's Rust crates on crates.io, at 0.1.0. Each
unchecked item in phases P1 to P5 is sized to be one commit, except P5's first four, which are the
release's own steps, done by hand. Phase P0's decisions were made with the maintainer on
2026-11-10, and the release tool on 2026-11-11.

## Goals

1. **Publish the Rust crates:** the sparse matrices, the interior-point solvers, the power flow
   and the experimental GPU path, under one `ybus_` family of names.
2. **Releases that are checked before they happen.** A CI job packages and builds every crate as
   crates.io would, and semver checks run before each release.
3. **Out of scope:** `ybus_ffi`. A C library has little use on crates.io; it will be
   distributed as a built library and header when its ABI policy reaches 1.0
   ([`FfiChecklist.md`](FfiChecklist.md)). Python wheels wait for phase F7 there too.

**Conventions carried over from the earlier roadmaps:** one commit per checklist item, on `dev/`
branches, with the PR's description as its report.

## What crates.io requires

- **A free name.** Names are first come, first served, and `-` and `_` count as the same
  character. crates.io's policy discourages publishing placeholder crates just to hold a name.
- **A `description` and a `license`**, as an SPDX expression, or a `license-file`.
- **Every dependency on crates.io**, optional ones included. A path dependency must also give a
  `version`, which is what the published manifest keeps, so crates are published in dependency
  order.
- **A package of at most 10 MB**, compressed.
- **Published versions are permanent.** A version can be yanked, which stops new lockfiles from
  choosing it, but not replaced. Deletion is possible only shortly after a first publish, and
  only for crates no one depends on.
- **Conventions, not rules:** a README, keywords and categories, docs that build on docs.rs, the
  license texts in the package, and semver. Under Cargo's rules a 0.x crate's breaking change
  bumps the minor version (0.1 to 0.2), and anything else bumps the patch.

## Where the crates stand

As of `main` at `1570570`:

| Crate | Name on crates.io | State |
|---|---|---|
| `sparse_mat` | free | Publishable now: `cargo publish --dry-run` passes, 48 files, 156 KiB compressed. It has a README, keywords, categories and links to the license files. Its README is also its crate documentation, and three of its links are relative (`../docs/…`, `../README.md`), which break on docs.rs and crates.io. |
| `interior_point` | **taken**, by an unrelated `interior-point` 1.0.0 | `publish = false`, a description only. Its MIPS port is under MIPS's BSD 3-clause license, which its `license` doesn't say yet, and `LICENSE-MIPS` isn't in its package. |
| `powerflow` | free, but an unrelated `power_flow` exists | `publish = false`, a description only. Its package is 182 files, 1.1 MiB (322 KiB compressed, far below the limit), with 704 KB of MATPOWER reference results in `tests/data`, plus benchmarks and survey examples that need the repository's fetched cases. It re-exports `interior_point`'s types in `powerflow::opf`, and its `cuda` feature needs `ybus_cuda` published. |
| `ybus_cuda` | free | `publish = false`, a description only. It builds anywhere, since CUDA is loaded at run time, so docs.rs can build it. |
| `ybus_ffi` | free | Stays unpublished. |

Path dependencies have no `version` yet, and no crate has docs.rs metadata.

## Phase P0: Decisions

- [x] **Decide** which crates to publish. Decided (2026-11-10): every crate but `ybus_ffi`, which
      keeps `publish = false`.
- [x] **Decide** the names. Decided (2026-11-10): one family, prefixed `ybus_`, with the
      underscores the workspace's crates already use:

      | Now | Published as |
      |---|---|
      | `sparse_mat` | `ybus_sparse` |
      | `interior_point` | `ybus_ipm` |
      | `powerflow` | `ybus_powerflow` |
      | `ybus_cuda` | `ybus_cuda` |

      The prefix marks the crates as one project, as `bevy_` and `polars-` do. It also avoids
      the taken `interior-point` and confusion with `power_flow`. Each crate's directory is
      renamed to match, so a path always names its crate. The plain name `ybus` stays
      unclaimed, since crates.io discourages holding names. A crate that re-exports the family
      can take it later.
- [x] **Decide** whether the experimental GPU path is published. Decided (2026-11-10): yes.
      `ybus_cuda` and `ybus_powerflow`'s `cuda` feature go out at 0.1.0, described as
      experimental, with breaking changes expected.
- [x] **Decide** what `ybus_powerflow`'s package holds. Decided (2026-11-10): the reference
      CSVs and the benchmark material stay out.
      - **Out:** `tests/data`, MATPOWER's results, which are also its data's second channel of
        redistribution; `benches/`; and the survey examples that produce `docs/benchmarks/`'s
        records.
      - **Also out:** the tests, which compare against those results and the fetched cases, so
        they couldn't run from the package anyway.
- [x] **Decide** the first version. Decided (2026-11-10): 0.1.0, with breaking changes expected
      in the 0.x releases that follow, each a minor bump.
- [x] **Decide** how releases are checked. Decided (2026-11-10): a CI job that dry-runs the
      publish of every crate, and semver checks before each release.
- [x] **Decide** the release tool, and whether the crates share one version. Decided
      (2026-11-11): `release-plz`, with one version for all four crates. The comparison is in
      [Release tooling](#release-tooling) below.
      - **The first release is published by hand** with plain Cargo. The maintainer works through
        that first-time workflow before `release-plz` is set up.
      - **`release-plz` runs every release after it,** once the crates exist and trusted
        publishing can be linked to them.

## Phase P1: Names

Each rename covers:
- the package and its directory
- every `use` and dependency on it in the workspace, `ybus_ffi` included
- the scripts, `check.sh`, CI, `deny.toml` and the live docs: the READMEs, `LongTermGoals.md`
  and the roadmaps still open

The completed roadmaps and the dated benchmark records keep the names they were written with,
and the root README notes the renames once.

- [x] Rename `sparse_mat` to `ybus_sparse`.
      - **`scripts/profile_summary.py`** matches each crate's frames under its old name too, so a
        profile recorded before the renames gives the same stage table.
      - **`rustfmt` re-sorted the imports** that the new name moved.
      - **`ybus_sparse/README.md`'s relative links** still resolve, since the directory depth
        is unchanged. P2 makes them absolute.
- [x] Rename `interior_point` to `ybus_ipm`, with `LICENSE-MIPS`'s list of files.
      - **Also renamed:** `AI-DISCLAIMER.md`'s account of the MIPS port, and the MIPS
        reference generator's output path (`scripts/generate-mips-reference.sh` and
        `generate_mips_reference.m`).
      - **`profile_summary.py`** finds the solvers' frames for its KKT assembly stage under
        either name. An old AC profile's stage table is unchanged.
- [x] Rename `powerflow` to `ybus_powerflow`. This includes `ybus_ffi`'s docs, header comment
      and `-p` flags. The C API's names don't change.
      - **The C library:** `ybus.h`, regenerated by cbindgen 0.29.4, changes in three comments
        only, and `symbols.txt` not at all.
      - **Kept:** `powerflow-ffi` in `LongTermGoals.md` and `FfiChecklist.md`, the name the
        long-term goal once suggested, which the checklist quotes.
      - **`profile_summary.py`** strips `ybus_powerflow::` before `powerflow::`, since the old
        name is a suffix of the new one.

## Phase P2: Manifests and licenses

- [ ] Shared dependencies through `[workspace.dependencies]`: each workspace crate once, with
      both `path` and `version`, so the published manifests get versions. Members use
      `ybus_sparse.workspace = true`. The one version all four share moves to
      `[workspace.package]`, as `version.workspace = true` in each. Shared `homepage` and
      `documentation` fields go there too, where they apply to every crate.
- [ ] Metadata for each published crate: its description reviewed for crates.io's search
      results, up to five keywords, and categories from crates.io's list of slugs.
      - `ybus_cuda`'s description says it's experimental.
      - `publish = false` comes off the four crates, so the dry run in P4 can cover them.
        Publishing itself stays a manual step until P5.
- [ ] A README for each crate. `ybus_ipm`, `ybus_powerflow` and `ybus_cuda` don't have one.
      Each says what the crate does, its status, a short example, and links back to the
      repository. Every link is absolute, `ybus_sparse`'s three relative ones included, since a
      README is shown on crates.io and, for `ybus_sparse`, as its docs.rs front page.
- [ ] License files in each package, the way `ybus_sparse` links `LICENSE-MIT` and
      `LICENSE-APACHE` to the root's.
      - **`ybus_ipm`'s license** becomes `(MIT OR Apache-2.0) AND BSD-3-Clause`, since its MIPS
        port is under MIPS's license. It ships `LICENSE-MIPS`.
      - **`AI-DISCLAIMER.md`** goes in each package too, linked the same way, since the root
        README cites it.

## Phase P3: Package contents and docs.rs

- [ ] An `include` list for each crate, so a package holds what builds and documents the
      library and nothing else: `src/`, `Cargo.toml`, the README and the license files.
      - **`ybus_powerflow`:** without `tests/`, `benches/` or `examples/`, as P0 decided.
      - **`ybus_cuda`:** `src/kernels.cu`, which `include_str!` compiles in, is under `src/`.
      - **Targets whose files are left out:** Cargo drops the `[[bench]]` and `[[example]]`
        declarations with a warning. The item makes the package warning-free.

      Each package's file count and size go in the PR.
- [ ] docs.rs: `[package.metadata.docs.rs]` with `all-features = true`, and labels on the
      feature-gated items (`rayon`, `serde`, `cuda`), so a reader sees what each one needs.
      These need nightly's `doc_cfg`, behind `--cfg docsrs` so stable builds don't see it. The
      item builds the docs as docs.rs would, `ybus_cuda` included, and fixes what breaks.
- [ ] A front page for each crate: crate-level docs with a short example that runs as a doc
      test, and for `ybus_powerflow` and `ybus_ipm` a guide to their modules.

## Phase P4: Release checks in CI

- [ ] A **Package** CI job: `cargo publish --workspace --dry-run`.
      - What it does: packages every publishable crate, then builds each from its package, in
        dependency order, as crates.io would.
      - **Before the first release,** unpublished crates resolve against a local registry of
        the ones packaged before them. Tried on a scratch copy of `main`, with versions added
        to the path dependencies: all four crates packaged and built, in 26 s on the
        maintainer's laptop.
      - **`ybus_ffi`** is skipped through its `publish = false`, with no `--exclude` needed.
      - **Packaging warnings** fail the job.
- [ ] **Semver checks** with `cargo-semver-checks`, against each crate's last published version.
      - **Required before each release:** a release isn't published until the check passes
        for the version it bumps to.
      - **A per-PR job,** reporting without failing, is decided with the item. A failing one
        would fail every intended break between releases.
      - **What it misses:** it catches breaks in the API's signatures and types, not in
        behaviour. A behaviour change, such as F4's to `DcModel::study`, still needs a person to
        mark it in the version and the changelog.
      - **When it starts:** it needs a published baseline, so it starts with the first
        release after 0.1.0.
- [ ] A changelog, `CHANGELOG.md` in Keep a Changelog's form, with 0.1.0's entry written by hand.
      `release-plz` writes the entries after it; how they're laid out across the four crates is
      settled with its configuration in P5.
- [ ] `docs/Releasing.md`, the guide to P5's first release by hand, step by step:
      - the account and token
      - the dry run, the publish, the tag and the GitHub release
      - the checks afterwards
      - linking trusted publishing

      The maintainer corrects it while working through those steps. It gains `release-plz`'s
      routine once that's set up.

## Phase P5: The first release

The maintainer does the first four items by hand, following `docs/Releasing.md`, as the
first-time workflow to learn before anything automates it. `release-plz` is set up after them.

- [ ] A crates.io account through GitHub, and an API token scoped to publishing new crates,
      with an expiry. The first publish of a new crate can't use trusted publishing.
- [ ] Publish 0.1.0 with `cargo publish --workspace`, in dependency order: `ybus_sparse`,
      `ybus_ipm`, `ybus_cuda`, `ybus_powerflow`. Tag the release and make a GitHub release
      from the changelog.
- [ ] Trusted publishing for each crate on crates.io, linked to the repository's release
      workflow, so later releases need no stored token. Then revoke the first token.
- [ ] Check what was published:
      - docs.rs built each crate, with the feature labels
      - each crates.io page shows its README with working links
      - a new project that depends on `ybus_powerflow` with `features = ["rayon", "cuda"]`
        builds

      Then the root README's crate table gets crates.io and docs.rs links.
- [ ] `release-plz`, for every release after 0.1.0:
      - **Its configuration:** `release-plz.toml`, with the four crates in one version group and
        the changelog's layout.
      - **A workflow with its two jobs:**
        - `release-pr` opens and updates the release PR.
        - `release` publishes, tags and makes the GitHub release when that PR merges, through
          trusted publishing.
      - **Its `cargo-semver-checks` run** is P4's required check before each release.
      - **The repository setting** that lets GitHub Actions open PRs.

      `docs/Releasing.md` gains the routine: review the release PR, correct its version for
      behaviour breaks the check can't see, edit its changelog, and merge. The first release
      it runs is the test of the setup.

## Phase P6: Later

- **The C library:** built libraries and the header as GitHub release assets, when the ABI
  policy reaches 1.0. `cargo-c` is one way to build them with a pkg-config file.
- **Python wheels** on PyPI, after FFI phase F7's bindings.
- **A `ybus` crate** that re-exports the family, if one name to depend on proves useful.
- **1.0:** what each crate needs before its API is promised stable.

## Release tooling

What each tool automates, as of 2026-11. The project's PRs are the unit of review, and its
commit messages are plain sentences, not Conventional Commits.

| | Plain Cargo | `cargo-release` | `release-plz` | `cargo-smart-release` | `release-please` |
|---|---|---|---|---|---|
| Runs | locally or in CI | locally | in CI (GitHub Action) | locally | in CI (GitHub Action) |
| Version bumps | by hand | one command, updating dependents' requirements | proposed in a release PR, editable before merging | from commits | proposed in a release PR |
| Changelog | by hand | by hand, with search-and-replace of the "Unreleased" heading | generated from commits, in the release PR | generated from commits | generated from commits |
| Semver checks | run separately | as a hook | built in: a detected break raises the bump | none | none |
| Publish order | `--workspace` orders it | orders it | orders it | orders it | doesn't publish; needs a separate step |
| Tags and GitHub releases | by hand | tags | both | both | both |
| Trusted publishing | through `crates-io-auth-action` | no; uses a local token | built in | no | through the separate step |
| Needs Conventional Commits | no | no | for `feat`/`fix` bumps; others are patch bumps | yes | yes |

**What this project's choice turns on:**
- **The release PR model matches the workflow.** `release-plz` puts every release through a PR
  with its version bumps and changelog, reviewed like any change, and publishes when it merges.
- **Conventional Commits barely matter at 0.x.** `release-plz` bumps the patch for a commit it
  can't classify. A break it detects bumps the minor, and a behaviour break the check misses
  can be fixed by editing the PR. So the project keeps its commit style. The changelog lists
  commit subjects; this project's are full sentences, which read well, and can be edited too.
- **`cargo-release` keeps the token on the maintainer's machine**, and its changelog is written
  by hand. That suits rare releases, but each still needs `cargo-semver-checks` run as a hook.
- **Plain Cargo needs no new tool:** `cargo publish --workspace`, `cargo semver-checks`, a
  hand-kept changelog, and `docs/Releasing.md` listing the steps. It's the fallback for any of
  the others, and the way the first publish goes regardless.
- **`cargo-smart-release` and `release-please`** fit worse. The first isn't meant for CI; the
  second needs Conventional Commits and doesn't publish to crates.io itself.

**Decided (2026-11-11), as recommended:**
- **The tool:** `release-plz` from the second release. The first publish uses plain Cargo, since
  trusted publishing can't create a crate. The maintainer works through that first by hand.
- **One shared version for the four crates**, `version.workspace = true`. `ybus_powerflow`
  re-exports `ybus_ipm`'s types, and both use `ybus_sparse`'s, so their breaks travel
  together anyway, as with `bevy`'s and `polars`' crates. One version is simpler
  to read than four. `release-plz` and `cargo-release` both support that.
