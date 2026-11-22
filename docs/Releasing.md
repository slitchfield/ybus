# Releasing

How the four published crates (`ybus_sparse`, `ybus_ipm`, `ybus_cuda` and `ybus_powerflow`)
are released to crates.io. They share one version, set in `[workspace.package]` in the root
`Cargo.toml`, and are always released together.

This guide covers the first release, 0.1.0, which is published by hand with plain Cargo
([`PublishingChecklist.md`](PublishingChecklist.md), P5). It's a first draft: correct it while
working through it, wherever a step turns out to be wrong or missing. Once `release-plz` is set
up, this guide also gains the routine for every release after 0.1.0.

## The rules for every release

- **The checks pass on `main`:** CI's Check, MSRV, Dependencies and Package jobs. The Docs job
  is advisory, but read its result: a failure there is what docs.rs would hit.
- **The semver check passes for the new version.** No release goes out until
  `cargo semver-checks` passes for the version it bumps to. The check needs a published
  baseline, so it starts with the release after 0.1.0. Until `release-plz` runs it, run it by
  hand with `cargo install cargo-semver-checks --locked`, then `cargo semver-checks --workspace`.
  - **At 0.x,** a minor bump (0.1 to 0.2) is the breaking one, and a patch bump must break
    nothing. With one shared version, a break in any crate bumps all four.
  - **The check sees only the API's signatures and types.** A change in behaviour, such as a
    function that now returns different results for the same input, is marked by hand: in the
    version, and in the changelog.
- **`CHANGELOG.md`** has the release's entry, with its version and date.

## The first release, 0.1.0

### 1. Before

- `main` is green, and your clone is at its tip, with no changes (`git status`).
- In `CHANGELOG.md`, replace `YYYY-MM-DD` in the 0.1.0 heading with the release date, then
  commit that through a PR as usual.
- **Check that each name is still free:** `https://crates.io/crates/<name>` for each of
  `ybus_sparse`, `ybus_ipm`, `ybus_cuda` and `ybus_powerflow` should show no crate. None was
  taken on 2026-11-22, checked through crates.io's API.

### 2. The account and the token

1. **Log in to [crates.io](https://crates.io) with GitHub.** It makes the account from the
   GitHub account.
2. **Verify your email address** under Account Settings. crates.io refuses to publish
   without a verified one.
3. **Make an API token** under Account Settings, then API Tokens:
   - **Scope:** `publish-new` only. The first release creates the crates, and later releases
     use trusted publishing (step 7), so the token needs nothing else.
   - **Crates:** limit it to the pattern `ybus_*`.
   - **Expiry:** a short one, such as 7 days. It's revoked in step 7 anyway.
4. **Run `cargo login`** and paste the token when it asks. Cargo stores it in
   `~/.cargo/credentials.toml`.

### 3. The dry run

```sh
scripts/check-package.sh
```

The Package CI job runs the same script. It packages each crate, builds it from its package as
crates.io will, and fails on any warning but the accepted ones: the tests, examples and benches
left out of the packages ([`PublishingChecklist.md`](PublishingChecklist.md), P3), and the dry
run's own notice.

### 4. The publish

```sh
cargo publish --workspace --locked
```

- **The order:** Cargo publishes the crates in dependency order (`ybus_sparse`, then
  `ybus_ipm` and `ybus_cuda`, then `ybus_powerflow`). It waits for each to reach the index
  before publishing the crates that need it. `ybus_ffi` is skipped, since it has
  `publish = false`.
- **A publish can't be undone.** A published version can be yanked, so new projects don't pick
  it up, but it can't be deleted or published again. A mistake is fixed by publishing the next
  version.
- **If it stops partway,** for example on a network error, the crates already published stay
  published. Find which ones are on crates.io. Then publish the rest, in the same order, with
  `cargo publish -p <crate> --locked` for each.

### 5. The tag and the GitHub release

```sh
git tag -a v0.1.0 -m "ybus 0.1.0"
git push origin v0.1.0
```

Copy the 0.1.0 section of `CHANGELOG.md`, from below its heading to the next one, into a file
outside the repository. Then make the release from it:

```sh
gh release create v0.1.0 --verify-tag --title "0.1.0" --notes-file <that file>
```

The changelog's links at its foot point to this tag, so they resolve from now on.

### 6. The checks afterwards

- **docs.rs** has built each crate (`https://docs.rs/<name>`), with the feature labels: for
  example, `rayon` on `DcStudy::par_solve_outages` and `cuda` on `DcStudy::gpu`. A build can
  take a while to start after the publish.
- **Each crates.io page** shows its README, and the README's links work.
- **A new project** builds with the published crates and every feature:

  ```sh
  cargo new --bin ybus-smoke && cd ybus-smoke
  cargo add ybus_powerflow@0.1.0 --features rayon,cuda
  cargo build
  ```

  It should build without the CUDA toolkit, since CUDA is loaded at run time.
- **The root README's crate table** then gets crates.io and docs.rs links, in a PR.

### 7. Trusted publishing, then revoking the token

Trusted publishing lets a GitHub Actions workflow publish without a stored token. crates.io
trusts the workflow's short-lived OpenID Connect token instead. It can't create a crate, which is
why 0.1.0 needed the token.

1. **On each crate's crates.io page,** under Settings, then Trusted Publishing, add a GitHub
   configuration:
   - the repository owner and name, `slitchfield` and `ybus`
   - the release workflow's file name, `release-plz.yml` (decided 2026-11-22), the bare name of
     a file in `.github/workflows/`. The file comes later, in P5's `release-plz` PR.
     - **Before the file exists:** crates.io's RFC (3691) describes matching the workflow's
       token against the configuration at publish time, which suggests a configuration can be
       saved first, but it doesn't say so, and this wasn't confirmed.
     - **If the form refuses it:** keep the token, merge the `release-plz` PR, then configure
       trusted publishing and revoke the token.
   - optionally, a GitHub environment the workflow must run in
2. **Revoke the token** from step 2 under Account Settings, then API Tokens, and run
   `cargo logout`.

## Every release after 0.1.0

[release-plz](https://release-plz.dev) prepares each release as a pull request, and publishes it
when that PR merges. Its configuration is `release-plz.toml`, and its workflow
`.github/workflows/release-plz.yml`.

### Once: the GitHub App

The workflow opens the release PR with a GitHub App's token, since a PR opened with the default
`GITHUB_TOKEN` doesn't trigger CI. Then the release PR runs every CI job, like any other.

1. **Make the App:** in GitHub, under Settings, Developer settings, GitHub Apps, then New
   GitHub App.
   - **Name:** anything, such as `ybus-release-plz`.
   - **Webhook:** off, since the workflow only uses the App's tokens.
   - **Repository permissions:** Contents and Pull requests, read and write. Add Administration,
     read and write, only if tags are protected.
   - **Where it can be installed:** only on this account.
2. **Make a private key** on the App's page, and note its **client ID**.
3. **Install the App** on the `ybus` repository only.
4. **In the repository's settings,** under Secrets and variables, then Actions:
   - **a variable** `RELEASE_PLZ_CLIENT_ID`, the client ID
   - **a secret** `RELEASE_PLZ_PRIVATE_KEY`, the private key's whole `.pem` file
5. **The trusted publishing environment:** if a crate's trusted publishing on crates.io was
   configured with a GitHub environment, the workflow's `release` job needs the same
   `environment:`. Without one, it needs nothing more.
6. **Revoke the first release's API token,** if it isn't already (step 7 above).

### Each release

1. **release-plz opens the release PR** on a push to `main` that changes a crate. It updates the
   PR as more changes merge.
   - **The version** bumps all four crates together (their version group), to the highest next
     version any of them needs. The PR reports cargo-semver-checks' result for each crate, and
     a breaking change raises the bump (at 0.x, the minor version).
   - **`CHANGELOG.md`** gets one generated section for the release, from the commits of all
     four crates.
2. **Review the PR.** CI runs on it like any other.
   - **The version:** correct it by hand for a behaviour break the semver check can't see, such
     as a function that returns different results for the same input. Remember that one break
     in any crate bumps all four.
   - **The changelog:** rewrite the generated section into the usual form: one subsection per
     crate, in Keep a Changelog's headings, describing changes rather than listing commits.
3. **Merge it.** The workflow's `release` job then:
   - publishes the four crates through trusted publishing, in dependency order
   - tags the release `v<version>`
   - makes the GitHub release
4. **Check, as for 0.1.0 (step 6):** docs.rs's builds, the crates.io pages, and a new project
   that depends on the new version.

**The first release PR is the test of this setup.** Some of release-plz's behaviour wasn't
confirmed by its documentation:
- **Versions:** whether it updates the version the crates inherit from `[workspace.package]`.
- **The changelog:** where it inserts the new section in the hand-written `CHANGELOG.md`.
- **The tokens:** whether `persist-credentials: false`, as in its quickstart, suits the App's
  token.

Read that PR's diff for each before merging, and fix `release-plz.toml` if any is wrong.
Nothing is published until it merges.
