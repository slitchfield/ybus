#!/usr/bin/env bash
# Builds each published crate's docs as docs.rs does: on nightly, with every feature, the
#   docsrs cfg that turns on the feature labels (doc_cfg), and DOCS_RS=1 set. Warnings are
#   errors, so a broken intra-doc link or label fails here before it reaches docs.rs
#   (docs/PublishingChecklist.md, P3 and P4).
#
# Usage: scripts/check-docsrs.sh
#   Needs a nightly toolchain: rustup toolchain install nightly. The docs go to target/docsrs,
#   so the nightly build doesn't evict the stable one from target/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

# The crates published, as in docs/PublishingChecklist.md; ybus_ffi isn't
for crate in ybus_sparse ybus_ipm ybus_cuda ybus_powerflow; do
    DOCS_RS=1 cargo +nightly rustdoc -p "${crate}" --lib --all-features \
        --target-dir target/docsrs -- --cfg docsrs -D warnings
done
echo "Documented every published crate as docs.rs would, with no warnings"
