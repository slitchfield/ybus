#!/usr/bin/env bash
# Downloads two of ARPA-E's GO Competition Challenge 1 dataset bundles into testdata/go-c1/, for
#   the GO reader's tests and survey (docs/GoParserChecklist.md, phase G2). The datasets are on
#   OEDI under CC BY 4.0 (https://doi.org/10.25984/2437761); they're not committed to this
#   repository, so this fetches them instead, as scripts/fetch-matpower.sh fetches MATPOWER.
#
# Usage: scripts/fetch-go-c1.sh
#   Prints the data directory on success: testdata/go-c1, with one subdirectory per bundle. Safe
#   to rerun: an existing, verified extraction is reused. If GO_C1_DIR is set, that existing
#   directory is checked and used instead, and nothing is downloaded.
set -euo pipefail

# Each bundle's archive and its SHA-256. OEDI publishes no checksums, so these were computed
#   from the first download (2026-11-26; 14,691,396 and 10,444,056 bytes) and pinned there: a
#   later download that differs is refused rather than trusted. Bump a pin only after checking
#   why the archive changed.
BASE_URL="https://data.openei.org/files/6153"
BUNDLES=(
    "Challenge_1_Original_Dataset_1_Scenarios ce4e56a23226ca929427544443b94bbeeadd86ff62875ccd6c71f1987ec0e9d1"
    "Challenge_1_Original_Dataset_2_Scenarios b438dab6d51492ff501cea1e9acb5ed3591873b188f8b844acaad0ef75749852"
)

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${ROOT}/testdata/go-c1"

# A directory is usable if it holds GO Challenge 1 data: a network folder's inputfiles.ini and
#   case.rop, and a scenario's case.raw and case.con, somewhere below it
check_dir() {
    local dir="$1"
    if [[ ! -d "${dir}" ]]; then
        echo "error: ${dir} isn't a directory" >&2
        return 1
    fi
    for name in inputfiles.ini case.rop case.raw case.con; do
        if [[ -z "$(find "${dir}" -name "${name}" -type f -print -quit)" ]]; then
            echo "error: ${dir} has no ${name}; is it GO Challenge 1 data?" >&2
            return 1
        fi
    done
}

sha256() {
    if command -v sha256sum > /dev/null; then
        sha256sum "$1" | cut -d' ' -f1
    else
        shasum -a 256 "$1" | cut -d' ' -f1
    fi
}

if [[ -n "${GO_C1_DIR:-}" ]]; then
    check_dir "${GO_C1_DIR}"
    echo "${GO_C1_DIR}"
    exit 0
fi

mkdir -p "${DEST}"
for bundle in "${BUNDLES[@]}"; do
    read -r name expected <<< "${bundle}"
    dir="${DEST}/${name}"

    # The marker is written only after a verified, complete extraction
    if [[ -f "${dir}/.ybus-verified" ]]; then
        continue
    fi

    archive="${DEST}/${name}.zip"
    if [[ ! -f "${archive}" ]] || [[ "$(sha256 "${archive}")" != "${expected}" ]]; then
        echo "Downloading ${name}.zip from ${BASE_URL}" >&2
        curl --fail --location --silent --show-error --output "${archive}.part" \
            "${BASE_URL}/${name}.zip"
        mv "${archive}.part" "${archive}"
    fi

    actual="$(sha256 "${archive}")"
    if [[ "${actual}" != "${expected}" ]]; then
        echo "error: checksum mismatch for ${archive}" >&2
        echo "  expected ${expected}" >&2
        echo "  actual   ${actual}" >&2
        echo "The download is corrupt or the dataset changed; not extracting it." >&2
        exit 1
    fi

    # Each bundle into a directory of its own, so the two can't collide
    rm -rf "${dir}"
    mkdir -p "${dir}"
    unzip -q "${archive}" -d "${dir}"
    check_dir "${dir}"
    touch "${dir}/.ybus-verified"
done
echo "${DEST}"
